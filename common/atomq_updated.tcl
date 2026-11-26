# Native atomic-contact fraction Q (VMD)
# Usage: vmd -dispdev text -e atomq_updated.tcl -args PSF DCD [TIMESTEP DCDFREQ OUTPUT_DIR STEP TIME_UNIT]
# Defaults: 2.0 fs, 5000 integration steps, analysis, 1 frame, ps.
# Requires VMDextensions.tcl (prepareNativeContacts / measureNativeContacts).

# Keep the analysis inside a procedure so Tcl errors stop the entire analysis,
# rather than allowing VMD's -e command reader to continue with later commands.
proc run_atomq {arguments} {
    set nargs [llength $arguments]
    if {$nargs < 2 || $nargs > 7} {
        error {Usage: vmd -dispdev text -e atomq_updated.tcl -args PSF DCD [TIMESTEP DCDFREQ OUTPUT_DIR STEP TIME_UNIT]}
    }

    set psf_file [lindex $arguments 0]
    set dcd_file [lindex $arguments 1]
    set timestep 2.0
    set dcdfreq 5000
    set output_dir analysis
    set step 1
    set time_unit ps

    foreach var {timestep dcdfreq output_dir step time_unit} idx {2 3 4 5 6} {
        if {$nargs > $idx} {set $var [lindex $arguments $idx]}
    }

    foreach f [list $psf_file $dcd_file] {
        if {![file isfile $f] || ![file readable $f]} {
            error "Cannot read input file: $f"
        }
    }
    if {![string is double -strict $timestep] || !($timestep > 0 && $timestep < Inf)} {
        error "TIMESTEP must be a finite positive number (fs)"
    }
    foreach var {dcdfreq step} {
        set value [set $var]
        if {![string is integer -strict $value] || $value <= 0} {
            error "$var must be a positive integer"
        }
    }
    if {$time_unit ni {ps ns}} {error "TIME_UNIT must be ps or ns"}
    if {[file exists $output_dir] && ![file isdirectory $output_dir]} {
        error "Output path exists but is not a directory: $output_dir"
    }

    # Original atomq.tcl needs helper procedures from this extension.
    # Search the working directory first, then next to this Tcl script.
    set ext "../../../../common/VMDextensions.tcl"
    if {![file readable $ext]} {
        # This fallback works when Tcl exposes the current script path.
        set current_script [info script]
        if {$current_script ne ""} {
            set ext [file join [file dirname [file normalize $current_script]] VMDextensions.tcl]
        }
    }
    if {![file isfile $ext] || ![file readable $ext]} {
        error "Cannot read VMDextensions.tcl (place it in the working directory or beside this script)"
    }
    source $ext
    foreach helper {prepareNativeContacts measureNativeContacts} {
        if {[llength [info commands $helper]] == 0} {
            error "VMDextensions.tcl does not define $helper"
        }
    }

    puts "Info) Arguments passed: $arguments"
    puts "Info) Timestep: $timestep fs; DCDfreq: $dcdfreq; stride: $step; output: $output_dir; time unit: $time_unit"

    set molid [mol new $psf_file type psf waitfor all]
    mol addfile $dcd_file type dcd waitfor all molid $molid
    set nf [molinfo $molid get numframes]
    if {$nf < 1} {error "No trajectory frames loaded"}

    # Preserve the original Q definition: heavy protein atoms, 4.5 A cutoff,
    # residue sequence separation strictly greater than 3, frame 0 reference.
    set cutoff 4.5
    set minSeqSep 3
    set sel [atomselect $molid "protein and noh"]
    if {[$sel num] == 0} {error "Protein heavy-atom selection is empty"}
    $sel frame 0

    set nativeContacts [prepareNativeContacts $cutoff $sel]
    array set indexToResid {}
    foreach idx [$sel get index] resid [$sel get resid] {
        set indexToResid($idx) $resid
    }
    set filteredContacts {}
    foreach pair $nativeContacts {
        lassign $pair i j
        set ri $indexToResid($i)
        set rj $indexToResid($j)
        if {abs($ri - $rj) > $minSeqSep} {
            lappend filteredContacts $pair
        }
    }
    set nNative [llength $filteredContacts]
    if {$nNative == 0} {
        error "No native contacts found. Check selection and reference structure."
    }
    puts "Info) Loaded $nf frames; number of native contacts: $nNative"
    puts "Warning: assumes the protein is whole under periodic boundary conditions."

    file mkdir $output_dir
    set output_file [file join $output_dir qvalue.dat]
    set out [open $output_file w]
    fconfigure $out -encoding utf-8 -translation lf
    puts $out "Time ($time_unit)\tQ\tNative contacts"

    # Preserve original starting frame 1; frame 0 defines the reference.
    # Step is the stride of sampled trajectory frames.
    for {set f 1} {$f < $nf} {incr f $step} {
        $sel frame $f
        set nRetained [measureNativeContacts $filteredContacts $cutoff $sel]
        set q [expr {double($nRetained) / $nNative}]
        set time_ps [expr {$f * $dcdfreq * $timestep / 1000.0}]
        if {$time_unit eq "ns"} {
            set time [expr {$time_ps / 1000.0}]
        } else {
            set time $time_ps
        }
        puts $out [join [list $time $q $nRetained] "\t"]
    }
    close $out
    $sel delete
    puts "Q-value calculation completed: [file normalize $output_file]"
}

if {[catch {run_atomq $argv} err options]} {
    puts stderr "Error) $err"
    puts stderr [dict get $options -errorinfo]
    exit
}
exit
