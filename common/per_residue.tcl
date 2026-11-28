# Per-residue native atomic-contact Q (VMD)
# Usage: vmd -dispdev text -e per_residue_q_pipeline.tcl -args PSF DCD [TIMESTEP DCDFREQ OUTPUT_DIR STEP TIME_UNIT]
# Defaults: 2.0 fs, 5000 integration steps, analysis, 1 frame, ps.
# Requires VMDextensions.tcl (prepareNativeContacts).

proc run_per_residue_q {arguments} {
    set nargs [llength $arguments]
    if {$nargs < 2 || $nargs > 7} {
        error {Usage: vmd -dispdev text -e per_residue_q_pipeline.tcl -args PSF DCD [TIMESTEP DCDFREQ OUTPUT_DIR STEP TIME_UNIT]}
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

    set ext "../../../../common/VMDextensions.tcl"
    if {![file readable $ext]} {
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

    # parameters

    set cutoff 4.5
    set minSeqSep 3
    ### as in paper
    set minNativeAtomContacts 30
    set sel [atomselect $molid "protein and noh"]
    if {[$sel num] == 0} {error "Protein heavy-atom selection is empty"}
    set residues [lsort -unique -integer [$sel get resid]]

    array set indexToResid {}
    foreach idx [$sel get index] resid [$sel get resid] {
        set indexToResid($idx) $resid
    }

    # native contacts section

    $sel frame 0
    set nativeAtomPairs [prepareNativeContacts $cutoff $sel]
    if {[llength $nativeAtomPairs] == 0} {
        error "No native atomic contacts found. Check the reference structure and cutoff."
    }

    array set nativeSet {}
    array set nativeCount {}
    foreach resid $residues {set nativeCount($resid) 0}

    foreach pair $nativeAtomPairs {
        lassign $pair i j
        set ri $indexToResid($i)
        set rj $indexToResid($j)
        if {abs($ri - $rj) <= $minSeqSep} {continue}

        if {$i < $j} {
            set key "$i,$j"
        } else {
            set key "$j,$i"
        }

        # Avoid counting duplicate atom pairs
        if {![info exists nativeSet($key)]} {
            set nativeSet($key) 1

            # Each atom contact contributes to both residues
            incr nativeCount($ri)
            incr nativeCount($rj)
        }
    }

    set totalNative [array size nativeSet]
    if {$totalNative == 0} {
        error "No native atomic contacts remain after applying the sequence-separation filter."
    }
    puts "Info) Loaded $nf frames; native atomic contacts: $totalNative"

    # one file per residue
    #
    file mkdir $output_dir
    set residue_output_dir [file join $output_dir qvalue_per_residue]
    file mkdir $residue_output_dir

    array set outFile {}
    foreach resid $residues {
        set output_file [file join $residue_output_dir "residue_${resid}.dat"]
        set outFile($resid) [open $output_file w]
        fconfigure $outFile($resid) -encoding utf-8 -translation lf
        puts $outFile($resid) [join [list Frame "Time ($time_unit)" Q Retained_native_atomic_contacts Native_atomic_contacts] "\t"]
    }

    # Q VAULE PER EACH RESIDUE
    # consider step defined as argument (or default)

    for {set f 1} {$f < $nf} {incr f $step} {
        $sel frame $f
        set currentAtomPairs [prepareNativeContacts $cutoff $sel]

        array unset retainedCount
        array set retainedCount {}
        foreach resid $residues {set retainedCount($resid) 0}

        array unset currentSeen
        array set currentSeen {}

        foreach pair $currentAtomPairs {
            lassign $pair i j
            set ri $indexToResid($i)
            set rj $indexToResid($j)
            # ignore less than 3 atom distance by default
            if {abs($ri - $rj) <= $minSeqSep} {continue}

            if {$i < $j} {
                set key "$i,$j"
            } else {
                set key "$j,$i"
            }

            if {[info exists nativeSet($key)] && ![info exists currentSeen($key)]} {
                set currentSeen($key) 1
                incr retainedCount($ri)
                incr retainedCount($rj)
            }
        }

        set time_ps [expr {$f * $dcdfreq * $timestep / 1000.0}]
        if {$time_unit eq "ns"} {
            set time [expr {$time_ps / 1000.0}]
        } else {
            set time $time_ps
        }

        foreach resid $residues {
            set nNative $nativeCount($resid)
            set nRetained $retainedCount($resid)
            # as in paper definition
            if {$nNative >= $minNativeAtomContacts} {
                set q_text [format "%.6f" [expr {double($nRetained) / $nNative}]]
            } else {
                set q_text NA
            }
            # N.B. ignored atomic contacts get NA, so they get ignored in this analysis
            puts $outFile($resid) [join [list $f $time $q_text $nRetained $nNative] "\t"]
        }

        if {$f == 1 || $f % 100 == 0} {
            puts "Info) Processed frame $f / [expr {$nf - 1}]"
        }
    }

    foreach resid $residues {close $outFile($resid)}
    $sel delete
    puts "Per-residue atomic Q calculation completed: [file normalize $residue_output_dir]"
}

if {[catch {run_per_residue_q $argv} err options]} {
    puts stderr "Error) $err"
    puts stderr [dict get $options -errorinfo]
    exit
}
exit
