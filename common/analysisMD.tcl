# Cofilin trajectory analysis (VMD)
# Usage: vmd -dispdev text -e analysisMD.tcl -args PSF DCD [TIMESTEP DCDFREQ OUTPUT_DIR STEP TIME_UNIT]
# Defaults: 2.0 fs, 5000 integration steps, analysis, 1 frame, ps.

# =================================================
# INPUT PARAMETERS
# =================================================

set nargs [llength $argv]

proc fatal {message} {
    puts stderr "Error) $message"
    exit
}

if {$nargs < 2 || $nargs > 7} {
    puts stderr {Usage: vmd -dispdev text -e analysisMD.tcl -args PSF DCD [TIMESTEP DCDFREQ OUTPUT_DIR STEP TIME_UNIT]}
    exit 1
}

set psf_file [lindex $argv 0]
set dcd_file [lindex $argv 1]
set timestep 2.0
set dcdfreq 5000
set output_dir analysis
set step 1
set time_unit ps

foreach var {timestep dcdfreq output_dir step time_unit} idx {2 3 4 5 6} {
    if {$nargs > $idx} {set $var [lindex $argv $idx]}
}

foreach f [list $psf_file $dcd_file] {
    if {![file isfile $f] || ![file readable $f]} {
        fatal "Cannot read input file: $f"
    }
}

if {![string is double -strict $timestep] || !($timestep > 0 && $timestep < Inf)} {
    fatal "TIMESTEP must be a finite positive number (fs)"
}

foreach var {dcdfreq step} {
    set value [set $var]
    if {![string is integer -strict $value] || $value <= 0} {
        fatal "$var must be a positive integer"
    }
}

if {$time_unit ni {ps ns}} {fatal "TIME_UNIT must be ps or ns"}
file mkdir $output_dir

puts "Info) Arguments passed: $argv"
puts "Info) Timestep: $timestep fs; DCDfreq: $dcdfreq; stride: $step; output: $output_dir; time unit: $time_unit"

# =================================================
# UTILITY FUNCTIONS
# =================================================

proc write_header {output_file args} {
    puts $output_file [join $args "\t"]
}

proc write_row {output_file args} {
    puts $output_file [join $args "\t"]
}

proc open_table {filename args} {
    global output_dir
    set out [open [file join $output_dir $filename] w]
    fconfigure $out -encoding utf-8 -translation lf
    write_header $out {*}$args
    return $out
}

proc time_ps {frame} {
    global timestep dcdfreq
    return [expr {$frame * $dcdfreq * $timestep / 1000.0}]
}

proc time_ns {frame} {
    return [expr {[time_ps $frame] / 1000.0}]
}

proc frame_time {frame} {
    global time_unit
    if {$time_unit eq "ns"} {return [time_ns $frame]}
    return [time_ps $frame]
}

# =================================================
# ANALYSIS SWITCHES (all enabled; set to 0 to disable)
# =================================================

set do_rmsd 1
set do_rmsf 1
set do_hbonds 1
set do_sasa 0
set do_rgyr 0
set do_beta 1
set do_alpha 1

# =================================================
# LOAD MOLECULE
# =================================================

set molid [mol new $psf_file type psf waitfor all]
mol addfile $dcd_file type dcd waitfor all molid $molid
set nf [molinfo $molid get numframes]
if {$nf < 1} {error "No trajectory frames loaded"}
set last [expr {$nf - 1}]
set time_label "Time ($time_unit)"

# =================================================
# ATOM SELECTIONS
# =================================================

set all [atomselect $molid all]
set bb [atomselect $molid "protein and backbone"]
set bb0 [atomselect $molid "protein and backbone" frame 0]
set protein [atomselect $molid protein]
set ca [atomselect $molid "protein and name CA"]
set numres [$ca num]
if {$numres == 0 || [$bb num] == 0} {error "Protein selection is empty"}

puts "Loaded $nf frames; selected $numres C-alpha atoms"
puts "Warning: assumes the protein is whole under periodic boundary conditions."

# =================================================
# ALIGNMENT: only frames that will be analyzed, to frame 0.
# =================================================

puts "Align"
for {set i 0} {$i < $nf} {incr i $step} {
    $all frame $i
    $bb frame $i
    set trans_mat [measure fit $bb $bb0]
    $all move $trans_mat
}

# =================================================
# RMSD
# =================================================

if {$do_rmsd} {
    puts "Calculate RMSD"
    set out [open_table rmsd.dat $time_label "RMSD (\u00c5)"]
    for {set i 0} {$i < $nf} {incr i $step} {
        $bb frame $i
        set rmsd [measure rmsd $bb0 $bb]
        write_row $out [frame_time $i] $rmsd
    }
    close $out
}

# =================================================
# RMSF (one row per C-alpha, not per frame)
# =================================================

if {$do_rmsf} {
    puts "Calculate RMSF"
    set out [open_table rmsfCA.dat Residue "RMSF (\u00c5)"]
    set rmsfCA [measure rmsf $ca first 0 last $last step $step]
    foreach resid [$ca get resid] rmsf $rmsfCA {
        write_row $out $resid $rmsf
    }
    close $out
}

# =================================================
# Intraprotein hydrogen bonds (backbone and side chains)
# =================================================

if {$do_hbonds} {
    puts "Calculate Hbonds"
    set out [open_table hbonds.dat $time_label "Hydrogen bonds"]
    for {set i 0} {$i < $nf} {incr i $step} {
        $protein frame $i
        set hb [measure hbonds 3.0 20 $protein]
        write_row $out [frame_time $i] [llength [lindex $hb 0]]
    }
    close $out
}

# =================================================
# SASA (water probe radius 1.4 Angstrom)
# =================================================

if {$do_sasa} {
    puts "Calculate SASA"
    set out [open_table sasa.dat $time_label "SASA (\u00c5\u00b2)"]
    for {set i 0} {$i < $nf} {incr i $step} {
        $protein frame $i
        write_row $out [frame_time $i] [measure sasa 1.4 $protein]
    }
    close $out
}

# =================================================
# Mass-weighted radius of gyration
# =================================================

if {$do_rgyr} {
    puts "Calculate Rgyr"
    set out [open_table rgyr.dat $time_label "Rgyr (\u00c5)"]
    for {set i 0} {$i < $nf} {incr i $step} {
        $protein frame $i
        write_row $out [frame_time $i] [measure rgyr $protein weight mass]
    }
    close $out
}

# =================================================
# Secondary structure: recalculate once per sampled frame for both metrics.
# =================================================

if {$do_beta || $do_alpha} {
    puts "Calculate secondary structure"
    if {$do_beta} {set beta_out [open_table beta.dat $time_label "Beta fraction" "Beta count"]}
    if {$do_alpha} {set alpha_out [open_table alpha.dat $time_label "Alpha fraction" "Alpha count"]}

    for {set i 0} {$i < $nf} {incr i $step} {
        molinfo $molid set frame $i
        $ca frame $i

        mol ssrecalc $molid
        set st [$ca get structure]
        if {$do_beta} {
            set Beta [regexp -all -- E $st]
            write_row $beta_out [frame_time $i] [expr {double($Beta)/$numres}] $Beta
        }
        if {$do_alpha} {
            set Alpha [regexp -all -- H $st]
            write_row $alpha_out [frame_time $i] [expr {double($Alpha)/$numres}] $Alpha
        }
    }
    if {$do_beta} {close $beta_out}
    if {$do_alpha} {close $alpha_out}
}

foreach sel [list $all $bb $bb0 $protein $ca] {$sel delete}
puts "Analysis complete: [file normalize $output_dir]"

exit
