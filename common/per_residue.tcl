# N.B. ignored atomic contacts get NA, so they get ignored in this analysis

# source VMDextensions.tcl

set molid [molinfo top]
set nf [molinfo $molid get numframes]

# parameters

set cutoff 4.5
set minSeqSep 3
### as in paper
set minNativeAtomContacts 30

# Heavy protein atoms
set sel [atomselect $molid "protein and noh"]

# Check required functions
if {[llength [info commands prepareNativeContacts]] == 0} {
    error "prepareNativeContacts not found. Source VMDextensions.tcl first."
}

# Get residue IDs
set residues [lsort -unique -integer [$sel get resid]]

# Map atom indices to residue IDs
array set indexToResid {}

foreach idx [$sel get index] r [$sel get resid] {
    set indexToResid($idx) $r
}

# native contacts section

$sel frame 0

set nativeAtomPairs [prepareNativeContacts $cutoff $sel]

# Native atomic contact lookup
array set nativeSet {}

# Number of native atomic contacts per residue
array set nativeCount {}

foreach r $residues {
    set nativeCount($r) 0
}

foreach pair $nativeAtomPairs {

    lassign $pair i j

    set ri $indexToResid($i)
    set rj $indexToResid($j)

    # Exclude nearby sequence neighbours
    if {abs($ri - $rj) <= $minSeqSep} {
        continue
    }

    # Canonical atom pair representation
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
set nExcluded 0

puts "Total native atomic contacts: $totalNative"
puts "Native atomic contacts per residue:"

foreach r $residues {

    set n $nativeCount($r)

    if {$n < $minNativeAtomContacts} {
        puts "Residue $r: $n contacts -- EXCLUDED"
        incr nExcluded
    } else {
        puts "Residue $r: $n contacts -- INCLUDED"
    }
}

puts "Excluded residues: $nExcluded"

# one file per residue

set outputDir "qvalue_per_residue"
file mkdir $outputDir

array set outFile {}

foreach r $residues {

    set filename [file join $outputDir "residue_${r}.dat"]

    set outFile($r) [open $filename w]

    puts $outFile($r) \
        "frame\tQ\tretained_native_atomic_contacts\tnative_atomic_contacts"
}

# Q VALUE PER EACH RESIDUE

for {set f 0} {$f < $nf} {incr f} {

    $sel frame $f

    # reset array for each frame
    array unset retainedCount
    array set retainedCount {}

    foreach r $residues {
        set retainedCount($r) 0
    }

    set currentAtomPairs [prepareNativeContacts $cutoff $sel]

    array unset currentSeen
    array set currentSeen {}

    # for each pair, get unique residue
    foreach pair $currentAtomPairs {

        lassign $pair i j

        set ri $indexToResid($i)
        set rj $indexToResid($j)

        # ignore less than 3 atom distance by default
        if {abs($ri - $rj) <= $minSeqSep} {
            continue
        }

        if {$i < $j} {
            set key "$i,$j"
        } else {
            set key "$j,$i"
        }

        if {[info exists nativeSet($key)] &&
            ![info exists currentSeen($key)]} {

            set currentSeen($key) 1

            incr retainedCount($ri)
            incr retainedCount($rj)
        }
    }

    foreach r $residues {

        set nNative $nativeCount($r)
        set nRetained $retainedCount($r)

        # paper definition
        if {$nNative >= $minNativeAtomContacts} {

            set q [expr {double($nRetained) / $nNative}]

            puts $outFile($r) [format "%d\t%.6f\t%d\t%d" \
                $f $q $nRetained $nNative]

        } else {
            # N.B. ignored atomic contacts get NA, so they get ignored in this analysis
            puts $outFile($r) "$f\tNA\t$nRetained\t$nNative"
        }
    }

    if {$f % 100 == 0} {
        puts "Processed frame $f / [expr {$nf - 1}]"
    }
}

foreach r $residues {
    close $outFile($r)
}

$sel delete

puts "Per-residue atomic Q calculation completed."
puts "Output directory: $outputDir"
