set molid [molinfo top]
set nf [molinfo $molid get numframes]

# Analysis parameters
set cutoff 4.5
set minSeqSep 3

# Heavy protein atoms
set sel [atomselect $molid "protein and noh"]

# Check that the required extension function is available
if {[llength [info commands residueContactPairs]] == 0} {
    error "residueContactPairs not found. Source VMDextensions.tcl first."
}

# --------------------------------------------------
# 1. Generate native residue contacts (frame 0)
# --------------------------------------------------

$sel frame 0

set nativePairs [residueContactPairs $cutoff $sel $sel]

# Store native contacts in an array for fast lookup.
# Each pair is represented canonically as r1,r2,
# with r1 < r2, to avoid orientation issues.

catch {array unset nativeSet}
array set nativeSet {}

foreach pair $nativePairs {

    lassign $pair r1 r2

    # Exclude nearby sequence neighbours
    # if {abs($r1 - $r2) > $minSeqSep} {

        if {$r1 < $r2} {
            set key "$r1,$r2"
        } else {
            set key "$r2,$r1"
        }

        set nativeSet($key) 1
    #}
}

set nNative [array size nativeSet]

if {$nNative == 0} {
    error "No native residue contacts found. Check selection and cutoff."
}

puts "Number of native residue contacts: $nNative"

# --------------------------------------------------
# 2. Calculate Q for every trajectory frame
# --------------------------------------------------

set out [open "qvalue_residue.dat" w]

puts $out "# frame Q retained_native_residue_contacts"

for {set f 0} {$f < $nf} {incr f} {

    $sel frame $f

    # Obtain current residue-residue contacts
    set currentPairs [residueContactPairs $cutoff $sel $sel]

    set nRetained 0

    foreach pair $currentPairs {

        lassign $pair r1 r2

        if {$r1 < $r2} {
            set key "$r1,$r2"
        } else {
            set key "$r2,$r1"
        }

        # Count only contacts present in the native reference
        if {[info exists nativeSet($key)]} {
            incr nRetained
        }
    }

    # Fraction of native residue contacts retained
    set q [expr {double($nRetained) / $nNative}]

    puts $out [format "%d %.6f %d" $f $q $nRetained]
}

close $out
$sel delete

puts "Residue-based Q-value calculation completed."
