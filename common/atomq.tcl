source VMDextensions.tcl

set molid [molinfo top]
set nf [molinfo $molid get numframes]

# Analysis parameters: minseqsep to not fall in same residue (3 atoms)
set cutoff 4.5
set minSeqSep 3 

# Heavy protein atoms
set sel [atomselect $molid "protein and noh"]

# starting frame
$sel frame 0

# generate all native atom-atom contacts
set nativeContacts [prepareNativeContacts $cutoff $sel]

# map atom indices to residue IDs
array set indexToResid {}

foreach idx [$sel get index] resid [$sel get resid] {
    set indexToResid($idx) $resid
}

# consider only out residue contacts
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

puts "Number of native contacts: $nNative"

# Calculate Q for every trajectory frame
set out [open "qvalue.dat" w]

puts $out "frame\tQ\tnative_contacts"

for {set f 1} {$f < $nf} {incr f} {

    $sel frame $f

    set nRetained [measureNativeContacts \
        $filteredContacts $cutoff $sel]

    set q [expr {double($nRetained) / $nNative}]

    puts $out "$f\t$q\t$nRetained"
}

close $out
$sel delete

puts "Q-value calculation completed."
