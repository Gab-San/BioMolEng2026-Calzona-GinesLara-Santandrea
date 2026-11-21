#!/bin/bash

export PREPROCESS_INPUT_PDB="../common/1COF.pdb" PREPROCESS_OUTPUT_PDB="cofp.pdb"
vmd -dispdev text -e /dev/stdin <<'EOF'
if {[catch {
    mol new $env(PREPROCESS_INPUT_PDB) type pdb waitfor all
    set protein [atomselect top protein]
    if {[$protein num] == 0} {error "No protein atoms found in input PDB"}
    $protein writepdb $env(PREPROCESS_OUTPUT_PDB)
    $protein delete
} message]} {
    puts stderr "Preprocessing failed: $message"
    exit 1
}
exit 0
EOF

