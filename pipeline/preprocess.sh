#!/bin/bash

pipeline_dir=$PWD

# Assigning k: either the first argument or the first 15 available (default)
k=${1:-"15"}
idx=0

# We can only preprocess as many models have been built by build_model
folder_count=$(find output/build_model -mindepth 1 -maxdepth 1 -type d | wc -l)

while read -r mutation rest; do
    if [[ $idx -ge $k || $idx -ge $folder_count ]]; then
        break
    fi
    
    # Setting input file from build_model outputs
    export PREPROCESS_INPUT_PDB="output/build_model/${mutation}/1COF_Repair_1.pdb"

    # Setting output directory and file
    output_dir="output/preprocessed/${mutation}"
    mkdir -p "$output_dir" || exit 1
    export PREPROCESS_OUTPUT_PDB="${output_dir}/cofp.pdb"

    # Exporting protein and removing noise atoms
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
    # Creating psf
    (
        cd $output_dir
        # Without /dev/null, the rest is interpreted as additional input
        vmd -dispdev text -e "${pipeline_dir}/cof.pgn" < /dev/null
    )

    idx=$((idx + 1))
done < "sorted_positions.txt"



