#!/bin/bash

# Create dictionary of 1-letter amino acids
declare -A amino_acid_codes

# Note: the read command reads from the < argument at the end
# of the loop
while read -r three_letter one_letter; do
    amino_acid_codes["$three_letter"]=$one_letter
done < "amino_acid_codes.txt"

# For every amino acid, try all the possible mutations
current_pos=0
while read -r record atom_number atom_name residue_name chain position rest; do
    # The .pdb lists every atom in every amino acid, but 
    # the position number writes the amino acid number in the chain.
    # Thus only when the position changes we need to run PositionScan
#     if [[ $record == ATOM && $position != "$current_pos" ]]; then
    if [[ $record == ATOM && $position != "$current_pos" && ($position == 19 || $position == 20) ]]; then
        # Selecting correct letter
        current_pos=$position
        amino_acid=${amino_acid_codes[$residue_name]-}

        # Output for keeping additional foldx outputs
        output_dir="output/position_scan/runs/${position}"
        mkdir -p "$output_dir" || exit 1

        pdb_dir=$PWD
        (   
            cd "$output_dir" || exit 1

            foldx --command=PositionScan \
                --pdb=1COF_Repair.pdb \
                --pdb-dir="$pdb_dir/" \
                --out-pdb=false \
                --positions="${amino_acid}A${position}a"

            # Moving \Delta \Delta g file to position_scan
            mv PS_1COF_Repair_scanning_output.txt ../../${position}.txt
        )
    fi
done < "1COF_Repair.pdb"
