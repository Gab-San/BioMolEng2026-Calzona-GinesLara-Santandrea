#!/bin/bash

pipeline_dir=$PWD

# Assigning k: either the first argument or the first 15 available (default)
k=${1:-"15"}
idx=0

# We can only preprocess as many models have been built by the simulation
folder_count=$(find output/simulated -mindepth 1 -maxdepth 1 -type d | wc -l)

while read -r mutation rest; do
    if [[ $idx -ge $k || $idx -ge $folder_count ]]; then
        break
    fi

    mkdir -p output/postprocessed/${mutation} || exit 1

    (
        cd output/postprocessed/${mutation}

        # Postprocessing analysis
        vmd -dispdev text -e ${pipeline_dir}/../common/analysisMD.tcl \
         -args ../../preprocessed/${mutation}/cof_ionized.psf ../../simulated/${mutation}/cof_eq.pdb.dcd < /dev/null|| exit 1

        # Q-value87 analysis
 #       vmd -dispdev text -e ${pipeline_dir}/../common/residueq.tcl \
 #        -args ../../preprocessed/${mutation}/cof_ionized.psf ../../simulated/${mutation}/cof_eq.pdb.dcd < /dev/null|| exit 1
        
        # Postprocessing Q-value calculation
        vmd -dispdev text -e ${pipeline_dir}/../common/atomq_updated.tcl \
         -args ../../preprocessed/${mutation}/cof_ionized.psf ../../simulated/${mutation}/cof_eq.pdb.dcd < /dev/null|| exit 1

        # Plotting
        poetry run python3 ${pipeline_dir}/../common/plot.py analysis/*.dat -o plots || exit 1
    )

    idx=$((idx + 1))
done < "sorted_positions.txt"
