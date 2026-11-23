#!/bin/bash

pipeline_dir=$PWD

# Assigning k: either the first argument or the first 15 available (default)
k=${1:-"15"}
idx=0

# Setting the correct processor
case ${3:-0} in
    0) config=cpu.conf ;;
    1) config=gpu.conf ;;
    *) echo "Third argument must be 0 (CPU) or 1 (GPU)." >&2; exit 1 ;;
esac

# We can only preprocess as many models have been built by build_model
folder_count=$(find output/preprocessed -mindepth 1 -maxdepth 1 -type d | wc -l)

while read -r mutation rest; do
    if [[ $idx -ge $k || $idx -ge $folder_count ]]; then
        break
    fi

    # Creating input and output environment variables
    export NAMD_STRUCTURE="${PWD}/output/preprocessed/${mutation}/cof_ionized.psf"
    export NAMD_COORDINATES="${PWD}/output/preprocessed/${mutation}/cof_ionized.pdb"

    # Creating output folder and environment variable
    output_dir="output/simulated/$mutation"
    mkdir -p "$output_dir" || exit 1
    export NAMD_OUTPUT="${PWD}/output/simulated/${mutation}/cof_eq.pdb"

    (
        cd "$output_dir" || exit 1
        # Actual namd3 simulation
        namd3 "+p$2" "${pipeline_dir}/${config}"

    ) || exit 1

    idx=$((idx + 1))
done < "sorted_positions.txt"
