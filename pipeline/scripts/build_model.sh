#!/bin/bash

pdb_dir=$PWD

# Assigning k: either the first argument or 15 (default)
k=${1:-15}
idx=0

while read -r mutation rest; do
    # Considering only the top-k mutations
    if [[ $idx -ge $k ]]; then
        break
    fi

    output_dir="output/build_model/$mutation"
    mkdir -p "$output_dir" || exit 1
    printf '%s;\n' "${mutation//_/,}" > "$output_dir/individual_list.txt"

    (
        cd "$output_dir" || exit 1
        foldx --command=BuildModel \
            --pdb=1COF_Repair.pdb \
            --pdb-dir="$pdb_dir/" \
            --mutant-file=individual_list.txt < /dev/null
    ) || exit 1

    idx=$((idx + 1))
done < "sorted_positions.txt"
