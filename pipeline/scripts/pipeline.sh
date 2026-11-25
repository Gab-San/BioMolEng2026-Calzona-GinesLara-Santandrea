#!/bin/bash

if [[ ${1-} == --help ]]; then
    cat <<'EOF'
Usage: bash ./scripts/pipeline.sh [-k COUNT] [-p THREADS] [-g] [-m]
Run from the pipeline/ directory.
  -k COUNT    Process up to COUNT mutation entries (default: 15).
  -p THREADS  Use THREADS for NAMD (default: 2).
  -g          Use the GPU configuration (default: CPU).
  -m          Skip the position scan and use sorted_positions.txt.
  -w          Run ONLY the wildtype simulation.
  --help      Show this help and exit.
EOF
    exit 0
fi

# Getting parameters
k=15
processors=2
gpu=0
manual_positions=0
wildtype=0

while getopts ":k:p:gmw" option; do
    case "$option" in
        k) k=$OPTARG ;;
        p) processors=$OPTARG ;;
        g) gpu=1 ;;
        m) manual_positions=1 ;;
        w) wildtype=1 ;;
        :) echo "Option -$OPTARG needs a value." >&2; exit 1 ;;
        \?) echo "Unknown option: -$OPTARG" >&2; exit 1 ;;
    esac
done


# Running pipeline.
bash scripts/repair.sh

# Wildtype simulation case: creating folder and files for wildtype
if [[ $wildtype -eq 1 ]]; then
    touch sorted_positions.txt
    printf 'WT\n' > sorted_positions.txt
    mkdir -p "output/build_model/WT" || exit 1
    cp 1COF_Repair.pdb output/build_model/WT/1COF_Repair_1.pdb
    k=1
elif [[ $manual_positions -eq 0 ]]; then
    # Scanning positions if mutations are not chosen manually.
    bash scripts/position_scan.sh
    bash scripts/build_model.sh $k
else
    echo "Manual positions chosen, position scan skipped"
    bash scripts/build_model.sh $k
fi

bash scripts/preprocess.sh $k
bash scripts/simulate.sh $k $processors $gpu
bash scripts/postprocess.sh $k
