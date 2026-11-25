#!/bin/bash

if [[ ${1-} == --help ]]; then
    cat <<'EOF'
Usage: bash ./scripts/pipeline.sh [-k COUNT] [-p THREADS] [-g] [-m]
Run from the pipeline/ directory.
  -k COUNT    Process up to COUNT mutation entries (default: 15).
  -p THREADS  Use THREADS for NAMD (default: 2).
  -g          Use the GPU configuration (default: CPU).
  -m          Skip the position scan and use sorted_positions.txt.
  --help      Show this help and exit.
EOF
    exit 0
fi

# Getting parameters
k=15
processors=2
gpu=0
manual_positions=0

while getopts ":k:p:gm" option; do
    case "$option" in
        k) k=$OPTARG ;;
        p) processors=$OPTARG ;;
        g) gpu=1 ;;
        m) manual_positions=1 ;;
        :) echo "Option -$OPTARG needs a value." >&2; exit 1 ;;
        \?) echo "Unknown option: -$OPTARG" >&2; exit 1 ;;
    esac
done

# Running pipeline
bash scripts/repair.sh

if [[ $manual_positions -eq 0 ]]; then
    bash scripts/position_scan.sh
else
    echo "Manual positions chosen, position scan skipped"
fi

bash scripts/build_model.sh $k
bash scripts/preprocess.sh $k
bash scripts/simulate.sh $k $processors $gpu
bash scripts/postprocess.sh $k
