#!/bin/bash

# Getting parameters
k=15
processors=2
gpu=0

while getopts ":k:p:g" option; do
    case "$option" in
        k) k=$OPTARG ;;
        p) processors=$OPTARG ;;
        g) gpu=1 ;;
        :) echo "Option -$OPTARG needs a value." >&2; exit 1 ;;
        \?) echo "Unknown option: -$OPTARG" >&2; exit 1 ;;
    esac
done

bash ./repair.sh
bash ./position_scan.sh
bash ./build_model.sh $k
bash ./preprocess.sh $k
bash ./simulate.sh $k $processors $gpu