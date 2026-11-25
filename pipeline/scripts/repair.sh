#!/bin/bash

# Checking if the file is already present
if [[ -f "1COF_Repair.pdb" ]]; then
    echo "Repaired file exists"
    exit 0
fi

# Output for keeping additional foldx outputs
output_dir="output/repair_pdb"
mkdir -p "$output_dir" || exit 1

# Switching to output folder
pdb_dir=$PWD 
cd "$output_dir" || exit 1

# Actual RepairPDB command
foldx --command=RepairPDB \
  --pdb=1COF.pdb \
  --pdb-dir=${pdb_dir}/../common/

# Extracting useful output and leaving the rest in runs
mv 1COF_Repair.pdb ${pdb_dir}/1COF_Repair.pdb