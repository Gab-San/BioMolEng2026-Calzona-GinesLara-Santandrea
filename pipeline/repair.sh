#!/bin/bash

# Running the script from the pipeline folder
cd -- "$(dirname -- "${BASH_SOURCE[0]}")" || exit 1

# Actual RepairPDB command
foldx --command=RepairPDB \
  --pdb=1COF.pdb \
  --pdb-dir=../common/

# Deleting byproducts of the command
rm -r molecules
rm 1COF_Repair.fxout