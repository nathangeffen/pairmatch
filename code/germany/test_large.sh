#!/bin/bash

# To skip recreation of all agent input files, run this script as:
# $ ./run_tests skip

# Stop on error
set -v
set -x
set -e

echo "Starting time: "
date

./germany_partners-rel -i agent_input_80620000.csv -s $RANDOM -k 30000 -c 10000 -a CD >output_fixed_80620000.csv


echo "Ending at:"
date
