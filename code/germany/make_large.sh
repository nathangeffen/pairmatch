#!/bin/bash

# To skip recreation of all agent input files, run this script as:
# $ ./run_tests skip

# Stop on error
set -v
set -x
set -e

echo "Starting time: "
date

Rscript step2.R 80620000 1
mv agent_list.csv agent_input_80620000.csv

echo "Ending at:"
date
