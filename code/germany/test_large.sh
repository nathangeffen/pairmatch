#!/bin/bash

# To skip recreation of all agent input files, run this script as:
# $ ./run_tests skip

# Stop on error
set -v
set -x
set -e

echo "Starting time: "
date

./germany_partners-rel -i agent_input_80620000.csv -s $RANDOM -k 30000 -c 10000 -a C >output_fixed_C_80620000_k_30000_c_10000.csv &

./germany_partners-rel -i agent_input_80620000.csv -s $RANDOM -k 35000 -c 10000 -a C >output_fixed_C_80620000_k_35000_c_10000.csv &

./germany_partners-rel -i agent_input_80620000.csv -s $RANDOM -k 40000 -c 10000 -a C >output_fixed_C_80620000_k_40000_c_10000.csv &

./germany_partners-rel -i agent_input_80620000.csv -s $RANDOM -k 45000 -c 6000 -a C >output_fixed_C_80620000_k_45000_c_6000.csv &

./germany_partners-rel -i agent_input_80620000.csv -s $RANDOM -k 45000 -c 8000 -a C >output_fixed_C_80620000_k_45000_c_8000.csv &

./germany_partners-rel -i agent_input_80620000.csv -s $RANDOM -k 45000 -c 10000 -a C >output_fixed_C_80620000_k_45000_c_10000.csv &

./germany_partners-rel -i agent_input_80620000.csv -s $RANDOM -k 50000 -c 10000 -a C >output_fixed_C_80620000_k_50000_c_10000.csv &

./germany_partners-rel -i agent_input_80620000.csv -s $RANDOM -k 50000 -c 20000 -a C >output_fixed_C_80620000_k_50000_c_20000.csv &

# ./germany_partners-rel -i agent_input_80620000.csv -s $RANDOM -k 5000 -a D >output_fixed_D_80620000_k_5000.csv &

# ./germany_partners-rel -i agent_input_80620000.csv -s $RANDOM -k 10000 -a D >output_fixed_D_80620000_k_10000.csv &

wait

echo "Ending at:"
date
