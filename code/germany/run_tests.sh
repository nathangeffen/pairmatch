#!/bin/bash

# Stop on error
set -e

# Create input files

Rscript step2b_simmatch_nathan_edits.R 10000
mv agent_list.csv agent_input_10000.csv

Rscript step2b_simmatch_nathan_edits.R 50000
mv agent_list.csv agent_input_50000.csv

Rscript step2b_simmatch_nathan_edits.R 100000
mv agent_list.csv agent_input_100000.csv

Rscript step2b_simmatch_nathan_edits.R 500000
mv agent_list.csv agent_input_500000.csv

Rscript step2b_simmatch_nathan_edits.R 1000000
mv agent_list.csv agent_input_1000000.csv

./germany_partners-rel -i agent_input_10000.csv -s $RANDOM -r 30 -k 350 -c 100 >output_fixed_10000.csv &

./germany_partners-rel -i agent_input_50000.csv -s $RANDOM -r 30 -k 700 -c 100 >output_fixed_50000.csv &

./germany_partners-rel -i agent_input_100000.csv -s $RANDOM -r 30 -k 1400 -c 100 >output_fixed_100000.csv &

./germany_partners-rel -i agent_input_500000.csv -s $RANDOM -r 30 -k 2000 -c 100 >output_fixed_500000.csv &

./germany_partners-rel -i agent_input_1000000.csv -s $RANDOM -r 30 -k 3500 -c 100 >output_fixed_1000000.csv &

wait

./germany_partners-rel -i agent_input_10000.csv -s $RANDOM -k 5 -c 5 -vk 10 -vc 10 -vt 30 -r 900 -a C >output_var_C_10000.csv &

./germany_partners-rel -i agent_input_50000.csv -s $RANDOM -k 5 -c 5 -vk 20 -vc 20 -vt 30 -r 900 -a C >output_var_C_50000.csv &

./germany_partners-rel -i agent_input_100000.csv -s $RANDOM -k 50 -c 50 -vk 40 -vc 40 -vt 30 -r 900 -a C >output_var_C_100000.csv &

./germany_partners-rel -i agent_input_500000.csv -s $RANDOM -k 100 -c 100 -vk 80 -vc 80 -vt 30 -r 900 -a C >output_var_C_500000.csv &

./germany_partners-rel -i agent_input_1000000.csv -s $RANDOM -k 100 -c 100 -vk 100 -vc 100 -vt 40 -r 1600 -a C >output_var_C_1000000.csv &

wait

./germany_partners-rel -i agent_input_10000.csv -s $RANDOM -k 5 -vk 10 -r 50 -a D >output_var_D_10000.csv &

./germany_partners-rel -i agent_input_50000.csv -s $RANDOM -k 5 -vk 20 -r 50 -a D >output_var_D_50000.csv &

./germany_partners-rel -i agent_input_100000.csv -s $RANDOM -k 50 -vk 40 -r 50 -a D >output_var_D_100000.csv &

./germany_partners-rel -i agent_input_500000.csv -s $RANDOM -k 100 -vk 80 -r 50 -a D >output_var_D_500000.csv &

./germany_partners-rel -i agent_input_1000000.csv -s $RANDOM -k 100 -vk 100 -r 50 -a D >output_var_D_1000000.csv &

wait
