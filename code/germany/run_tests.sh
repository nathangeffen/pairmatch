#!/bin/bash

# To skip recreation of all agent input files, run this script as:
# $ ./run_tests skip

# Stop on error
set -v
set -x
set -e

echo "Arguments:" $0 $1

echo "Starting time: "
date

if [ "$1" != "skip" ]
    then
    ./germany_partners-rel -s $RANDOM -a I -n 10000 -o
    mv Initialpop_0_0.csv agent_input_10000.csv

    ./germany_partners-rel -s $RANDOM -a I -n 50000 -o
    mv Initialpop_0_0.csv agent_input_50000.csv

    ./germany_partners-rel -s $RANDOM -a I -n 100000 -o
    mv Initialpop_0_0.csv agent_input_100000.csv

    ./germany_partners-rel -s $RANDOM -a I -n 250000 -o
    mv Initialpop_0_0.csv agent_input_250000.csv

    ./germany_partners-rel -s $RANDOM -a I -n 500000 -o
    mv Initialpop_0_0.csv agent_input_500000.csv

    ./germany_partners-rel -s $RANDOM -a I -n 750000 -o
    mv Initialpop_0_0.csv agent_input_750000.csv

    ./germany_partners-rel -s $RANDOM -a I -n 1000000 -o
    mv Initialpop_0_0.csv agent_input_1000000.csv
fi


# Fixed k,c runs for timing

./germany_partners-rel -s $RANDOM -a 30,IRKWCD -n 10000 -k 75 -c 215 >output_fixed_10000.csv &

./germany_partners-rel -s $RANDOM -a 30,IRKWCD -n 50000 -k 225 -c 565 >output_fixed_50000.csv &

./germany_partners-rel -s $RANDOM -a 30,IRKWCD -n 100000 -k 225 -c 565 >output_fixed_100000.csv &

./germany_partners-rel -s $RANDOM -a 30,IRKWCD -n 250000 -k 740 -c 1280 >output_fixed_250000.csv &

./germany_partners-rel -s $RANDOM -a 30,IRKWCD -n 500000 -k 1620 -c 1700 >output_fixed_500000.csv &

./germany_partners-rel -s $RANDOM -a 30,IRKWCD -n 750000 -k 2170 -c 810 >output_fixed_750000.csv &

./germany_partners-rel -s $RANDOM -a 30,IRKWCD -n 1000000 -k 2800 -c 3900 >output_fixed_1000000.csv &

wait

./germany_partners-rel -i agent_input_10000.csv -s $RANDOM -k 5 -c 5 -vk 10 -vc 10 -vt 30 -r 900 -a C >output_var_C_10000.csv &

./germany_partners-rel -i agent_input_50000.csv -s $RANDOM -k 5 -c 5 -vk 20 -vc 20 -vt 30 -r 900 -a C >output_var_C_50000.csv &

./germany_partners-rel -i agent_input_100000.csv -s $RANDOM -k 50 -c 50 -vk 40 -vc 40 -vt 30 -r 900 -a C >output_var_C_100000.csv &

./germany_partners-rel -i agent_input_250000.csv -s $RANDOM -k 80 -c 80 -vk 60 -vc 60 -vt 30 -r 900 -a C >output_var_C_250000.csv &

./germany_partners-rel -i agent_input_500000.csv -s $RANDOM -k 100 -c 100 -vk 80 -vc 80 -vt 30 -r 900 -a C >output_var_C_500000.csv &

./germany_partners-rel -i agent_input_750000.csv -s $RANDOM -k 90 -c 90 -vk 90 -vc 90 -vt 35 -r 1225 -a C >output_var_C_750000.csv &

./germany_partners-rel -i agent_input_1000000.csv -s $RANDOM -k 100 -c 100 -vk 100 -vc 100 -vt 40 -r 1600 -a C >output_var_C_1000000.csv &

wait

./germany_partners-rel -i agent_input_10000.csv -s $RANDOM -k 5 -vk 10 -r 50 -a D >output_var_D_10000.csv &

./germany_partners-rel -i agent_input_50000.csv -s $RANDOM -k 5 -vk 20 -r 50 -a D >output_var_D_50000.csv &

./germany_partners-rel -i agent_input_100000.csv -s $RANDOM -k 50 -vk 40 -r 50 -a D >output_var_D_100000.csv &

./germany_partners-rel -i agent_input_250000.csv -s $RANDOM -k 80 -vk 60 -r 50 -a D >output_var_D_250000.csv &

./germany_partners-rel -i agent_input_500000.csv -s $RANDOM -k 100 -vk 80 -r 50 -a D >output_var_D_500000.csv &

./germany_partners-rel -i agent_input_750000.csv -s $RANDOM -k 100 -vk 90 -r 50 -a D >output_var_D_750000.csv &

./germany_partners-rel -i agent_input_1000000.csv -s $RANDOM -k 100 -vk 100 -r 50 -a D >output_var_D_1000000.csv &

wait

Rscript analyse_output.R

echo "Ending at:"
date
