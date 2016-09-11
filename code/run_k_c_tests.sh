#!/bin/bash

set -x

make release

# ./partnermatch -p 20000 -a C -r 3 -i 10 -n 50 -c 100 -s $RANDOM >output20k_k50.csv &
# ./partnermatch -p 20000 -a C -r 3 -i 10 -n 100 -c 100 -s $RANDOM >output20k_k100.csv &
# ./partnermatch -p 20000 -a C -r 3 -i 10 -n 150 -c 100 -s $RANDOM >output20k_k150.csv &
# ./partnermatch -p 20000 -a C -r 3 -i 10 -n 200 -c 100 -s $RANDOM >output20k_k200.csv &
# ./partnermatch -p 20000 -a C -r 3 -i 10 -n 250 -c 100 -s $RANDOM >output20k_k250.csv &
# ./partnermatch -p 20000 -a C -r 3 -i 10 -n 300 -c 100 -s $RANDOM >output20k_k300.csv &
# ./partnermatch -p 20000 -a C -r 3 -i 10 -n 350 -c 100 -s $RANDOM >output20k_k350.csv &
# ./partnermatch -p 20000 -a C -r 3 -i 10 -n 400 -c 100 -s $RANDOM >output20k_k400.csv &
# ./partnermatch -p 20000 -a C -r 3 -i 10 -n 450 -c 100 -s $RANDOM >output20k_k450.csv &
# ./partnermatch -p 20000 -a C -r 3 -i 10 -n 500 -c 100 -s $RANDOM >output20k_k500.csv &
# ./partnermatch -p 20000 -a C -r 3 -i 10 -n 550 -c 100 -s $RANDOM >output20k_k550.csv &
# ./partnermatch -p 20000 -a C -r 3 -i 10 -n 600 -c 100 -s $RANDOM >output20k_k600.csv &
# ./partnermatch -p 20000 -a C -r 3 -i 10 -n 650 -c 100 -s $RANDOM >output20k_k650.csv &
# ./partnermatch -p 20000 -a C -r 3 -i 10 -n 700 -c 100 -s $RANDOM >output20k_k700.csv &
# ./partnermatch -p 20000 -a C -r 3 -i 10 -n 750 -c 100 -s $RANDOM >output20k_k750.csv &
# ./partnermatch -p 20000 -a C -r 3 -i 10 -n 800 -c 100 -s $RANDOM >output20k_k800.csv &
# ./partnermatch -p 20000 -a C -r 3 -i 10 -n 850 -c 100 -s $RANDOM >output20k_k850.csv &
# ./partnermatch -p 20000 -a C -r 3 -i 10 -n 900 -c 100 -s $RANDOM >output20k_k900.csv &
# ./partnermatch -p 20000 -a C -r 3 -i 10 -n 950 -c 100 -s $RANDOM >output20k_k950.csv &
# ./partnermatch -p 20000 -a C -r 3 -i 10 -n 1000 -c 100 -s $RANDOM >output20k_k1000.csv &

# wait


./partnermatch -p 20000 -a C -r 3 -i 10 -n 200 -c 1 -s $RANDOM >output20k_C1.csv &
./partnermatch -p 20000 -a C -r 3 -i 10 -n 200 -c 50 -s $RANDOM >output20k_C50.csv &
./partnermatch -p 20000 -a C -r 3 -i 10 -n 200 -c 100 -s $RANDOM >output20k_C100.csv &
./partnermatch -p 20000 -a C -r 3 -i 10 -n 200 -c 150 -s $RANDOM >output20k_C150.csv &
./partnermatch -p 20000 -a C -r 3 -i 10 -n 200 -c 200 -s $RANDOM >output20k_C200.csv &
./partnermatch -p 20000 -a C -r 3 -i 10 -n 200 -c 250 -s $RANDOM >output20k_C250.csv &
./partnermatch -p 20000 -a C -r 3 -i 10 -n 200 -c 300 -s $RANDOM >output20k_C300.csv &
./partnermatch -p 20000 -a C -r 3 -i 10 -n 200 -c 350 -s $RANDOM >output20k_C350.csv &
./partnermatch -p 20000 -a C -r 3 -i 10 -n 200 -c 400 -s $RANDOM >output20k_C400.csv &
./partnermatch -p 20000 -a C -r 3 -i 10 -n 200 -c 450 -s $RANDOM >output20k_C450.csv &
./partnermatch -p 20000 -a C -r 3 -i 10 -n 200 -c 500 -s $RANDOM >output20k_C500.csv &
./partnermatch -p 20000 -a C -r 3 -i 10 -n 200 -c 550 -s $RANDOM >output20k_C550.csv &
./partnermatch -p 20000 -a C -r 3 -i 10 -n 200 -c 600 -s $RANDOM >output20k_C600.csv &
./partnermatch -p 20000 -a C -r 3 -i 10 -n 200 -c 650 -s $RANDOM >output20k_C650.csv &
./partnermatch -p 20000 -a C -r 3 -i 10 -n 200 -c 700 -s $RANDOM >output20k_C700.csv &
./partnermatch -p 20000 -a C -r 3 -i 10 -n 200 -c 750 -s $RANDOM >output20k_C750.csv &
./partnermatch -p 20000 -a C -r 3 -i 10 -n 200 -c 800 -s $RANDOM >output20k_C800.csv &
./partnermatch -p 20000 -a C -r 3 -i 10 -n 200 -c 850 -s $RANDOM >output20k_C850.csv &
./partnermatch -p 20000 -a C -r 3 -i 10 -n 200 -c 900 -s $RANDOM >output20k_C900.csv &
./partnermatch -p 20000 -a C -r 3 -i 10 -n 200 -c 950 -s $RANDOM >output20k_C950.csv &

wait # There are only 20 processors! Don't want the processes to contest for the CPUs.

./partnermatch -p 20000 -a C -r 3 -i 10 -n 200 -c 1000 -s $RANDOM >output20k_C1000.csv

wait
