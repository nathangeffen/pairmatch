#!/bin/bash

time ./partnermatch -p 50000 -a C -r 10 -c 64 -n 128 -t >timing_50000.csv
time ./partnermatch -p 100000 -a C -r 10 -c 64 -n 128 -t  >timing_100000.csv
time ./partnermatch -p 500000 -a C -r 10 -c 64 -n 128 -t  >timing_500000.csv
time ./partnermatch -p 1000000 -a C -r 10 -c 64 -n 128 -t  >timing_1000000.csv
