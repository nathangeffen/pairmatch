#!/bin/bash

time ./partnermatch -p 500000 -a C -r 10 -c 50 -n 128 -t  >timing_C_500000.csv
time ./partnermatch -p 500000 -a C -r 10 -c 100 -n 128 -t  >>timing_C_500000.csv
time ./partnermatch -p 500000 -a C -r 10 -c 150 -n 128 -t  >>timing_C_500000.csv
time ./partnermatch -p 500000 -a C -r 10 -c 200 -n 128 -t  >>timing_C_500000.csv
time ./partnermatch -p 500000 -a C -r 10 -c 250 -n 128 -t  >>timing_C_500000.csv
time ./partnermatch -p 500000 -a C -r 10 -c 300 -n 128 -t  >>timing_C_500000.csv
time ./partnermatch -p 500000 -a C -r 10 -c 350 -n 128 -t  >>timing_C_500000.csv
time ./partnermatch -p 500000 -a C -r 10 -c 400 -n 128 -t  >>timing_C_500000.csv
time ./partnermatch -p 500000 -a C -r 10 -c 450 -n 128 -t  >>timing_C_500000.csv
time ./partnermatch -p 500000 -a C -r 10 -c 500 -n 128 -t  >>timing_C_500000.csv
