#!/bin/bash

time ./partnermatch -p 500000 -a C -r 10 -c 64 -n 50 -t  >timing_K_500000.csv
time ./partnermatch -p 500000 -a C -r 10 -c 64 -n 100 -t  >>timing_K_500000.csv
time ./partnermatch -p 500000 -a C -r 10 -c 64 -n 150 -t  >>timing_K_500000.csv
time ./partnermatch -p 500000 -a C -r 10 -c 64 -n 200 -t  >>timing_K_500000.csv
time ./partnermatch -p 500000 -a C -r 10 -c 64 -n 250 -t  >>timing_K_500000.csv
time ./partnermatch -p 500000 -a C -r 10 -c 64 -n 300 -t  >>timing_K_500000.csv
time ./partnermatch -p 500000 -a C -r 10 -c 64 -n 350 -t  >>timing_K_500000.csv
time ./partnermatch -p 500000 -a C -r 10 -c 64 -n 400 -t  >>timing_K_500000.csv
time ./partnermatch -p 500000 -a C -r 10 -c 64 -n 450 -t  >>timing_K_500000.csv
time ./partnermatch -p 500000 -a C -r 10 -c 64 -n 500 -t  >>timing_K_500000.csv
