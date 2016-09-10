#!/bin/bash

make release
time ./partnermatch -p 20000 -r 10 -i 20 -c 100 -n 200 -t >timing_20000.csv

