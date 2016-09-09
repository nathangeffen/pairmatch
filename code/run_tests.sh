#!/bin/bash

set -x

make release

./partnermatch -p 5000 -r 1 -i 20 -n 200 -c 100 -s $RANDOM -b >output5k_1.csv &

./partnermatch -p 5000 -r 1 -i 20 -n 200 -c 100 -s $RANDOM -b >output5k_2.csv &

./partnermatch -p 5000 -r 1 -i 20 -n 200 -c 100 -s $RANDOM -b >output5k_3.csv &

./partnermatch -p 5000 -r 1 -i 20 -n 200 -c 100 -s $RANDOM -b >output5k_4.csv &

./partnermatch -p 5000 -r 1 -i 20 -n 200 -c 100 -s $RANDOM -b >output5k_5.csv &

./partnermatch -p 5000 -r 1 -i 20 -n 200 -c 100 -s $RANDOM -b >output5k_6.csv &

./partnermatch -p 5000 -r 1 -i 20 -n 200 -c 100 -s $RANDOM -b >output5k_7.csv &

./partnermatch -p 5000 -r 1 -i 20 -n 200 -c 100 -s $RANDOM -b >output5k_8.csv &

./partnermatch -p 5000 -r 1 -i 20 -n 200 -c 100 -s $RANDOM -b >output5k_9.csv &

./partnermatch -p 5000 -r 1 -i 20 -n 200 -c 100 -s $RANDOM -b >output5k_10.csv &

./partnermatch -p 5000 -r 1 -i 20 -n 200 -c 100 -s $RANDOM -b >output5k_11.csv &

./partnermatch -p 5000 -r 1 -i 20 -n 200 -c 100 -s $RANDOM -b >output5k_12.csv &

./partnermatch -p 5000 -r 1 -i 20 -n 200 -c 100 -s $RANDOM -b >output5k_13.csv &

./partnermatch -p 5000 -r 1 -i 20 -n 200 -c 100 -s $RANDOM -b >output5k_14.csv &

./partnermatch -p 5000 -r 1 -i 20 -n 200 -c 100 -s $RANDOM -b >output5k_15.csv &

./partnermatch -p 5000 -r 1 -i 20 -n 200 -c 100 -s $RANDOM -b >output5k_16.csv &



wait

make release_attract

for i in `seq 0 0.25 1`;
do
      j=$(echo 1 - $i|bc)
      echo "# Running partnermatch with attract=$i and reject=$j"
      ./partnermatch -p 5000 -r 4 -i 1 -n 200 -c 100 -A $i -R $j -b > output5k_ra_$i.csv &
done

wait
