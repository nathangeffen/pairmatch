#!/bin/bash

set -x

make release

./partnermatch -p 20000 -r 1 -i 20 -n 200 -c 100 -s $RANDOM >output20k_1.csv &

./partnermatch -p 20000 -r 1 -i 20 -n 200 -c 100 -s $RANDOM >output20k_2.csv &

./partnermatch -p 20000 -r 1 -i 20 -n 200 -c 100 -s $RANDOM >output20k_3.csv &

./partnermatch -p 20000 -r 1 -i 20 -n 200 -c 100 -s $RANDOM >output20k_4.csv &

./partnermatch -p 20000 -r 1 -i 20 -n 200 -c 100 -s $RANDOM >output20k_5.csv &

./partnermatch -p 20000 -r 1 -i 20 -n 200 -c 100 -s $RANDOM >output20k_6.csv &

./partnermatch -p 20000 -r 1 -i 20 -n 200 -c 100 -s $RANDOM >output20k_7.csv &

./partnermatch -p 20000 -r 1 -i 20 -n 200 -c 100 -s $RANDOM >output20k_8.csv &

./partnermatch -p 20000 -r 1 -i 20 -n 200 -c 100 -s $RANDOM >output20k_9.csv &

./partnermatch -p 20000 -r 1 -i 20 -n 200 -c 100 -s $RANDOM >output20k_10.csv &

./partnermatch -p 20000 -r 1 -i 20 -n 200 -c 100 -s $RANDOM >output20k_11.csv &

./partnermatch -p 20000 -r 1 -i 20 -n 200 -c 100 -s $RANDOM >output20k_12.csv &

./partnermatch -p 20000 -r 1 -i 20 -n 200 -c 100 -s $RANDOM >output20k_13.csv &

./partnermatch -p 20000 -r 1 -i 20 -n 200 -c 100 -s $RANDOM >output20k_14.csv &

./partnermatch -p 20000 -r 1 -i 20 -n 200 -c 100 -s $RANDOM >output20k_15.csv &

./partnermatch -p 20000 -r 1 -i 20 -n 200 -c 100 -s $RANDOM >output20k_16.csv &



wait

make release_attract

for i in `seq 0 0.25 1`;
do
      j=$(echo 1 - $i|bc)
      echo "# Running partnermatch with attract=$i and reject=$j"
      ./partnermatch -p 20000 -r 4 -i 1 -n 200 -c 100 -A $i -R $j > output20k_ra_$i.csv &
done

wait
