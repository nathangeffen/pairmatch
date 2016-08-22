#!/bin/bash

make release_attract

for i in `seq 0 0.1 1`;
do
    j=$(echo 1 - $i|bc)
    echo Running partnermatch with attract=$i and reject=$j
    ./partnermatch -p 16384 -A $i -R $j -a "RNWCB"
done
