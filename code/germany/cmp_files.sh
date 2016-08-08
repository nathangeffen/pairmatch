#!/bin/bash
cat $1 | sort --field-separator="," --key=1 -g | awk --field-separator="," '$5 == 1 { print $0 }' >f1.tmp
cat $2 | sort --field-separator="," --key=1 -g | awk --field-separator="," '$5 == 1 { print $0 }' >f2.tmp
meld f1.tmp f2.tmp
