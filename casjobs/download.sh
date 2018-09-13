#!/bin/bash

#Downloads extracted files from casjobs Database

OUTDIR="./output"
NAME_BASE="part"

for c in {0..48};
do
     HTTP="ps1images.stsci.edu/datadelivery/outgoing/casjobs/csv/${NAME_BASE}_${c}.csv"
     wget $HTTP -w 2 --random-wait -O "${OUTDIR}/${NAME_BASE}_$c"
done
