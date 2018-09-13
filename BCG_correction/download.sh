#! /bin/bash

HTTP_FILE = 'PS_https.dat'
NAME_FILE = 'PS_names.dat'
OUTDIR = './orignial_picuters'


LINES=$(wc -l < ../${NAME_FILE})

for (( c=1; c<=$LINES; c++ ))
do 
    HTTP=`cat ../${HTTP_FILE} | head -$c | tail -1`
    NAME=`cat ../${NAME_FILE} | head -$c | tail -1`
    wget $HTTP -w 2 --random-wait -O ${OUTDIR}/${NAME} 
done
