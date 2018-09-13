#!/bin/bash

#Submits SQL jobs to the casjobs server to extract data from the catalog to downloadable files.

NAME_BASE="part"
Q_DIR="./queries"

for i in {0..48};
do
    java -jar casjobs.jar run -t 'PANSTARRS_DR1/720' -n "${NAME_BASE}_$i" -f "${Q_DIR}/query_$i.sql"
    java -jar casjobs.jar extract -b "${NAME_BASE}_$i" -type CSV 
    sleep 15m        
    java -jar casjobs.jar execute -n 'Dropper' "drop table ${NAME_BASE}_$i"
done
