#! /bin/bash

function combine_data {
    OUTDIR="/work/dominik.zuercher/Output/Mest"
    OutFileName="${OUTDIR}/$1/dev2_full.txt"
    rm $OutFileName
    base="${OUTDIR}/$1/dev2_parts/MCMC_dev2"
    for i in {0..979}
    do
        filename=${base}_${i}.txt
        echo $filename
        if [ "$filename"  != "$OutFileName" ]
        then
            cat $filename >>  $OutFileName
        fi
    done
}

##########################
#Main
##########################
prior='best'
types='Planck_PS_21_no_mc Planck_PS_21.5_no_mc Planck_PS_22_no_mc'
for type in $types
do
    combine_data $type 
done


