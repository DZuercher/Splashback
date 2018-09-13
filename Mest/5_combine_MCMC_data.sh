#! /bin/bash

function combine_data {
    OUTDIR="/work/dominik.zuercher/Output/Mest"
    OutFileName="${OUTDIR}/$1/Data_full.txt"
    rm $OutFileName
    base="${OUTDIR}/$1/data_parts/MCMC_data"
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
types='Planck_PS_21.5_red_spline Planck_PS_21.5_blue_spline'
for type in $types
do
    combine_data $type 
done


