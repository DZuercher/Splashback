#! /bin/bash

function combine_chains {
    OUTDIR="/work/dominik.zuercher/Output/Mest"
    OutFileName="${OUTDIR}/$1/chainstate_full.txt"
    rm $OutFileName
    base="${OUTDIR}/$1/chainstates/chainstate_${1}_${prior}"
    if [ $mc -eq 1 ]
    then
        echo "bla"
        base=$base'_mc'
    fi
    for i in {1..100}
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

types='Planck_PS_21.5_blue_spline Planck_PS_21.5_red_spline'
prior='best'
mc=1

for type in $types
do
    combine_chains $type 
done


