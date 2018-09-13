#! /bin/bash

OUTDIR="/work/dominik.zuercher/DataStore/Pan-Starrs/Randoms/"
INDIR="/work/dominik.zuercher/DataStore/Pan-Starrs/Randoms/randoms3/"

OutFileName="${OUTDIR}PS_randoms_new.dat"


rm $OutFileName
i=0 # Reset a counter
for filename in "${INDIR}*";
do
    if [ "$filename"  != "$OutFileName" ] ;      # Avoid recursion 
    then
        cat $filename >>  $OutFileName # Append from the 2nd line each file
        i=$(( $i + 1 ))                        # Increase the counter
	echo "$i"    
    fi
done

