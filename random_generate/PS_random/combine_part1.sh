#! /bin/bash

OFFSET=635
SIZE=28
INDIR="/work/dominik.zuercher/DataStore/Pan-Starrs/Randoms/randoms/finished/"
OUTDIR="/work/dominik.zuercher/DataStore/Pan-Starrs/Randoms/randoms/parts/"

TOTAL=$(((2643 - $OFFSET) / $SIZE))
MIN=$(( $OFFSET + $OMPI_COMM_WORLD_RANK * $TOTAL ))
MAX=$(( $OFFSET + ($OMPI_COMM_WORLD_RANK + 1) * $TOTAL - 1 ))

rm "${OUTDIR}*"
if [ "$OMPI_COMM_WORLD_RANK" == 27 ];
        then
        MAX=2643
        fi
echo "$MIN"
echo "$MAX"
for cell in $(seq $MIN $MAX);
do

OutFileName="${OUTDIR}random_cell_${cell}.dat"                 

rm $OutFileName
i=0 # Reset a counter
for filename in "${INDIR}random_cell_${cell}_rank*.dat";
        do
        if [ "$filename"  != "$OutFileName" ] ;      # Avoid recursion 
                then
                cat $filename >>  $OutFileName # Append from the 2nd line each file
                i=$(( $i + 1 ))                        # Increase the counter
        fi
        done

echo "${cell} done"

done

