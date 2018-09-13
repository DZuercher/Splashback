#! /bin/bash
INDIR="/work/dominik.zuercher/DataStore/Pan-Starrs/PS_parts2/parts"
OUTDIR="/work/dominik.zuercher/DataStore/Pan-Starrs"

OutFileName="${OUTDIR}/PS_star_catalog.csv" # Fix the output name
rm $OutFileName # Reset a counter
for filename in "${INDIR}/*"
do        
    if [ "$filename"  != "$OutFileName" ] ;      # Avoid recursion 
    then
        tail -n+2 $filename >>  $OutFileName # Append from the 2nd line each file
    fi
done

