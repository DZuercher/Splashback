#!/bin/bash
#Generates configuration files for the SQL queries to access the casjobs Database
#Splits the query into 49 pieces
#Requires dust data as dust_data object in MyDB
OUTDIR="./queries"

for i in {0..48};
do
    min=$((69+$i*3))
    max=$((69+($i+1)*3))
    
    if [ "$i" -eq "48" ]
    then
	min=213
	max=215
    fi

    echo "select objid,projectionID,skyCellID,primaryDetection,bestDetection,ira,iraErr,idec,idecErr,rPSFMag,rPSFMagErr,rKronMag,rKronMagErr,gPSFMag,gPSFMagErr,gKronMag,gKronMagErr,iPSFMag,iPSFMagErr,iKronMag,iKronMagErr,iinfoflag, iinfoflag2, iinfoflag3, SDSS_i_avg
    from PanSTARRS_DR1.StackObjectThin o,
    Mydb.dust_data d
    where o.objid between 1E15*$min and 1E15*$max
    and o.primaryDetection=1
    and o.bestDetection=1
    and o.projectionID=d.procell
    and d.skycell=o.skyCellID
    and o.iKronMag between 0.0 and 22.0+d.SDSS_i_avg
    into mydb.deep_$i" > ${OUTDIR}/query_$i.sql

done
