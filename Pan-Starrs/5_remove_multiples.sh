#! /bin/bash
DIR="/work/dominik.zuercher/DataStore/Pan-Starrs"
uniq -f5 "${DIR}/PS_catalog.csv" > "${DIR}/PS_catalog_deep_fin_stripped.csv"
echo "multiples printed"
