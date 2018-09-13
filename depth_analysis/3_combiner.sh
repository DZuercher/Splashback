OutFileName="bad_pixels.dat"
INDIR="./pix_parts/"

rm $OutFileName
i=0                           # Reset a counter
for filename in "${INDIR}*";
        do
        echo $i;
        if [ "$filename"  != "$OutFileName" ] ;      # Avoid recursion 
        then
                cat $filename >> $OutFileName
                #tail -n+2 $filename >>  $OutFileName # Append from the 2nd line each file
                i=$(( $i + 1 ))                        # Increase the counter
        fi
        done

