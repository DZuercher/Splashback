OutFileName="bad_mask_pixels.dat"
INDIR="./mask_pix_parts/"
rm $OutFileName
i=0                           # Reset a counter
for filename in "${INDIR}*";
        do
        echo $i;
        if [ "$filename"  != "$OutFileName" ] ;      # Avoid recursion 
        then
                cat $filename >> $OutFileName
                i=$(( $i + 1 ))                        # Increase the counter
        fi
        done

