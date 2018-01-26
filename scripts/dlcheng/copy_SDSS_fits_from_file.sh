#!/bin/bash
#input="/home/dlcheng/sparcfire/scripts/dlcheng/P_SP.8.uniq.txt"
input="/home/dlcheng/sparcfire/scripts/dlcheng/P_SP.8-.5.uniq.txt"
START=$1
END=$2
COPYDIR="/home/dlcheng/SDSS/G.in"
COUNTER=0
while IFS= read line
do
    FILE="/extra/wayne1/research/drdavis/SDSS/FITS/hash/${line:(-2)}/${line}/${line}_g.fits.gz"
    if [ -f ${FILE} ]; then
        if [ $COUNTER -ge $START ]; then
            echo "zcat /extra/wayne1/research/drdavis/SDSS/FITS/hash/${line:(-2)}/${line}/${line}_g.fits.gz > $COPYDIR/${line}_g.fits"
            zcat /extra/wayne1/research/drdavis/SDSS/FITS/hash/${line:(-2)}/${line}/${line}_g.fits.gz > $COPYDIR/${line}_g.fits
            COUNTER=$[COUNTER + 1]
            if [ $COUNTER -eq $END ]; #Break out when set amount of galaxies have been read out. 
            then
                break 
            fi
        fi
    else
        echo "File not found: ${FILE}"
    fi
done < "$input"
