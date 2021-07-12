#!/bin/bash
for i in $(seq $1 $2); do
    if [[ $i -lt 10 ]]; then
        ./runShiftGal.sh /extra/wayne1/preserve/antholn1/SDSS_DR12/00$i /extra/wayne1/preserve/antholn1/SDSS_DR12_Aligned/00$i
    elif [[ $i -lt 100 ]]; then
        ./runShiftGal.sh /extra/wayne1/preserve/antholn1/SDSS_DR12/0$i /extra/wayne1/preserve/antholn1/SDSS_DR12_Aligned/0$i
    else
        ./runShiftGal.sh /extra/wayne1/preserve/antholn1/SDSS_DR12/$i /extra/wayne1/preserve/antholn1/SDSS_DR12_Aligned/$i
    fi
done;
