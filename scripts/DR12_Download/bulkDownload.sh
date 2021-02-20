#!/bin/bash

USAGE=$(cat <<-END
    USAGE: bulkdownlaod.sh outDir inputFile inputFileType
    
    outDir: output directory location, will be created if it does not exist

    inputFile: file that contains the download information either by object (where each
        each line contains the "name", "ra", and "dec" of each object), or by frame (where
        line containes the "frame", "camcol", and "field" values for each frame).  Examples
        of these files are located in /inputExamples. 

    intputFileType: {object, frame} if set to "object", then each line in the input file 
        should contain information for an object, if "frame" then each line should contain 
        information for a frame.
    
END
)

if [[ $# -lt 3 ]] || [[ $1 == "-h" ]]; then
    echo "$USAGE"
    exit
fi

if [[ ! -d $1 ]]; then
    mkdir $1
fi

if [[ $3 == "object" ]]; then
    awk -v out="$1" '{system("python2 -W ignore getGal.py " $2 " " $3 " " $1 " -out_dir " out " -min_num_stars 15")}' $2

elif [[ $3 == "frame" ]]; then
    while read -r line; do
        names=`grep "$line" frame_object_list.tsv | cut -f 4-`
        for objid in $(echo $names); do
            info=`grep $objid "SDSS_name_ra_dec_clean.tsv"`
            name=`echo $info | awk '{print $1}'`
            ra=`echo $info | awk '{print $2}'`
            dec=`echo $info | awk '{print $3}'`
            python2 -W ignore getGal.py $ra $dec $name -out_dir $1 -min_num_stars 15
        done
    done < $2
else
    echo "$3 is not a valid input file type.  Run bulkdownload.sh -h for help."
fi




