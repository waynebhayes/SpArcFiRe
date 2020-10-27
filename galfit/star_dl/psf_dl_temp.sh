#!bin/bashuthor: Vivek Gite under GPL v2.0+
# ------------------------------------------
INPUT=10_star_dl.csv
mkdir psfield_files
OLDIFS=$IFS
IFS=','
url="http://das.sdss.org/raw"
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read name objid ra dec run rerun camcol field rowc colc type magu magg magr magi magz pmagu pmagg pmagr pmagi pmagz
do
	input_str="${url}/${run}/${rerun}/objcs/${camcol}/psField-00" 
	if [[ ${#run} -eq 3 ]]; then run="0${run}"; fi
	if [[ ${#field} -eq 2 ]]; then field="0${field}"; fi
	
	wget "${input_str}${run}-${camcol}-0${field}.fit" -O "psfield_files/${name}_${rowc%.*}_${colc%.*}_psfield.fit"

done < $INPUT
IFS=$OLDIFS
rm "psfield_files/name_rowC_colC_psfield.fit"


