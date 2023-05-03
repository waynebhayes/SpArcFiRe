#!bin/bash/

long_lat="name,ra,dec\n" # "name,ra,dec\n"
output_file="long_lat.csv"

rm "$output_file"

# Relative paths for now
files=(../sparcfire-in/*.fits) #($( ls ./sparcfire-in/*.fits ))

for galaxy in ${files[@]}
do
		
	galaxy_num=(${galaxy/.fits/})
	long_lat="$long_lat${galaxy_num##*/},"

	# dumping the header into a variable
	crvals=$(hexdump -e '80/1 "%_p" "\n"' $galaxy | grep -i "crval")
	crvals=${crvals//[[:space:]]/}
	# echo $crvals

	long=${crvals%%/*}
	# long=${long//[[:space:]]/}
	long=${long//CRVAL1=/}
	long_lat=$long_lat"${long},"

	lat=${crvals##*CRVAL2=}
	# lat=${lat//[[:space:]]/}
	lat=${lat%%/*}"\n" #,'
	#echo $lat
	long_lat=$long_lat$lat
done

echo -e "${long_lat%'\n'}" >> $output_file

#iconv -f ascii -t utf8 $output_file > ${output_file/.txt/.csv}

