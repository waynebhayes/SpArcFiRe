#Sparcfire to Galfit
#	Description: Takes in Sparcfire *_galaxy.csv and *_galaxy_arcs.csv files, and outputs Galfit *.feedme files.
#	Usage: python main.py 1237666338116993064r_galaxy.csv 1237666338116993064r_galaxy_arcs.csv
#	Example: python main.py 1237666338116993064r_galaxy.csv 1237666338116993064r_galaxy_arcs.csv

import csv
import sys
import math

if __name__ == "__main__":
	with open(sys.argv[1], "rt") as fin:
		reader = csv.reader(fin)
		sfdict = {rows[0].strip():rows[1].strip() for rows in reader}
		headers = ["name", "bulgeAvgBrt", "standardizedCenterR", "standardizedCenterC", "diskMinAxsLen", "diskAxisRatio", "diskMajAxsAngleRadians"]
		
	with open(sys.argv[2], "rt") as fin:
		reader = csv.reader(fin)
		sfarcsdict = dict()
		for rows in reader:
			sfarcsdict[rows[1]] = rows

	with open(sfdict[headers[0]] + ".feedme", "wt") as fout:
		with open("template.feedme", "rt") as fin:
			for line in fin:
				fout.write(line.replace('$name', sfdict[headers[0]]).replace('$integratedmagnitude', str(float(sfdict[headers[1]]) * 10)).replace('$centerx', sfdict[headers[2]]).replace('$centery', sfdict[headers[3]]).replace('$radius', str(float(sfdict[headers[4]]) / 2)).replace('$axisratio', sfdict[headers[5]]).replace('$positionangle', str(float(sfdict[headers[6]]) * 180 / 3.14592 + 90)))
		with open("template_arcs.feedme", "rt") as fin_arcs:
			template_arcs = fin_arcs.readlines();
			object_number = 3 #1: disk, 2: sky
			for key in sfarcsdict:
				arc_headers = sfarcsdict[key]
				if key.isdigit() and (float(arc_headers[6]) * 180 / 3.14592 + 90) > 90: #if header is not column labels and arc has rotation greater than 90, add to feedme
					for line in template_arcs:
						fout.write(line.replace('$name', sfdict[headers[0]]).replace('$integratedmagnitude', str(float(sfdict[headers[1]]) * 10)).replace('$centerx', sfdict[headers[2]]).replace('$centery', sfdict[headers[3]]).replace('$radius', str(float(sfdict[headers[4]]) / 2)).replace('$axisratio', sfdict[headers[5]]).replace('$positionangle', str(float(sfdict[headers[6]]) * 180 / 3.14592 + 90)).replace('$alenRank', str(object_number)).replace('$pitch_angle', str(float(sfdict[headers[6]]) * 180 / 3.14592 + 90 + float(arc_headers[3]))).replace('$inclination', str(math.atan(math.log1p(float(sfdict[headers[5]]))) * 180 / 3.14592)).replace('$relative_theta_end', str(float(arc_headers[6]) * 180 / 3.14592 + 90)).replace('$r_start', arc_headers[7]).replace('$r_end', arc_headers[8]))
				object_number += 1
				break
