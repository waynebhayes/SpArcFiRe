import csv
import math
from collections import defaultdict

#Parse galaxy.csv file
galaxies = defaultdict(list) # each value in each column is appended to a list
with open('out/galaxy.csv') as f:
	header = [h.strip() for h in f.next().split(',')] #strip whitespace in headers
	reader = csv.DictReader(f, fieldnames=header) # read rows into a dictionary format
	for row in reader: # read a row as {column1: value1, column2: value2,...}
		for (k,v) in row.items(): # go over each column name and value 
			galaxies[k].append(v) # append the value into the appropriate list based on column name k

#Parse galaxy_arcs file
arcs = defaultdict(list) # each value in each column is appended to a list
with open('out/galaxy_arcs.csv') as f:
	header = [h.strip() for h in f.next().split(',')] #strip whitespace in headers
	reader = csv.DictReader(f, fieldnames=header) # read rows into a dictionary format
	for row in reader: # read a row as {column1: value1, column2: value2,...}
		for (k,v) in row.items(): # go over each column name and value 
			arcs[k].append(v) # append the value into the appropriate list based on column name k

#Create list of arcs indicies corresponding to each galaxy
arcsIndex = 0 # Use arcsIndex to traverse arcs["name"] 
arcsKeys = defaultdict(list)
for i in range(len(galaxies["name"])): 
	while arcsIndex < len(arcs["gxyName"]) and arcs["gxyName"][arcsIndex] == galaxies["name"][i]:
		arcsKeys[galaxies["name"][i]].append(arcsIndex)
		arcsIndex += 1

#Using the galfit.feedme templates, subsitute values from galaxy.csv into a new feeedme.
#i  =galaxyIndex, j = arcIndex
for i in range(len(galaxies["name"])): 
	with open("galfit_in/" + galaxies["name"][i] + ".feedme", "wt") as fout, open("template.feedme", "rt") as fin, open("template_arcs.feedme", "rt") as finArcs:
		print "Writing to galfit_in/" + galaxies["name"][i] + ".feedme"
		iptSz = galaxies["iptSz"][i][1:-1].split() #Transform string '[xxx yyy]' into array ['xxx', 'yyy']
		for line in fin:
			line = line.replace('$name', galaxies["name"][i])
			line = line.replace('$integrated_magnitude_disk', str(float(galaxies["bulgeAvgBrt"][i]) * 10))
			line = line.replace('$integrated_magnitude_bulge', str(float(galaxies["bulgeAvgBrt"][i]) * 10))
			line = line.replace('$x_center', galaxies["inputCenterC"][i])
			line = line.replace('$y_center', galaxies["inputCenterR"][i])
			line = line.replace('$radius_disk', str(float(galaxies["diskMajAxsLen"][i]) / 4))
			line = line.replace('$radius_bulge', str(float(galaxies["bulgeMajAxsLen"][i]) / 4))
			line = line.replace('$sersic_index', "1.2")
			line = line.replace('$axis_ratio', galaxies["diskAxisRatio"][i])
			line = line.replace('$position_angle', str(90 - float(galaxies["diskMajAxsAngleRadians"][i]) * 180 / 3.14592)) #Transform Sparcfire output into degrees, then correct for flipped x/y axes by subtracting from 90
			line = line.replace('$x_max', iptSz[0])
			line = line.replace('$y_max', iptSz[1])
			fout.write(line)

		componentNumber = 4 # 1 is the disk, 2 is the bulge, 3 is the sky
		continue
		for j in arcsKeys[galaxies["name"][i]]:
			#if (arcs["arc_length"][arcIndex] > 100): #only take in arcs larger than a certain size
			for line in finArcs:
				line = line.replace('$component_number', str(componentNumber))
				line = line.replace('$x_center', galaxies["inputCenterC"][i])
				line = line.replace('$y_center', galaxies["inputCenterR"][i])
				line = line.replace('$radius', str(float(galaxies["diskMajAxsLen"][i]) / 4))
				line = line.replace('$sersic_index', "0.7")
				line = line.replace('$axis_ratio', galaxies["diskAxisRatio"][i])
				line = line.replace('$position_angle', str(90 - float(galaxies["diskMajAxsAngleRadians"][i]) * 180 / 3.14592)) #Transform Sparcfire output into degrees, then correct for flipped x/y axes by subtracting from 90
				line = line.replace('$inclination', str(math.atan(math.log1p(float(galaxies["diskAxisRatio"][i]))) * 180 / 3.14592))
				line = line.replace('$integrated_magnitude', str(float(galaxies["bulgeAvgBrt"][i]) * 10))
				#line = line.replace('$integrated_magnitude', str(float(arcs["mean_intensity"][j]) * 10))
				line = line.replace('$r_start', str(float(arcs["r_start"][j]) / 2))
				line = line.replace('$r_end', str(float(arcs["r_end"][j]) / 2))
				line = line.replace('$relative_theta_end', str(-float(arcs["relative_theta_end"][j]) * 180 / 3.14592))
				#line = line.replace('$pitch_angle', str(float(galaxies["diskMajAxsAngleRadians"][i]) * 180 / 3.14592 + 90 + float(arcs["pitch_angle"][j])))
				line = line.replace('$pitch_angle', arcs["pitch_angle"][j])
				fout.write(line)
			componentNumber += 1
			break #Only get the first arc for now
