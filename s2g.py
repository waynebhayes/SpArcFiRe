import csv
from collections import defaultdict

galaxies = defaultdict(list) # each value in each column is appended to a list
arcs = defaultdict(list) # each value in each column is appended to a list

with open('out/galaxy.csv') as f:
	header = [h.strip() for h in f.next().split(',')] #strip whitespace in headers
	reader = csv.DictReader(f, fieldnames=header) # read rows into a dictionary format
	for row in reader: # read a row as {column1: value1, column2: value2,...}
		for (k,v) in row.items(): # go over each column name and value 
			galaxies[k].append(v) # append the value into the appropriate list based on column name k
print(galaxies["name"])

#Using the galfit.feedme templates, subsitute values from galaxy.csv into a new feeedme.
for i in range(len(galaxies["name"])): 
	with open("galfit_in/" + galaxies["name"][i] + ".feedme", "wt") as fout: #For each input, generate a unique feedme file
		print("Writing to galfit_in/" + galaxies["name"][i] + ".feedme")
		with open("template.feedme", "rt") as fin:
			for line in fin:
				iptSz = galaxies["iptSz"][i][1:-1].split() #Transform string '[xxx yyy]' into array ['xxx', 'yyy']
				line = line.replace('$name', galaxies["name"][i])
				line = line.replace('$integratedmagnitude', str(float(galaxies["bulgeAvgBrt"][i]) * 10))
				line = line.replace('$centerx', galaxies["inputCenterR"][i])
				line = line.replace('$centery', galaxies["inputCenterC"][i])
				line = line.replace('$xmax', iptSz[0])
				line = line.replace('$ymax', iptSz[1])
				line = line.replace('$radius', str(float(galaxies["diskMinAxsLen"][i]) / 2))
				line = line.replace('$axisratio', galaxies["diskAxisRatio"][i])
				line = line.replace('$positionangle', str(float(galaxies["diskMajAxsAngleRadians"][i]) * 180 / 3.14592 + 90))
				fout.write(line)

with open('out/galaxy_arcs.csv') as f:
	reader = csv.DictReader(f) # read rows into a dictionary format
	for row in reader: # read a row as {column1: value1, column2: value2,...}
		for (k,v) in row.items(): # go over each column name and value 
			arcs[k].append(v) # append the value into the appropriate list based on column name k
#print(arcs["gxyName"])
