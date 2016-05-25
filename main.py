import csv
import sys

if __name__ == "__main__":
	with open(sys.argv[1], "rt") as fin:
		reader = csv.reader(fin)
		sfdict = {rows[0].strip():rows[1].strip() for rows in reader}
		headers = ["name", "bulgeAvgBrt", "standardizedCenterR", "standardizedCenterC", "diskMinAxsLen", "diskAxisRatio", "diskMajAxsAngleRadians"]

	with open(sfdict[headers[0]] + ".feedme", "wt") as fout:
		with open("template.feedme", "rt") as fin:
			for line in fin:
				fout.write(line.replace('$name', sfdict[headers[0]]).replace('$integratedmagnitude', str(float(sfdict[headers[1]]) * 10)).replace('$centerx', sfdict[headers[2]]).replace('$centery', sfdict[headers[3]]).replace('$radius', str(float(sfdict[headers[4]]) / 2)).replace('$axisratio', sfdict[headers[5]]).replace('$positionangle', str(float(sfdict[headers[6]]) * 180 / 3.14592)))