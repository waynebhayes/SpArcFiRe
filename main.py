from sys import argv
from spiralgalaxy import SpiralGalaxy
from sparcfire2galfit import Sparcfire2Galfit


def file_input():
	with open(argv[1]) as f:
		tsv_headers = f.readline().strip().split('\t')
		needed_headers = ["name","bulgeAxisRatio","bulgeAvgBrightness","diskAxisRatio","diskMajAxisAngle","totalNumArcs"]
		col_num = [tsv_headers.index(i) for i in needed_headers]
		for line in f:
			line = line.strip().split('\t')
			params = [line(i) for i in len(line) if i in col_num]
			S2G = Sparcfire2Galfit(SpiralGalaxy(params(0),float(params(1)),float(params(2)),float(params(3)),float(params(4))))
			S2G.output_galfit()

def single_input():
	Sparcfire2Galfit(SpiralGalaxy("Default",1.0,1.0,1.0,45.0,0)).output_galfit()

if __name__ == "__main__":
	if len(argv) == 2:
		tsv_input()
	else:
		single_input()
