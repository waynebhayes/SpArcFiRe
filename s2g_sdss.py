#!/usr/bin/env python

import os
import csv
import math
from collections import defaultdict
from s2g import s2g
import subprocess
import compare

def read_csv(path):
	dictionary = defaultdict(list) # each value in each column is appended to a list
	with open(path) as f:
		header = [h.strip() for h in f.next().split(',')] #strip whitespace in headers
		reader = csv.DictReader(f, fieldnames=header) # read rows into a dictionary format
		for row in reader: # read a row as {column1: value1, column2: value2,...}
			for (k,v) in row.items(): # go over each column name and value 
				dictionary[k].append(v) # append the value into the appropriate list based on column name k
	return dictionary

def s2g_sdss(sdssOutDir, sdssIdFile, galfitInDir, quantity):
	sdssIds = read_csv(sdssIdFile)
	#for i in range(len(sdssIds["name"])): 

	for i in range(quantity): 
		objId = sdssIds["objID"][i]
		objIdHash = objId[-3:]
		objIdPath = os.path.join(sdssOutDir, objIdHash, objId)
		s2g(objIdPath, galfitInDir, objId)
		convertFITSCommand = 'convert '  + os.path.join(objIdPath, objId + '-A_input.png') + ' fits/' + objId + '.fits'
		subprocess.call(convertFITSCommand, shell=True)  

	for i in range(quantity): 
		objId = sdssIds["objID"][i]
		objIdHash = objId[-3:]
		objIdPath = os.path.join(sdssOutDir, objIdHash, objId)
		galfitCommand = 'galfit galfit_in/' + objId + '.feedme'
		subprocess.call(galfitCommand, shell=True)  

	for i in range(quantity): 
		objId = sdssIds["objID"][i]
		objIdHash = objId[-3:]
		objIdPath = os.path.join(sdssOutDir, objIdHash, objId)
		try:
			compare.compare_images('galfit_out/' + objId + '-out.fits', os.path.join(objIdPath, objId + '-K_clusMask-reprojected.png'), 'compare_out', objId) 
		except:
			print objId + " failed - galfit ouput is likely to be missing."
			

if __name__ == "__main__":
	sdssOutDir = '/extra/wayne1/research/drdavis/SDSS/SpArcFiRe/r'
	sdssIdFile = '/extra/wayne1/research/drdavis/SDSS/SpiralsIDS.csv'
	galfitInDir = 'galfit_in'
	quantity = 25
	s2g_sdss(sdssOutDir, sdssIdFile, galfitInDir, quantity)
