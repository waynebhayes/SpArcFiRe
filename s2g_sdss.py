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

def s2g_sdss(sdssOutDir, sdssIds, galfitInDir, useHash, quantity, steps):
	if (steps[0]):
		for i in range(quantity): 
			objId = sdssIds[i]
			objIdHash = objId[-3:]
			if (useHash):
				objIdPath = os.path.join(sdssOutDir, objIdHash, objId)
			else:
				objIdPath = os.path.join(sdssOutDir, objId)
			s2g(objIdPath, galfitInDir, objId)
			convertFITSCommand = 'convert '  + os.path.join(objIdPath, objId + '-A_input.png') + ' fits/' + objId + '.fits'
			print 'Executing ' + convertFITSCommand
			subprocess.call(convertFITSCommand, shell=True)  

	if (steps[1]):
		for i in range(quantity): 
			objId = sdssIds[i]
			objIdHash = objId[-3:]
			if (useHash):
				objIdPath = os.path.join(sdssOutDir, objIdHash, objId)
			else:
				objIdPath = os.path.join(sdssOutDir, objId)
			galfitCommand = 'galfit galfit_in/' + objId + '.feedme'
			print 'Executing ' + galfitCommand 
			subprocess.call(galfitCommand, shell=True)  

	if (steps[2]):
		for i in range(quantity): 
			objId = sdssIds[i]
			objIdHash = objId[-3:]
			if (useHash):
				objIdPath = os.path.join(sdssOutDir, objIdHash, objId)
			else:
				objIdPath = os.path.join(sdssOutDir, objId)
			print 'Comparing ' + objId 
			compare.compare_images('galfit_out/' + objId + '-out.fits', os.path.join(objIdPath, objId + '-K_clusMask-reprojected.png'), os.path.join(objIdPath, objId + '.csv'), 'compare_out', objId) 
			
def read_SpiralsIDS():
	sdssOutDir = '/extra/wayne1/research/drdavis/SDSS/SpArcFiRe/r'
	sdssIdFile = '/extra/wayne1/research/drdavis/SDSS/SpiralsIDS.csv'
	sdssIds = read_csv(sdssIdFile)['objID']
	galfitInDir = 'galfit_in'
	useHash = True 
	quantity = 5
	steps = [True, True, True]
	s2g_sdss(sdssOutDir, sdssIds, galfitInDir, useHash, quantity, steps)
	

def process_directory():
	sdssOutDir = '/home/wayne/public_html/SDSS/G.out'
	sdssIdFile = '/home/wayne/public_html/SDSS/G.out/galaxy.csv'
	sdssIds = read_csv(sdssIdFile)['name']
	print str(sdssIds)
	galfitInDir = 'galfit_in'
	useHash = False
	quantity = len(sdssIds)
	steps = [False, False, True]
	s2g_sdss(sdssOutDir, sdssIds, galfitInDir, useHash, quantity, steps)

if __name__ == "__main__":
	process_directory()

