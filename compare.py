#!/usr/bin/env python
import sys
import os
import csv
import math
from collections import defaultdict

import numpy as np
from scipy.misc import imread, imsave
from scipy.linalg import norm
from scipy import sum, average
import pyfits

def compare():
	#Step through all files in galfit output directory
	galfitDir = 'galfit_out'
	sparcfireDir = 'sparcfire_out'
	outputDir = 'compare_out'
	for fileFits in os.listdir(galfitDir):
		sdssId = os.path.splitext(fileFits)[0][0:-4] #Remove '-out' from filename
		pathMask = os.path.join(sparcfireDir, sdssId, sdssId + "-K_clusMask-reprojected.png")
		pathContrast = os.path.join(sparcfireDir, sdssId, sdssId + "-K_clusMask-reprojected.png")
		pathFits = os.path.join(galfitDir, fileFits)
		pathCSV = os.path.join(sparcfireDir, sdssId, sdssId + '.csv')
		compare_images(pathFits, pathMask, pathCSV, outputDir, sdssId)

def compare_images(pathFits, pathMask, pathCSV, outputDir, sdssId):
	thresholds = [0.7, 0.725, 0.75, 0.775, 0.8]
	print 'SDSS id: ' + sdssId
	bestThreshold = thresholds[0]
	bestConfidence = 0
	for threshold in thresholds: 
		print '	Threshold: ' + str(threshold)
		# read images as 2D arrays (convert to grayscale for simplicity)
		imgInput = normalize(grayscale(pyfits.open(pathFits)[1].data))
		imgModel = normalize(grayscale(pyfits.open(pathFits)[2].data))
		imgResidual = normalize(grayscale(pyfits.open(pathFits)[3].data))
		imgDiskMask = mask_disk(imgInput.shape, pathCSV)
		imgClusterMask = np.flipud(normalize(grayscale(imread(pathMask).astype(float))))
		
		(confidence, imgInputMasked, imgClusterMask, relevantElements, selectedElements, truePositives, imgMaskedResidual) = calculate_confidence(imgInput, imgResidual, imgClusterMask, imgDiskMask, threshold)
		print '		Confidence: ' + str(confidence)
		if confidence > bestConfidence:
			bestConfidence = confidence
			bestThreshold = threshold

		img = np.concatenate((imgInput, imgModel, imgResidual, relevantElements, selectedElements, truePositives, imgMaskedResidual), axis=1)
		try:
			imgOutput = np.concatenate((imgOutput, img), axis=0)
		except:
			imgOutput = img
		
	print '		Best threshold: ' + str(bestThreshold)
	print '		Best confidence: ' + str(bestConfidence)
	pathOutput= os.path.join(outputDir, str(bestThreshold) + "t-" + str(bestConfidence) + "c-" + sdssId + ".png")
	imsave(pathOutput, imgOutput)
	print '		Saved at: ' + pathOutput
	
def calculate_confidence(imgInput, imgResidual, imgClusterMask, imgDiskMask, threshold):
	#Step 1: Apply disk mask - only interested in fit within disk of galaxy
	imgInputMasked = np.multiply(imgInput, imgDiskMask)
	imgClusterMaskMasked = np.multiply(imgClusterMask, imgDiskMask)
	
	#Step 2: Threhold and binarize the images to determine relevant and selected elements.
	relevantElements = np.where(imgInputMasked > 255 * threshold, 1, 0)	
	selectedElements = np.where(imgClusterMaskMasked > 127, 1, 0)
	print '		sum(selectedElements): ' + str(sum(selectedElements))
	print '		sum(relevantElements): ' + str(sum(relevantElements))
	
	#Step 3: Calculate number of true positives.
	truePositives = float(sum(np.multiply(selectedElements, relevantElements)))
	print '		truePositives: ' + str(truePositives)

	#Step 4: Calculate F1 confidence.
	precision = truePositives / float(sum(selectedElements))
	recall = truePositives / float(sum(relevantElements))
	print '		precision: ' + str(precision)
	print '		recall: ' + str(recall)

	truePositives = np.multiply(selectedElements, relevantElements)
	truePositives = np.where(truePositives > 0, 255, 0)	
	imgMaskedResidual = np.multiply(np.where(truePositives > 127, 0, 1), imgResidual)
	relevantElements = np.where(relevantElements > 0, 255, 0)	
	selectedElements = np.where(selectedElements > 0, 255, 0)

	return ((2 * precision * recall) / (precision + recall), imgInputMasked, imgClusterMaskMasked, relevantElements, selectedElements, truePositives, imgMaskedResidual)

def mask_disk(shape, pathCSV):
	# Parse galaxy.csv file
	mask = np.zeros(shape)
	galaxies = defaultdict(list) # each value in each column is appended to a list
	with open(pathCSV) as f:
		header = [h.strip() for h in f.next().split(',')] #strip whitespace in headers
		reader = csv.DictReader(f, fieldnames=header) # read rows into a dictionary format
		for row in reader: # read a row as {column1: value1, column2: value2,...}
			for (k,v) in row.items(): # go over each column name and value 
				galaxies[k].append(v) # append the value into the appropriate list based on column name k
	# Define disk ellipse	
	center = (float(galaxies["inputCenterC"][0]), float(galaxies["inputCenterR"][0]))
	width = float(galaxies["diskMajAxsLen"][0])
	height = float(galaxies["diskMajAxsLen"][0]) * float(galaxies["diskAxisRatio"][0])
	angle = float(galaxies["diskAxisRatio"][0])
	#print center, width, height, angle
	
	# Mask pixels within disk ellipse
	for index, value in np.ndenumerate(mask):
		mask[index[0], index[1]] = in_ellipse(index, center, width, height, angle) 	

	return mask

def in_ellipse(index, center, width, height, angle):
	cos_angle = np.cos(angle)
	sin_angle = np.sin(angle)

	xc = index[0] - center[0]
	yc = index[1] - center[1]

	xca = xc * cos_angle - yc * sin_angle
	yca = xc * sin_angle + yc * cos_angle 

	distance = (xca**2 / (width / 2.0)**2) + (yca**2 / (height / 2.0)**2)
	# If point is within ellipse, return True. 
	return abs(distance) <= 1

def grayscale(arr):
	"If arr is a color image (3D array), convert it to grayscale (2D array)."
	if len(arr.shape) == 3:
		return np.amax(arr, axis=-1)  # take max value over the last axis (color channels)
	else:
		return arr

def normalize(arr):
	arr = np.log(np.clip(arr, 1.0, float("inf")))
	rng = arr.max() - arr.min()
	amin = arr.min()
	return (arr - amin) / rng  * 255

if __name__ == "__main__":
	compare()
