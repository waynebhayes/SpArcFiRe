#!/usr/bin/env python
"""Compare two aligned images of the same size.

Usage: python compare.py first-image second-image
Original author: astanin @ http://stackoverflow.com/a/3935002/1905613
"""

import sys
import os
import csv
import math
from collections import defaultdict

import numpy
from scipy.misc import imread, imsave
from scipy.linalg import norm
from scipy import sum, average

def main():
	#Step through all subdirs in SpArcFiRe directory
	d = sys.argv[1]
	subdirs = [s for s in os.listdir(d) if os.path.isdir(os.path.join(d, s)) and s != "matout"]
	for s in subdirs:
		#Generate image residuals for each input
		fileA = os.path.join(d, s, s + "-A_input.png")
		fileK = os.path.join(d, s, s + "-K_clusMask-reprojected.png")
		fileL = os.path.join(d, s, s + "-L_input-blurred.png")
		fileM = os.path.join(d, s, s + "-M_residual.png")
		fileN = os.path.join(d, s, s + "-L_blurred-residual.png")
		fileCSV = os.path.join(d, s, s + ".csv")
		compare_images(fileA, fileK, fileL, fileM, fileN, fileCSV)
	
def compare_images(fileA, fileK, fileL, fileM, fileN, fileCSV):
	# read images as 2D arrays (convert to grayscale for simplicity)
	imgA = to_grayscale(imread(fileA).astype(float))
	imgK = to_grayscale(imread(fileK).astype(float))
	# Compare difference without bulge masking
	n_m, n_0, chi2nu, bright_ratio, diff = compare(imgA, imgK)
	print "Stats for residual at ", fileM
	print "	Without bulge mask:"
	print "		Manhattan norm/Nu:	", n_m/imgA.size
	print "		Zero norm/Nu:		", n_0*1.0/imgA.size
	print "		Chi2/Nu:		", chi2nu
	print "		Brightness ratio:	", bright_ratio

	# Compare difference with bulge masking
	imgK = mask_bulge(imgK, fileCSV)
	n_m_b, n_0_b, chi2nu_b, bright_ratio_b, diff_b = compare(imgA, imgK)
	imsave(fileM, diff_b)
	print "	With bulge mask:"
	print "		Manhattan norm/Nu:	", n_m_b/imgA.size
	print "		Zero norm/Nu:		", n_0_b*1.0/imgA.size
	print "		Chi2/Nu:		", chi2nu_b
	print "		Brightness ratio:	", bright_ratio_b
	print "	Mask/no mask differences:"
	print "		Manhattan norm/Nu:	", n_m_b/imgA.size - n_m/imgA.size
	print "		Zero norm/Nu:		", n_0_b*1.0/imgA.size - n_0*1.0/imgA.size
	print "		Chi2/Nu:		", chi2nu_b - chi2nu
	print "		Brightness ratio:	", bright_ratio_b - bright_ratio

def compare(imgA, imgK):
	# normalize to compensate for exposure difference
	imgA = normalize(imgA)
	imgK = normalize(imgK)
	# calculate the difference and its norms
	diff = imgA - imgK
	diff = numpy.clip(imgA - imgK, 0, 255)  # elementwise for scipy arrays
	#print "diff min/max", diff.min(), " ", diff.max()
	#print "A min/max", imgA.min(), " ", imgA.max()
	#print "K min/max", imgK.min(), " ", imgK.max()
	m_norm = sum(abs(diff))  # Manhattan norm - the sum of the absolute values
	z_norm = norm(diff.ravel(), 0)  # Zero norm - the number of elements not equal to zero
	chi2nu_value = chi2nu(diff) # Chi2Nu - sum of chi squared values divided by number of degrees of freedom
	b_ratio = bright_ratio(imgA, diff) # brightness ratio - the ratio of ratios of upper quartile of pixels over lower quartile of pixels between input and masked images 
	return (m_norm, z_norm, chi2nu_value, b_ratio, diff)

def bright_ratio(imgA, diff):
	# Calculate ratio of pixels between top and bottom quartile of pixels
	# We expect a higher brightness ratio value to correspond to better masking of bright pixels in the input image by the cluster mask.
	histogramA = numpy.histogram(imgA, bins=[0, 63, 127, 191, 255])
	#print "	histogramA = ", histogramA
	bright_ratioA = float((imgA >= numpy.percentile(imgA, 75)).sum()) / float((imgA <= numpy.percentile(imgA, 25)).sum()) 
	bright_ratioA = float(histogramA[0][3]) / float(histogramA[0][0])
	#print "	bright_ratioA = ", bright_ratioA
 	histogramDiff = numpy.histogram(diff, bins=[0, 63, 127, 191, 255])
	#print "	histogramDiff = ", histogramDiff
	bright_ratioDiff = float((diff >= numpy.percentile(diff, 75)).sum()) / float((diff <= numpy.percentile(diff, 25)).sum()) 
	bright_ratioDiff = float(histogramDiff[0][3]) / float(histogramDiff[0][0])
	#print "	bright_ratioDiff = ", bright_ratioDiff
	bright_ratio = bright_ratioA / bright_ratioDiff
	return bright_ratio

def chi2nu(arr):
	# sum of chi squared values divided by number of degrees of freedom
	expected2 = average(arr) ** 2
	return sum(arr ** 2 / expected2) / arr.size

def mask_bulge(img, fileCSV):
	# Parse galaxy.csv file
	galaxies = defaultdict(list) # each value in each column is appended to a list
	with open(fileCSV) as f:
		header = [h.strip() for h in f.next().split(',')] #strip whitespace in headers
		reader = csv.DictReader(f, fieldnames=header) # read rows into a dictionary format
		for row in reader: # read a row as {column1: value1, column2: value2,...}
			for (k,v) in row.items(): # go over each column name and value 
				galaxies[k].append(v) # append the value into the appropriate list based on column name k
	# Define bulge ellipse	
	center = (float(galaxies["inputCenterC"][0]), float(galaxies["inputCenterR"][0]))
	width = float(galaxies["bulgeMajAxsLen"][0])
	height = float(galaxies["bulgeMajAxsLen"][0]) * float(galaxies["bulgeAxisRatio"][0])
	angle = float(galaxies["bulgeAxisRatio"][0])
	#print center, width, height, angle
	
	# Mask pixels within bulge ellipse
	for index, value in numpy.ndenumerate(img):
		img[index[0], index[1]] = in_ellipse(index, value, center, width, height, angle) 	

	return img

def in_ellipse(index, value, center, width, height, angle):
	cos_angle = numpy.cos(angle)
	sin_angle = numpy.sin(angle)

	xc = index[0] - center[0]
	yc = index[1] - center[1]

	xca = xc * cos_angle - yc * sin_angle
	yca = xc * sin_angle + yc * cos_angle 

	distance = (xca**2 / (width / 2.0)**2) + (yca**2 / (height / 2.0)**2)
	#print distance
	# If point is within ellipse, mask to max value of image. Else, return existing value. 
	if (abs(distance) <= 1):
		#print value, distance
		return 255.0
	else:
		return value

def to_grayscale(arr):
	"If arr is a color image (3D array), convert it to grayscale (2D array)."
	if len(arr.shape) == 3:
		return numpy.amax(arr, axis=-1)  # average over the last axis (color channels)
	else:
		return arr

def normalize(arr):
	rng = arr.max()-arr.min()
	amin = arr.min()
	return (arr-amin)*255/rng

if __name__ == "__main__":
	main()
