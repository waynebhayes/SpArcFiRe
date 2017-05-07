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
import pyfits

def main():
	#Step through all files in galfit output directory
	galfitDir = 'galfit_out'
	sparcfireDir = 'sparcfire_out'
	outputDir = 'compare_out'
	for fileFits in os.listdir(galfitDir):
		sdssId = os.path.splitext(fileFits)[0][0:-4] #Remove '-out' from filename
		pathMask = os.path.join(sparcfireDir, sdssId, sdssId + "-K_clusMask-reprojected.png")
		pathFits = os.path.join(galfitDir, fileFits)
		compare_galfit_images(pathFits, pathMask, outputDir, sdssId)

def compare_galfit_images(pathFits, pathMask, outputDir, sdssId):
	# read images as 2D arrays (convert to grayscale for simplicity)
	imgInput = normalize(grayscale(pyfits.open(pathFits)[1].data))
	imgModel = normalize(grayscale(pyfits.open(pathFits)[2].data))
	imgResidual = normalize(grayscale(pyfits.open(pathFits)[3].data))
	imgMask = normalize(grayscale(imread(pathMask).astype(float)))
	imgMaskedResidual = numpy.multiply(mask(imgMask), imgResidual)

	# Compare difference without arc masking
	n_m, n_0, chi2nu, bright_ratio, diff = compare_galfit(imgInput, imgResidual)
	print "Stats for residual at ", pathFits
	print "		Manhattan norm/Nu:	", n_m/imgInput.size
	print "		Zero norm/Nu:		", n_0*1.0/imgInput.size
	print "		Chi2/Nu:		", chi2nu
	print "		Brightness ratio:	", bright_ratio

	# Compare difference with arc masking
	n_m_b, n_0_b, chi2nu_b, bright_ratio_b, diff_b = compare_galfit(imgInput, imgMaskedResidual)
	print "	Arc masked residual"
	print "		Manhattan norm/Nu:	", n_m_b/imgInput.size
	print "		Zero norm/Nu:		", n_0_b*1.0/imgInput.size
	print "		Chi2/Nu:		", chi2nu_b
	print "		Brightness ratio:	", bright_ratio_b
	print "	Mask/no mask differences:"
	print "		Manhattan norm/Nu:	", n_m_b/imgInput.size - n_m/imgInput.size
	print "		Zero norm/Nu:		", n_0_b*1.0/imgInput.size - n_0*1.0/imgInput.size
	print "		Chi2/Nu:		", chi2nu_b - chi2nu
	print "		Brightness ratio:	", bright_ratio_b - bright_ratio

	confidence_ratio = numpy.clip(0.034403 + 0.252975 * numpy.log(numpy.clip(bright_ratio_b, 0.873, 50)), 0.0, 1.0) #equation through (1, 0), (5, 0.5), (50, 1.0) 
	print "	Confidence Ratio:		", confidence_ratio
	
	pathOutput= os.path.join(outputDir, sdssId + ".png")
	imgOutput = numpy.concatenate((imgInput, imgModel, imgResidual, imgMask, imgMaskedResidual), axis=1)
	imsave(pathOutput, imgOutput)
	
def compare_galfit(imgInput, diff):
	m_norm = sum(abs(diff))  # Manhattan norm - the sum of the absolute values
	z_norm = norm(diff.ravel(), 0)  # Zero norm - the number of elements not equal to zero
	chi2nu_value = chi2nu(diff) # Chi2Nu - sum of chi squared values divided by number of degrees of freedom
	b_ratio = bright_ratio(imgInput, diff) # brightness ratio - the ratio of ratios of upper quartile of pixels over lower quartile of pixels between input and masked images 
	return (m_norm, z_norm, chi2nu_value, b_ratio, diff)

def bright_ratio(imgInput, diff):
	# Calculate ratio of pixels between top and bottom quartile of pixels
	# We expect a higher brightness ratio value to correspond to better masking of bright pixels in the input image by the cluster mask.
	try:
		histogramA = numpy.histogram(imgInput, bins=[0, 63, 127, 191, 255])
		bright_ratioA = float(histogramA[0][3]) / float(histogramA[0][0])
		histogramDiff = numpy.histogram(diff, bins=[0, 63, 127, 191, 255])
		bright_ratioDiff = float(histogramDiff[0][3]) / float(histogramDiff[0][0])
		bright_ratio = bright_ratioA / bright_ratioDiff
	except:
		print "Likely the bulge mask is too large or too small - check fit"
		bright_ratio = float("inf")
	
	return bright_ratio

def chi2nu(arr):
	# sum of chi squared values divided by number of degrees of freedom
	expected2 = average(arr) ** 2
	return sum(arr ** 2 / expected2) / arr.size

def grayscale(arr):
	"If arr is a color image (3D array), convert it to grayscale (2D array)."
	if len(arr.shape) == 3:
		return numpy.amax(arr, axis=-1)  # average over the last axis (color channels)
	else:
		return arr

def normalize(arr):
	arr = numpy.log(numpy.clip(arr, 1.0, float("inf")))
	rng = arr.max() - arr.min()
	amin = arr.min()
	return (arr - amin) / rng  * 255

def mask(arr):
	return numpy.where(arr > 127, 0, 1)

if __name__ == "__main__":
	main()
