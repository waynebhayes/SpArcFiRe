#!/bin/csh
########################################################################
#
# Sample C-Shell script for fv 2.6
#
# This script will open 2 sample files (available online and distributed
# with the fv executables) and displays their contents in several forms...
# header keyword list, a curve, and an image.  It uses the XPA program
# xpaset to control fv.  XPA is not distributed with fv, but can be
# obtained from the SAO/HEAD R&D group at 
#          http://hea-www.harvard.edu/RD/xpa/
# You will need to build the TCL version as well as the executables.
#
# USAGE:
#
#   After starting fv, executed this script file on the command line:
#       ./XPAScript.csh
#
########################################################################

# Set an environment variable here which points to the fits files
# to be opened.  This is just to make it easier to specify files
# later.

setenv FitsDir ftp://heasarc.gsfc.nasa.gov/software/ftools/release/other/pdw

#
# Open 2 sample files
#

xpaset -p fv open $FitsDir/ngc1316r.fit $FitsDir/rate.fit

#
# Select one of the files and open a header window of extension #1
#

xpaset -p fv select rate.fit
xpaset -p fv display header 1 

#
# Plot a curve of Time vs Rate in POW and alter the graph's appearance
#

xpaset -p fv display curve 1 time rate
xpaset -p fv pow bounds 770 -30 1070 300
xpaset -p fv pow curve pDisp No lDisp Yes lColor Blue

#
# Select the other file and plot an image, setting its colormap to histogram
#

xpaset -p fv select ngc1316r.fit
xpaset -p fv display image 0
xpaset -p fv pow colormap scale histo
