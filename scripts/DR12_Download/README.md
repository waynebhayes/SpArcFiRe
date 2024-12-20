# SDSS DR12 Downloader by Anthony Navarrette (antholn1@uci.edu)

## About

A program that automatically downloads SDSS DR12 data using two methods:
1. a single object given its name, ra, and dec (name does not have to be objId, but should be for consistency)
2. all galaxies in a frame given its field, run, and camcol (requires the provided frame_object_list.tsv and SDSS_name_ra_dec_clean.tsv files to function)

## Requirements

Python 2.7
   - astropy        (1.4.9)
   - numpy          (1.16.6)

## Setup

If using method one to download data, have a file with three columns in the following order (name ra dec), with each line being a new object.
If using method two to download data, have a file with three columns in the following order (frame camcol field), with each line being a new frame.
In both cases, the delimiter of the input files do not matter as long as it is consistent throughout the whole file. Examples of how these files should be formatted are in ./inputExamples.

For documentation on command line options, run 
>./bulkDownload.sh -h
