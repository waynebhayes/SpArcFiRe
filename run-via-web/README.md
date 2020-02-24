# Run-via-web

For those who don't want to install SpArcFiRe but have too many galaxies to process by hand on the website,
the script in this directory will automatically upload images to our webserver, where they will be processed,
and then download the results.

WARNING: the web server is not fast.

This script runs on Python 2.7.17

## Setup
```
python2.7 -m pip install requests pathlib wget --user
```

## Usage
First, place your source images into a single directory, and run the following script.
```
python2.7 run-via-web.py <srcdir> <destdir> [-o <output>]
```
The first argument, `srcdir`, is a directory of the source images. The second argument, `destdir`, is a directory where the galaxy data (as zip files) will be downloaded to.
The optional argument, `output`, is a file where the query ID will be stored.
