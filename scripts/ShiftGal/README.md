# ShiftGal by Anthony Navarrette (antholn1@uci.edu)

## About

ShiftGal is a program that alignes inter-color-band images of the same galaxy to a 100th of a pixel accuracy.  A detailed analysis of this program was submitted to MNRAS for publication (link comming soon).  

Credit to Sergey Anatolyevich with the ALGLIB Project for their [data processing library.](https://www.alglib.net/)

## Requirements

Python 2.7
   - astropy        (1.4.9)
   - matplotlib     (2.2.5)
   - numpy          (1.16.6)
   - opencv-python  (2.4.5)
   - Pillow         (6.2.2)
   - psutil         (5.7.0)
   - scipy          (1.2.3)

It is recommended that you have at least 16GB of memory available, or more specifically, MaxImgDimensionX\*MaxImgDimensionY\*800 bytes of memory to run this program. If you have limited memeory, first set -runInParallel 0, if that is not enough, you must lower the -upscaleFactor until it runs sucessfully (note: the lower the upscale factor the less precison in the output).

## Setup

Before running, make sure to have your input images organized as follows:

>* InputDirectory
>    * GalaxyOneDirectory
>        * GalaxyOne_g.fits
>        * GalaxyOne_i.fits
>        * GalaxyOne_r.fits
>        * GalaxyOne_u.fits
>        * GalaxyOne_z.fits
>
>     * GalaxyTwoDirectory
>        * GalaxyTwo_g.fits
>        * GalaxyTwo_i.fits
>        * GalaxyTwo_r.fits
>        * GalaxyTwo_u.fits
>        * GalaxyTwo_z.fits
>...

The fits image names must contain an '_color' in the filename, otherwise ShiftGal will not work correctly.  The images can also be compressed as .fits.xz files and will automatically be decompressed at runtime.

For documentation on command line options, run 
>./runShiftGal.sh -h
