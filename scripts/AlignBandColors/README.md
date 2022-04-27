# AlignBandColors by Anthony Navarrette

## About

AlignBandColors is a program that aligns inter-waveband images of the same galaxy to a 1/30th of a pixel accuracy.  A detailed analysis of this program was published in the Journal of Astronomy and Computing [link tbd].  

## Requirements

- Using a x86 Linux based operating system

- Python 3.10
   - astropy        (5.0.4)
   - matplotlib     (3.5.1)
   - numpy          (1.22.3)
   - Pillow         (9.1.0)
   - psutil         (5.9.0)
   - scipy          (1.8.0)

- [Recommended] 16GB of memory or more
    - If running into memory limitations, turn `-runInParallel` off
    - If still running into issues then lower `-upscaleFactor` (this will lead to lower accuracy if less than 100)

## Setup

Install Python 3.10 and run `setup.sh` to install all of the above libraries.

Before running AlignBandColors on images, make sure to have your input images organized as follows:

InputDirectory/  
|-- GalaxyOne/  
|  |-- GalaxyOne_u.fits  
|  |-- GalaxyOne_g.fits  
|  |-- GalaxyOne_r.fits  
|-- GalaxyTwo/  
|  |-- GalaxyTwo_u.fits.xz  
|  |-- GalaxyTwo_g.fits.xz  
|  |-- GalaxyTwo_r.fits.xz  


The fits image names must contain a '_color' in the filename, otherwise the program will not work correctly.  The images can also be compressed with xz and will automatically be decompressed at runtime (no other compression types are currently supported).

For documentation on command line options, run 
>./abc.sh -h


## Version History
- July 2021
    - Initial release

- April 2022
    - Upgraded to Python 3.10 (package versions listed above)
    - Added option for waveband labels (-colorsToProcess [COLORS]) 
    - Fixed issue with updating non-SDSS FITS headers
    - Updated output data to be stored in a CSV

