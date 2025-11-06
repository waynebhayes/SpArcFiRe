# AlignBandColors by Anthony Navarrette

## About

AlignBandColors is a program that aligns inter-waveband images of the same galaxy to a 1/30th of a pixel accuracy.  A detailed analysis of this program was published in the [Journal of Astronomy and Computing](https://www.sciencedirect.com/science/article/pii/S2213133722000415).  

## Requirements

- Using a x86 Linux based operating system
- Python 3.10 environment outlined below
- (Recommended) 16GB of memory or more
    - If running into memory limitations, turn `-runInParallel` off, this will process only a single image at a time.
    - If still running into issues then lower `-upscaleFactor` (this will lead to lower accuracy if less than 100).

## Environment Setup
Install Python 3.10 and the included pip packages requirements.txt.  This can be done globally or use through pyenv as outlined below.

### Pyenv Setup
1. Follow [this link](https://realpython.com/intro-to-pyenv/) for instructions on setting up pyenv.  Be sure to install all system packages listed and add the required changes to your shell profile.
2. Install Python 3.10 with the following `
pyenv install 3.10`
3. Confirm it is available with `pyenv versions`
4. Navigate to the `AlignBandColors` directory and set the local python version using `pyenv local 3.10`, verify it is selected by `pyenv versions`
5. Install the required packages using `python -m pip install -r requirements.txt`

## Image Setup
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

- November 2025
    - Fix file loading path issue
    - Fix g-band needing to be present
    - Update README and environment setup

