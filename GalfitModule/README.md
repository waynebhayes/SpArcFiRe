<h2> README (work in progress) </h2>

The scripts in this folder are written to automatically generate GALFIT input files from SpArcFiRe output.

Please note, this module is still a work-in-progress and relies on some assumptions having to do with SpArcFiRe's conventions.
These will be fixed in a future update, likely after Matthew finishes his Ph.D.

Written by Matthew Portman, portmanm@uci.edu

<h3> INTRODUCTION </h3>

The previous version of the GALFIT processing was simply a few scripts, **s2g.py** and others, written to automatically generate the disk and bulge
parameters for GALFIT. **s2g.py** was intended to include the generation of spiral arc parameters but was halted
before it could be completed.

The current version of the code has been written as a Python module in order to use an OOP framework.
In the GalfitModule, OOP constitutes the underlying framework and much of the pre and post-processing relies on it. 
The creation of the module in this framework is twofold:
1) OOP allows us to utilize distributed computing with multi-step fitting (necessary for certain GALFIT practices) and
2) The OOP framework can be used by others to streamline their usage of GALFIT with Python.
Note, the GalfitModule in its current state _requires_ **Python3.7 or greater**. This choice has been made to utilize
the CaptureOutput feature of the Subprocess Python Module, which is crucial during multi-step fitting.

There are three main scripts which are run to perform the fitting, 
**sparc_to_galfit_feedme_gen.py**, which generates the galfit input file from the results from SpArcFiRe,
**go_go_galfit.py**, which more or less runs GALFIT and is the script that is sent to the compute nodes in distributed computing, and
**control_script.py**, which prepares everything for **go_go_galfit.py** and processes the output.

**sparc_to_galfit_feedme_gen.py** is written such that it can be run on its own from the command line to generate the GALFIT input
files without actually running them. These files are placed in the individual galaxy directories. We recommend doing this on
a first installation to verify that everything is working correctly before running **control_script.py**.

The **control_script.py** operates on the assumption that there are three directories, input, temporary, and output to coincide
with SpArcFiRe's convention of the same. The script is currently restricted to directories which end in _"-in", "-tmp", "-out"_
for convention reasons and some assumptions used in the module when processing. Again, this will likely disappear in a production version.
**control_script.py** defaults to using on-node parallelization via **parallel** which can be found in the ParallelDrivers directory,
an executable written by Wayne Hayes at UCI. It can be run in serial and has several options which can be accessed by running
**control_script.py --help** on the command line.
See the 'Current Directory Structure' section below for an example of what files are put where when **control_script.py** is run.

**go_go_galfit.py** need not be touched. In fact, it may be a little confusing since it uses many functionalities from the OOP framework
to do most of the leg work, so it may be more of a rabbit hole than anything else. It is currently only somewhat documented but feel
free to contact Matthew if something is unclear.

DEPRECATED ~~There is now an additional script, *in_out_comparison.py* which compares our automatically generated input with GALFIT's output to 
determine where it may be improved. A version of this will be used for error detection in the future but for now, it simply offers
the user a chance to see how much the two differ. *in_out_comparison.py* outputs two types of files: ~~~~
* ~~one file per galaxy and located in the galaxy's folder which contains the parameter information from the input, output, and the difference 
between the two as a text file called *galfit_io_compare.txt*~~
~~* one file which contains *all* of the difference values between the two for every galaxy but just the differences. This is called~~
~~*comparison_params.csv* and will be generated in the folder *all_galfit_out*.~~

---

<h3> Current Directory Structure </h3> 

```
control_script.sh

sparc_to_galfit_feedme_gen.py

*...-in  -- input FITS files

               *galfit_masks -- SExtractor Star Masks
*...-tmp -- <- *galfits      -- Output FITS models (and empty files to indicate failure)
               *galfit_png   -- PNG files generated from output FITS

*...-out -- <- *galfit_png   -- FITS models converted to PNG

                                Output 'combined' PNG (observation, model, residual)
               *galaxy_dirs  <- Output FITS model
                                File(s) used for input to GALFIT
                                (PSFs go here as galaxyname_psf.fits)

```

---

<h3> TO RUN: </h3>

These two scripts should be placed in the same directory as the containing directory for the in, tmp, and out folders
used with SpArcFiRe. Since it's still a WIP, ***please name these folders [name]-in, [name]-tmp, [name]-out or
the script(s) won't work. Preferably sparcfire-in, sparfire-tmp, sparcfire-out just in case...*** 

---

<h3> JUST THE FEEDME PARAMETERS </h3>

`python sparc_to_galfit_feedme_gen.py`
(it must be run in Python >3.6 to work with f-string implementation)

We recommend running this alone upon first installation to ensure it works correctly.

---

<h3> THE FEEDMEs + GALFIT </h3>

`python3 control_script.py`

---
<h3> OUTPUT FILE STRUCTURE </h3>

The folder ~~*all_galfit_out*~~ *sparcfire-tmp/galfits* contains a summary of sorts of all pertinent output: models as FITS files, the models converted to PNGs,
and the *comparison_params.csv* file comparing the values between the input and output parameter files for GALFIT. This folder and
its subfolders are generated by **control_script.py**. **control_script.py** also creates folders inside *sparcfire-tmp* to hold
the masks, PSF files, and galfits (just in case GALFIT 'explodes') while the control script is still running. This alleviates the 
error catching burden when tidying up! 

You can also find all of the above in the individual galaxy folders - with the exception of PNGs and the overall comparison between
galaxies.~~; there is instead a single text file, *galfit_io_compare.txt* which contains the input, output, and difference between the two.~~

---

It is possible to perform the feedme generation with just these two but there are other optional
inputs which we recommend and have also included already built into the control script - if you do not wish to include
these, the scripts should run anyway albeit with many messages:

* Star Masks: Requires the use of SExtractor. Parameters and python script found in star_removal. **control_script.sh**
automatically generates a folder to hold the masks and checks if they've already been generated. Placed in a folder in
sparcfire-tmp. 

* PSF Generation (SDSS only): Requires the use of the *Read_PSF* program from SDSS's website. **sparc_to_galfit_feedme_gen.py**
takes information from the csv included in star_dl according to our own naming convention - as controlled by a function in
**sparc_to_galfit_feedme_gen.py** to be easily modified to follow SDSS' naming convention - and
outputs the information necessary to find and download the psfield file (work in progress) into the generated feedme. The 
psfield should then be processed by the Read_PSF executable.
** NOTE, once you compile the Read_PSF executable, you must replace the file, phConsts.h with the file included here or 
otherwise change the value of the SOFT_BIAS to 0. If you don't, the PSF generated will not work with GALFIT and will 
produce bad output due to the existence of a pedestal. 

* Output conversion to PNG: Requires the use of the *fitspng* utility from http://integral.physics.muni.cz/fitspng/. 
fitspng provides an easy and efficient way to convert the output from GALFIT to pngs for ease of viewing. This is also
already included in the control script. If the pngs seem to be blown out, it's possible the conversion utility is having
trouble with some flags in control_script.sh. If this is the case, try removing the flag -fr "1,150" and then run. This has
only occurred on one machine for reasons I'm not familiar with. control_script.sh then uses the *ImageMagick* utility to 
automatically panel the PNG converted output with some of SpArcFiRe's to reduce bloat and make comparison easy. The converted
PNGs are then deleted and all that is left is the combined/paneled image in a separate folder within sparcfire-out.


------------------------------------------------------------------------------------
Newest scripts written by Matthew Portman under the tutelage of Dr. Wayne Hayes, UCI.
Others written by the SDSS team and Darren Davis who also made SpArcFiRe. 
