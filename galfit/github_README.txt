--------------------------------------------- README (work in progress)------------------------------------------

The scripts in this folder are written to automatically generate GALFIT input files from SpArcFiRe output.
Please note, this script is still a work-in-progress and relies on the hard-coding of some things including
directory structure. 

The previous version, s2g.py and its related scripts, are written to automatically generate the disk and bulge
parameters for GALFIT. This version intended to include the generation of spiral arc parameters but was halted
before it could be completed.

The new version, sparc_to_galfit_feedme_gen.py, includes the generation of spiral arc parameters. Due to the nature
of the work however, one cannot run sparc_to_galfit_feedme_gen.py and then immediately GALFIT; sparc_to_galfit_feedme_gen.py
places each feedme into the individual output folder for each galaxy. Furthermore, the feedme generated allots for the inclusion
of star masking, (and soon) PSFs, which are stored separately. For these reasons and more, we use a bash script aptly named
control_script.sh to handle all the nuts and bolts. To reiterate, for the components and parameters, you may run
sparc_to_galfit_feedme_gen.py alone. 

--------------------------------------------- TO RUN: ---------------------------------------------------------

These two scripts should be placed in the same directory as the containing directory for the in, tmp, and out folders
used with SpArcFiRe. Since it's still a WIP, please name these folders sparcfire-in, sparcfire-tmp, sparcfire-out or
the script(s) won't work. 

------------------------------------------ JUST THE FEEDME PARAMETERS ---------------------------------------

python sparc_to_galfit_feedme_gen.py (it was written in python3 but should work in python2 as well)

----------------------------------------- THE FEEDMEs + GALFIT --------------------------------------

bash control_script.sh

------------------------------------------------------------------------------------------------------

It is possible to perform the feedme generation with just these two but there are other optional
inputs which we recommend and have also included already built into the control script - if you do not wish to include
these, the scripts should run anyway albeit with many messages:

* Star Masks: Requires the use of SExtractor. Parameters and python script found in star_removal. control_script.sh
automatically generates a folder to hold the masks and checks if they've already been generated. Placed in a folder in
sparcfire-tmp. 

* PSF Generation (SDSS only): Requires the use of the Read_PSF program from SDSS's website. sparc_to_galfit_feedme_gen.py
takes information from the csv included in star_dl according to our own naming convention - as controlled by a function in
sparc_to_galfit_feedme_gen.py to be easily modified to follow SDSS' naming convention - and
outputs the information necessary to find and download the psfield file (work in progress) into the generated feedme. The 
psfield should then be processed by the Read_PSF executable.
** NOTE, once you compile the Read_PSF executable, you must replace the file, phConsts.h with the file included here or 
otherwise change the value of the SOFT_BIAS to 0. If you don't, the PSF generated will not work with GALFIT and will 
produce bad output due to the existence of a pedestal. 

* Output conversion to PNG: Requires the use of the fitspng utility from http://integral.physics.muni.cz/fitspng/. 
fitspng provides an easy and efficient way to convert the output from GALFIT to pngs for ease of viewing. This is also
already included in the control script. If the pngs seem to be blown out, it's possible the conversion utility is having
trouble with some flags in control_script.sh. If this is the case, try removing the flag -fr "1,150" and then run. This has
only occurred on one machine for reasons I'm not familiar with. control_script.sh then uses the ImageMagick utility to 
automatically panel the PNG converted output with some of SpArcFiRe's to reduce bloat and make comparison easy. The converted
PNGs are then deleted and all that is left is the combined/paneled image in a separate folder within sparcfire-out.


------------------------------------------------------------------------------------
Some scripts written by Matthew Portman under the tutelage of Dr. Wayne Hayes, UCI.
Others written by the SDSS team and Darren Davis who also made SpArcFiRe. 
