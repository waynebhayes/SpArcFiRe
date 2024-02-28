import numpy as np
import astropy as ap
import pandas as pd
from astropy.io import fits
from glob import glob

import os, sys
from os.path import join as pj

_HOME_DIR = os.path.expanduser("~")
try:
    _SPARCFIRE_DIR = os.environ["SPARCFIRE_HOME"]
    _MODULE_DIR = pj(_SPARCFIRE_DIR, "GalfitModule")
except KeyError:
    # print("SPARCFIRE_HOME is not set. Please run 'setup.bash' inside SpArcFiRe directory if not done so already.")
    # print("Running on the assumption that GalfitModule is in your home directory... (if not this will fail and quit!)") 
    _MODULE_DIR = pj(_HOME_DIR, "GalfitModule")
    
sys.path.append(_MODULE_DIR)
from Classes.Components import *
from Classes.Containers import *
from Classes.FitsHandlers import *
from Functions.helper_functions import *

in_dir  = "sparcfire-in"
tmp_dir = "sparcfire-tmp"
out_dir = "sparcfire-out"
model_dir = "." #"no_opt_models"

for fits_file in glob(pj(tmp_dir, "galfits", "*_galfit_out.fits")):
#for fits_file in glob(pj(model_dir, "*_galfit_out.fits")):
    gname = os.path.basename(fits_file).rstrip("_galfit_out.fits")
    
    observation_obj = FitsFile(pj(in_dir, f"{gname}.fits"))
    model_obj       = FitsFile(fits_file, name = "model")
    
    feedme = FeedmeContainer()
    feedme.from_file(pj(out_dir, gname, f"{gname}.in"))
    crop_box = feedme.header.region_to_fit.value
    
    xbox_min = crop_box[0] - 1
    xbox_max = crop_box[1]
    ybox_min = crop_box[2] - 1
    ybox_max = crop_box[3]
    
    observation = observation_obj.data[xbox_min:xbox_max, ybox_min:ybox_max]
    residual = observation_obj.data[xbox_min:xbox_max, ybox_min:ybox_max] - model_obj.data

    #with fits.open(pj(model_dir, os.path.basename(fits_file)), mode = 'update') as f:
    #hdu = fits.PrimaryHDU(np.zeros(1)) 
    observation_hdu = fits.ImageHDU(observation, 
                                    header = fits.getheader(pj(in_dir, f"{gname}.fits"))
                                   )
    model_hdu       = fits.ImageHDU(model_obj.data, 
                                    header = fits.getheader(fits_file)
                                   )
    
    hdr = fits.Header()
    hdr['OBSERVER'] = 'Edwin Hubble'
    hdr['COMMENT'] = "Here's some commentary about this FITS file."
    
    residual_hdu  = fits.ImageHDU(residual, header = hdr)
    empty_primary = fits.PrimaryHDU(header = hdr)
    
    hdul = fits.HDUList([empty_primary, observation_hdu, model_hdu, residual_hdu])
    
    new_file = pj(tmp_dir, os.path.basename(fits_file))
    if exists(new_file):
        sp(f"rm -f {new_file}")
        
    hdul.writeto(new_file) #pj(tmp_dir, os.path.basename(fits_file)))
    
    #new_fits = OutputFits(new_file)#pj(model_dir, os.path.basename(fits_file)))
    new_fits = FitsFile(new_file, from_galfit = True, names = ["observation", "model", "residual"])#pj(model_dir, os.path.basename(fits_file)))
    new_fits.all_hdu.pop("residual")
    new_fits.num_imgs = 2
    new_fits.to_png(out_png_dir = model_dir)

