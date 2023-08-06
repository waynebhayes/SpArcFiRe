import numpy as np
import astropy as ap
from astropy.io import fits

import os

from os.path import join as pj
from os.path import exists

import pickle
from joblib import Parallel, delayed

import sys

# _HOME_DIR = os.path.expanduser("~")
# try:
#     _SPARCFIRE_DIR = os.environ["SPARCFIRE_HOME"]
#     _MODULE_DIR = pj(_SPARCFIRE_DIR, "GalfitModule")
# except KeyError:
#     # print("SPARCFIRE_HOME is not set. Please run 'setup.bash' inside SpArcFiRe directory if not done so already.")
#     # print("Running on the assumption that GalfitModule is in your home directory... (if not this will fail and quit!)") 
#     _MODULE_DIR = pj(_HOME_DIR, "GalfitModule")
    
# sys.path.append(_MODULE_DIR)
# from Classes.Components import *
# from Classes.Containers import *
# from Classes.FitsHandlers import *
# from Functions.helper_functions import *

#from parallel_residual_calc import parallel_wrapper

if __name__ == "__main__":
    basename         = sys.argv[1]
    out_dir          = os.path.dirname(basename)
    galaxy_names     = sys.argv[2].split(",")
    
    #out_nmr = parallel_wrapper(galfit_tmp_path, galfit_mask_path, out_png_dir, all_gname_tmp_out)
    output_dict = {}
    for gname in galaxy_names:
        output_file = pj(out_dir, gname, f"{gname}_galfit_out.fits")
        if exists(output_file):
            with fits.open(output_file) as hdul: 
                output_dict[gname] = (hdul[2].header.get("NMR", None), 
                                      hdul[2].header.get("ks_p", None),
                                      hdul[2].header.get("ks_stat", None)
                                     )
    
    # In the future, drop this in out_dir
    pickle_filename = f'{basename}_output_nmr.pkl'
    pickle.dump(output_dict, open(pickle_filename, 'wb'))