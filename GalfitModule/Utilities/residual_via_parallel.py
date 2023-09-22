import numpy as np
import astropy as ap
from astropy.io import fits

import os

from os.path import join as pj

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

from parallel_residual_calc import parallel_wrapper

if __name__ == "__main__":
    basename         = sys.argv[1]
    pkl_end_str      = sys.argv[2]
    galfit_tmp_path  = sys.argv[3]
    galfit_mask_path = sys.argv[4]
    out_dir          = sys.argv[5]
    to_png           = sys.argv[6]
    
    if str.lower(to_png) == "false":
        out_png_dir = False
    elif str.lower(to_png) == "true":
        out_png_dir = True
        
    all_gname_tmp_out = sys.argv[7].split(",")
    
    out_df = parallel_wrapper(galfit_tmp_path, galfit_mask_path, out_dir, to_png, all_gname_tmp_out)
    
    # In the future, drop this in out_dir
    pickle_filename = f'{basename}_{pkl_end_str}.pkl'
    #pickle.dump(out_nmr, open(pickle_filename, 'wb'))
    out_df.set_index("gname", inplace = True)
    out_df.to_pickle(pickle_filename)