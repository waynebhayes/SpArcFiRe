import numpy as np
import astropy as ap
# import pandas as pd
from astropy.io import fits

# from math import ceil

#import glob
import os
import argparse

from os.path import join as pj
from os.path import abspath as absp

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
    name              = sys.argv[1]
    galfit_tmp_path   = sys.argv[2]
    galfit_mask_path  = sys.argv[3]
    out_png_dir       = sys.argv[4]
    all_gname_tmp_out = sys.argv[5].split(",")
    
    out_nmr = parallel_wrapper(galfit_tmp_path, galfit_mask_path, out_png_dir, all_gname_tmp_out)
    
    # In the future, drop this in out_dir
    pickle_filename = f'{name}_output_nmr.pkl'
    pickle.dump(out_nmr, open(pickle_filename, 'wb'))