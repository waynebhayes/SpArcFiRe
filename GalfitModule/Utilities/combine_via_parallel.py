import numpy as np
import astropy as ap
from astropy.io import fits
import pandas as pd

import os

from os.path import join as pj
from os.path import exists

import pickle
from joblib import Parallel, delayed

import sys

_HOME_DIR = os.path.expanduser("~")

try:
    _SPARCFIRE_DIR = os.environ["SPARCFIRE_HOME"]
    _MODULE_DIR = pj(_SPARCFIRE_DIR, "GalfitModule")
except KeyError:

    if __name__ == "__main__":
        print("SPARCFIRE_HOME is not set. Please run 'setup.bash' inside SpArcFiRe directory if not done so already.")
        print("Checking the current directory for GalfitModule, otherwise quitting.")

    _MODULE_DIR = pj(os.getcwd(), "GalfitModule")

    if not exists(_MODULE_DIR):
        raise Exception("Could not find GalfitModule!")
    
sys.path.append(_MODULE_DIR)
from Classes.Components import *
from Classes.Containers import *
from Classes.FitsHandlers import *
from Functions.helper_functions import *

#from parallel_residual_calc import parallel_wrapper

def main(*args):
    basename         = args[0]
    out_dir          = args[1] #os.path.dirname(basename)
    galaxy_names     = args[2].split(",")
    
    #out_nmr = parallel_wrapper(galfit_tmp_path, galfit_mask_path, out_png_dir, all_gname_tmp_out)
    all_df = pd.DataFrame()
    DNE_dict = {"gname" : []}
    for gname in galaxy_names:
        output_file = pj(out_dir, gname, f"{gname}_galfit_out.fits")
        if exists(output_file):
            output_fits = OutputFits(output_file)
            output_df   = output_fits.feedme.to_pandas()
            output_df["gname"]   = gname
            output_df["NMR"]     = output_fits.model.header.get("NMR", None) 
            output_df["KS_P"]    = output_fits.model.header.get("KS_P", None)
            output_df["KS_STAT"] = output_fits.model.header.get("KS_STAT", None)
            
            all_df = pd.concat([all_df, output_df])
            # with fits.open(output_file) as hdul: 
            #     output_dict[gname] = (hdul[2].header.get("NMR", None), 
            #                           hdul[2].header.get("ks_p", None),
            #                           hdul[2].header.get("ks_stat", None)
            #                          )
            
        else:
            DNE_dict["gname"].append(gname)
            #output_dict[gname] = (None, None, None)
            
    DNE_dict = {col : [None]*len(DNE_dict["gname"]) for col in output_df.columns if col != "gname"}
    all_df = pd.concat([all_df, pd.DataFrame(DNE_dict)])
    all_df.set_index("gname", inplace = True)
    
    if basename:
        pickle_filename = f'{basename}_output_nmr.pkl'
        all_df.to_pickle(pickle_filename)
    #pickle.dump(output_dict, open(pickle_filename, 'wb'))
    
    return all_df

if __name__ == "__main__":
    basename         = sys.argv[1]
    out_dir          = sys.argv[2] #os.path.dirname(basename)
    galaxy_names     = sys.argv[3] #.split(",")
    
    main(basename, out_dir, galaxy_names)