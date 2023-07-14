import numpy as np
import astropy as ap
import pandas as pd
from astropy.io import fits
# # import scipy.linalg as slg
# # from scipy.stats import norm
from math import ceil

# import plotly.express as px
# import plotly.graph_objects as go
# from plotly.subplots import make_subplots

# import matplotlib.pyplot as plt
# import matplotlib.colors as colors
# from matplotlib.colors import LinearSegmentedColormap
import glob
import os
import argparse

from os.path import join as pj
from os.path import abspath as absp

import PIL
import pickle
from joblib import Parallel, delayed

import sys

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

if __name__ == "__main__":
    run_dir = os.getcwd()

    #in_dir = pj(_HOME_DIR, "29k_galaxies_obs")
    #run_dir = pj(_HOME_DIR, "29k_galaxies")
    # in_dir  = "sparcfire-in"
    # tmp_dir = "sparcfire-tmp"
    # out_dir = "sparcfire-out"
    #galfits_tmp = "galfits"
    #galfit_masks = "galfit_masks"
    #galfit_out = "all_galfit_out"
    #nmr = "norm_masked_residual"
    
    parser = argparse.ArgumentParser()#description = USAGE)
    
    # TODO: invert this to running without slurm for when it comes time for the big runs
    parser.add_argument('--no-slurm',
                        dest     = 'slurm',
                        action   = 'store_const',
                        const    = False,
                        # Soon to be the other way around
                        default  = True,
                        help     = 'Run GALFITs using Slurm.')

    parser.add_argument('-drS', '--dont-remove-slurm',
                        dest     = 'dont_remove_slurm',
                        action   = 'store_const',
                        const    = True,
                        default  = False,
                        help     = 'Choose NOT to remove all old slurm files (they may contain basic info about each fit but there will be a bunch!)')
    
    parser.add_argument('-r', '--restart',
                        dest     = 'restart', 
                        action   = 'store_const',
                        const    = True,
                        default  = False,
                        help     = 'Restart script on the premise that some have already run (likely with SLURM).')

    parser.add_argument(dest     = 'paths',
                        nargs    = "*",
                        type     = str,
                        help     = "IN-DIRECTORY TMP-DIRECTORY OUT-DIRECTORY from SpArcFiRe. \
                                    Must follow -in, -tmp, out or this won't work.")
    
    args = parser.parse_args() # Using vars(args) will call produce the args as a dict
    slurm = args.slurm
    dont_remove_slurm = args.dont_remove_slurm
    restart = args.restart
    
    if len(args.paths) == 3:
            in_dir, tmp_dir, out_dir = args.paths[0], args.paths[1], args.paths[2]
            #print(f"Paths are, {in_dir}, {tmp_dir}, {out_dir}")
    else:
        in_dir = pj(run_dir, "sparcfire-in")
        tmp_dir = pj(run_dir, "sparcfire-tmp")
        out_dir = pj(run_dir, "sparcfire-out")

        print(f"Paths incorrectly specified, defaulting to (in, tmp, out)...")
        print(f"{in_dir}\n{tmp_dir}\n{out_dir}")
        print()
        
    #in_dir = pj(_HOME_DIR, "29k_galaxies_obs")
    in_dir = absp(in_dir)
    tmp_dir = absp(tmp_dir)
    out_dir = absp(out_dir)
    
# ==================================================================================================================

def fill_objects(in_tuple, galfit_mask_path, out_png_dir):
    gname  = in_tuple[0]
    fitspath = in_tuple[1]
    count = in_tuple[2]
    
    #galfit_mask_path = "/home/portmanm/29k_galaxies/sparcfire-tmp/galfit_masks"
    # if gname in output_fits_dict and gname in mask_dict:
    #     continue

    if not count % 1000:
        print(count, gname)

    star_mask_name = f"{gname}_star-rm.fits"
    mask_fits_name = pj(galfit_mask_path, star_mask_name)

    try:
        fits_file = OutputFits(fitspath)
        mask_fits_file = FitsFile(mask_fits_name)
    except Exception as e:
        print(f"There was an issue opening galaxy {gname}. Continuing...")
        print(e)
        return None, None
    
    masked_residual_normalized = fits_file.generate_masked_residual(mask_fits_file)
    if masked_residual_normalized is None:
        #print(f"Could not calculate nmr") # for galaxy {gname}. Continuing...")
        return None, None
    
    fits_file.to_png(out_png_dir = out_png_dir) #"/home/portmanm/29k_galaxies/sparcfire-out/galfit_png")
    
    return gname, fits_file.nmr#, fits_file.nmrr

    # output_fits_dict[gname] = fits_file
    # mask_dict[gname] = mask_fits_file
# ==================================================================================================================

# TODO chunk this up via slurm
if __name__ == "__main__":
    
    outpath = pj(run_dir, out_dir)
    total_gal = len(glob.glob(pj(outpath, "/123*/")))

    galfit_tmp_path = pj(run_dir, tmp_dir, "galfits")
    galfit_mask_path  = pj(run_dir, tmp_dir, "galfit_masks")
    out_png_dir = pj(outpath, "galfit_png")
    
    out_nmr = Parallel(n_jobs = -2, timeout = 30)(
                       delayed(fill_objects)(
                                             (os.path.basename(i).rstrip("_galfit_out.fits") , i, count),
                                             galfit_mask_path,
                                             out_png_dir
                                            )
                       for count, i in enumerate(glob.glob(pj(galfit_tmp_path, "*_galfit_out.fits"))) if not os.path.basename(i).startswith("failed")
                                                  )
    
    pickle_filename = 'output_nmr.pkl'
    pickle.dump(out_nmr, open(pickle_filename, 'wb'))
    
    #output_fits_dict = dict(zip(out_nmr))

       
