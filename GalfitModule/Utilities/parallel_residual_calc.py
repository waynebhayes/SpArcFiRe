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

#from go_go_galfit import check_programs
    
# ==================================================================================================================

def fill_objects(gname, count, galfit_tmp_path, galfit_mask_path, out_png_dir = "", slurm = True):

    report_num = 1000
    if slurm:
        report_num = 100
        
    if not count % report_num:
        print(count, gname)

    star_mask_name = f"{gname}_star-rm.fits"
    mask_fits_name = pj(galfit_mask_path, star_mask_name)

    try:
        fits_filename = f"{gname}_galfit_out.fits"
        fits_file = OutputFits(pj(galfit_tmp_path, fits_filename))
    except Exception as e:
        print(f"There was an issue opening galaxy {gname}. Continuing...")
        print(e)
        return None, None
    
    try:
        mask_fits_file = FitsFile(mask_fits_name)
    except Exception as e:
        print(f"There was an issue opening star mask for galaxy {gname}. Proceeding without mask...")
        # Logic implemented to handle None
        mask_fits_file = None #np.zeros((500,500))
    
    masked_residual_normalized = fits_file.generate_masked_residual(mask_fits_file)
    if masked_residual_normalized is None:
        print(f"Could not calculate nmr for galaxy {gname}. Continuing...")
        return None, None
    
    # Doesn't work on Openlab (sadly)
    # Keep this in for actual slurming since it's hard to read booleans
    # from command line (use the default value to our advantage)
    if not slurm and out_png_dir:
        fits_file.to_png(out_png_dir = out_png_dir)
    
    return gname, fits_file.nmr#, fits_file.nmrr

    # output_fits_dict[gname] = fits_file
    # mask_dict[gname] = mask_fits_file
    
# ==================================================================================================================

def parallel_wrapper(galfit_tmp_path, galfit_mask_path, out_png_dir, all_gname_tmp_out, slurm = True):
    out_nmr = Parallel(n_jobs = -2)( #, timeout = 30)(
                   delayed(fill_objects)(
                                         gname,
                                         count,
                                         galfit_tmp_path,
                                         galfit_mask_path,
                                         out_png_dir,
                                         slurm = slurm
                                        )
                   for count, gname in enumerate(all_gname_tmp_out) if not gname.startswith("failed")
                                    )
    return out_nmr

# ==================================================================================================================

def main(**kwargs):
    
    # Dirs and name of the run for which the calculation was done
    run_dir       = kwargs.get("run_dir", os.getcwd())
    #in_dir        = kwargs.get("in_dir", pj(run_dir, "sparcfire-in"))
    tmp_dir       = kwargs.get("tmp_dir", pj(run_dir, "sparcfire-tmp"))
    out_dir       = kwargs.get("out_dir", pj(run_dir, "sparcfire-out"))
    basename      = kwargs.get("basename", "GALFIT")
    
    # Of course the important things
    slurm             = kwargs.get("slurm", False)
    dont_remove_slurm = kwargs.get("dont_remove_slurm", False)
    restart           = kwargs.get("restart", False)
    
    # For verbosity, default to capturing output
    # Keep both for clarity
    verbose        = kwargs.get("verbose", False)
    capture_output = kwargs.get("capture_output", True)
    
    final_pkl_file = f'{pj(run_dir, basename)}_output_nmr.pkl'
    if not slurm and exists(final_pkl_file):
        ans = input(f"Do you wish to delete the current final output nmr pickle file? Y/N\n{final_pkl_file}\n")
        if ans.upper() == "Y":
            _ = sp(f"rm -f {final_pkl_file}")
        else:
            print("You selected something other than y/Y. Not deleting. Please backup/rename/delete and try again.")
            print("Quitting.")
            sys.exit()
    
    outpath = pj(run_dir, out_dir)

    galfit_tmp_path = pj(run_dir, tmp_dir, "galfits")
    galfit_mask_path  = pj(run_dir, tmp_dir, "galfit_masks")
    out_png_dir = pj(outpath, "galfit_png")
    
    all_gname_tmp_out = [os.path.basename(i).replace("_galfit_out.fits","") for i in glob.glob(pj(galfit_tmp_path, "*_galfit_out.fits"))]
    
    out_nmr = []
    
    chunk_size = 1000
    if slurm and chunk_size > len(all_gname_tmp_out)*0.5:
        print("No need to slurm!")
        out_nmr = parallel_wrapper(galfit_tmp_path, galfit_mask_path, None, all_gname_tmp_out, slurm = False)
        slurm = False
        
    elif slurm:
        _, _, run_python = check_programs()
        python_slurm   = pj(_MODULE_DIR, "Utilities", "residual_via_slurm.py")
        slurm_file     = "parallel_residual_slurm"
        run_slurm      = "~wayne/bin/distrib_slurm"
        slurm_run_name = "CALCULATE_NMR"
        slurm_options  = "-M all"
        slurm_verbose  = "-v" if verbose else ""
        
        finished_pkl_num = 0
        if restart:
            finished_pkl_num = max(
                                   [int(os.path.basename(i).replace(basename, "")) 
                                    for i in glob.glob(pj(run_dir, f'{basename}*_output_nmr.pkl'))
                                   ]
                                  )
        
        with open(slurm_file, "w") as sf:
            for i, chunk in enumerate(range(chunk_size, len(all_gname_tmp_out) +  chunk_size, chunk_size)):
                if i < finished_pkl_num:
                    continue
                    
                gal_to_slurm = all_gname_tmp_out[chunk - chunk_size:][:chunk_size]
                #num_str = f"{i:0>3}"
                sf.write(f"{run_python} {python_slurm} {pj(run_dir, basename + str(i))} {galfit_tmp_path} {galfit_mask_path} {out_png_dir} {','.join(gal_to_slurm)}\n")
                
        print("Slurming")
        slurm_run_cmd = f"cat {slurm_file} | {run_slurm} {slurm_run_name} {slurm_options} {slurm_verbose}"
        completed = sp(slurm_run_cmd, capture_output = capture_output)
        
        # TODO: Add check if not on Slurm capable machine
        if completed.stderr:
            print(completed.stderr)
            # Rerun out_nmr
        
        else:
            all_output_pkl = glob.glob(pj(run_dir, f'{basename}*_output_nmr.pkl'))
            for file in all_output_pkl:
                out_nmr.extend(pickle.load(open(file, 'rb')))

    else:
        out_nmr = parallel_wrapper(galfit_tmp_path, galfit_mask_path, out_png_dir, all_gname_tmp_out, slurm = False)
            
    # out_nmr = Parallel(n_jobs = -2, timeout = 30)(
    #                    delayed(fill_objects)(
    #                                          (os.path.basename(i).rstrip("_galfit_out.fits") , i, count),
    #                                          galfit_mask_path,
    #                                          out_png_dir
    #                                         )
    #                    for count, i in enumerate(all_tmp_out) if not os.path.basename(i).startswith("failed")
    #                                               )
    
    # In the future, drop this in out_dir
    pickle_filename_temp = f'{pj(run_dir, basename)}_output_nmr_final.pkl'
    pickle.dump(out_nmr, open(pickle_filename_temp, 'wb'))
    
    if not dont_remove_slurm and slurm:
        _ = sp(f"rm -r \"$HOME/SLURM_turds/{slurm_run_name}\"", capture_output = capture_output)
        _ = sp(f"rm -f {pj(run_dir, basename)}*_output_nmr.pkl", capture_output = capture_output)
        _ = sp(f"rm -f {slurm_file}", capture_output = capture_output)
        
    pickle_filename = f'{pj(run_dir, basename)}_output_nmr.pkl'
    _ = sp(f"mv {pickle_filename_temp} {pickle_filename}", capture_output = capture_output)

        #output_fits_dict = dict(zip(out_nmr))

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
    
    parser.add_argument('-n', '--name',
                        dest     = 'basename',
                        action   = 'store',
                        default  = 'GALFIT',
                        #required = True, 
                        help     = 'Name of the run (for save file).')
    
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
    
    parser.add_argument('-v', '--verbose',
                        dest     = 'verbose', 
                        action   = 'store_const',
                        const    = True,
                        default  = False,
                        help     = 'Verbosity.')

    parser.add_argument(dest     = 'paths',
                        nargs    = "*",
                        type     = str,
                        help     = "IN-DIRECTORY TMP-DIRECTORY OUT-DIRECTORY from SpArcFiRe. \
                                    Must follow -in, -tmp, out or this won't work.")
    
    args = parser.parse_args() # Using vars(args) will call produce the args as a dict
    basename = args.basename
    slurm = args.slurm
    dont_remove_slurm = args.dont_remove_slurm
    restart = args.restart
    verbose = args.verbose
    
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
    # Don't actually need in_dir but to follow the structure of all the others...
    #in_dir  = absp(in_dir)
    tmp_dir = absp(tmp_dir)
    out_dir = absp(out_dir)
    
    kwargs = {
              #"in_dir"            : in_dir,
              "tmp_dir"           : tmp_dir,
              "out_dir"           : out_dir,
              "basename"          : basename,
              "slurm"             : slurm,
              "dont_remove_slurm" : dont_remove_slurm,
              "restart"           : restart,
              "verbose"           : verbose,
              "capture_output"    : not verbose
             }
    
    main(**kwargs)
       
