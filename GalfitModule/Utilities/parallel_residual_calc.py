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
import joblib

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

def fill_objects(gname, count, galfit_tmp_path, galfit_mask_path, out_dir = "", to_png = True, parallel = True):
    
    use_bulge_mask = False
    output_df = pd.DataFrame()
    
    report_num = 1000
    if parallel:
        report_num = 100
        
    if not count % report_num:
        print(count, gname)
    
    try:
        fits_filename = f"{gname}_galfit_out.fits"
        fits_file = OutputFits(pj(galfit_tmp_path, fits_filename))
    except Exception as e:
        print(f"There was an issue opening galaxy {gname}. Continuing...")
        print(e)
        #return gname, None, None, None
        return None
    
    #last_feedme_file = pj(out_dir, gname, f"{gname}.in")
    #if not exists(last_feedme_file):
    #feedme_files = glob.glob(pj(out_dir, gname, f"*.in"))
    mask_fits_name = pj(galfit_mask_path, f"{gname}_star-rm.fits")
    
    # For simultaneous fitting
#     if feedme_files: 
#         feedme_files.sort(key = os.path.getmtime)
#         last_feedme_file = feedme_files[-1]

#         feedme = FeedmeContainer()
#         feedme.from_file(last_feedme_file)
        
        # if feedme.bulge.position == fits_file.feedme.bulge.position:       
        #     fits_file.feedme.header.region_to_fit = feedme.header.region_to_fit
        #     mask_fits_name = feedme.header.pixel_mask
        
       # fits_file.feedme.header.update_param_values()
        
           
    #if "sim_fitting" in fits_file.feedme.header.input_image:
        #mask_fits_name = pj(sf_mask_path, star_mask_name)
    
    try:
        mask_fits_file = FitsFile(mask_fits_name)
    except Exception as e:
        print(f"There was an issue opening star mask for galaxy {gname}. Proceeding without mask...")
        # Logic implemented to handle None
        mask_fits_file = None #np.zeros((500,500))
    
    # If skip is enabled then arms are not fit so bulge masking doesn't make sense
    c_types = [comp.component_type for comp in fits_file.feedme.components.values()]
    if out_dir and "power" in c_types:
        _ = fits_file.generate_bulge_mask(pj(out_dir, gname, f"{gname}.csv"))
        use_bulge_mask = True
        
    masked_residual_normalized = fits_file.generate_masked_residual(mask_fits_file, use_bulge_mask = use_bulge_mask)
    if masked_residual_normalized is None:
        print(f"Could not calculate nmr for galaxy {gname}. Continuing...")
        #return gname, None, None, None
        return None
    
    # Doesn't work on Openlab (sadly)
    # Keep this in for actual parallelizing since it's a PITA to read booleans
    # from command line (use the default value to our advantage)
    if sp(f"hostname").stdout.split(".")[0] == "bayonet-09" and to_png:
        if isinstance(to_png, str):
            fits_file.to_png(out_png_dir = to_png)
        elif out_dir:
            fits_file.to_png(out_png_dir = pj(out_dir, "galfit_png"))
    
    #return gname, fits_file.nmr, fits_file.kstest.pvalue, fits_file.kstest.statistic#, fits_file.nmrr
    reloaded  = OutputFits(pj(galfit_tmp_path, fits_filename))
    output_df = reloaded.feedme.to_pandas()
    output_df["gname"]   = gname
    output_df["NMR"]     = reloaded.model.header.get("NMR", None) 
    output_df["KS_P"]    = reloaded.model.header.get("KS_P", None)
    output_df["KS_STAT"] = reloaded.model.header.get("KS_STAT", None)
    
    return output_df
    
# ==================================================================================================================

def parallel_wrapper(galfit_tmp_path, galfit_mask_path, out_dir, to_png, all_gname_tmp_out, parallel = True):
    out_df = Parallel(n_jobs = -2)( #, timeout = 30)(
                   delayed(fill_objects)(
                                         gname,
                                         count,
                                         galfit_tmp_path,
                                         galfit_mask_path,
                                         out_dir,
                                         to_png,
                                         parallel = parallel
                                        )
                   for count, gname in enumerate(all_gname_tmp_out) if not gname.startswith("failed")
                                    )
        
    #return {i[0] : tuple(i[1:]) for i in out_nmr}
    return pd.concat(out_df)

# ==================================================================================================================

def main(**kwargs):
    
    # Dirs and name of the run for which the calculation was done
    #run_dir       = kwargs.get("run_dir", os.getcwd())
    #in_dir        = kwargs.get("in_dir", pj(run_dir, "sparcfire-in"))
    tmp_dir       = kwargs.get("tmp_dir", pj(run_dir, "sparcfire-tmp"))
    out_dir       = kwargs.get("out_dir", pj(run_dir, "sparcfire-out"))
    basename      = kwargs.get("basename", "GALFIT")
    #pkl_end_str   = kwargs.get("pkl_end_str", "output_results")
    # No need to make this an option... for now(?)
    pkl_end_str = "output_results"
    
    # Of course the important things
    parallel          = kwargs.get("parallel", 1)
    dont_remove_slurm = kwargs.get("dont_remove_slurm", False)
    restart           = kwargs.get("restart", False)
    #simultaneous_fitting = kwargs.get("simultaneous_fitting", True)
    
    # For verbosity, default to capturing output
    # Keep both for clarity
    verbose        = kwargs.get("verbose", False)
    capture_output = kwargs.get("capture_output", True)
    
    final_pkl_file = f'{pj(out_dir, basename)}_{pkl_end_str}.pkl'
    if not parallel and exists(final_pkl_file):
        ans = input(f"Do you wish to delete the current final output nmr pickle file? Y/N\n{final_pkl_file}\n")
        if ans.upper() == "Y":
            _ = sp(f"rm -f {final_pkl_file}")
        else:
            print("You selected something other than y/Y. Not deleting. Please backup/rename/delete and try again.")
            print("Quitting.")
            sys.exit()

    galfit_tmp_path = pj(tmp_dir, "galfits")
    galfit_mask_path  = pj(tmp_dir, "galfit_masks")
    #sf_mask_path = pj(tmp_dir, "sim_fitting", "sparcfire-tmp", "galfit_masks")
    out_png_dir = pj(out_dir, "galfit_png")
    
    all_gname_tmp_out = [os.path.basename(i).replace("_galfit_out.fits","") for i in glob.glob(pj(galfit_tmp_path, "*_galfit_out.fits"))]
    
    #out_nmr = {}
    out_df = pd.DataFrame()
    
    chunk_size = 1000
    if parallel and chunk_size > len(all_gname_tmp_out)*0.5:
        print("No need to (massively) parallelize!")
        out_df = parallel_wrapper(galfit_tmp_path, 
                                  galfit_mask_path,
                                  out_dir, 
                                  True, 
                                  all_gname_tmp_out,
                                  parallel = False
                                 )
        parallel = False
        
    elif parallel:
        _, _, run_python = check_programs()
        python_parallel   = pj(_MODULE_DIR, "Utilities", "residual_via_parallel.py")
        parallel_file     = "parallel_combine_residual"
        
        if exists(parallel_file):
            _ = sp(f"rm -f {parallel_file}", capture_output = capture_output)
        
        if parallel == 1:
            run_parallel      = "/home/sana/bin/parallel"
            parallel_run_name = ""
            parallel_options  = joblib.cpu_count()
            
        elif parallel == 2:
            # TODO: Add check if not on Slurm capable machine
            run_parallel      = "~wayne/bin/distrib_slurm"
            parallel_run_name = "CALCULATE_COMBINE_GALFIT_RESULTS"
            parallel_options  = "-M all"
            
        parallel_verbose  = "-v" if verbose else ""
        
        finished_pkl_num = 0
        if restart:
            check_output_pkl = [int(os.path.basename(i).split("_")[0].replace(basename, "")) 
                                for i in find_files(tmp_dir, f'{basename}*_{pkl_end_str}.pkl', "f")
                                if os.path.basename(i).split("_")[0].replace(basename, "")
                                ]
            if check_output_pkl:
                finished_pkl_num = max(check_output_pkl)
        
        # Use count for restarting purposes i.e. if all pkl files have been generated 
        # but just need to be combined
        count = 0
        with open(parallel_file, "w") as sf:
            for i, chunk in enumerate(range(chunk_size, len(all_gname_tmp_out) +  chunk_size, chunk_size)):
                if i < finished_pkl_num:
                    continue
                
                gal_to_parallel = all_gname_tmp_out[chunk - chunk_size:][:chunk_size]
                #num_str = f"{i:0>3}"
                sf.write(f"{run_python} {python_parallel} {pj(tmp_dir, basename + str(i))} {galfit_tmp_path} {galfit_mask_path} {out_dir} {out_png_dir} {','.join(gal_to_parallel)}\n")
                count += 1
        
        if count:
            print("parallelizing")
            parallel_run_cmd = f"cat {parallel_file} | {run_parallel} {parallel_run_name} {parallel_options} {parallel_verbose}"
            completed = sp(parallel_run_cmd, capture_output = capture_output)

            if completed.stderr:
                print(completed.stderr)
                # Rerun out_nmr
        
        #else:
        all_output_pkl = glob.glob(pj(tmp_dir, f'{basename}*_{pkl_end_str}.pkl'))
        #_ = [out_nmr.update(pickle.load(open(file, 'rb'))) for file in all_output_pkl if os.path.basename(file) != f"{basename}_output_nmr.pkl"]
        out_df = pd.concat(
                           [pd.read_pickle(file) for file in all_output_pkl 
                            if os.path.basename(file) != f"{basename}_{pkl_end_str}.pkl"
                           ]
                          ) 

    else:
        out_df = parallel_wrapper(galfit_tmp_path, 
                                  galfit_mask_path,
                                  out_dir,
                                  out_png_dir, 
                                  all_gname_tmp_out, 
                                  parallel = False
                                 )
        
    # In the future, drop this in out_dir
    #pickle_filename_temp = f'{pj(out_dir, basename)}_output_nmr_final.pkl'
    print(f"Outputting results to {final_pkl_file}")
    #pickle.dump(out_nmr, open(final_pkl_file, 'wb'))
    
    # For when I do it in one of the other scripts
    if "gname" in out_df.columns:
        out_df.set_index("gname", inplace = True)
    
    out_df.to_pickle(final_pkl_file)
    
    if not dont_remove_slurm and parallel:
        _ = sp(f"rm -r \"$HOME/SLURM_turds/{parallel_run_name}\"", capture_output = capture_output)
        _ = sp(f"rm -f {pj(tmp_dir, basename)}*_{pkl_end_str}.pkl", capture_output = capture_output)
        _ = sp(f"rm -f {parallel_file}", capture_output = capture_output)
        
    #pickle_filename = f'{pj(out_dir, basename)}_output_nmr.pkl'
    
    #_ = sp(f"mv {pickle_filename_temp} {pickle_filename}", capture_output = capture_output)

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
    
    parser.add_argument('-p', '--parallel',
                        dest     = 'parallel',
                        action   = 'store',
                        type     = int,
                        choices  = range(0,3),
                        default  = 1,
                        help     = 'Run algorithm with/without intensive parallelization. Defaults to on machine parallel.\nOptions are:\n\t\
                                    0: in serial,\n\t\
                                    1: on machine parallel,\n\t\
                                    2: cluster computing via SLURM'
                       )

    parser.add_argument('-drS', '--dont-remove-slurm',
                        dest     = 'dont_remove_slurm',
                        action   = 'store_const',
                        const    = True,
                        default  = False,
                        help     = 'Choose NOT to remove all old slurm files (they may contain basic info about each fit but there will be a bunch!)')
    
    # parser.add_argument('-nsf', '--no-simultaneous-fitting',
    #                     dest     = 'simultaneous_fitting',
    #                     action   = 'store_const',
    #                     const    = False,
    #                     default  = True,
    #                     help     = 'Turn off simultaneous fitting.'
    #                    )
    
    parser.add_argument('-r', '--restart',
                        dest     = 'restart', 
                        action   = 'store_const',
                        const    = True,
                        default  = False,
                        help     = 'Restart script on the premise that some have already run (likely in parallel).')
    
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
    parallel = args.parallel
    dont_remove_slurm = args.dont_remove_slurm
    restart = args.restart
    verbose = args.verbose
    #simultaneous_fitting = args.simultaneous_fitting
    
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
              "tmp_dir"              : tmp_dir,
              "out_dir"              : out_dir,
              "basename"             : basename,
              "parallel"             : parallel,
              "dont_remove_slurm"    : dont_remove_slurm,
              "restart"              : restart,
              #"simultaneous_fitting" : simultaneous_fitting,
              #"sf_tmp_dir"           : sf_tmp_dir
              "verbose"              : verbose,
              "capture_output"       : not verbose
             }
    
    main(**kwargs)
       
