#!/usr/bin/env python
# coding: utf-8

# **Author: Matthew Portman**
# 
# **Date (Github date will likely be more accurate): 4/17/23**

# In[19]:


import sys
import os
from os.path import join as pj

from astropy.io import fits

import argparse
import shutil
import subprocess

import time
import pickle
import joblib


# In[3]:


# For debugging purposes
from IPython import get_ipython
def in_notebook():
    ip = get_ipython()
    
    if ip:
        return True
    else:
        return False


# In[4]:


_HOME_DIR = os.path.expanduser("~")
if in_notebook():
    # We hardcode the directory name because the notebooks 
    # are used for prototyping and debugging. 
    _SPARCFIRE_DIR = pj(_HOME_DIR, "sparcfire_matt") 
    _MODULE_DIR    = pj(_SPARCFIRE_DIR, "GalfitModule")
else:
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
from Functions.helper_functions import *
import Utilities.parallel_residual_calc as parallel_residual_calc
import Utilities.combine_via_parallel as combine_via_parallel

# This should give me numpy and pandas and whatnot
# also gives this, from os.path import join as pj
from sparc_to_galfit_feedme_gen import *
import go_go_galfit


# In[13]:


if __name__ == "__main__":
    # TODO: Add run_dir? and default to cwd if not specified
    
    # Force >python 3.7 for various compatabilities
    out_str = "\t Python3.7 or greater required! Exitting without generating feedmes..."
    assert sys.version_info >= (3, 7), out_str
    
    cwd = absp(os.getcwd()) # Doesn't work *in* notebook
    old_cwd = absp(cwd) # Strings are immutable
    
    username = os.environ["USER"]
    
    USAGE = f"""USAGE:

    python3 ./{sys.argv[0]} [OPTION] [[RUN-DIRECTORY] IN-DIRECTORY TMP-DIRECTORY OUT-DIRECTORY]
    
    OPTIONS =>[-p | --parallel]
              [-drs | --dont-remove-slurm]
              [-t  | --tmp]
              [-ac | --aggressive-clean]
              [-NS | --num-steps] 
              [-r | --restart]
              [-RrG | --rerun-galfit]
              [-v | --verbose]
              [-n | --name]

    This script is the wrapping script for running GALFIT using SpArcFiRe to inform 
    the input. By default, it runs from the RUN (or current) directory and uses the
    '-in' '-tmp' and '-out' directories as specified or otherwise defaults to 
    'sparcfire-in', 'sparcfire-tmp', 'sparcfire-out'. 

    Please do not specify symlinks for the above, they discomfort the programmer.
    """
    
    parser = argparse.ArgumentParser(description = USAGE)
    
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

    parser.add_argument('-drs', '--dont-remove-slurm',
                        dest     = 'dont_remove_slurm',
                        action   = 'store_const',
                        const    = True,
                        default  = False,
                        help     = 'Choose NOT to remove all old slurm files (they may contain basic info about each fit but there will be a bunch!)'
                       )
    
    parser.add_argument('-t', '--tmp',
                        dest     = 'run_from_tmp',
                        action   = 'store_const',
                        const    = True,
                        default  = False,
                        help     = 'Indicates to the program a run from a directory in /tmp/ local to each cluster machine.\n\
                                    WARNING: either output to a different location or copy from tmp to said location \
                                    under the assumption that tmp will be wiped at some point in the near future.'
                       )
    
    parser.add_argument('-ac', '--aggressive-clean',
                        dest     = 'aggressive_clean',
                        action   = 'store_const',
                        const    = True,
                        default  = False,
                        help     = 'Aggressively clean-up directories, removing -in, temp output, psf, and mask files after galfit runs'
                       )
    
    parser.add_argument('-NS', '--num-steps',
                        dest     = 'steps', 
                        action   = 'store',
                        type     = int,
                        choices  = range(1,4),
                        default  = 2,
                        help     = 'Run GALFIT using step-by-step component selection (up to 3), i.e.\n\t\
                                    1: Bulge + Disk + Arms,\n\t\
                                    2: Bulge -> Bulge + Disk + Arms,\n\t\
                                    3: Bulge -> Bulge + Disk -> Bulge + Disk + Arms'
                       )
    
    parser.add_argument('-r', '--restart',
                        dest     = 'restart',
                        action   = 'store_const',
                        const    = True,
                        default  = False,
                        help     = 'Restart control script on the premise that some have already run (likely in parallel).'
                       )
    
    parser.add_argument('-RrG', '--rerun-galfit',
                        dest     = 'rerun', 
                        action   = 'store_const',
                        const    = True,
                        default  = False,
                        help     = 'Run GALFIT again after the final fit to hopefully refine said fit.'
                       )
    
    parser.add_argument('-n', '--basename',
                        dest     = 'basename', 
                        action   = 'store',
                        type     = str,
                        default  = "GALFIT",
                        help     = 'Basename of the output results pkl file ([name]_output_results.pkl).'
                       )
    
    parser.add_argument('-v', '--verbose',
                        dest     = 'verbose', 
                        action   = 'store_const',
                        const    = True,
                        default  = False,
                        help     = 'Verbose output for all bash commands in control script.'
                       )
    
    parser.add_argument(dest     = 'paths',
                        nargs    = "*",
                        type     = str,
                        help     = "RUN-DIRECTORY [IN-DIRECTORY TMP-DIRECTORY OUT-DIRECTORY] from SpArcFiRe. \
                                    SpArcFiRe directories should follow -in, -tmp, -out."
                       )
    
    if not in_notebook():
        args              = parser.parse_args() # Using vars(args) will call produce the args as a dict
        num_steps         = args.steps
        parallel          = args.parallel
        dont_remove_slurm = args.dont_remove_slurm
        run_from_tmp      = args.run_from_tmp
        aggressive_clean  = args.aggressive_clean
        
        rerun             = args.rerun
        restart           = args.restart
        basename          = args.basename
        
        verbose           = args.verbose
        capture_output    = not args.verbose
        
        # if num_steps not in range(1,4):
        #     print("The number of steps you selected cannot be used!")
        #     print("Using two.")
        #     print()
        #     num_steps = 2

        if len(args.paths) == 1:
            cwd = args.paths[0]
            in_dir = pj(cwd, "sparcfire-in")
            tmp_dir = pj(cwd, "sparcfire-tmp")
            out_dir = pj(cwd, "sparcfire-out")
            
        elif len(args.paths) == 3:
            in_dir, tmp_dir, out_dir = args.paths[0], args.paths[1], args.paths[2]
            #print(f"Paths are, {in_dir}, {tmp_dir}, {out_dir}")
            
        elif len(args.paths) == 4:
            cwd, in_dir, tmp_dir, out_dir = args.paths[0], args.paths[1], args.paths[2], args.paths[3]
            #print(f"Paths are, {in_dir}, {tmp_dir}, {out_dir}")
            
        else:
            in_dir = pj(cwd, "sparcfire-in")
            tmp_dir = pj(cwd, "sparcfire-tmp")
            out_dir = pj(cwd, "sparcfire-out")
            print(f"Paths incorrectly specified, defaulting to {cwd} (-in, -tmp, -out)...")
            print(f"{in_dir}\n{tmp_dir}\n{out_dir}")
            print()
            
        check_dir_names = [1 for i in (in_dir, tmp_dir, out_dir) if "-" not in i ]
        if check_dir_names:
            raise Exception("Directory paths must end in '-in' '-tmp' and '-out'")
            
    else:
        parallel = 0
        rerun = ""
        num_steps = 2
        # Avoid some... nasty surprises for when debugging
        restart = True
        verbose = False
        capture_output = True
        
        cwd = cwd.replace("ics-home", username)
        in_dir = pj(cwd, "sparcfire-in")
        tmp_dir = pj(cwd, "sparcfire-tmp")
        out_dir = pj(cwd, "sparcfire-out")
        
        sys.path.append(pj(_HOME_DIR, ".local", "bin"))
        
    # Making these absolute paths
    cwd     = absp(cwd)
    in_dir  = absp(in_dir)
    tmp_dir = absp(tmp_dir)
    out_dir = absp(out_dir)
    
    # Changing to specified working dir
    os.chdir(cwd)


# ## Checking file directories and installed programs

# In[ ]:


if __name__ == "__main__":
    # Checking dirs
    
    assert all((exists(in_dir), exists(tmp_dir), exists(out_dir))), \
           f"Cannot find one of: in, tmp, or out directories:\n{in_dir}\n{tmp_dir}\n{out_dir}\n" \
           "Do those look right?"
    
#     # to redirect stderr to /dev/null as well:
#     #subprocess.run(['ls', '-l'], stderr=subprocess.DEVNULL)

    if not in_notebook():
        # This is in helper_functions
        run_galfit, run_fitspng, run_python = check_programs()
#         # This seems to work in Python directly so I'm leaving it as-is
#         # Checking galfit
#         run_galfit = shutil.which("galfit")
#         #run_galfit = response.stdout.strip()
        
#         # Checking fitspng
#         run_fitspng   = shutil.which("fitspng")
#         fitspng_param = "0.25,1" #1,150"
        
#         # Checking exact python3 call
#         run_python = shutil.which("python3")

    else:
        run_galfit  = pj(_HOME_DIR, ".local/bin/galfit")
        run_fitspng = pj(_HOME_DIR, ".local/bin/fitspng")
        run_python  = "/opt/conda/bin/python3"


# ## Setting up directories and handy variables

# In[ ]:


if __name__ == "__main__":
    # Setting up paths and variables
    tmp_fits_dir    = pj(tmp_dir, "galfits")
    tmp_masks_dir   = pj(tmp_dir, "galfit_masks")
    # Now dropping these in the individual galaxy folders
    #tmp_psf_dir     = pj(tmp_dir, "psf_files")
    tmp_png_dir     = pj(tmp_dir, "galfit_png")
    #need_masks_dir  = pj(tmp_dir, "need_masks")
    
    #all_galfit_out = pj(out_dir, "all_galfit_out")
    out_png_dir     = pj(out_dir, "galfit_png")
    
    # Should be same as SpArcFiRe's but I'm packaging it with just in case
    star_removal_path = pj(_MODULE_DIR, "star_removal")
    
    if parallel == 1:
        # CPU Parallel
        pipe_to_parallel_cmd = "/home/sana/bin/parallel"
        
    elif parallel == 2:
        # SLURM/Cluster Computing
        pipe_to_parallel_cmd = "~wayne/bin/distrib_slurm"
    
    if not restart:
        # Remove old
        print("Removing old results in tmp folder (if they exist).")
        _ = sp(f"rm -rf {tmp_fits_dir}", capture_output = capture_output)
        # try:
        #     shutil.rmtree(tmp_fits_dir)
        # except OSError as e:
        #     pass
    else:
        print("Restarting the run (still have to find everything).")

    # Making sub-directories
    _ = [os.mkdir(i) for i in (tmp_fits_dir, 
                               tmp_masks_dir, 
                               #tmp_psf_dir, 
                               tmp_png_dir, 
                               out_png_dir
                              )#, 
                               #need_masks_dir) 
         if not exists(i)
        ]
    
    # makedirs will make both at once, handy!
    #if not exists(out_png): os.makedirs(out_png)


# ## Running sextractor (if necessary)

# In[ ]:


if __name__ == "__main__":    
    # Grabbing list of file names and masks with bash variable expansion
    print("Setting up, finding files...")
    # input_filenames = glob.glob(pj(in_dir, "*.fits"))
    # output_folders  = glob.glob(pj(out_dir, "123*/"))
    # star_masks      = glob.glob(pj(tmp_dir, "*_star-rm.fits"))
    input_filenames = find_files(in_dir, "*.fits", "f")
    if not input_filenames:
        print(f"No input files found in {in_dir}. Is something wrong?")
        sys.exit()
        
    output_folders  = [i for i in find_files(out_dir, "*", "d") 
                       if os.path.basename(i) != os.path.basename(out_dir) and
                          os.path.basename(i) != "galfit_png"
                      ]
    
    star_masks      = find_files(tmp_dir, "*_star-rm.fits", "f")
    
    # The ONLY reason we need this is because of how remove_stars_with... works
    # If we change that, the whole mess of code and operations that use can go away
    # the need_masks_dir can go away
    #in_need_masks   = glob.glob(pj(need_masks_dir, "*.fits"))
#     in_need_masks    = find_files(need_masks_dir, "*.fits", "f")
    
#     check_need_masks = set(output_folders).intersection(
#                        set([os.path.basename(i) for i in in_need_masks])
#                                                         )
#     if not check_need_masks:
#         _ = sp(f"rm -rf {need_masks_dir}", capture_output = capture_output)
#         os.mkdir(need_masks_dir)
#         in_need_masks = []
    
    # Sparcfire plops them into tmp so I just move them to where I need them
    if star_masks and not restart:
        sp(f"mv {pj(tmp_dir,'*_star-rm.fits')} {tmp_masks_dir}", capture_output = capture_output)
        
    #star_masks = glob.glob(pj(tmp_masks_dir, "*_star-rm.fits"))
    star_masks      = find_files(tmp_masks_dir, "*_star-rm.fits", "f")


# In[ ]:


if __name__ == "__main__":
    generate_starmasks = True
    if not restart:
        #generate_starmasks = True
        if not exists(pj(tmp_masks_dir, 'remove_stars_with_sextractor.log')):
            # remove_stars_with_sextractor needs this available before it can log
            _ = sp(f"touch {pj(tmp_masks_dir, 'remove_stars_with_sextractor.log')}", capture_output = capture_output)
        
#         try:
#             star_removal_path = pj(os.environ["SPARCFIRE_HOME"], "star_removal")
            
#         except KeyError as k:
#             print("SPARCFIRE_HOME environment variable is not set.")
#             print("Will use 'remove_stars_with_sextractor.py' from GalfitModule")
#             print("If they have already been generated, ignore this message.")
#             print()
            
        # else:
        
        # Compare against output folders because extra observations may be in input directory
        # GALFIT can only run on what's there... I mean there are defaults for SpArcFiRe
        # but this is the more appropriate choice. 
        if len(output_folders) != len(star_masks):
            print("The temp directory has a different number of star masks than the number of output directories.") 
        else:
            generate_starmasks = False
            
# DEPRECATED
                           
#             print("Getting things ready to generate star masks...")
#             #galaxy_folder_names = [os.path.basename(i.rstrip("/")) for i in output_folders]
#             #star_mask_names = [os.path.basename(i) for i in star_masks]
            
#             parallel_copy_input = "parallel_copy_inputs"
#             if parallel:
#                 if exists(parallel_copy_input):
#                     _ = sp(f"rm {parallel_copy_input}", capture_output = capture_output)

#                 _ = sp(f"touch {parallel_copy_input}", capture_output = capture_output)

#                 sci = open(parallel_copy_input, "a")
                
#             # TODO: Could also tar then transfer(?)
#             if len(in_need_masks) + len(star_mask_names) != len(output_folders):
#                 print("Copying input galaxies without masks to 'need_masks' in tmp folder")
#                 for gname in galaxy_folder_names:
#                     star_mask_filename = f"{gname}_star-rm.fits"
#                     in_file = f"{gname}.fits"
#                     cp_cmd = f"cp -uv {pj(in_dir, in_file)} {need_masks_dir}"

#                     if in_file in in_need_masks:
#                         continue
                    
#                     # This is probably unecessary but keeping just in case
#                     elif in_file.split(".fits")[0] not in galaxy_folder_names:
#                         _ = sp(f"rm -f {pj(need_masks_dir, in_file)}", capture_output = capture_output)

#                     elif star_mask_filename not in star_mask_names:

#                         if parallel:
#                             cp_cmd += "\n"
#                             sci.write(cp_cmd)

#                         else:
#                             try:
#                                 shutil.copy2(pj(in_dir, in_file), need_masks_dir)
#                             except FileNotFoundError:
#                                 print(f"Could not find {gname} in {in_dir}. Continuing...")
#                             #result = sp(f"cp -u {pj(need_masks_dir, in_file)}", capture_output = True)
#                             # if result.stderr:
#                             #     print(f"Could not find {gname} in {in_dir}, {result.stderr}. Continuing...")
#                     else:
#                         # To remove extras, print out the rest of this list
#                         star_mask_names.remove(star_mask_filename)
#                         _ = sp(f"rm -f {pj(need_masks_dir, in_file)}", capture_output = capture_output)

#                 if parallel:
#                     sci.close()

#                     extra_parallel = ""
#                     if verbose:
#                         extra_parallel = "-v"

#                     parallel_run_name = "COPYING_INPUT_FILES"
#                     parallel_run_cmd = f"cat {parallel_copy_input} | {pipe_to_parallel_cmd} {parallel_run_name} -M all {extra_parallel}"
#                     # Running without timeout for now
#                     print("Performing the copy with parallel...")
#                     _ = sp(f"{parallel_run_cmd}", capture_output = capture_output)
                
#             print("Generating starmasks...")
#             os.chdir(star_removal_path)
#             out_text = sp(f"python3 {pj(star_removal_path, 'remove_stars_with_sextractor.py3')} {need_masks_dir} {tmp_masks_dir}", capture_output = capture_output)
#             os.chdir(cwd)
            
#             #star_masks = glob.glob(pj(tmp_masks_dir, "*_star-rm.fits"))
#             star_masks = find_files(tmp_masks_dir, "*_star-rm.fits", "f")
#             issue_masking_gnames = [os.path.basename(sm).rstrip("_star-rm.fits") 
#                                     for sm in star_masks
#                                     if sm not in galaxy_folder_names
#                                    ]
            
#             # GALFIT can still run on these galaxies without a mask 
#             # so I don't have to worry about letting it know
#             if issue_masking_gnames:
#                 issue_file = pj(tmp_dir, "issue_masking.txt")
                
#                 # This is so that I don't have to generate all star masks again...
#                 with open(issue_file, "w") as f:
#                     for im in issue_masking_gnames:
#                         _ = sp(f"touch {pj(tmp_masks_dir, im)}-failed_star-rm.fits")
#                         f.write(im)
                    
#                 print("Could not produce a star mask for all galaxies")
#                 print(f"See {issue_file} for more details, creating empties...")
            
#             if out_text.stderr:
#                 print(f"Something went wrong running 'remove_stars_with_sextractor.py'! Printing debug info...")
#                 print(out_text)
#                 #print(type(out_text.stderr))
#             else:
#                 print("Star masks have been generated successfully.")
#                 print()
                
    # else:
    #     generate_starmasks = False
    #     print("Star masks have already been generated, proceeding.")


# ## Galfitting!

# In[ ]:


def write_to_parallel(cwd, 
                   kwargs_main, 
                   galfit_script_name = pj(_MODULE_DIR, "go_go_galfit.py"), 
                   parallel_file = pj(cwd, "parallel_cmd_file"),
                   # determined by parallel proc limit 
                   chunk_size = 20
                  ):
    _, _, run_python = check_programs()
    kwargs_in        = deepcopy(kwargs_main)
    
    print(f"Generating input file for parallelization in {cwd}: {parallel_file}")
    with open(pj(cwd, parallel_file), "w") as scf:
        #for gname in kwargs_main["galaxy_names"]:
        for i, chunk in enumerate(range(chunk_size, len(kwargs_main["galaxy_names"]) +  chunk_size, chunk_size)):
            chunk_o_galaxies = kwargs_main["galaxy_names"][chunk - chunk_size:][:chunk_size]
            kwargs_in["galaxy_names"] = ",".join(chunk_o_galaxies)
            
            # kwargs_in["petromags"] = ",".join(kwargs_main["petromags"][chunk - chunk_size:][:chunk_size])
            # kwargs_in["bulge_axis_ratios"] = ",".join(kwargs_main["bulge_axis_ratios"][chunk - chunk_size:][:chunk_size])
            

            cmd_str = ""
            for k,v in kwargs_in.items():
                cmd_str += f"{k}={v} "
                
            # Good thing dictionaries retain order now, *whistles innocently*
            scf.write(f"{run_python} {galfit_script_name} {cmd_str}\n")

    sp(f"chmod a+x {parallel_file}", capture_output = capture_output)


# In[ ]:


def check_galfit_out_hangups(tmp_fits_dir, out_dir, kwargs_main):
    # For hang-ups, check if final copy has occurred
    # Check for file which indicates galfit outputs nothing at all
    # This avoids conflating ones which haven't run and ones which have nothing to show for it
    kwargs_main["galaxy_names"] = [gname for gname in kwargs_main["galaxy_names"]
                                   if not exists(f"{pj(out_dir, gname, gname)}_galfit_out.fits")
                                   and not exists(f"{pj(tmp_fits_dir, 'failed_' + gname)}_galfit_out.fits")]
    # TODO: CHECK FOR FINAL .IN FILE TO ENSURE IT TRIED TO RUN THE FULL FIT
    
    # because of mutability I don't need to do this but
    # because it's good practice...
    return kwargs_main


# In[ ]:


def write_failed(failed_dir = cwd, failures = []):
    if failures:
        fail_filename = "galfit_failed.txt"
        fail_filepath = pj(failed_dir, fail_filename)
        print(f"{len(failures)} galax(y)ies completely failed. Writing the list of these to {fail_filepath}")
        with open(fail_filepath, "w") as ff:
            ff.writelines("\n".join(failures))
            ff.write("\n")


# In[ ]:


if __name__ == "__main__":    
    #print("Finding all galaxies...")
    # galaxy_names = [os.path.basename(i).rstrip(".fits") 
    #                 for i in input_filenames]
    
    # Replace is safer than rstrip because rstrip could remove those characters
    # from the end of the filename
    galaxy_names = [i.replace(".fits", "") for i in input_filenames]
    
    if not restart:
        # Backing up old fits
        # So we can check for new ones
        # and also backup previous runs if there's an accidental overwrite
        print("Backing up previous output (if found)")
        for gname in galaxy_names:
            new_out = pj(out_dir, gname, f"{gname}_galfit_out.fits")
            old_out = pj(out_dir, gname, f"{gname}_galfit_out_old.fits")
            if exists(new_out):
                shutil.move(new_out, old_out)
                
    #print("Reading in SDSS information")
    # TODO: don't hardcode this
    # ASSUME ALL ARE IN SAME COLOR BAND (FOR NOW)
    # with fits.open(pj(in_dir, galaxy_names[0] + ".fits")) as f:
    #     #SURVEY = SDSS-r  DR7
    #     color = f[0].header["SURVEY"].split()[0][-1]
        
#     gzoo_file = pj(_HOME_DIR, "kelly_stuff", "Kelly-Final-GZ-all-psfield-incl.csv")
#     try:
#         gzoo_data = pd.read_csv(gzoo_file, 
#                                 #sep = "\t", 
#                                 usecols = ["GZ_dr8objid", "petroMag_r", "deVAB_r"],
#                                 index_col = "GZ_dr8objid", 
#                                 dtype = {"GZ_dr8objid" : str}
#                                )
        
#         petromags = [str(gzoo_data.loc[gname, f"petroMag_{color}"]) if gname in gzoo_data.index else "16" 
#                      for gname in galaxy_names]
        
#         bulge_axis_ratios = [str(gzoo_data.loc[gname, f"deVAB_{color}"]) if gname in gzoo_data.index else "1"
#                              for gname in galaxy_names]
        
#     except FileNotFoundError:
#         print(f"Could not find {gzoo_file}. Proceeding.")
#         gzoo_data = None

    kwargs_main = {"cwd"                : cwd,
                   "in_dir"             : in_dir,
                   "tmp_dir"            : tmp_dir,
                   "out_dir"            : out_dir,
                   "num_steps"          : num_steps,
                   "rerun"              : rerun,
                   "parallel"           : parallel,
                   "verbose"            : verbose,
                   "capture_output"     : capture_output,
                   "generate_starmasks" : generate_starmasks,
                   "run_from_tmp"       : run_from_tmp,
                   "aggressive_clean"   : aggressive_clean,
                   # "petromags"          : petromags,
                   # "bulge_axis_ratios"  : bulge_axis_ratios,
                   # Keep this last just in case
                   "galaxy_names"       : galaxy_names
                  }
    
    # In case we're running back to back, this will reduce galaxy_names appropriately
    kwargs_main = check_galfit_out_hangups(tmp_fits_dir, out_dir, kwargs_main)
    parallel_file = "parallel_cmd_file"
    
#    raise(AssertionError())
    # One at a time
    if parallel:
        if not kwargs_main["galaxy_names"] and not restart:
            print("No galaxies to fit, exitting.")
            sys.exit()
            
        #print("Piping to parallel")
        print(f"{len(kwargs_main['galaxy_names'])} galaxies")
        
        chunk_size = 5
        if parallel == 1:
            # For CPU parallel
            parallel_run_name = ""#"GALFITTING"
            parallel_options  = joblib.cpu_count()
            parallel_verbose  = ""
            chunk_size = len(kwargs_main["galaxy_names"])//joblib.cpu_count() + 1
            # Two whole days for big runs
            timeout = 2880 # Minutes
            
        elif parallel == 2:
            # For SLURM/Cluster Computing
            parallel_run_name = "GALFITTING"
            # Slurm needs different timeout limits
            timeout = 60 # Minutes
            parallel_options  = f"-M all -t {timeout}"
            parallel_verbose  = "-v" if verbose else ""
            chunk_size = 20

            
        parallel_run_cmd = f"cat {parallel_file} | {pipe_to_parallel_cmd} {parallel_run_name} {parallel_options} {parallel_verbose}"
        
        if not restart:
            write_to_parallel(cwd, kwargs_main, parallel_file = parallel_file, chunk_size = chunk_size)
            print("Galfitting via parallelization...")
            try:
                sp(f"{parallel_run_cmd}", capture_output = capture_output, timeout = 60*(timeout + 1))
            except subprocess.TimeoutExpired:
                print("Timed out.")
                pass

            # Python needs a moment to catch-up it seems
            time.sleep(60)
            kwargs_main = check_galfit_out_hangups(tmp_fits_dir, out_dir, kwargs_main)
        else:
            print("Restart fitting commencing")
        
        count = 2
        while kwargs_main["galaxy_names"] and count < 10:
            print("Did not finish all galaxies, parallelizing again...\n")
            print(f"{len(kwargs_main['galaxy_names'])} galaxies to go.")
            write_to_parallel(cwd, kwargs_main, parallel_file = parallel_file, chunk_size = chunk_size)
            
            try:
                print("Piping to parallel")
                sp(f"{parallel_run_cmd}", capture_output = capture_output, timeout = 60*(timeout + 1))
            except subprocess.TimeoutExpired:
                pass
            
            time.sleep(60)
            kwargs_main = check_galfit_out_hangups(tmp_fits_dir, out_dir, kwargs_main)
            count += 1
            
        if count == 10 or kwargs_main["galaxy_names"]:
            print("Did not finish all galaxies. There is likely an error. Check parallel output files if possible. Restart with -r option.")
            print("Quitting.")
            sys.exit()
            
#             run_galfit.main(cwd          = cwd,
#                             in_dir       = in_dir,
#                             tmp_dir      = tmp_dir,
#                             out_dir      = out_dir,
#                             num_steps    = num_steps,
#                             rerun        = rerun,
#                             galaxy_names = gname
#                            )
        
    else:
        failures = go_go_galfit.main(**kwargs_main, 
                                     run_galfit = run_galfit, 
                                     run_fitspng = run_fitspng, 
                                     run_python = run_python)
        
        write_failed(out_dir, failures)

    # Unused but here for a good time
    boom = """Numerical Recipes run-time error...
gaussj: Singular Matrix-1
...now exiting to system...

                 __/~*##$%@@@******~\-__
               /f=r/~_-~ _-_ --_.^-~--\=b\
             4fF / */  .o  ._-__.__/~-. \*R\
            /fF./  . /- /' /|/|  \_  * *\ *\R\
           (iC.I+ '| - *-/00  |-  \  )  ) )|RB
           (I| (  [  / -|/^^\ |   )  /_/ | *)B
           (I(. \ `` \   \m_m_|~__/ )_ .-~ F/
            \b\\=_.\_b`-+-~x-_/ .. ,._/ , F/
             ~\_\= =  =-*###%#x==-#  *=- =/
                ~\**U/~  | i i | ~~~\===~
                        | I I \\
                       / // i\ \\
                  (   [ (( I@) )))  )
                       \_\_VYVU_/
                         || * |
                        /* /I\ *~~\
                      /~-/*  / \ \ ~~M~\
            ____----=~ // /WVW\* \|\ ***===--___

   Doh!  GALFIT crashed because at least one of the model parameters
   is bad.  The most common causes are: effective radius too small/big,
   component is too far outside of fitting region (also check fitting
   region), model mag too faint, axis ratio too small, Sersic index
   too small/big, Nuker powerlaw too small/big.  If frustrated or
   problem should persist, email for help or report problem to:
                     Chien.Y.Peng@gmail.com"""
# ## Tidying Up in case anything is leftover

# In[ ]:


if __name__ == "__main__":
    print("Cleaning up...")
    _ = sp("rm galfit.* fit.log", capture_output = capture_output)
    _ = sp("rm *.png", capture_output = capture_output)
    
    # We use the negative of remove slurm because we want cleanup to be the default
    if parallel and not dont_remove_slurm:
        _ = sp(f"rm -r \"$HOME/SLURM_turds/{parallel_run_name}\"", capture_output = capture_output)
        _ = sp(f"rm {parallel_file}", capture_output = capture_output)
        #_ = sp(f"rm {parallel_copy_input}", capture_output = capture_output)


# ## Combining Residuals

# In[ ]:


if __name__ == "__main__":
    #basename = "GALFIT"
    pkl_end_str = "output_results"
    final_pkl_file = pj(out_dir, f"{basename}_{pkl_end_str}.pkl")
    print(f"Combining all the residual calculations into {final_pkl_file}")

    #all_nmr = {}
    output_df = pd.DataFrame()
    
    if parallel:
        python_parallel   = pj(_MODULE_DIR, "Utilities", "combine_via_parallel.py")
        parallel_file     = "parallel_combine_residual"
        if exists(parallel_file):
            _ = sp(f"rm -f {parallel_file}", capture_output = capture_output)
        
        if parallel == 1:
            # For CPU parallel
            parallel_run_name = ""#"GALFITTING"
            parallel_options  = joblib.cpu_count() #"-M all"
            
        elif parallel == 2:
            # For SLURM/Cluster Computing
            parallel_run_name = "COMBINE_RESULTS"
            parallel_options  = "-M all"

        finished_pkl_num = 0
        if restart:
            check_output_pkl = [int(os.path.basename(i).split("_")[0].replace(basename, "")) 
                                for i in find_files(tmp_dir, f'{basename}*_{pkl_end_str}.pkl', "f")
                                if os.path.basename(i).split("_")[0].replace(basename, "")
                                ]
            if check_output_pkl:
                finished_pkl_num = max(check_output_pkl)

        chunk_size = 20
        # Use count for restarting purposes i.e. if all pkl files have been generated 
        # but just need to be combined
        count = 0
        with open(parallel_file, "w") as sf:
            for i, chunk in enumerate(range(chunk_size, len(galaxy_names) + chunk_size, chunk_size)):
                if i < finished_pkl_num:
                    continue

                gal_to_parallel = galaxy_names[chunk - chunk_size:][:chunk_size]
                #num_str = f"{i:0>3}"
                sf.write(f"{run_python} {python_parallel} {pj(tmp_dir, basename + str(i))} {out_dir} {','.join(gal_to_parallel)}\n")
                count += 1

        if count:
            print("parallelizing to combine residuals")
            parallel_run_cmd = f"cat {parallel_file} | {pipe_to_parallel_cmd} {parallel_run_name} {parallel_options} {parallel_verbose}"
            _ = sp(parallel_run_cmd, capture_output = capture_output)

        all_output_pkl = [pj(tmp_dir, fname) 
                          for fname in find_files(tmp_dir, f'{basename}*_{pkl_end_str}.pkl', "f")
                          if fname != f"{basename}_{pkl_end_str}.pkl"
                         ]
        #_ = [all_nmr.update(pickle.load(open(file, 'rb'))) for file in all_output_pkl]
        out_df = pd.concat(
                           [pd.read_pickle(file) for file in all_output_pkl 
                            if os.path.basename(file) != f"{basename}_{pkl_end_str}.pkl"
                           ]
                          ) 
        
    else:
        # for gname in galaxy_names:
        #     output_file = pj(out_dir, gname, f"{gname}_galfit_out.fits")
        #     if exists(output_file):
        #         with fits.open(output_file) as hdul: 
        #             all_nmr[gname] = (hdul[2].header.get("NMR", None), 
        #                               hdul[2].header.get("ks_p", None),
        #                               hdul[2].header.get("ks_stat", None)
        #                              )
        #basename         = args[0]
        #out_dir          = args[1] #os.path.dirname(basename)
        #galaxy_names     = args[2].split(",")
        
        # In this case it's not parallel but I'm just saving some hassle here
        out_df = combine_via_parallel.main("", out_dir, ",".join(galaxy_names))

    # Could split this into the above if/else but this keeps everything output
    # related in one place
    print(f"Outputting results to {final_pkl_file}")
    
    # For when I do it in one of the other scripts
    if "gname" in out_df.columns:
        out_df.set_index("gname", inplace = True)
        
    #pickle.dump(all_nmr, open(final_pkl_file, 'wb'))
    out_df.to_pickle(final_pkl_file)
    
    if not dont_remove_slurm and parallel:
        _ = sp(f"rm -r \"$HOME/SLURM_turds/{parallel_run_name}\"", capture_output = capture_output)
        _ = sp(f"rm -f {pj(tmp_dir, basename)}*_{pkl_end_str}.pkl", capture_output = capture_output)
        _ = sp(f"rm -f {parallel_file}", capture_output = capture_output)
        
    #_ = sp(f"mv {pickle_filename_temp} {pkl_file}", capture_output = capture_output)


# In[ ]:


# Deprecated
# Combine the residuals now that they're in the headers
    # print("Calculating residuals")
    # parallel_residual_calc.main(run_dir           = cwd,
    #                             # Don't actually need in_dir
    #                             #in_dir            = in_dir,
    #                             tmp_dir           = tmp_dir,
    #                             out_dir           = out_dir,
    #                             basename          = "GALFIT",
    #                             parallel             = parallel,
    #                             dont_remove_parallel = dont_remove_parallel,
    #                             restart           = restart,
    #                             verbose           = verbose,
    #                             capture_output    = not verbose
    #                            )
    # if aggressive_clean:
    #     # May need to use find and delete
    #     print("Aggressively cleaning, removing star masks.")
    #     _ = sp(f"rm -rf {pj(tmp_masks_dir)} {pj(tmp_png_dir)}", capture_output = capture_output)


# In[ ]:


if __name__ == "__main__":
    if aggressive_clean:
        print("Final tidying...")
        _ = sp(f"rm -rf {out_dir}", capture_output = capture_output)
        _ = sp(f"mkdir -p {out_dir}", capture_output = capture_output)
        
    print("All done!")
    # Moving back to original directory
    os.chdir(old_cwd)


# In[43]:


if __name__ == "__main__":
    # in_notebook() is checked in the export function
    export_to_py("control_script_debug", pj(_MODULE_DIR, "control_script"))

