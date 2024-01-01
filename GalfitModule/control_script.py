#!/usr/bin/env python
# coding: utf-8

# **Author: Matthew Portman**
# 
# **Date (Github date will likely be more accurate): 11/29/23**


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


# For debugging purposes
from IPython import get_ipython
def in_notebook():
    ip = get_ipython()
    
    if ip:
        return True
    else:
        return False

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

# This should give me numpy and pandas and whatnot
# also gives this, from os.path import join as pj
from sparc_to_galfit_feedme_gen import *
import go_go_galfit

#from Classes.Components import *
from Classes.Containers import *
from Functions.helper_functions import *
import Utilities.parallel_residual_calc as parallel_residual_calc
import Utilities.combine_via_parallel as combine_via_parallel

# Separations were previously code blocks in jupyter notebook
# ==========================================================================================================
# Parse input
# ==========================================================================================================
if __name__ == "__main__":
    
    # out_str = "\t Python3.7 or greater required! Exitting without generating feedmes..."
    # assert sys.version_info >= (3, 7), out_str
    
    cwd = absp(os.getcwd())
    old_cwd = absp(cwd)
    
    username = os.environ["USER"]
    
    USAGE = f"""USAGE:

    python3 ./{sys.argv[0]} [OPTION] [[RUN-DIRECTORY] IN-DIRECTORY TMP-DIRECTORY OUT-DIRECTORY]
    
    OPTIONS =>[-p | --parallel]
              [-drs | --dont-remove-slurm]
              [-t  | --tmp]
              [-ac | --aggressive-clean]
              [-NS | --num-steps] 
              [-r | --restart]
              [-nsf | --no-simultaneous-fitting]
              [-v | --verbose]
              [-n | --basename]

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
                        #help     = 'Aggressively clean-up directories, removing -in, temp output, psf, and mask files after galfit runs'
                        help     = 'Aggressively clean-up directories, removing temp output and mask files after galfit runs'
                       )
    
    parser.add_argument('-NS', '--num-steps',
                        dest     = 'steps', 
                        action   = 'store',
                        type     = int,
                        choices  = range(1,4),
                        default  = 2,
                        help     = 'Run GALFIT using step-by-step component selection (up to 3), i.e.\n\t\
                                    1: Bulge + Disk + Arms,\n\t\
                                    2: Disk -> Bulge + Disk + Arms,\n\t\
                                    3: Disk -> Bulge + Disk -> Bulge + Disk + Arms'
                       )
    
    parser.add_argument('-r', '--restart',
                        dest     = 'restart',
                        action   = 'store_const',
                        const    = True,
                        default  = False,
                        help     = 'Restart control script on the premise that some have already run (likely in parallel).'
                       )
    
    parser.add_argument('-nsf', '--no-simultaneous-fitting',
                        dest     = 'simultaneous_fitting',
                        action   = 'store_const',
                        const    = False,
                        default  = True,
                        help     = 'Turn off simultaneous fitting.'
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
                        help     = 'Verbose output. Includes stdout.'
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
        
        restart           = args.restart
        basename          = args.basename
        
        # TODO: This isn't working, forcing false for now
        simultaneous_fitting = False #args.simultaneous_fitting
        
        verbose           = args.verbose
        capture_output    = not args.verbose

        print(f"Number of fitting steps: {num_steps}")
        print(f"Parallel choice (0 - serial, 1 - CPU, 2 - SLURM): {parallel}")
        print(f"Run from temporary directory? {run_from_tmp}")
        print(f"Aggressively clean? {aggressive_clean}")
        print(f"Verbosity? {verbose}")
        if parallel == 2:
            print(f"Don't remove SLURM dump files? {dont_remove_slurm}")
        
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
        #rerun = ""
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

# ==========================================================================================================
# Checking file directories and installed programs
# ==========================================================================================================

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

    else:
        run_galfit  = pj(_HOME_DIR, ".local/bin/galfit")
        run_fitspng = pj(_HOME_DIR, ".local/bin/fitspng")
        run_python  = "/opt/conda/bin/python3"

# ==========================================================================================================
# Setting up directories and handy variables
# ==========================================================================================================

if __name__ == "__main__":
    # Setting up paths and variables
    tmp_fits_dir    = pj(tmp_dir, "galfits")
    tmp_masks_dir   = pj(tmp_dir, "galfit_masks")
    # Now dropping these in the individual galaxy folders
    #tmp_psf_dir     = pj(tmp_dir, "psf_files")
    tmp_png_dir     = pj(tmp_dir, "galfit_png")
    
    sim_fitting_dir = pj(tmp_dir, "sim_fitting")
    if simultaneous_fitting:
        sf_in_dir       = pj(sim_fitting_dir, "sparcfire-in")
        sf_tmp_dir      = pj(sim_fitting_dir, "sparcfire-tmp")
        sf_out_dir      = pj(sim_fitting_dir, "sparcfire-out")
        sf_masks_dir    = pj(sf_tmp_dir, "galfit_masks")
        
        # These really the necessary ones, tmp and masks not so much but kept for posterity
        for path in [sim_fitting_dir, sf_in_dir, sf_out_dir]:
            if not exists(path):
                simultaenous_fitting = False
                print(f"Simultaneous Fitting cannot be performed, {path} does not exist!")
                print("Continuing...")
                break
        
    out_png_dir     = pj(out_dir, "galfit_png")
    
    # Minor changes to SpArcFiRe's version
    star_removal_path = pj(_MODULE_DIR, "star_removal")
    
    if parallel == 1:
        # CPU Parallel
        pipe_to_parallel_cmd = pj(_MODULE_DIR, "ParallelDrivers", "parallel")
        
    elif parallel == 2:
        # SLURM/Cluster Computing
        pipe_to_parallel_cmd = pj(_MODULE_DIR, "ParallelDrivers", "distrib_slurm")
    
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
                              )
         if not exists(i)
        ]
    
    if simultaneous_fitting:
        _ = [os.mkdir(i) for i in (sf_tmp_dir, 
                                   sf_masks_dir
                                  )
         if not exists(i)
        ]
    
# ==========================================================================================================
# Finding input FITS files (observations), shuffling things around from SpArcFiRe
# ==========================================================================================================

if __name__ == "__main__":    
    # Grabbing list of file names and masks with bash variable expansion
    print("Setting up, finding files...")

    input_filenames = find_files(in_dir, "*.fits", "f")
    if not input_filenames:
        print(f"No input files found in {in_dir}. Is something wrong?")
        sys.exit()
        
    output_folders  = [i for i in find_files(out_dir, "*", "d") 
                       if os.path.basename(i) != os.path.basename(out_dir) and
                          os.path.basename(i) != "galfit_png"
                      ]
    
    star_masks      = find_files(tmp_dir, "*_star-rm.fits", "f")
       
    # Sparcfire plops them into tmp so I just move them to where I need them
    if star_masks and not restart:
        sp(f"mv {pj(tmp_dir,'*_star-rm.fits')} {tmp_masks_dir}", capture_output = capture_output)

# ==========================================================================================================
# Getting setup to generate starmasks
# ==========================================================================================================

if __name__ == "__main__":
    generate_starmasks = True
    if not restart:
        #generate_starmasks = True
        if not exists(pj(tmp_masks_dir, 'remove_stars_with_sextractor.log')):
            # remove_stars_with_sextractor needs this available before it can log
            _ = sp(f"touch {pj(tmp_masks_dir, 'remove_stars_with_sextractor.log')}", capture_output = capture_output)
            
        if simultaneous_fitting and not exists(pj(sf_masks_dir, 'remove_stars_with_sextractor.log')):
            # remove_stars_with_sextractor needs this available before it can log
            _ = sp(f"touch {pj(sf_masks_dir, 'remove_stars_with_sextractor.log')}", capture_output = capture_output)

# ==========================================================================================================
# write_to_parallel function sets up a text file for piping to parallel scripts
# ==========================================================================================================

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


# ==========================================================================================================
# check_galfit_out_hangups function checks for output being successful/insuccessful
# ==========================================================================================================


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

# ==========================================================================================================
# write_failed function writes a list of failures to file :(
# ==========================================================================================================

def write_failed(failed_dir = cwd, failures = []):
    if failures:
        fail_filename = "galfit_failed.txt"
        fail_filepath = pj(failed_dir, fail_filename)
        print(f"{len(failures)} galax(y)ies completely failed. Writing the list of these to {fail_filepath}")
        with open(fail_filepath, "w") as ff:
            ff.writelines("\n".join(failures))
            ff.write("\n")

# ==========================================================================================================
# GALFITting! Here it all comes together.
# ==========================================================================================================

if __name__ == "__main__":    
    
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

    kwargs_main = {"cwd"                  : cwd,
                   "in_dir"               : in_dir,
                   "tmp_dir"              : tmp_dir,
                   "out_dir"              : out_dir,
                   "num_steps"            : num_steps,
                   #"rerun"                : rerun,
                   "parallel"             : parallel,
                   "verbose"              : verbose,
                   "capture_output"       : capture_output,
                   "generate_starmasks"   : generate_starmasks,
                   "run_from_tmp"         : run_from_tmp,
                   "aggressive_clean"     : aggressive_clean,
                   "simultaneous_fitting" : simultaneous_fitting,
                   "sim_fitting_dir"      : sim_fitting_dir,
                   # "petromags"          : petromags,
                   # "bulge_axis_ratios"  : bulge_axis_ratios,
                   # Keep this last just in case
                   "galaxy_names"       : galaxy_names
                  }
    
    # In case we're running back to back, this will reduce galaxy_names appropriately
    kwargs_main = check_galfit_out_hangups(tmp_fits_dir, out_dir, kwargs_main)
    parallel_file = "parallel_cmd_file"
    
    # Getting read for parallelizing
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
            timeout = 480 # Minutes
            # TODO: Consider SLURM + CPU parallel
            # --ntasks-per-node=1 and --ntasks=1 ensures processes will stay
            # on the same node which is crucial for asyncio
            parallel_options  = f"-M all --ntasks=1 --ntasks-per-node=1 -t {timeout}"
            parallel_verbose  = "-v" if verbose else ""
            chunk_size = 20

            
        # Running things via distributed computing           
        parallel_run_cmd = f"cat {parallel_file} | nice -19 {pipe_to_parallel_cmd} {parallel_run_name} {parallel_options} {parallel_verbose}"
        
        if not restart:
            write_to_parallel(cwd, kwargs_main, parallel_file = parallel_file, chunk_size = chunk_size)
            print("Galfitting via parallelization...")
            try:
                sp(f"{parallel_run_cmd}", capture_output = capture_output, timeout = (timeout + 1) * 60)
            except subprocess.TimeoutExpired:
                print("Timed out.")
                pass

            # Python needs a moment to catch-up it seems
            time.sleep(60)
            kwargs_main = check_galfit_out_hangups(tmp_fits_dir, out_dir, kwargs_main)
        else:
            print("Restart fitting commencing")
        
        # Just in case something happened on one of the nodes i.e. too slow/frozen/etc.
        count = 2
        while kwargs_main["galaxy_names"] and count < 10:
            if not restart:
                print("Did not finish all galaxies, parallelizing again...\n")
                print(f"{len(kwargs_main['galaxy_names'])} galaxies to go.")

            write_to_parallel(cwd, kwargs_main, parallel_file = parallel_file, chunk_size = chunk_size)
            
            try:
                print("Piping to parallel")
                sp(f"{parallel_run_cmd}", capture_output = capture_output, timeout = (timeout + 1)* 60)
            except subprocess.TimeoutExpired:
                pass
            
            time.sleep(60)
            kwargs_main = check_galfit_out_hangups(tmp_fits_dir, out_dir, kwargs_main)
            count += 1
            
        if count == 10 or kwargs_main["galaxy_names"]:
            print("Did not finish all galaxies. There is likely an error. Check parallel output files if possible. Restart with -r option.")
            print("Quitting.")
            sys.exit()
        
    # Running in serial
    else:
        failures = go_go_galfit.main(**kwargs_main, 
                                     run_galfit  = run_galfit, 
                                     run_fitspng = run_fitspng, 
                                     run_python  = run_python
                                    )
        
        write_failed(out_dir, failures)

# ==========================================================================================================
# Unused but here for a good time
# ==========================================================================================================
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
    
# ==========================================================================================================
# Tidying Up in case anything is leftover
# ==========================================================================================================

if __name__ == "__main__":
    print("Cleaning up...")
    _ = sp("rm galfit.* fit.log", capture_output = capture_output)
    _ = sp("rm *.png", capture_output = capture_output)
    
    # We use the negative of remove slurm because we want cleanup to be the default
    if parallel and not dont_remove_slurm:
        _ = sp(f"rm -r \"$HOME/SLURM_turds/{parallel_run_name}\"", capture_output = capture_output)
        _ = sp(f"rm {parallel_file}", capture_output = capture_output)
        #_ = sp(f"rm {parallel_copy_input}", capture_output = capture_output)


# ==========================================================================================================
# Combining Residuals
# ==========================================================================================================
if __name__ == "__main__":
    pkl_end_str = "output_results"
    final_pkl_file = pj(out_dir, f"{basename}_{pkl_end_str}.pkl")
    num = 1
    
    # Iterative naming for subsequent runs
    while exists(final_pkl_file) and num < 100:
        final_pkl_file = pj(out_dir, f"{basename}_{pkl_end_str}{num}.pkl")
        num += 1
        
    if num == 100:
        print("Stopping the iterative naming just in case...")
        print("Use the -n option to specify a new run name. Exitting.")
        sys.exit()
        
    print(f"Combining all the residual calculations into {final_pkl_file}")

    output_df = pd.DataFrame()
    
    # Piping to parallel again because we're already setup to do so
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
        # Performing a similar check to that above for failures/success
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
                sf.write(f"{run_python} {python_parallel} {pj(tmp_dir, basename + str(i))} {pkl_end_str} {out_dir} {','.join(gal_to_parallel)}\n")
                count += 1

        if count:
            print("parallelizing to combine residuals")
            parallel_run_cmd = f"cat {parallel_file} | nice -19 {pipe_to_parallel_cmd} {parallel_run_name} {parallel_options} {parallel_verbose}"
            _ = sp(parallel_run_cmd, capture_output = capture_output)

        all_output_pkl = [pj(tmp_dir, fname) 
                          for fname in find_files(tmp_dir, f'{basename}*_{pkl_end_str}.pkl', "f")
                          if fname != f"{basename}_{pkl_end_str}.pkl"
                         ]
        #_ = [all_nmr.update(pickle.load(open(file, 'rb'))) for file in all_output_pkl]
        # Joining everything together
        out_df = pd.concat(
                           [pd.read_pickle(file) for file in all_output_pkl 
                            if os.path.basename(file) != f"{basename}_{pkl_end_str}.pkl"
                           ]
                          ) 
    # Serial
    else:
        
        # In this case it's not parallel but I'm just saving some hassle here
        out_df = combine_via_parallel.main("", pkl_end_str, out_dir, ",".join(galaxy_names))

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
        
        
# ==========================================================================================================
# Aggressively tidying up if necessary and moving back to original directory
# ==========================================================================================================
if __name__ == "__main__":
    if aggressive_clean and run_from_tmp:
        print("Final tidying...")
        _ = sp(f"rm -rf {out_dir}", capture_output = capture_output)
        _ = sp(f"mkdir -p {out_dir}", capture_output = capture_output)
        
    print("All done!")
    # Moving back to original directory
    os.chdir(old_cwd)
