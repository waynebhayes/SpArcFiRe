#!/usr/bin/env python
# coding: utf-8

# **Author: Matthew Portman**
# 
# **Date (Github date will likely be more accurate): 4/17/23**

# In[3]:


import sys
import os
from os.path import join as pj

import argparse
import shutil
import subprocess

import time


# In[4]:


# For debugging purposes
from IPython import get_ipython
def in_notebook():
    ip = get_ipython()
    
    if ip:
        return True
    else:
        return False


# In[13]:


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
    
    OPTIONS =>[-s | --serial] 
              [-drs | --dont-remove-slurm] 
              [-NS | --num-steps] 
              [-r | --restart]
              [-RrG | --rerun-galfit]
              [-v | --verbose]

    This script is the wrapping script for running GALFIT using SpArcFiRe to inform 
    the input. By default, it runs from the RUN (or current) directory and uses the
    '-in' '-tmp' and '-out' directories as specified or otherwise defaults to 
    'sparcfire-in', 'sparcfire-tmp', 'sparcfire-out'. 

    Please do not specify symlinks for the above, they discomfort the programmer.
    """
    
    parser = argparse.ArgumentParser(description = USAGE)
    
    parser.add_argument('-s', '--serial',
                        dest     = 'slurm',
                        action   = 'store_const',
                        const    = False,
                        default  = True,
                        help     = 'Run GALFITs without using Slurm.'
                       )

    parser.add_argument('-drs', '--dont-remove-slurm',
                        dest     = 'dont_remove_slurm',
                        action   = 'store_const',
                        const    = True,
                        default  = False,
                        help     = 'Choose NOT to remove all old slurm files (they may contain basic info about each fit but there will be a bunch!)'
                       )
    
    parser.add_argument('-NS', '--num-steps',
                        dest     = 'steps', 
                        action   = 'store',
                        type     = int,
                        choices  = range(1,4),
                        default  = 2,
                        help     = 'Run GALFIT using step-by-step component selection (up to 3), i.e. \n \
                                    1: Bulge + Disk + Arms,\n \
                                    2: Bulge -> Bulge + Disk + Arms,\n \
                                    3: Bulge -> Bulge + Disk -> Bulge + Disk + Arms'
                       )
    
    parser.add_argument('-r', '--restart',
                        dest     = 'restart',
                        action   = 'store_const',
                        const    = True,
                        default  = False,
                        help     = 'Restart control script on the premise that some have already run (likely with SLURM).'
                       )
    
    parser.add_argument('-RrG', '--rerun-galfit',
                        dest     = 'rerun', 
                        action   = 'store_const',
                        const    = True,
                        # Python cannot convert 'False' to a boolean... it's true *cries*
                        default  = "",
                        help     = 'Run GALFIT again after the final fit to hopefully refine said fit.'
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
        slurm             = args.slurm
        dont_remove_slurm = args.dont_remove_slurm
        rerun             = args.rerun
        restart           = args.restart
        
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
        slurm = False
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
        run_galfit, run_fitspng, run_python = go_go_galfit.check_programs()
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
    tmp_psf_dir     = pj(tmp_dir, "psf_files")
    tmp_png_dir     = pj(tmp_dir, "galfit_png")
    #need_masks_dir  = pj(tmp_dir, "need_masks")
    
    #all_galfit_out = pj(out_dir, "all_galfit_out")
    out_png_dir     = pj(out_dir, "galfit_png")
    
    # Should be same as SpArcFiRe's but I'm packaging it with just in case
    star_removal_path = pj(_MODULE_DIR, "star_removal")
    
    pipe_to_slurm_cmd = "~wayne/bin/distrib_slurm"
    
    if not restart:
        # Remove old
        try:
            shutil.rmtree(tmp_fits_dir)
        except OSError as e:
            pass

    # Making sub-directories
    _ = [os.mkdir(i) for i in (tmp_fits_dir, 
                               tmp_masks_dir, 
                               tmp_psf_dir, 
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
    output_folders  = find_files(out_dir, "123*", "d")
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
    
    if not restart:
        generate_starmasks = True
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
            
# DEPRECATED
                           
#             print("Getting things ready to generate star masks...")
#             #galaxy_folder_names = [os.path.basename(i.rstrip("/")) for i in output_folders]
#             #star_mask_names = [os.path.basename(i) for i in star_masks]
            
#             slurm_copy_input = "slurm_copy_inputs"
#             if slurm:
#                 if exists(slurm_copy_input):
#                     _ = sp(f"rm {slurm_copy_input}", capture_output = capture_output)

#                 _ = sp(f"touch {slurm_copy_input}", capture_output = capture_output)

#                 sci = open(slurm_copy_input, "a")
                
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

#                         if slurm:
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

#                 if slurm:
#                     sci.close()

#                     extra_slurm = ""
#                     if verbose:
#                         extra_slurm = "-v"

#                     slurm_run_name = "COPYING_INPUT_FILES"
#                     slurm_run_cmd = f"cat {slurm_copy_input} | {pipe_to_slurm_cmd} {slurm_run_name} -M all {extra_slurm}"
#                     # Running without timeout for now
#                     print("Performing the copy with slurm...")
#                     _ = sp(f"{slurm_run_cmd}", capture_output = capture_output)
                
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
                
        else:
            generate_starmasks = False
            print("Star masks have already been generated, proceeding.")
            print()


# ## Galfitting/Slurming!

# In[ ]:


def write_to_slurm(cwd, 
                   kwargs_main, 
                   galfit_script_name = pj(_MODULE_DIR, "go_go_galfit.py"), 
                   slurm_file = pj(cwd, "slurm_cmd_file"),
                   chunk_size = 10
                  ):
    _, _, run_python = go_go_galfit.check_programs()
    kwargs_in        = deepcopy(kwargs_main)
    
    print(f"Generating distrib-slurm input file in {cwd}: {slurm_file}")
    with open(pj(cwd, slurm_file), "w") as scf:
        #for gname in kwargs_main["galaxy_names"]:
        for i, chunk in enumerate(range(chunk_size, len(kwargs_main["galaxy_names"]) +  chunk_size, chunk_size)):
            chunk_o_galaxies = kwargs_main["galaxy_names"][chunk - chunk_size:][:chunk_size]
            kwargs_in["galaxy_names"] = ",".join(chunk_o_galaxies)

            cmd_str = ""
            for k,v in kwargs_in.items():
                cmd_str += f"{k}={v} "
                
            # Good thing dictionaries retain order now, *whistles innocently*
            scf.write(f"{run_python} {galfit_script_name} {cmd_str}\n")

    sp(f"chmod a+x {slurm_file}", capture_output = capture_output)


# In[ ]:


def check_galfit_out_hangups(tmp_fits_dir, out_dir, kwargs_main):
    # For hang-ups, check if final copy has occurred
    # Check for file which indicates galfit outputs nothing at all
    # This avoids conflating ones which haven't run and ones which have nothing to show for it
    kwargs_main["galaxy_names"] = [gname for gname in kwargs_main["galaxy_names"]
                                   if not exists(f"{pj(out_dir, gname, gname)}_galfit_out.fits")
                                   and not exists(f"{pj(tmp_fits_dir, 'failed_' + gname)}_galfit_out.fits")]
    
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
    print("Finding all galaxies...")
    # galaxy_names = [os.path.basename(i).rstrip(".fits") 
    #                 for i in input_filenames]
    
    # Replace is safer than rstrip because rstrip could remove those characters
    # from the end of the filename
    galaxy_names = [i.replace(".fits", "") for i in input_filenames]
    
    if not restart:
        # Backing up old fits
        # So we can check for new ones
        # and also backup previous runs if there's an accidental overwrite
        for gname in galaxy_names:
            new_out = pj(out_dir, gname, f"{gname}_galfit_out.fits")
            old_out = pj(out_dir, gname, f"{gname}_galfit_out_old.fits")
            if exists(pj(out_dir, gname, new_out)):
                shutil.move(new_out, old_out)

    kwargs_main = {"cwd"                : cwd,
                   "in_dir"             : in_dir,
                   "tmp_dir"            : tmp_dir,
                   "out_dir"            : out_dir,
                   "num_steps"          : num_steps,
                   "rerun"              : rerun,
                   "slurm"              : slurm,
                   "verbose"            : verbose,
                   "capture_output"     : capture_output,
                   "generate_starmasks" : generate_starmasks,
                   # THIS MUST BE LAST FOR SENDING SEVERAL TO SLURM
                   "galaxy_names"       : galaxy_names
                  }
    
    # In case we're running back to back, this will reduce galaxy_names appropriately
    kwargs_main = check_galfit_out_hangups(tmp_fits_dir, out_dir, kwargs_main)
    slurm_file = "slurm_cmd_file"
    
#    raise(AssertionError())
    # One at a time
    if slurm:
        if not kwargs_main["galaxy_names"]:
            print("No galaxies to fit, exitting.")
            sys.exit()
            
        #print("Piping to slurm")
        print(f"{len(kwargs_main['galaxy_names'])} galaxies")
        slurm_run_name = "GALFITTING"
        timeout = 29 # Minutes
        
        extra_slurm = ""
        if verbose:
            extra_slurm = "-v"
            
        slurm_run_cmd = f"cat {slurm_file} | {pipe_to_slurm_cmd} {slurm_run_name} -M all {extra_slurm}"
        
        if not restart:
            write_to_slurm(cwd, kwargs_main, slurm_file = slurm_file)
            print("Galfitting via SLURM...")
            try:
                sp(f"{slurm_run_cmd} -t {timeout}", capture_output = capture_output, timeout = 60*(timeout + 1))
            except subprocess.TimeoutExpired:
                pass

        # Python needs a moment to catch-up it seems
        time.sleep(60)
        kwargs_main = check_galfit_out_hangups(tmp_fits_dir, out_dir, kwargs_main)
        
        count = 2
        while kwargs_main["galaxy_names"] and count < 10:
            print("Did not finish all galaxies, slurming again, increasing timeout...\n")
            print(f"{len(kwargs_main['galaxy_names'])} galaxies to go.")
            write_to_slurm(cwd, kwargs_main, slurm_file = slurm_file)
            
            timeout *= count
            try:
                print("Piping to slurm")
                sp(f"{slurm_run_cmd} -t {timeout}", capture_output = capture_output, timeout = 60*(timeout + 1))
            except subprocess.TimeoutExpired:
                pass
            
            kwargs_main = check_galfit_out_hangups(tmp_fits_dir, out_dir, kwargs_main)
            count += 1
            
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
        
        write_failed(tmp_dir, failures)

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
    print("Done! Cleaning up...")
    _ = sp("rm galfit.* fit.log", capture_output = capture_output)
    _ = sp("rm *.png", capture_output = capture_output)
    
    # We use the negative of remove slurm because we want cleanup to be the default
    if slurm and not dont_remove_slurm:
        _ = sp(f"rm -r \"$HOME/SLURM_turds/{slurm_run_name}\"", capture_output = capture_output)
        _ = sp(f"rm {slurm_file}", capture_output = capture_output)
        #_ = sp(f"rm {slurm_copy_input}", capture_output = capture_output)
        
    # Moving back to original directory
    os.chdir(old_cwd)


# In[14]:


if __name__ == "__main__":
    # in_notebook() is checked in the export function
    export_to_py("control_script_debug", pj(_MODULE_DIR, "control_script"))

