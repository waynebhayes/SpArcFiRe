#!/usr/bin/env python
# coding: utf-8

# **Author: Matthew Portman**
# 
# **Date (Github date will likely be more accurate): 4/17/23**
# 1/26/21
# New usage ---
# Will default to running on top of sparc directory 
# unless otherwise specified
# 
# Command line args
#
# -r, --rerun-galfit
#		Re-run GALFIT using the galfit.01 file
#
# -p, --path
#		Specify path to Sparcfire's in/tmp/out directories in that order;
#		e.g. '-p /home/sparcfire_in /home/sparcfire_tmp /home/sparcfire_out'
#

# Change log - 1/19/21 
# Updated silent to catch everything
# Updated call to feedme gen to include paths
# Changed 'which' to 'type -P' which seems better
# See note for files= and masks=

# For controlling galfitting via sparcfire

# ***************************************************************
# In[23]:


from galfit_objects import *
# This should give me numpy and pandas and whatnot
# also gives this, from os.path import join as pj
from sparc_to_galfit_feedme_gen import *
import go_go_galfit

import argparse
import shutil
from IPython import get_ipython
import subprocess

# For debugging purposes
def in_notebook():
    ip = get_ipython()
    
    if ip:
        return True
    else:
        return False
# In[ ]:


def sp(cmd_str, capture_output = True):
    # Because it is a pain in the butt to call subprocess with all those commands every time
    return subprocess.run(cmd_str, capture_output = capture_output, text = True, shell = True, executable="/bin/bash")


# In[ ]:


if __name__ == "__main__":
    
    # Force >python 3.6 for various compatabilities
    out_str = "\t Python3.6 or greater required! Exitting without generating feedmes..."
    assert sys.version_info >= (3, 6), out_str
    
    user_home = os.environ["HOME"]
    cwd = os.getcwd() # Doesn't work *in* notebook
    username = os.environ["USER"]
    
    USAGE = f"""USAGE:

    python3 ./{sys.argv[0]} [-S|--slurm] [-R|--rerun-galfit] IN-DIRECTORY TMP-DIRECTORY OUT-DIRECTORY

    This script is the wrapping script for running GALFIT using SpArcFiRe to inform 
    the input. By default, it runs from the directory it is called and checks for the
    existence (verbatim) of the 'sparcfire-in' 'sparcfire-tmp' and 'sparcfire-out' 
    directories in the current directory. If this is not the case, it will fail and exit. 

    You may specify a different in, tmp, and out directory using the '-p' option. Please
    do not specify a symlink, it discomforts the programmer.
    """
    
    parser = argparse.ArgumentParser(description = USAGE)
    
    # TODO: invert this to running without slurm for when it comes time for the big runs
    parser.add_argument('-S', '--slurm',
                        dest     = 'slurm',
                        action   = 'store_const',
                        const    = True,
                        # Soon to be the other way around
                        default  = False,
                        help     = 'Run GALFITs using Slurm.')

    parser.add_argument('-drS', '--dont-remove-slurm',
                        dest     = 'dont_remove_slurm',
                        action   = 'store_const',
                        const    = True,
                        default  = False,
                        help     = 'Choose NOT to remove all old slurm files (they may contain basic info about each fit but there will be a bunch!)')
    
    
    parser.add_argument('-ns', '--num-steps',
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
    
    parser.add_argument('-R', '--rerun-galfit',
                        dest     = 'rerun', 
                        action   = 'store_const',
                        const    = True,
                        # Python cannot convert 'False' to a boolean... it's true *cries*
                        default  = "",
                        help     = 'Run GALFIT again after the final fit to hopefully refine said fit.')
    
    parser.add_argument(dest     = 'paths',
                        nargs    = "*",
                        type     = str,
                        help     = "IN-DIRECTORY TMP-DIRECTORY OUT-DIRECTORY from SpArcFiRe. \
                                    Must follow -in, -tmp, out or this won't work.")
    
    if not in_notebook():
        args = parser.parse_args() # Using vars(args) will call produce the args as a dict
        num_steps = args.steps
        slurm = args.slurm
        dont_remove_slurm = args.dont_remove_slurm
        rerun = args.rerun
        
        # if num_steps not in range(1,4):
        #     print("The number of steps you selected cannot be used!")
        #     print("Using two.")
        #     print()
        #     num_steps = 2

        if len(args.paths) == 3:
            in_dir, tmp_dir, out_dir = args.paths[0], args.paths[1], args.paths[2]
        else:
            in_dir = pj(cwd, "sparcfire-in")
            tmp_dir = pj(cwd, "sparcfire-tmp")
            out_dir = pj(cwd, "sparcfire-out")

            print(f"Paths incorrectly specified, defaulting to (in, tmp, out)...")
            print(f"{in_dir}\n{tmp_dir}\n{out_dir}")
            print()
            
    else:
        slurm = False
        rerun = ""
        num_steps = 2
        
        cwd = cwd.replace("ics-home", username)
        in_dir = pj(cwd, "sparcfire-in")
        tmp_dir = pj(cwd, "sparcfire-tmp")
        out_dir = pj(cwd, "sparcfire-out")
        
        sys.path.append(pj(user_home, ".local", "bin"))


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
        run_galfit  = pj(user_home, ".local/bin/galfit")
        run_fitspng = pj(user_home, ".local/bin/fitspng")
        run_python  = "/opt/conda/bin/python3"


# ## Setting up directories and handy variables

# In[ ]:


if __name__ == "__main__":
    # Setting up paths and variables
    tmp_fits_dir = pj(tmp_dir, "galfits")
    tmp_masks_dir   = pj(tmp_dir, "galfit_masks")
    tmp_psf_dir     = pj(tmp_dir, "psf_files")
    tmp_png_dir     = pj(tmp_dir, "galfit_png")
    
    #all_galfit_out = pj(out_dir, "all_galfit_out")
    out_png_dir     = pj(out_dir, "galfit_png")
    
    # Remove old
    try:
        shutil.rmtree(tmp_fits_dir)
    except OSError as e:
        pass
    
    # Making sub-directories
    _ = [os.mkdir(i) for i in (tmp_fits_dir, tmp_masks_dir, tmp_psf_dir, tmp_png_dir, out_png_dir) if not exists(i)]
    
    # makedirs will make both at once, handy!
    #if not exists(out_png): os.makedirs(out_png)


# ## Running sextractor (if necessary)

# In[ ]:


if __name__ == "__main__":    
    # Grabbing list of file names and masks with bash variable expansion
    input_filenames = glob.glob(pj(in_dir, "*.fits"))
    star_masks      = glob.glob(pj(tmp_masks_dir, "*_star-rm.fits"))
    
    if star_masks:
        sp(f"mv {pj(tmp_dir,'*_star-rm.fits')} {tmp_masks_dir}")
        
    star_masks      = glob.glob(pj(tmp_masks_dir, "*_star-rm.fits"))        
    
    try:
        if exists(pj(cwd, "star_removal")):
            star_removal_path = pj(cwd, "star_removal")
        else:
            star_removal_path = pj(os.environ["SPARCFIRE_HOME"], "star_removal")
    except KeyError as k:
        print("SPARCFIRE_HOME environment variable is not set.")
        print("Will not be able to create star masks using 'remove_stars_with_sextractor.py'")
        print("If they have already been generated, ignore this message.")
        print()
        
    else:
        os.chdir(star_removal_path)
        if len(input_filenames) != len(star_masks):
            print("Generating starmasks...")
            out_text = sp(f"python3 {pj(star_removal_path, 'remove_stars_with_sextractor.py')} {in_dir} {tmp_masks_dir}")

            if out_text.stderr.strip():
                print(f"Something went wrong running 'remove_stars_with_sextractor.py'! Printing debug info...")
                print(out_text)
                #print(type(out_text.stderr))
            else:
                print("Star masks have been generated successfully.")
                print()
        else:
            print("Star masks have already been generated, proceeding.")
            print()
            
        os.chdir(cwd)


# ## Galfitting/Slurming!

# In[ ]:


def write_to_slurm(cwd, kwargs_main, galfit_script_name = "go_go_galfit.py", slurm_file = "slurm_cmd_file"):
    _, _, run_python = go_go_galfit.check_programs()
    kwargs_in    = deepcopy(kwargs_main)
    
    print(f"Generating distrib-slurm input file in {cwd}: {slurm_file}")
    with open(pj(cwd, slurm_file), "w") as scf:
        for gname in kwargs_main["galaxy_names"]:
            kwargs_in["galaxy_names"] = gname

            cmd_str = ""
            for k,v in kwargs_in.items():
                cmd_str += f"{k}={v} "
            # Good thing dictionaries retain order now, *whistles innocently*
            scf.write(f"{run_python} {galfit_script_name} {cmd_str}\n")

    sp(f"chmod a+x {slurm_file}", capture_output = False)


# In[41]:


def check_galfit_out_hangups(galaxy_names, tmp_fits_dir, out_dir, kwargs_main):
    # For hang-ups, check if final copy has occurred
    # Check for file which indicates galfit outputs nothing at all
    # This avoids conflating ones which han'vet run and ones which have nothing to show for it
    kwargs_main["galaxy_names"] = [gname for gname in galaxy_names 
                                   if not exists(f"{pj(out_dir, gname, gname)}_galfit_out.fits")
                                   and not exists(f"{pj(tmp_fits_dir, 'failed_' + gname)}_galfit_out.fits")]
    
    # because of mutability I don't need to do this but
    # because it's good practice...
    return kwargs_main


# In[ ]:


if __name__ == "__main__":
    galaxy_names = [os.path.basename(i).rstrip(".fits") 
                for i in input_filenames]
    
    # Backing up old fits
    # So we can check for new ones
    # and also backup previous runs if there's an accidental overwrite
    for gname in galaxy_names:
        new_out = pj(out_dir, gname, f"{gname}_galfit_out.fits")
        old_out = pj(out_dir, gname, f"{gname}_galfit_out_old.fits")
        if exists(pj(out_dir, gname, new_out)):
            shutil.move(new_out, old_out)

    kwargs_main = {"cwd"          : cwd,
                   "in_dir"       : in_dir,
                   "tmp_dir"      : tmp_dir,
                   "out_dir"      : out_dir,
                   "num_steps"    : num_steps,
                   "rerun"        : rerun,
                   "galaxy_names" : galaxy_names
                  }
    
#    raise(AssertionError())
    # One at a time
    if slurm:
        slurm_file = "slurm_cmd_file"
        write_to_slurm(cwd, kwargs_main, slurm_file = slurm_file)
#         print(f"Generating distrib-slurm input file in {cwd}: {slurm_file}")
#         with open(pj(cwd, slurm_file), "w") as scf:
#             for gname in galaxy_names:
#                 kwargs_main["galaxy_names"] = gname
                
#                 cmd_str = ""
#                 for k,v in kwargs_main.items():
#                     cmd_str += f"{k}={v} "
#                 # Good thing dictionaries retain order now, *whistles innocently*
#                 scf.write(f"{run_python} go_go_galfit.py {cmd_str}\n")
        
#         sp(f"chmod a+x {slurm_file}", capture_output = False)
        
        print("Running with slurm")
        slurm_run_name = "GALFITTING"
        timeout = 5 # Minutes
        slurm_run_cmd = f"cat {slurm_file} | ~wayne/bin/distrib_slurm {slurm_run_name} -M all"
        sp(f"{slurm_run_cmd} -t {timeout}", capture_output = False)
        
        kwargs_main = check_galfit_out_hangups(galaxy_names, tmp_fits_dir, out_dir, kwargs_main)
        
        count = 2
        while kwargs_main["galaxy_names"] and count < 10:
            print("Did not finish all galaxies, slurming again, increasing timeout...")
            write_to_slurm(cwd, kwargs_main, slurm_file = slurm_file)
            
            timeout *= count
            sp(f"{slurm_run_cmd} -t {timeout}", capture_output = False)    
            kwargs_main = check_galfit_out_hangups(galaxy_names, tmp_fits_dir, out_dir, kwargs_main)
            count += 1
            
        #os.stat("file").st_size == 0
        #sp(f"")
            
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
        # TODO: Use failures for something later? Maybe write to a file?

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
    print("\nDone! Cleaning up...")
    _ = sp("rm galfit.* fit.log") # , capture_output = False)
    _ = sp("rm *.png") # , capture_output = False)
    
    # We use the negative of remove slurm because we want cleanup to be the default
    if slurm and not dont_remove_slurm:
        _ = sp(f"rm -r \"$HOME/SLURM_turds/{slurm_run_name}\"", capture_output = False)


# In[ ]:


if __name__ == "__main__":
    if in_notebook():
        export_to_py("control_script_debug", output_filename = "control_script.py")

