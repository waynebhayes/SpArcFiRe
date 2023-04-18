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
# In[ ]:


from galfit_objects import *
# This should give me numpy and pandas and whatnot
# also gives this, from os.path import join as pj
from sparc_to_galfit_feedme_gen import *
import argparse
import shutil

from os.path import exists
from IPython import get_ipython
import subprocess

# # For debugging purposes
# def in_notebook():
#     ip = get_ipython()
    
#     if ip:
#         return True
#     else:
#         return False
# # In[123]:


def sp(cmd_str, capture_output = True):
    # Because it is a pain in the butt to call subprocess with all those commands every time
    return subprocess.run(cmd_str, capture_output = capture_output, text = True, shell = True, executable="/bin/bash")


# In[153]:


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
                        type     = bool, 
                        help     = 'Run GALFITs using Slurm.')
    
    parser.add_argument('-ns', '--num-steps',
                        dest     = 'steps', 
                        action   = 'store', 
                        help     = 'Run GALFIT using step-by-step component selection (up to 3), i.e. \n \
                                    1: Bulge + Disk + Arms,\n \
                                    2: Bulge -> Bulge + Disk + Arms,\n \
                                    3: Bulge -> Bulge + Disk -> Bulge + Disk + Arms'
                       )
    
    parser.add_argument('-R', '--rerun-galfit',
                        dest     = 'rerun', 
                        type     = bool, 
                        help     = 'Run GALFIT again after a successful fit to refine said fit.')
    
    parser.add_argument(dest     = 'paths',
                        nargs    = "*",
                        type     = str,
                        help     = "IN-DIRECTORY TMP-DIRECTORY OUT-DIRECTORY from SpArcFiRe. \
                                    Must follow -in, -tmp, out or this won't work.")
    
    if not in_notebook():
        args = parser.parse_args() # Using vars(args) will call produce the args as a dict
        num_steps = int(args.steps)
        if num_steps not in range(1,4):
            print("The number of steps you selected cannot be used!")
            print("Using two.")
            num_steps = 2

        if len(args.paths) == 3:
            in_dir, tmp_dir, out_dir = args.paths[0], args.paths[1], args.paths[2]
        else:
            in_dir = pj(cwd, "sparcfire-in")
            tmp_dir = pj(cwd, "sparcfire-tmp")
            out_dir = pj(cwd, "sparcfire-out")

            print(f"Paths incorrectly specified, defaulting to (in, tmp, out)...")
            print(f"{in_dir}\n{tmp_dir}\n{out_dir}")
            
    else:
        cwd = cwd.replace("ics-home", username)
        in_dir = pj(cwd, "sparcfire-in")
        tmp_dir = pj(cwd, "sparcfire-tmp")
        out_dir = pj(cwd, "sparcfire-out")
        
        sys.path.append(pj(user_home, ".local", "bin"))


# ## Checking file directories and installed programs

# In[135]:


if __name__ == "__main__":
    # Checking dirs
    assert all((exists(in_dir), exists(tmp_dir), exists(out_dir))), \
           f"Cannot find one of: in, tmp, or out directories:\n{in_dir}\n{tmp_dir}\n{out_dir}\n" \
           "Do those look right?"
    
    # to redirect stderr to /dev/null as well:
    #subprocess.run(['ls', '-l'], stderr=subprocess.DEVNULL)
    
    if not in_notebook():
        # This seems to work in Python directly so I'm leaving it as-is
        # Checking galfit
        run_galfit = shutil.which("galfit")
        #run_galfit = response.stdout.strip()
        
        # Checking fitspng
        run_fitspng = shutil.which("fitspng")
        fitspng_param = "0.25,1" #1,150"

    else:
        run_galfit = pj(home_dir, ".local/bin/galfit")
        run_fitspng = pj(home_dir, ".local/bin/fitspng")


# ## Setting up directories and handy variables

# In[1]:


if __name__ == "__main__":
    # Setting up paths and variables
    tmp_galfits = pj(tmp_dir, "galfits")
    tmp_masks   = pj(tmp_dir, "galfit_masks")
    tmp_psf     = pj(tmp_dir, "psf_files")
    tmp_png     = pj(tmp_dir, "galfit_png")
    
    # I was going to take this out but
    #all_galfit_out = pj(out_dir, "all_galfit_out")
    out_png     = pj(out_dir, "galfit_png")
    
    # Remove old
    try:
        shutil.rmtree(tmp_galfits)
    except OSError as e:
        pass
    
    # Making sub-directories
    _ = [os.mkdir(i) for i in (tmp_galfits, tmp_masks, tmp_psf, out_png) if not exists(i)]
    
    # makedirs will make both at once, handy!
    #if not exists(out_png): os.makedirs(out_png)


# ## Running feedme gen and sextractor (if necessary)

# In[137]:


if __name__ == "__main__":
    print("Running feedme generator...")
    feedme_info = write_to_feedmes(top_dir = cwd)
    
    # Grabbing list of file names and masks with bash variable expansion
    input_filenames = glob.glob(pj(in_dir, "*.fits"))
    star_masks      = glob.glob(pj(tmp_dir, "*_star-rm.fits"))
    
    if star_masks:
        sp(f"mv {pj(tmp_dir,'*_star-rm.fits')} {tmp_masks}")
        
    star_masks      = glob.glob(pj(tmp_masks, "*_star-rm.fits"))        
    
    try:
        if exists(pj(cwd, "star_removal")):
            star_removal_path = pj(cwd, "star_removal")
        else:
            star_removal_path = pj(os.environ["SPARCFIRE_HOME"], "star_removal")
    except KeyError as k:
        print("SPARCFIRE_HOME environment variable is not set.")
        print("Will not be able to create star masks using 'remove_stars_with_sextractor.py'")
        print("If they have already been generated, ignore this message.")
        
    else:
        os.chdir(star_removal_path)
        if len(input_filenames) != len(star_masks):
            print("Generating starmasks...")
            out_text = sp(f"python3 {pj(star_removal_path, 'remove_stars_with_sextractor.py')} {in_dir} {tmp_masks}")

            if out_text.stderr.strip():
                print(f"Something went wrong running 'remove_stars_with_sextractor.py'! Printing debug info...")
                print(out_text)
                #print(type(out_text.stderr))
            else:
                print("Star masks have been generated successfully.")
        else:
            print("Star masks have already been generated, proceeding.")
            
        os.chdir(cwd)

# ## Galfitting/preparing for slurm!

# In[161]:


def check_update(log_out, galfit_num, components, rerun = False):
    
    # Put this here (for now) so it's easier to examine other components
    # of CompletedProcess from subprocess should I desire
    # i.e. I don't have to find it somewhere below
    log_out = log_out.stdout
    if exists(f"galfit.{galfit_num}"):
        final_model = log_out.split("Iteration")[-1]

        # bulge, disk, arms, fourier, sky
        components = update_components(final_model, 
                                       *components)
        
        galfit_num = f"{int(galfit_num) + 1:0>2}"
        
        if rerun:
            run_galfit_cmd = f"{run_galfit} -imax {max_it} galfit.{galfit_num}"
            galfit_num = f"{int(galfit_num) + 1:0>2}"
            out_text = sp(run_galfit_cmd)
            # Recursion... but forcing false to avoid rerunning over and over
            components, galfit_num = check_update(out_text, galfit_num, rerun = False)
            
    else:
        print(f"Galfit failed this run!")
        #print(log_out)
        #print("Expect output .# files to be less one albeit not necessarily missing this one.")
        
    
    return components, galfit_num


# In[155]:


if __name__ == "__main__":
    galaxy_names = [os.path.basename(i).rstrip(".fits") 
                    for i in input_filenames]
    #out_dirs = [os.path.abspath(i).replace("-in", "-out").rstrip(".fits") 
    #            for i in input_filenames]
    
    max_it = 150
    failed = []
    
    # Working on non-slurm for now
    print()
    print(f"Galfitting with {num_steps} component steps... (stdout is captured):")
    for gname in galaxy_names:
        # For debugging
        # if gname != "1237667783896924233":
        #     continue
        print(gname)
        
        # This doesn't change
        header  = feedme_info[gname]["header"]
        
        # These get *really* updated
        bulge   = feedme_info[gname]["bulge"]
        disk    = feedme_info[gname]["disk"]
        arms    = feedme_info[gname]["arms"]
        fourier = feedme_info[gname]["fourier"]
        sky     = feedme_info[gname]["sky"]
        components = (bulge, disk, arms, fourier, sky)
        
        if num_steps == 1:
            run_galfit_cmd = f"{run_galfit} -imax {max_it} {feedme_info[gname]['path']}"
            _ = sp(run_galfit_cmd)
            
        # Todo: how to get this to play nice with slurm?
        elif num_steps >= 2:
            bulge_in = pj(out_dir, gname, f"{gname}_bulge.in")
            header.to_file(bulge_in, bulge, sky)
            
            run_galfit_cmd = f"{run_galfit} -imax {max_it} {bulge_in}"
            print("Bulge")
            out_text = sp(run_galfit_cmd)
            # For now(?) force rerun false
            galfit_num = "01"
            components, galfit_num = check_update(out_text, galfit_num, components, rerun = False)
            
            if num_steps == 3:
                disk_in = pj(out_dir, gname, f"{gname}_disk.in")
                header.to_file(disk_in, bulge, disk, sky)
                
                run_galfit_cmd = f"{run_galfit} -imax {max_it} {disk_in}"
                print("Bulge + Disk")
                out_text = sp(run_galfit_cmd)
                
                components, galfit_num = check_update(out_text, galfit_num, components, rerun = False)
                
            # Overwrite original... for now 
            header.to_file(f"{feedme_info[gname]['path']}", *components)
            run_galfit_cmd = f"{run_galfit} -imax {max_it} {feedme_info[gname]['path']}"
            print("Bulge + Disk + Arms")
            _ = sp(run_galfit_cmd)
        
        tmp_png_path  = pj(tmp_png, gname)
        tmp_fits_path = pj(tmp_galfits, gname)
        
        if not exists(f"{tmp_fits_path}_out.fits"):
            print(f"{gname} completely failed!!!")
            failed.append(gname)
            continue
        
        fitspng_cmd1 = f"{run_fitspng} -fr \"{fitspng_param}\" -o \
                        {tmp_png_path}.png {tmp_fits_path}.fits[1]"
        fitspng_cmd2 = f"{run_fitspng} -fr \"{fitspng_param}\" -o \
                        {tmp_png_path}_out.png {tmp_fits_path}.fits[2]"
        fitspng_cmd3 = f"{run_fitspng} -fr \"{fitspng_param}\" -o \
                        {tmp_png_path}_residual.png {tmp_fits_path}.fits[3]"
        
        sp(fitspng_cmd1)
        sp(fitspng_cmd2)
        sp(fitspng_cmd3)
        
        # Combining the images using ImageMagick
        montage_cmd = f"montage {tmp_png_path}.png \
                                {tmp_png_path}_out.png \
                                {tmp_png_path}_residual.png \
                                -tile 3x1 -geometry \"175x175+2+0<\" \
                                {pj(out_png, gname)}_combined.png"
        
        for galfit_out in glob.glob(pj(cwd, "galfit.*")):
            ext_num = galfit_out.split(".")[-1]
            shutil.move(galfit_out, pj(out_dir, gname, f"{gname}_galfit.{ext_num}"))
        
        shutil.copy2(pj(tmp_galfits, f"{gname}_out.fits"), pj(out_dir, gname))
        
        print()


# In[122]:


## Tidying Up in case anything is leftover


# In[ ]:


if __name__ == "__main__":
    sp("rm galfit.* fit.log", capture_output = False)
    sp("rm *.png", capture_output = False)


# In[ ]:


if __name__ == "__main__":
    export_to_py("notebook_feedme_gen", output_filename = "sparc_to_galfit_feedme_gen.py")

