import sys
import os
from os.path import join as pj

import argparse
import shutil
import subprocess

import glob

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

from Classes.Components import *
from Classes.Containers import *
from Functions.helper_functions import *

if __name__ == "__main__":

    # The file or base directory from which to gather a list of fits files
    OLD_SPARC_IN = sys.argv[1]
    gal_sparc_in = find_files(arg1, "123*.fits", "f")
    assert gal_sparc_in, f"No fits files found in {OLD_SPARC_IN}."
    
    NEW_BASE_DIR = sys.argv[2] #/extra/wayne_scratch0/preserve/portman/29k_galaxies_gband
    #SDSS_filepath = "/home/wayne/research/drdavis/SDSS/FITS/2015/color/r"
    COLOR    = sys.argv[3]
    DATA_DIR = f"/extra/wayne1/research/drdavis/SDSS/SpArcFiRe/2016-09/{COLOR}/"

    NEW_SPARC_IN  = pj(NEW_BASE_DIR, "sparcfire-in") 
    NEW_SPARC_TMP = pj(NEW_BASE_DIR, "sparcfire-tmp") 
    NEW_SPARC_OUT = pj(NEW_BASE_DIR, "sparcfire-out") #"/extra/wayne1/preserve/portmanm/sparcfire-out"

    # Make sparcfire-in if not found
    # No need to do this for sparcfire-out, I use makedirs to cover that later
    if not exists(NEW_SPARC_IN):
        os.mkdir(NEW_SPARC_IN)

    if not exists(NEW_SPARC_TMP):
        os.mkdir(NEW_SPARC_TMP)
    
    copy_filename = f"{os.path.basename(NEW_BASE_DIR)}_copy_galaxy_folders.sh" #pj(NEW_BASE_DIR, f"slurm_copy_galaxy_folders")
    
    write2 = True
    if write2:
        gnames = [os.path.basename(gal).replace(".fits", "") for gal in gal_sparc_in]
        bad_gnames = []
        
        if exists(copy_filename):
            _ = sp(f"rm {copy_filename}", capture_output = False)
        
        with open(copy_filename, "w") as f:
            #for gname in gnames:
            for gname in gnames:
                folder_path = pj(DATA_DIR, gname[-3:], gname)
                new_folder  = pj(NEW_SPARC_OUT, gname)
                if exists(folder_path):
                    if not exists(new_folder):
                        os.makedirs(new_folder)

                    #_ = sp(f"rsync {folder_path}/*.csv {new_folder}")
                    #f.write(f"rsync -u {folder_path}/*.csv* {new_folder}\n")
                    f.write(f"ln -sf {OLD_SPARC_IN}/{gname}.fits {pj(NEW_SPARC_IN, gname + '.fits')} && ")
                    f.write(f"ln -sf {folder_path}/{gname}.csv {pj(new_folder, gname + '.csv')} && ")
                    f.write(f"ln -sf {folder_path}/{gname}_arcs.csv {pj(new_folder, gname + '_arcs.csv')} && ")
                    f.write(f"ln -sf {folder_path}/{gname}.csv+CR {pj(new_folder, gname + '.csv+CR')}\n")
                    # Include tsvs just in case no csvs are found
                    f.write(f"ln -sf {folder_path}/{gname}.tsv {pj(new_folder, gname + '.tsv')} && ")
                    f.write(f"ln -sf {folder_path}/{gname}_arcs.tsv {pj(new_folder, gname + '_arcs.tsv')} && ")
                    f.write(f"ln -sf {folder_path}/{gname}.tsv+CR {pj(new_folder, gname + '.tsv+CR')}\n")
                    # Include elps just in case no crop rad is in the csv or tsv file
                    f.write(f"ln -sf {folder_path}/{gname}-elps-fit-params.txt {pj(new_folder, gname + '-elps-fit-params.txt')} && ")
                else:
                    bad_gnames.append(gname)

        if bad_gnames:
            print("This many bad galaxies...", len(bad_gnames))
            #print("No output folder found for:")
            #print(*bad_gnames, sep = "\n")
    
    DUMP_NAME = "SYMLINKING_FOLDERS"

    slurm2 = False
    if slurm2:
        _ = sp(f"chmod a+rx {copy_filename}")
        _ = sp(f"cat {copy_filename} | ~wayne/bin/distrib_slurm {DUMP_NAME} -M all", capture_output = False)

    cleanup = False
    if cleanup:
        _ = sp(f"rm -rf ~/SLURM_turds/{DUMP_NAME}", capture_output = False)

