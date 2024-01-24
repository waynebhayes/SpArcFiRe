from joblib import Parallel, delayed
import joblib
import os
import argparse

from os.path import join as pj
from os.path import abspath as absp

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

def wrapper(gfits, out_png_dir, vertical = False):
    output = OutputFits(gfits)
    output.to_png(out_png_dir = out_png_dir, vertical = vertical)


if __name__ == "__main__":
    
    USAGE = f"""USAGE:

    python3 ./{sys.argv[0]} [OPTION] PNG-OUT-DIRECTORY [FITS-DIRECTORY | SPACE SEPARATED FITS FILES]
    
    OPTIONS =>[-vert | --vertical]
              [-v | --verbose]

    This script is a utility to generate images from output GALFIT FITS files using the to_png 
    functionality of the OutputFits class.
    """
    
    parser = argparse.ArgumentParser(description = USAGE)
    
    parser.add_argument('-vert', '--vertical',
                        dest     = 'vertical',
                        action   = 'store_const',
                        const    = True,
                        default  = False,
                        help     = 'Choose to orient the output images from GALFIT vertically as opposed to horizontally.'
                       )
    
    parser.add_argument('-v', '--verbose',
                        dest     = 'verbose', 
                        action   = 'store_const',
                        const    = True,
                        default  = False,
                        help     = 'Verbose output. Includes stdout.'
                       )
    
    parser.add_argument(dest     = 'out_png_dir',
                        help     = "Directory designated to hold the output PNGs."
                       )
    
    parser.add_argument(dest     = 'list_o_fits',
                        nargs    = "*",
                        help     = "Either a list of fits files or a directory containing the files to be output to PNG."
                       )
    
    args           = parser.parse_args() # Using vars(args) will call produce the args as a dict
    out_png_dir    = args.out_png_dir
    list_o_fits    = args.list_o_fits
    
    vertical       = args.vertical
    verbose        = args.verbose
    capture_output = not args.verbose
    
    # Check if dir or file
    if len(list_o_fits) == 1:
        file_or_dirname = list_o_fits[0]
        
        # If dir, prep list, if not do nothing
        if os.path.isdir(file_or_dirname):
            list_o_fits = [pj(file_or_dirname, gfit) for gfit in find_files(file_or_dirname, "*.fits")]
    
    if "bayonet" in sp(f"hostname").stdout.split(".")[0]:
        if len(list_o_fits) > 500:
            Parallel(n_jobs = -2)(delayed(wrapper)(gfits, out_png_dir, vertical) for gfits in list_o_fits)
        else:
            _ = [wrapper(gfits, out_png_dir, vertical) for gfits in list_o_fits]
    else:
        print("Not on bayonet, can't generate galaxies.")
