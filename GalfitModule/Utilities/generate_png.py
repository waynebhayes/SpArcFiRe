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


def wrapper(
    gfits, 
    out_png_dir, 
    single_image = False, 
    hoirzontal = False, 
    cleanup = True
):
    
    if os.stat(gfits).st_size == 0:
        print(f"{os.path.basename(gfits)} is empty!!!")
        return

    if not single_image:
        output = OutputFits(gfits)
    else:
        output = FitsFile(gfits)
        output.gname = os.path.basename(gfits).split(".fits")[0]
        
    output.to_png(
        out_png_dir = out_png_dir, 
        horizontal = horizontal, 
        cleanup = cleanup
    )


if __name__ == "__main__":
    
    USAGE = f"""USAGE:

    python3 ./{sys.argv[0]} [OPTION] PNG-OUT-DIRECTORY [FITS-DIRECTORY | SPACE SEPARATED FITS FILES]
    
    OPTIONS =>[-h  | --help]
              [-s  | --single-image]
              [-nc | --no-cleanup]
              [-H  | --horizontal]
              [-v  | --verbose]
              [-ff | --from-file]

    This script is a utility to generate images from output GALFIT FITS files using the to_png 
    functionality of the OutputFits class.
    """
    
    parser = argparse.ArgumentParser(description = USAGE)
    
    parser.add_argument('-s', '--single-image',
                        dest     = 'single_image',
                        action   = 'store_const',
                        const    = True,
                        default  = False,
                        help     = 'Specify that the FITS files being read in are only a single image data block, i.e. model or observation only.'
                       )
    
    parser.add_argument('-H', '--horizontal',
                        dest     = 'horizontal',
                        action   = 'store_const',
                        const    = True,
                        default  = False,
                        help     = 'Choose to orient the output images from GALFIT horizontally as opposed to vertically (better for papers).'
                       )
    
    parser.add_argument('-nc', '--no-cleanup',
                        dest     = 'cleanup',
                        action   = 'store_const',
                        const    = False,
                        default  = True,
                        help     = 'Do not cleanup the individual pngs which are combined into the final montage\'d image.'
                       )
    
    parser.add_argument('-ff', '--from-file',
                        dest     = 'from_file',
                        action   = 'store_const',
                        const    = True,
                        default  = False,
                        help     = 'Use a file containing the image names from GALFIT rather than a list of fits read in on the command line.'
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
    
    single_image   = args.single_image
    horizontal     = args.horizontal
    cleanup        = args.cleanup
    
    from_file      = args.from_file
    
    verbose        = args.verbose
    capture_output = not args.verbose
    
    # Check if dir or file
    if len(list_o_fits) == 1:
        file_or_dirname = list_o_fits[0]
        glob_attempt = glob(file_or_dirname)
        
        if from_file:
            with open(file_or_dirname, "r") as f:
                list_o_fits = [fname.strip() for fname in f.readlines()]
        
        # If dir, prep list, if not do nothing
        elif os.path.isdir(file_or_dirname):
            list_o_fits = [pj(file_or_dirname, gfit) for gfit in find_files(file_or_dirname, "*.fits")]
            
        elif glob_attempt:
            list_o_fits = glob_attempt
        
        elif not file_or_dirname.endswith((".fits", ".fit", ".FITS", ".FIT")):
            print("Did you mean to specify that we are reading from file?")
            print("Please do so using the '-ff' option. Quitting.")
            sys.exit()
            
    #if "bayonet" in sp(f"hostname").stdout.split(".")[0]:
    if len(list_o_fits) > 500:
        Parallel(n_jobs = -2)(delayed(wrapper)(gfits, out_png_dir, single_image, horizontal, cleanup) for gfits in list_o_fits)
    else:
        _ = [wrapper(gfits, out_png_dir, single_image, horizontal, cleanup) for gfits in list_o_fits]
    #else:
    #    print("Not on bayonet, can't generate galaxies.")
