import os
from os.path import basename
from os.path import join as pj
import argparse
import sys

from joblib import Parallel, delayed
import joblib

_HOME_DIR = os.path.expanduser("~")

try:
    _SPARCFIRE_DIR = os.environ["SPARCFIRE_HOME"]
    _MODULE_DIR = pj(_SPARCFIRE_DIR, "GalfitModule")
except KeyError:
    # print("SPARCFIRE_HOME is not set. Please run 'setup.bash' inside SpArcFiRe directory if not done so already.")
    # print("Running on the assumption that GalfitModule is in your home directory... (if not this will fail and quit!)") 
    _MODULE_DIR = pj(_HOME_DIR, "GalfitModule")
    
sys.path.append(_MODULE_DIR)
    
#from Classes.Components import *
#from Classes.Containers import *
# FitsHandlers imports all of the above
from go_go_galfit import generate_model_for_sparcfire as gmfs
from Classes.FitsHandlers import *
from Functions.helper_functions import *

from glob import glob

# NOTE this is also in calculate pitch-angle...
# TODO: Decide which to keep
def generate_model_for_pitch_angle(
    gfits, 
    base_galfit_cmd, 
    in_dir, 
    tmp_fits_dir = "",
    output_dir   = ""
):
    gname  = os.path.basename(gfits).replace("_galfit_out.fits", "")
    suffix = "for_pitch_angle"
    if not tmp_fits_dir:
        tmp_fits_dir = os.path.split(gfits)[0]
        
    fits_file = OutputFits(gfits, load_default = False)
    feedme = fits_file.feedme
    
    feedme.header.input_image.value   = pj(in_dir, f"{gname}.fits")
    feedme.header.output_image.value  = pj(tmp_fits_dir, f"{gname}_{suffix}.fits")
    if output_dir:
        feedme.header.output_image.value  = pj(output_dir, f"{gname}_{suffix}.fits")
        
    feedme.header.optimize.value      = 1
    
    if "power_0" in feedme.components: #or "arms" in feedme.components:
        # power_comp_number = feedme.power_0.component_number
        # Turning off inclination
        feedme.power_0.inclination.value = 0
        
    else:
        print("power_0 component not found in feedme header. Something went wrong!")
        return feedme.header.output_image.value
        
        
    filename = pj(tmp_fits_dir, f"{gname}_{suffix}.in")
    feedme.to_file(filename = filename)
    _ = sp(f"{base_galfit_cmd} {filename}", capture_output = True)
    
    return feedme.header.output_image.value

if __name__ == "__main__":
    
    USAGE = f"""USAGE:

    python3 ./{sys.argv[0]} [OPTIONS] IN-DIR TMP-DIR [GALFITS-DIR] [OUTPUT-DIR]
    
    OPTIONS =>[-h   | --help]
              [-fs  | --for-sparcfire] (default)
              [-noi | --no-inclination] (for pitch angle plots)
              [-v   | --verbose]

    This script is a utility to generate model-only FITS from output GALFIT FITS files.
    If GALFITS-DIR is not supplied, it is assumed to be TMP-DIR/galfits.
    If OUTPUT-DIR is not supplied, it is assumed to be GALFITS-DIR.
    """
    
    parser = argparse.ArgumentParser(description = USAGE)
    
    parser.add_argument('-fs', '--for-sparcfire',
                        dest     = 'for_sparcfire',
                        action   = 'store_const',
                        const    = True,
                        default  = False,
                        help     = 'Generate models to aid SpArcFiRe in its processing.'
                       )
    
    parser.add_argument('-noi', '--no-inclination',
                        dest     = 'for_pitch_angle',
                        action   = 'store_const',
                        const    = True,
                        default  = False,
                        help     = 'Generate models without inclination to create pitch angle overlay.'
                       )
    
    parser.add_argument('-nc', '--no-cleanup',
                        dest     = 'no_cleanup', 
                        action   = 'store_const',
                        const    = True,
                        default  = False,
                        help     = 'Do not clean-up generated .in files.'
                       )
    
    parser.add_argument('-v', '--verbose',
                        dest     = 'verbose', 
                        action   = 'store_const',
                        const    = True,
                        default  = False,
                        help     = 'Verbose output.'
                       )
    
    parser.add_argument(dest     = 'dirs',
                        nargs    = "*",
                        help     = "Dirs used for the generation of the models: the in-dir, out-dir, and optional directory for outputting the new models."
                       )
    
    args            = parser.parse_args() # Using vars(args) will call produce the args as a dict
    for_sparcfire   = args.for_sparcfire
    for_pitch_angle = args.for_pitch_angle
    
    cleanup        = not args.no_cleanup
    verbose        = args.verbose
    capture_output = not args.verbose
    
    dirs = args.dirs
    
    if len(dirs) >= 2:
        in_dir  = dirs[0]
        tmp_dir = dirs[1]
        galfits_dir = pj(tmp_dir, "galfits")
        output_dir  = galfits_dir
        
        if len(dirs) >= 3:
            galfits_dir = dirs[2]
            
        if len(dirs) >= 4:
            output_dir = dirs[3]
        
    else:
        cwd = os.getcwd()
        in_dir = pj(cwd, "sparcfire-in")
        tmp_dir = pj(cwd, "sparcfire-tmp")
        galfits_dir = pj(tmp_dir, "galfits")
        output_dir  = galfits_dir
        #out_dir = pj(cwd, "sparcfire-out")
    
    run_galfit, _, _ = check_programs()
    
    if for_sparcfire:
        func = gmfs
    
    elif for_pitch_angle:
        func = generate_model_for_pitch_angle
        
    else:
        func = gmfs
    
    if not exists(galfits_dir):
        os.makedirs(galfits_dir)
        
    if not exists(output_dir):
        os.makedirs(output_dir)
        
    list_o_fits = glob(pj(galfits_dir, "*_galfit_out.fits"))
                       
    if len(list_o_fits) > 500:
        _ = Parallel(n_jobs = -2)(
            delayed(func)(gfits, run_galfit, in_dir, galfits_dir, output_dir) 
            for gfits in list_o_fits
            if not basename(gfits).startswith("failed")
        )
    else:
        _ = [
            func(gfits, run_galfit, in_dir, galfits_dir, output_dir) 
            for gfits in list_o_fits
            if not basename(gfits).startswith("failed")
            ]
        
    if cleanup:
        rm_files(*glob(pj(galfits_dir, "*.in")))
