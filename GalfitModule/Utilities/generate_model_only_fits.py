import os
from os.path import basename
from os.path import join as pj
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
    
#from Classes.Components import *
#from Classes.Containers import *
# FitsHandlers imports all of the above
from go_go_galfit import generate_model_for_sparcfire as gmfs
from Classes.FitsHandlers import *
from Functions.helper_functions import *

from glob import glob

if __name__ == "__main__":
    
    cleanup = True
    
    if len(sys.argv) == 3:
        in_dir  = sys.argv[1]
        tmp_dir = sys.argv[2]
        #out_dir = sys.argv[3]
        
    else:
        cwd = os.getcwd()
        in_dir = pj(cwd, "sparcfire-in")
        tmp_dir = pj(cwd, "sparcfire-tmp")
        #out_dir = pj(cwd, "sparcfire-out")
    
    run_galfit, _, _ = check_programs()
    
    galfits_dir = pj(tmp_dir, "galfits")
    _ = [
        gmfs(gfits, run_galfit, in_dir, tmp_fits_dir = galfits_dir)
        for gfits in glob(pj(galfits_dir, "*_galfit_out.fits"))
        if not basename(gfits).startswith("failed")
    ]
    
    if cleanup:
        rm_files(*glob(pj(galfits_dir, "*.in")))
