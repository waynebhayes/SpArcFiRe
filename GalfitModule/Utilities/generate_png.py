from joblib import Parallel, delayed
import joblib
import os

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

def wrapper(gfits):
    output = OutputFits(gfits)
    output.to_png(out_png_dir = out_png_dir)


if __name__ == "__main__":
    
    out_png_dir = sys.argv[1]
    list_o_fits = sys.argv[2].split(" ")
    
    if sp(f"hostname").stdout.split(".")[0] == "bayonet-09":
        if len(list_o_fits) > 500:
            Parallel(n_jobs = -2)(delayed(wrapper)(gfits) for gfits in list_o_fits)
        else:
            _ = [wrapper(gfits) for gfits in list_o_fits]
    else:
        print("Not on bayonet, can't generate galaxies.")
