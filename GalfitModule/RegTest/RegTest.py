import os
import sys

from os.path import join as pj
#from os.path import exists

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
    _SPARCFIRE_DIR = pj(_HOME_DIR, "sparcfire_matt") 
    _MODULE_DIR    = pj(_SPARCFIRE_DIR, "GalfitModule")
else:
    try:
        _SPARCFIRE_DIR = os.environ["SPARCFIRE_HOME"]
        _MODULE_DIR = pj(_SPARCFIRE_DIR, "GalfitModule")
    except KeyError:
        # print("SPARCFIRE_HOME is not set. Please run 'setup.bash' inside SpArcFiRe directory if not done so already.")
        # print("Running on the assumption that GalfitModule is in your home directory... (if not this will fail and quit!)") 
        _MODULE_DIR = pj(_HOME_DIR, "GalfitModule")

sys.path.append(_MODULE_DIR)

# File paths to import
REG_TEST_DIR    = pj(_MODULE_DIR, "RegTest")
TEST_OUTPUT_DIR = pj(REG_TEST_DIR, "TestOutput")
TEST_DATA_DIR   = pj(REG_TEST_DIR, "TestData")

base_out = pj(TEST_OUTPUT_DIR, "UnitTest")