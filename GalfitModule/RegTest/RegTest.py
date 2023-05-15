import os
import sys
import subprocess

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

# Force >python 3.6 for f string compatibility
out_str = """\t Python3.6 or greater required! Exitting without generating feedmes... 
            if feedmes have already been generated, galfit will run with those.\n"""
assert sys.version_info >= (3, 6), out_str

if __name__ == "__main__":
    # Run all unit tests
    list_of_classes = ["Components", "Containers", "FitsHandlers"]
    list_of_helpers = ["HelperFunctions"]

    def sp(cmd_str, capture_output = True, timeout = None):
        # Because it is a pain in the butt to call subprocess with all those commands every time
        return subprocess.run(cmd_str, 
                              capture_output = capture_output, 
                              text = True, 
                              shell = True,
                              timeout = timeout,
                              executable="/bin/bash")

    def run_unit_tests(things_to_test, PyDir = ""):
        fail_count = 0
        error_list = []

        for thing_name in things_to_test:
            path_to_thing = pj(_MODULE_DIR, PyDir, thing_name)
            result = sp(f"python3 {path_to_thing}.py")
            if result.stderr:
                print(f"The unit test(s) for {thing_name} failed! Storing stderr.")
                error_list.append(result.stderr)
                fail_count += 1

        print(f"{fail_count} unit tests for {PyDir} failed.")
        return error_list

    all_error = {}
    all_error["Helpers"] = run_unit_tests(list_of_helpers, PyDir = "Functions")
    all_error["Classes"] = run_unit_tests(list_of_classes, PyDir = "Classes")

    error_file = "OutputError.txt"
    if all_error:
        with open(pj(TEST_OUTPUT_DIR, error_file), "w") as ef:
            for name, err_list in all_error.items():
                err_str = "\n".join(err_list)
                ef.write(f"{name}\n{err_str}")
            
# Reg tests come next