import os
import sys
import subprocess
import shutil
from glob import glob
import argparse

from os.path import join as pj
from os.path import exists

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
    
    USAGE = f"""USAGE:

    python3 ./{sys.argv[0]} [OPTION]
    
    OPTIONS =>[-v  | --verbose]

    This script is used for running the built-in unit and regression tests 
    for the GalfitModule.
    """
    # [--no-cleanup]
    
    parser = argparse.ArgumentParser(description = USAGE)
    
    # parser.add_argument('--no-cleanup',
    #                     dest     = 'cleanup',
    #                     action   = 'store_const',
    #                     const    = False,
    #                     # Soon to be the other way around
    #                     default  = True,
    #                     help     = 'Do not clean-up test output.')

    parser.add_argument('-v', '--verbose',
                        dest     = 'verbose', 
                        action   = 'store_const',
                        const    = True,
                        default  = False,
                        help     = 'Verbose output for as many bash commands as possible.')
    
    args    = parser.parse_args()
    # cleanup = args.cleanup
    verbose = args.verbose
    
    total_fail_count = 0
    
    # Cleaning up old reg test (if it was previously run)
    print("Cleaning up old unit/regression test files from TestOutput directory.")
    shutil.rmtree(TEST_OUTPUT_DIR)
    os.mkdir(TEST_OUTPUT_DIR)

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
            
            if verbose:
                print(result.stdout)
                print(result.stderr)
            
            if result.stderr:
                print(f"The unit test(s) for {thing_name} failed to run! Storing stderr.")
                error_list.append(result.stderr)
                fail_count += 1

        print(f"{fail_count} unit tests for {PyDir} failed to run.")
        return error_list, fail_count

    # Run all unit tests
    list_of_classes = ["Components", "Containers", "FitsHandlers"]
    list_of_helpers = ["HelperFunctions"]
    
    all_unit_error = {}
    all_unit_error["Helpers"], fail_helpers = run_unit_tests(list_of_helpers, PyDir = "Functions")
    all_unit_error["Classes"], fail_classes = run_unit_tests(list_of_classes, PyDir = "Classes")
    
    total_fail_count += fail_helpers + fail_classes

    error_file = "OutputError.txt"
    error_path = pj(TEST_OUTPUT_DIR, error_file)
    if exists(error_path):
        os.remove(error_path)
        
    with open(error_path, "w") as ef:
        for name, err_list in all_unit_error.items():
            err_str = "\n".join(err_list)
            ef.write(f"{name}\n{err_str}")
            
# Reg tests come next
if __name__ == "__main__":
    
    # Prepping for GALFIT run
    in_dir  = pj(TEST_DATA_DIR, "test-in")
    tmp_dir = pj(TEST_OUTPUT_DIR, "test-tmp")
    out_dir = pj(TEST_OUTPUT_DIR, "test-out")
    
    if not exists(tmp_dir):
        os.mkdir(tmp_dir)
        
    # if not exists(out_dir):
    #     os.mkdir(out_dir)
        
    data_out_dir  = pj(TEST_DATA_DIR, "test-out")
    shutil.copytree(data_out_dir, out_dir)
    
    # Run control script!
    ctrl_script = pj(_MODULE_DIR, "control_script.py")
    print("Running GALFIT (capturing output)...")
    result = sp(f"python3 {ctrl_script} {in_dir} {tmp_dir} {out_dir}")
    # try:
    #     shutil.move(pj(os.getcwd(), "galfit_failed.txt"), TEST_OUTPUT_DIR)
    # except FileNotFoundError:
    #     pass
    
    if verbose:
        print(result.stdout)
        print(result.stderr)
        
    # As output by HelperFunctions script per its own unit test
    with open(pj(TEST_OUTPUT_DIR, "UnitTestStdOuput.txt"), "w") as f:
        f.write(result.stdout)
    
    # DIFF CHECKS
    print("Performing diff check for unit tests.")
    
    # Starting with unit test output
    # Validate these are the same length
    unit_test_data = sorted(glob(pj(TEST_DATA_DIR, "*.txt")) + glob(pj(TEST_DATA_DIR, "*.in")))
    unit_test_out  = sorted(glob(pj(TEST_OUTPUT_DIR, "*.txt")) + glob(pj(TEST_OUTPUT_DIR, "*.in")))
    
    data_filenames   = [os.path.basename(i) for i in unit_test_data]
    output_filenames = [os.path.basename(i) for i in unit_test_out]
    
    compared = set(output_filenames).difference(set(data_filenames))
    
    if compared:
        print("Missing some output files from unit tests.")
        print("Diff check will fail so skipping that...")
        print("\n".join(compared))
        
    fail_count = 0
    all_diff_error = {}
    
    if not compared:
        for data, output, filename in zip(unit_test_data, unit_test_out, output_filenames):
            all_diff_error[filename] = ""

            result = sp(f"diff {data} {output}")

            if verbose:
                print(result.stdout)
                print(result.stderr)

            if result.stdout:
                print(f"Diff check for {filename} failed!")
                fail_count += 1
                all_diff_error[filename] = result.stdout
            
    # Diff checking GALFIT input/output
    # Diff does not work on FITS files
    # That and all the info we really need from those is in the galfit.# files anyway
    print("Performing diff check for Galfit input/output.")
    gnames = [os.path.basename(i) for i in glob(pj(in_dir, "*.fits"))]
    things_to_check = ["galfit.01", "galfit.02", ".in", "bulge.in"]
    
    for gname in gnames:
        data = pj(TEST_DATA_DIR, gname)
        output  = pj(TEST_OUTPUT_DIR, gname)
        for suffix in things_to_check:
            # This may produce a stderr if the file doesn't exist
            # say when galfit failes. Just a heads-up, shouldn't be a problem.
            result = sp(f"diff {data}_{suffix} {output}_{suffix}")
            
            if verbose:
                print(result.stdout)
                print(result.stderr)
        
            if result.stdout:
                print(f"Diff check for {gname}_{suffix} failed!")
                fail_count += 1
                all_diff_error[f"{gname}_{suffix}"] = result.stdout
                
    with open(error_path, "w") as ef:
        for name, err_str in all_diff_error.items():
            #err_str = "\n".join(err_list)
            ef.write(f"{name}\n{err_str}")
            
    print(f"{fail_count} regression tests failed.")
    if fail_count:
        print(f"See {error_path} for more information.")
        total_fail_count += fail_count
            
    print(f"Total number of tests failed: {total_fail_count}")