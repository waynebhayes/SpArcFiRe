import os
import sys
import subprocess
import shutil
from glob import glob
import argparse
import numpy as np
import pandas as pd

from os.path import join as pj
from os.path import exists

from astropy.io import fits

# For debugging purposes
# Keep for importing RegTest to Notebooks
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

# Force >python 3.7 for f string compatibility
out_str = """\t Python3.7 or greater required! Exitting without generating feedmes... 
            if feedmes have already been generated, galfit will run with those.\n"""
assert sys.version_info >= (3, 7), out_str

# ignore feedme filepaths
def iff(feedme_str):
    out_str = feedme_str.split("\n")
    G_idx = [idx for idx,i in enumerate(out_str) if i.startswith("G)")][0]
    return "\n".join(out_str[G_idx:])

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
    if exists(TEST_OUTPUT_DIR):
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
        fail_count  = 0
        error_list  = []
        stdout_list = []

        for thing_name in things_to_test:
            path_to_thing = pj(_MODULE_DIR, PyDir, thing_name)
            result = sp(f"python3 {path_to_thing}.py")
            stdout_list.append(result.stdout)
            
            if verbose:
                print(result.stdout)
                print(result.stderr)
            
            if result.stderr:
                print(f"The unit test(s) for {thing_name} failed to run! Storing stderr.")
                error_list.append(result.stderr)
                fail_count += 1

        print(f"{fail_count} unit tests for {PyDir} failed to run.")
        return stdout_list, error_list, fail_count

    # Run all unit tests
    list_of_classes = ["Components", "Containers", "FitsHandlers"]
    list_of_helpers = ["helper_functions"]
    
    all_stdout     = {}
    all_unit_error = {}
    
    all_stdout["Helpers"], all_unit_error["Helpers"], fail_helpers = run_unit_tests(list_of_helpers, PyDir = "Functions")
    all_stdout["Classes"], all_unit_error["Classes"], fail_classes = run_unit_tests(list_of_classes, PyDir = "Classes")
    
    total_fail_count += fail_helpers + fail_classes

    error_file = "OutputError.txt"
    error_path = pj(TEST_OUTPUT_DIR, error_file)
    # if exists(error_path):
    #     os.remove(error_path)
        
    with open(error_path, "a") as ef:
        for name, err_list in all_unit_error.items():
            err_str = "\n".join(err_list)
            ef.write(f"{name}\n{err_str}")

    stdout_path = pj(TEST_OUTPUT_DIR, "UnitTestStdOutput.txt")
    if not exists(stdout_path):
        print("Function helper must have failed! Adding to failure count and creating file.")
        total_fail_count += 1
    
    # As output by helper_functions script per its own unit test
    with open(stdout_path, "a") as f:
        for name, out_list in all_stdout.items():
            out_str = "\n".join(out_list)
            f.write(f"{name}\n{out_str}")
            f.write("="*80)
            f.write("="*80)
            
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
    base_run_name = "RegTest"
    steps = ("1-step", "2-step")
    
    print("Running GALFIT with single-step fitting (capturing output)...")
    result = sp(f"python3 {ctrl_script} -p 0 -NS 1 -n {base_run_name}_{steps[0]} {in_dir} {tmp_dir} {out_dir}")
    
    if verbose:
        print(result.stdout)
        print(result.stderr)
    
    print("Running GALFIT with default two-step fitting (capturing output)...")
    result = sp(f"python3 {ctrl_script} -p 0 -n {base_run_name}_{steps[1]} {in_dir} {tmp_dir} {out_dir}")
    
    # try:
    #     shutil.move(pj(os.getcwd(), "galfit_failed.txt"), TEST_OUTPUT_DIR)
    # except FileNotFoundError:
    #     pass
    
    if verbose:
        print(result.stdout)
        print(result.stderr)
        
    # As output by helper_functions script per its own unit test
    # Will have a filepath issue at the end of this as well as potentially an ordering
    # issue so we don't use this output after all.
    # with open(stdout_path, "a") as f:
    #     f.write(result.stdout)
    
    # DIFF CHECKS
    print("Performing diff check for unit tests.")
    
    # Starting with unit test output
    # Validate these are the same length
    unit_test_data = sorted(glob(pj(TEST_DATA_DIR, "*.txt")) + glob(pj(TEST_DATA_DIR, "*.in")))
    unit_test_out  = sorted(glob(pj(TEST_OUTPUT_DIR, "*.txt")) + glob(pj(TEST_OUTPUT_DIR, "*.in")))
    
    unit_test_data = [i for i in unit_test_data if os.path.basename(i) != "OutputError.txt"]
    unit_test_out  = [i for i in unit_test_out  if os.path.basename(i) != "OutputError.txt"]
    
    data_filenames   = [os.path.basename(i) for i in unit_test_data]
    output_filenames = [os.path.basename(i) for i in unit_test_out]
    
    compared = set(output_filenames).difference(set(data_filenames))
    
    if compared:
        print("Missing some output files from unit tests.")
        print("Diff check will fail so skipping the following...")
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
                #print(result.stdout)
                print(f"Diff check for {filename} failed!")
                if filename == "UnitTestStdOutput.txt":
                    print(f"Since {filename} contains the standard output for several tests")
                    print(f"the number of unit tests which actually failed may be higher than the final count.")
                
                fail_count += 1
                all_diff_error[filename] = f"{filename}\n{result.stdout}"
            
    # Diff checking GALFIT input 
    # Output may vary depending on compute differences b/t clusters (GALFIT issue)
    print("Performing diff check for Galfit input.")
    gnames = [os.path.basename(i).rstrip(".fits") for i in glob(pj(in_dir, "*.fits"))]
    things_to_check = [".in", "_disk.in"]
    
    for gname in gnames:
        data    = pj(TEST_DATA_DIR, "test-out", gname, gname)
        output  = pj(TEST_OUTPUT_DIR, "test-out", gname, gname)
        
        # Check galfit text output files
        for suffix in things_to_check:
            # This may produce a stderr if the file doesn't exist
            # say when galfit failes. Just a heads-up, shouldn't be a problem.
            # tail is to start at line 15 of each file to ignore the file path differences
            # for input/output/star mask/etc.
            result = sp(f"diff <(tail -n +15 {data}{suffix}) <(tail -n +15 {output}{suffix})")
            
            if verbose:
                print(result.stdout)
                print(result.stderr)
        
            if result.stdout:
                # Deprecated choice
                # Carve a special exception for 1237655463239155886 since we randomize the spin parity
                #if result.stdout.split("\n")[0] == "60c60" and gname == "1237655463239155886":
                #    continue
                    
                print(f"Diff check for {gname}{suffix} failed!")
                fail_count += 1
                all_diff_error[f"{gname}{suffix}"] = f"{data}{suffix}\n{result.stdout}"
                
        # Check _out_old FITS files using Astropy since this RegTest shouldn't call FitsHandler except to check it
#         suffixes = ["galfit_out.fits", "galfit_out_old.fits"]
#         for suffix in suffixes:
#             if exists(f"{data}_{suffix}"):
#                 data_fits   = fits.open(f"{data}_{suffix}")[2].data
                
#                 try:
#                     output_fits = fits.open(f"{output}_{suffix}")[2].data
#                 except FileNotFoundError as fnfe:
#                     fail_count += 1 
#                     all_diff_error[f"{gname}_{suffix}"] = f"{gname} did not successfully fit, something went wrong."
#                     continue

#                 # The tolerance can be quite high because when GALFIT is off, it will differ significantly
#                 # I set a tolerance in the first place for machine error/differences under the hood across clusters
#                 if not np.allclose(data_fits, output_fits, atol = 1e-2):
#                     fail_count += 1
#                     all_diff_error[f"{gname}_{suffix}"] = f"Single Component Fit differs for {gname}!"

#             elif exists(f"{output}_{suffix}"):
#                 fail_count += 1
#                 all_diff_error[f"{gname}_{suffix}"] = f"{gname} should not have successfully fit, something went... wrong (right?)."
          
    with open(error_path, "a") as ef:
        for name, err_str in all_diff_error.items():
            #err_str = "\n".join(err_list)
            ef.write(f"{name}\n{err_str}\n")
            
    print(f"{fail_count} unit tests failed.")
    if fail_count:
        print(f"See {error_path} for more information.")
        total_fail_count += fail_count
        
    print("Checking existence and readability of output pkl files.")
    for s in steps:
        filename = pj(out_dir, f"{base_run_name}_{s}_output_results.pkl")
        if not exists(filename):
            print(f"{filename} does not exist!")
            total_fail_count += 1
        else:
            try:
                _ = pd.read_pickle(filename)
                
            except Exception as e:
                print(f"Something went wrong reading {filename} via pandas.")
                total_fail_count += 1
                
                if verbose:
                    print(e)
            
    print(f"Total number of tests failed: {total_fail_count}")