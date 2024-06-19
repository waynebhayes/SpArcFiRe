#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import sys
from glob import glob

from os.path import join as pj
from os.path import abspath as absp
from os.path import exists

import shutil
import subprocess
from copy import deepcopy
from astropy.io import fits


# In[2]:


# For debugging purposes
from IPython import get_ipython
def in_notebook():
    ip = get_ipython()
    
    if ip:
        return True
    else:
        return False


# In[3]:


_HOME_DIR = os.path.expanduser("~")
if in_notebook():
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


# In[4]:


def export_to_py(notebook_name, output_filename = ""):
    from IPython import get_ipython
    
    if not notebook_name.endswith(".ipynb"):
        notebook_name += ".ipynb"
    
    if in_notebook():
        print(f"Converting {notebook_name}")
        
        result = get_ipython().getoutput('jupyter nbconvert --to script {notebook_name}')
        
        if output_filename:
            filename = result[1].split()[-1]
            try:
                if not output_filename.endswith(".py"):
                    output_filename += ".py"
                    
                os.rename(filename, output_filename)
                
            except FileNotFoundError as f:
                print(f"Could not find {filename} per error {f}...")
                print("Output from nbconvert: ", *result)


# In[5]:


def sp(cmd_str, capture_output = True, timeout = None, **kwargs):
    # Because it is a pain in the butt to call subprocess with all those commands every time
    return subprocess.run(
        cmd_str, 
        capture_output = capture_output, 
        text           = True, 
        shell          = True,
        timeout        = timeout,
        executable     = "/bin/bash",
        #stdin          = subprocess.PIPE,
        **kwargs
     )


# In[6]:


def check_programs():

    # This seems to work in Python directly so I'm leaving it as-is
    # Checking galfit
    hostname = sp(f"hostname").stdout.split(".")[0]
    
    run_galfit = shutil.which(f"galfit_{hostname}")
    if not run_galfit:
        run_galfit = shutil.which(f"galfit")

    # Checking fitspng
    run_fitspng   = shutil.which("fitspng")

    # Checking exact python3 call
    run_python = shutil.which("python3")

    return run_galfit, run_fitspng, run_python

global run_galfit
global run_fitspng
global run_python
run_galfit, run_fitspng, run_python = check_programs()


# In[7]:


def find_files(search_dir = ".", search_pattern = "*", filetype = "f"):
    
    if filetype in ("d", "folder", "dir", "directory"):
        type_cmd = "d"
        
    elif filetype in ("f", "file"):
        type_cmd = "f"
        
    result = sp(f"find -L {pj(search_dir)} -maxdepth 1 -type {filetype} -name \"{search_pattern}\"")
    
    return [os.path.basename(i) for i in result.stdout.split("\n") if i]


# In[8]:


# Writing this to do generic deletion without calling subprocess
def rm_files(*args):
    # Assume list/tuple given by accident
    if not args:
        print("Nothing was fed into rm_files...")
        return

    if isinstance(args[0], (tuple, list)):
        print("Please expand your iterable(s) before feeding into rm_files, thanks!")
        print("Assuming the first argument is the only thing that needs to be deleted.")
        args = args[0]
    # Thanks! https://stackoverflow.com/a/8915613
    # May the list comp live on
    def catch(func, *argss, handle=lambda e : e, **kwargs):
        try:
            return func(*argss, **kwargs)
        except FileNotFoundError as e:
            return handle(e)
        
    return [catch(os.remove, i) for i in args]


# In[9]:


# Writing this to replace os.path.exists since that's too slow
def exists(filename):
    result = sp(f"[ -e {filename} ] && echo 1 || echo 0")
    return bool(int(result.stdout))


# In[10]:


def generate_get_set(input_dict): #, exclude = []):
    exec_str = ""
    for key,v in input_dict.items():
        exec_str += f"""@property
def {key}(self):
    return self.{v}
            
@{key}.setter
def {key}(self, new_val):
    self.{v} = new_val

"""
    return exec_str


# In[11]:


if __name__ == "__main__":
    from RegTest.RegTest import *


# In[12]:


# Unit test for sp
# The components and things will overwrite files rather than append
# so the second touch is unnecessary
if __name__ == "__main__":
    stdout_file   = "UnitTestStdOutput.txt"
    # writeout_file = "UnitTestWriteOuput.txt"
    
    stdout_dest   = pj(TEST_OUTPUT_DIR, stdout_file)
    # writeout_dest = pj(_MODULE_DIR, "RegTest", "TestOutput", writeout_file)
    
    touch_stdout  = sp(f"touch {stdout_dest}")
    # touch_writeout = sp(f"touch {writeout_dest}")
    
    # if touch_stdout.stderr or touch_writeout.stderr:
    if touch_stdout.stderr:
        print("Touch failed in helper_functions unit test.")
        print(touch_stdout.stderr)
        # print(touch_writeout.stderr)
        raise(Exception())


# In[13]:


# Unit test for rm files
if __name__ == "__main__":
    fake_files = [pj(TEST_OUTPUT_DIR, f"fake_{i}.fake") for i in range(10)]
    
    for fake in fake_files:
        with open(fake, mode='a'): pass
    
    # Check try catch
    fake_files.append(pj(TEST_OUTPUT_DIR, "fakest_of_them_all.fake"))
    
    # Check warning message
    print(rm_files(fake_files)) #*fake_files
    assert not glob(pj(TEST_OUTPUT_DIR, f"fake_*")), "Files were not deleted, something went wrong!!!"


# In[14]:


# Unit test for list_files
if __name__ == "__main__":
    
    print(sorted(find_files(pj(TEST_DATA_DIR, "test-in"), "*.fits", "f")))
    print()
    print(sorted(find_files(pj(TEST_DATA_DIR, "test-out"), "123*", "d")))


# In[15]:


# Unit test for exists
if __name__ == "__main__":
    
    print("Does test-in exist?", exists(pj(TEST_DATA_DIR, "test-in")))
    print("Does test-spout exist?", exists(pj(TEST_DATA_DIR, "test-spout")))


# In[16]:


if __name__ == "__main__":
    print(generate_get_set({"x" : "_x1"}))
    
    class Fake():
        def __init__(self):
            self._x1 = "bye"

        exec(generate_get_set({"x" : "_x1"}))
        
    y = Fake()
    print(y.x)
    y.x = "hi"
    print(y._x1)


# In[17]:


if __name__ == "__main__":
    export_to_py("helper_functions", pj(_MODULE_DIR, "Functions", "helper_functions"))

