#!/usr/bin/env python
# coding: utf-8

# In[7]:


import os
import sys

from os.path import join as pj
from os.path import abspath as absp
from os.path import exists

import shutil
import subprocess
from copy import deepcopy
from astropy.io import fits


# In[6]:


# For debugging purposes
from IPython import get_ipython
def in_notebook():
    ip = get_ipython()
    
    if ip:
        return True
    else:
        return False


# In[5]:


_HOME_DIR = os.path.expanduser("~")
if in_notebook():
    _SPARCFIRE_DIR = pj(_HOME_DIR, "sparcfire_matt") 
    _MODULE_DIR    = pj(_SPARCFIRE_DIR, "GalfitModule")
else:
    try:
        _SPARCFIRE_DIR = os.environ["SPARCFIRE_HOME"]
        _MODULE_DIR = pj(_SPARCFIRE_DIR, "GalfitModule")
    except KeyError:
        print("SPARCFIRE_HOME is not set. Please run 'setup.bash' inside SpArcFiRe directory if not done so already.")
        print("Running on the assumption that GalfitModule is in your home directory... (if not this will fail and quit!)") 
        _MODULE_DIR = pj(_HOME_DIR, "GalfitModule")

sys.path.append(_MODULE_DIR)


# In[11]:


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


# In[12]:


def check_programs():

    # This seems to work in Python directly so I'm leaving it as-is
    # Checking galfit
    run_galfit = shutil.which("galfit")
    #run_galfit = response.stdout.strip()

    # Checking fitspng
    run_fitspng   = shutil.which("fitspng")

    # Checking exact python3 call
    run_python = shutil.which("python3")

    return run_galfit, run_fitspng, run_python

global run_galfit
global run_fitspng
global run_python
run_galfit, run_fitspng, run_python = check_programs()


# In[15]:


def sp(cmd_str, capture_output = True, timeout = None):
    # Because it is a pain in the butt to call subprocess with all those commands every time
    return subprocess.run(cmd_str, 
                          capture_output = capture_output, 
                          text = True, 
                          shell = True,
                          timeout = timeout,
                          executable="/bin/bash")


# In[18]:


if __name__ == "__main__":
    export_to_py("HelperFunctions", pj(_MODULE_DIR, "Functions", "HelperFunctions"))

