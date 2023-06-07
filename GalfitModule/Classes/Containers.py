#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import sys
from os.path import join as pj
from os.path import exists
import subprocess
from copy import deepcopy
from IPython import get_ipython
from astropy.io import fits

import numpy as np
import pandas as pd


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
        # print("SPARCFIRE_HOME is not set. Please run 'setup.bash' inside SpArcFiRe directory if not done so already.")
        # print("Running on the assumption that GalfitModule is in your home directory... (if not this will fail and quit!)") 
        _MODULE_DIR = pj(_HOME_DIR, "GalfitModule")

sys.path.append(_MODULE_DIR)

from Classes.Components import *
# Relies on relative file hierarchy
from Functions.helper_functions import *


# In[4]:


class ComponentContainer:
    def __init__(self, **kwargs):
        self.bulge   = kwargs.get("bulge", Sersic(1))
        self.disk    = kwargs.get("disk", Sersic(2))
        self.arms    = kwargs.get("arms", Power(2))
        self.fourier = kwargs.get("fourier", Fourier(2))
        self.sky     = kwargs.get("sky", Sky(3))
        
    def to_dict(self):
        # Order matters
        return {"bulge"   : self.bulge,
                "disk"    : self.disk,
                "arms"    : self.arms,
                "fourier" : self.fourier,
                "sky"     : self.sky}
    
    def to_tuple(self):
        #self.header,
        return (self.bulge,
                self.disk,
                self.arms,
                self.fourier,
                self.sky
               )
    
    def to_list(self):
        #self.header,
        return [self.bulge,
                self.disk,
                self.arms,
                self.fourier,
                self.sky
               ]
    
    def to_pandas(self):
        return pd.concat([comp.to_pandas().reset_index() 
                          for comp in ComponentContainer.to_list(self)], 
                         axis = 1).drop(columns = ["index"])
        
    def from_pandas(self, input_df):
        pass
    
    def update_components(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)
       
    # This is more for the daughter classes
    def extract_components(self):
        return ComponentContainer(**vars(self))
    
    def __str__(self):
        out_str = "\n".join([str(comp) for comp in ComponentContainer.to_list(self)])
        return out_str


# In[5]:


class FeedmeContainer(ComponentContainer):
    def __init__(self, **kwargs):
        ComponentContainer.__init__(self, **kwargs)
        # The path to the feedme that *generated* the components
        self.header  = kwargs.get("header", GalfitHeader())
        self.path_to_feedme = kwargs.get("path_to_feedme", "")
    
    def to_dict(self):
        return vars(self)
    
    def to_list(self):
        return [self.header] + ComponentContainer.to_list(self)
    
    def __str__(self):
        out_str = f"{str(self.header)}\n" + "\n".join(str(comp) for comp in ComponentContainer.to_list(self))
        return out_str
        
    def to_file(self, *args, filename = ""):
        if not filename:
            filename = self.path_to_feedme
            
        if args:
            self.header.to_file(filename, *args)
        else:
            self.header.to_file(filename, *ComponentContainer.to_list(self))
            
    def from_file(self, obj_in):
        # This function handles grabbing and storing the values from galfit files (input and output???)
        # It's written to generally handle both and stores everything in the respective component objects
        # obj_in could be a filename or a dict (per fitshandler)
        
        def open_file(self, open_type, filename = obj_in):
            try: 
                # Grabbing the filename
                #input_filename = glob_name(galaxy_path, '', filename) 
                input_file = open_type(filename)

            except FileNotFoundError:
                print(f"Can't open to read the file, {filename}. Check name/permissions/directory.")
                return None

            except OSError as ose:
                print(f"Something went wrong! {ose}")
                return None
            
            return input_file

        def from_fits(self, input_dict = obj_in):
            
            input_keys = list(input_dict.keys())
            
            component_list = self.to_list()
            component_num = 0
            send_to_helper = {}
            
            for idx, (key, value) in enumerate(input_dict.items()):
                component = component_list[component_num]
                
                try:
                    component_idx_start = input_keys.index(component.start_dict)
                    component_idx_end   = input_keys.index(component.end_dict)
                except ValueError as ve:
                    # Trying to recover...
                    # End will *always* (header excluded) be #_param                
                    component_end = [k for k in input_keys if k.endswith(self.end_dict[2:])][0]

                    if component_end[0].isnumeric():
                        component_idx_start = input_keys.index(f"{component_end[0]}_{self.start_dict[2:]}")
                        component_idx_end   = input_keys.index(component_end)

                    else:
                        print(f"Can't find start/end of {self.component_type} segment.")
                        print(f"Check the filename or start/end_dict variables.")
                        print(f"Filename: {filename}")
                        print(f"Start/End: {self.start_dict}/{self.end_dict}")
                        raise ValueError(ve)

                if component_idx_start <= idx <= component_idx_end:
                    send_to_helper[key] = value
                
                if idx == component_idx_end:
                    component.from_file_helper(send_to_helper)
                    send_to_helper = {}
                    component_num += 1
                    
                    if component_num == len(component_list): break
                    
            return
        
        def from_text(self, input_file):
            
            component_list = self.to_list()
            component_num = 0
            send_to_helper = []
            store = False
            
            for line in input_file:
                component = component_list[component_num]
                
                if line.strip().startswith(component.start_text):
                    store = True
                    
                if store:
                    send_to_helper.append(line)
                    
                if line.strip().startswith(component.end_text):
                    store = False
                    component.from_file_helper(send_to_helper)
                    send_to_helper = []
                    component_num += 1
                    
                    if component_num == len(component_list): break
        
            return
        
        # 
        if isinstance(obj_in, dict):
            from_fits(self)
            
        elif os.path.splitext(obj_in)[1] == ".fits":
            input_file = open_file(self, fits.open)
            input_dict = dict(input_file[2].header)
            from_fits(self, input_dict)
            input_file.close()
            
        else:
            input_file = open_file(self, open)
            from_text(self, input_file)
            input_file.close()
        
        # This is not needed now that we update
        # param_values in the file helpers
        #_ = [c.update_param_values() for c in self.to_list()]


# In[6]:


class OutputContainer(FeedmeContainer):
    def __init__(self, galfit_out_obj = subprocess.CompletedProcess("", 0), **kwargs):
        
        FeedmeContainer.__init__(self, **kwargs)
        
        galfit_out_text = galfit_out_obj.stdout        
        galfit_err_text = galfit_out_obj.stderr
        
        # Default to this so it doesn't break if no text is fed in
        self.success = False
        
        # This is more for cleaning up some functions elsewhere
        # self.galfit_num = kwargs.get("galfit_num", "01")
        
        def check_success(self, galfit_out_text) -> None:
            # I don't like checking for the full line because there are embedded quotes
            # and one is a backtick I think...: Fit summary is now being saved into `fit.log'.
            # To be safe I just check the first part
            success = "Fit summary is now being saved"
            failure = "...now exiting to system..."

            # Constrain to last 50 lines to save search time
            # The explosion is 23 lines
            last_out_lines = galfit_out_text.split("\n")[-50:]
            if any(line.strip().startswith(success) for line in last_out_lines):
                self.success = True
                
            elif any(line.strip().startswith(failure) for line in last_out_lines):
                print(f"Galfit failed this run!")
                last_out_lines = '\n'.join(last_out_lines)
                #print(f"{last_out_lines}")
                self.success = False
                # For debugging
                # print(galfit_out_text)

            else:
                print(f"Did not detect either '{success}' or '{failure}' in galfit output. Something must have gone terribly wrong! Printing output...")
                print(f"{err_text}")
                self.success = False
        
        if galfit_out_text:
            check_success(self, galfit_out_text)

        # For reading from galfit stdout to update classes
        def update_components(self, galfit_out_text) -> None: #, bulge, disk, arms, fourier, sky):
            
            last_it = galfit_out_text.split("Iteration")[-1]

            s_count = 0
            by_line = last_it.splitlines()
            for line in by_line:
                if line.strip().startswith("sersic"):
                    if s_count == 0:
                        comp = self.bulge
                        s_count += 1 
                    else:
                        comp = self.disk               

                elif line.strip().startswith("power"):
                    comp = self.arms

                elif line.strip().startswith("fourier"):
                    # Fourier never follows the rules...
                    # By how the class operates, we don't
                    # need to do this, but to be consistent...
                    comp = self.fourier

                elif line.strip().startswith("sky"):
                    comp = self.sky

                else:
                    continue

                comp.update_from_log(line)
                # This is now done in Components
                #comp.update_param_values()

        if self.success:
            update_components(self, galfit_out_text)
            
        if kwargs.get("store_text", False):
            self.galfit_out_text = galfit_out_text
            self.galfit_err_text = galfit_err_text
        else:
            #self.galfit_out_text = ""
            self.galfit_out_text = "Galfit output text not stored.\n" + \
                                   "Set keyword argument 'store' = True in OutputContainer call to store it " + \
                                   "or otherwise pass the text as an argument to the object when converting to string."
            
            self.galfit_err_text = ""
            
    def __str__(self, galfit_out_text = "") -> str:
        if galfit_out_text:
            return(galfit_out_text)
        elif self.galfit_out_text:
            return self.galfit_out_text
        # Shouldn't happen but just in case.
        else:
            return ""


# In[7]:


if __name__ == "__main__":
    # Testing basic functionality
    
    container = ComponentContainer()
    print(container)
    container_df = container.to_pandas()
    print()
    print(container_df)


# In[8]:


if __name__ == "__main__":
    from RegTest.RegTest import *


# In[9]:


# Testing FeedmeContainer kwargs and to_file
if __name__ == "__main__":
    
    header = GalfitHeader(galaxy_name = "tester")
    bulge = Sersic(1, position = (25,25))
    disk  = Sersic(2, position = (25,25))
    arms  = Power(2)
    fourier = Fourier(2)
    sky   = Sky(3)
    
    container = FeedmeContainer(**{"header"  : header,
                                      "bulge"   : bulge,
                                      "disk"    : disk,
                                      "arms"    : arms,
                                      "fourier" : fourier,
                                      "sky"     : sky}
                                )

    print()
    print(container.to_pandas())
    
    #print(container)
    container.path_to_feedme = pj(TEST_OUTPUT_DIR, f"{base_out}_FeedmeContainer.in")
    
    container.to_file()


# In[10]:


# Testing FeedmeContainer from_file
if __name__ == "__main__":
    
    header = GalfitHeader(galaxy_name = "fake_name")
    container.update_components(header = header)
    # bulge = Sersic(1)
    # disk  = Sersic(2)
    # arms  = Power()
    # fourier = Fourier()
    # sky   = Sky(3)

    example_fits = pj(TEST_DATA_DIR, "test-out", "1237667911674691747", "1237667911674691747_galfit_out.fits")
    example_feedme = pj(TEST_DATA_DIR, "test-out", "1237667911674691747", "1237667911674691747.in")
    
    print("These are feedme -> output")
    print("ignoring filepaths for reg tests...\n")
    
    container.from_file(example_feedme)
    print(iff(str(container)))
    
    print("*"*80)
    print("*"*80)
    
    container.from_file(example_fits)
    print(iff(str(container)))


# In[11]:


# if __name__ == "__main__":
#     # Testing additional functionality
    
#     out_str = """Iteration : 12    Chi2nu: 3.571e-01     dChi2/Chi2: -3.21e-08   alamda: 1e+02
#  sersic    : (  [62.90],  [62.90])  14.11     13.75    0.30    0.63    60.82
#  sky       : [ 63.00,  63.00]  1130.51  -4.92e-02  1.00e-02
#  sersic    : (  [62.90],  [62.90])  14.19     12.43    1.45    0.62   100.55
#    power   :     [0.00]   [23.51]  219.64     -0.16     ---  -44.95   -15.65
#    fourier : (1:  0.06,   -6.67)   (3:  0.05,    0.18)
# COUNTDOWN = 0"""
    
#     bulge1, disk1, arms1, fourier1, sky1 = update_components(out_str, 
#                                                              bulge, 
#                                                              disk, 
#                                                              arms, 
#                                                              fourier, 
#                                                              sky)
#     # Checking difference
#     # Use subtraction! When it's done
#     print(bulge)
#     print(disk)
#     print(arms)
#     print(fourier)
#     print(sky)
#     print("="*80)
#     print(bulge1)
#     print(disk1)
#     print(arms1)
#     print(fourier1)
#     print(sky1)
    


# In[12]:


# Testing extraction into FeedmeContainer attributes
if __name__ == "__main__":
    example_feedme = FeedmeContainer(path_to_feedme = "somewhere/out_there", 
                                     header         = header, 
                                     bulge          = bulge, 
                                     disk           = disk, 
                                     arms           = arms, 
                                     fourier        = fourier, 
                                     sky            = sky)
    feedme_components = example_feedme.extract_components()
    #print(feedme_components.to_list())
    #print()
    _ = [print("Key:", k) for k in example_feedme.to_dict().keys()]
    print()
    print(iff(str(example_feedme)))


# In[13]:


# Testing OutputContainer
if __name__ == "__main__":
    
    good_example = """Iteration : 6     Chi2nu: 3.205e-01     dChi2/Chi2: -2.24e-08   alamda: 1e+02
     sersic    : (  [67.38],  [67.77])  13.19     15.54    0.34    0.62   -19.23
     sersic    : (  [67.38],  [67.77])  14.58      8.89    1.56    0.68    30.23
       power   :     [0.00]   [22.01]   78.86     -2.39     ---   39.74    22.52
       fourier : (1:  0.15,   46.48)   (3:  0.10,  -33.06)
     sky       : [ 67.00,  68.00]  1133.43  1.16e-02  -1.35e-02
    COUNTDOWN = 6

    Iteration : 7     Chi2nu: 3.205e-01     dChi2/Chi2: -2.24e-08   alamda: 1e+03
     sersic    : (  [67.38],  [67.77])  13.19     15.54    0.34    0.62   -19.23
     sersic    : (  [67.38],  [67.77])  14.58      8.89    1.56    0.68    30.23
       power   :     [0.00]   [22.01]   78.86     -2.39     ---   39.74    22.52
       fourier : (1:  0.15,   46.48)   (3:  0.10,  -33.06)
     sky       : [ 67.00,  68.00]  1133.43  1.16e-02  -1.35e-02
    COUNTDOWN = 5

    Iteration : 8     Chi2nu: 3.205e-01     dChi2/Chi2: -2.24e-08   alamda: 1e+04
     sersic    : (  [67.38],  [67.77])  13.19     15.54    0.34    0.62   -19.23
     sersic    : (  [67.38],  [67.77])  14.58      8.89    1.56    0.68    30.23
       power   :     [0.00]   [22.01]   78.86     -2.39     ---   39.74    22.52
       fourier : (1:  0.15,   46.48)   (3:  0.10,  -33.06)
     sky       : [ 67.00,  68.00]  1133.43  1.16e-02  -1.35e-02
    COUNTDOWN = 4

    Iteration : 9     Chi2nu: 3.205e-01     dChi2/Chi2: -3.33e-08   alamda: 1e+03
     sersic    : (  [67.38],  [67.77])  13.19     15.54    0.34    0.62   -19.23
     sersic    : (  [67.38],  [67.77])  14.58      8.89    1.56    0.68    30.23
       power   :     [0.00]   [22.01]   78.86     -2.39     ---   39.74    22.52
       fourier : (1:  0.15,   46.48)   (3:  0.10,  -33.06)
     sky       : [ 67.00,  68.00]  1133.43  1.16e-02  -1.35e-02
    COUNTDOWN = 3

    Iteration : 10    Chi2nu: 3.205e-01     dChi2/Chi2: -3.33e-08   alamda: 1e+04
     sersic    : (  [67.38],  [67.77])  13.19     15.54    0.34    0.62   -19.23
     sersic    : (  [67.38],  [67.77])  14.58      8.89    1.56    0.68    30.23
       power   :     [0.00]   [22.01]   78.86     -2.39     ---   39.74    22.52
       fourier : (1:  0.15,   46.48)   (3:  0.10,  -33.06)
     sky       : [ 67.00,  68.00]  1133.43  1.16e-02  -1.35e-02
    COUNTDOWN = 2

    Iteration : 11    Chi2nu: 3.205e-01     dChi2/Chi2: -8.01e-08   alamda: 1e+03
     sersic    : (  [67.38],  [67.77])  13.19     15.54    0.34    0.62   -19.23
     sersic    : (  [67.38],  [67.77])  14.58      8.89    1.56    0.68    30.23
       power   :     [0.00]   [22.01]   78.86     -2.39     ---   39.74    22.52
       fourier : (1:  0.15,   46.48)   (3:  0.10,  -33.06)
     sky       : [ 67.00,  68.00]  1133.43  1.16e-02  -1.35e-02
    COUNTDOWN = 1

    Iteration : 12    Chi2nu: 3.205e-01     dChi2/Chi2: 1.75e-08    alamda: 1e+04
     sersic    : (  [67.38],  [67.77])  13.19     15.54    0.34    0.62   -19.23
     sersic    : (  [67.38],  [67.77])  14.58      8.89    1.56    0.68    30.23
       power   :     [0.00]   [22.01]   78.86     -2.39     ---   39.74    22.52
       fourier : (1:  0.15,   46.48)   (3:  0.10,  -33.06)
     sky       : [ 67.00,  68.00]  1133.43  1.16e-02  -1.35e-02
    COUNTDOWN = 0


    Fit summary is now being saved into `fit.log'.

    """

    bad_example = "\n".join(good_example.split("\n")[:-3] + ["...now exiting to system...\n"])
    
    dummy_obj = subprocess.CompletedProcess("", 0)
    dummy_obj.stdout = bad_example
    
    print("Checking the bad example")
    
    bad_output  = OutputContainer(dummy_obj)
    print(bad_output.bulge)
    print(bad_output.disk)
    print(bad_output.arms)
    print(bad_output.fourier)
    print(bad_output.sky)
    print(f"Did the bad example succeed? {bad_output.success}")
    
    print("*"*80)
    print("*"*80)
    
    print("And now checking the 'good' example (these should all be updated from the default values)\n")
    dummy_obj.stdout = good_example
    print("(this should produce failure text since we didn't store text)")
    good_output = OutputContainer(dummy_obj)
    print(good_output)
    print("\nNow it should succeed... Re-printing output text.\n")
    good_output = OutputContainer(dummy_obj, store_text = True)
    print(good_output)
    print()
    
    print(good_output.bulge)
    print(good_output.disk)
    print(good_output.arms)
    print(good_output.fourier)
    print(good_output.sky)
    print(f"Did the good example succeed? {good_output.success}")
    
    print("*"*80)
    print("*"*80)
    print("Testing extraction into ComponentContainer...")
    _ = [print(str(comp)) for comp in good_output.extract_components().to_list()]
    
    output_filename = pj(TEST_OUTPUT_DIR, f"{base_out}_OutputContainer.in")
    good_output.to_file(filename = output_filename)
    #good_output.header.to_file(output_filename, good_output.bulge, good_output.disk, good_output.arms, good_output.fourier, good_output.sky)


# In[14]:


if __name__ == "__main__":
    export_to_py("Containers", pj(_MODULE_DIR, "Classes", "Containers"))

