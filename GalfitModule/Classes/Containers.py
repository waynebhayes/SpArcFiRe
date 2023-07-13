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
        if __name__ == "__main__":
            print("SPARCFIRE_HOME is not set. Please run 'setup.bash' inside SpArcFiRe directory if not done so already.")
            print("Checking the current directory for GalfitModule, otherwise quitting.")
            
        _MODULE_DIR = pj(os.getcwd(), "GalfitModule")
        
        if not exists(_MODULE_DIR):
            raise Exception("Could not find GalfitModule!")

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
        # Skipping the output of power and fourier if 'skipped'
        # i.e. don't exist in this current implementation
        out_str = "\n".join(
            str(comp) for comp in ComponentContainer.to_list(self)
             if comp.param_values.get("skip", 0) != 1 or
                comp.component_type not in ("power", "fourier")
            )
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
        out_str = f"{str(self.header)}\n" + \
                    "\n".join(
                        str(comp) for comp in ComponentContainer.to_list(self)
                        if comp.param_values.get("skip", 0) != 1 or
                           comp.component_type not in ("power", "fourier")
                    )
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
            input_list = list(input_dict.items())
            
            # For now assume this will be the last one!
            try:
                final_idx  = input_keys.index("FLAGS")
            except ValueError:
                final_idx  = -1
            
            component_list = self.to_list()
            component_list_num = 1
            
            send_to_helper = {}
            
            component_nums = [(i, k) for i, k in enumerate(input_keys) if k.startswith("COMP")]
            
            for idx_component_nums, (i_begin, component_key) in enumerate(component_nums):
                i_end = final_idx
                
                if idx_component_nums + 1 < len(component_nums):
                    i_end = component_nums[idx_component_nums + 1][0]
                        
                chunk = input_list[i_begin:i_end]
                
                rotation_func  = [i for i, k in enumerate(chunk) if k[0].endswith("ROTF")]
                # Assume at least F1
                fourier_modes  = [i for i, k in enumerate(chunk) if k[0].endswith("F1")]

                component_name = input_dict[component_key.strip()].strip()
                component = component_list[component_list_num] #[c for c in component_list if c.component_type == component_name][0]
                
                # Aligning ourselves correctly
                count = 0
                while component.component_type != component_name or count == 100:
                    #print(component.component_type)
                    component.add_skip(skip_val = 1)
                    component_list_num += 1
                    component = component_list[component_list_num]
                    count += 1
                
                component_list_num += 1
                
                # Also assume only a single rotation function per component... which I think is fair
                if not rotation_func:
                    send_to_helper = dict(input_list[i_begin : i_end])
                    component.from_file_helper(send_to_helper)
                    
                else:
                    r_start = rotation_func[0]
                    send_to_helper = dict(chunk[:r_start])
                    component.from_file_helper(send_to_helper)
                    
                    component = component_list[component_list_num] #[c for c in component_list if c.component_type == component_name][0]
                    component_list_num += 1
                    
                    if not fourier_modes:
                        send_to_helper = dict(chunk[r_start:])
                        component.from_file_helper(send_to_helper)
                        
                    else:
                        f_start = fourier_modes[0]
                        send_to_helper = dict(chunk[r_start : f_start])
                        component.from_file_helper(send_to_helper)

                        component = component_list[component_list_num] #[c for c in component_list if c.component_type == "fourier"][0]
                        component_list_num += 1
                        
                        send_to_helper = dict(chunk[f_start:])
                        component.from_file_helper(send_to_helper)
                    
            
# DEPRECATED
# New implementation flips the logic, instead of starting from bulge + disk + spiral
# we start from what's actually in the file and go from there (working around my defaults that is)

            #components_to_pop = []
#             for idx, (key, value) in enumerate(input_dict.items()):
                
#                 if component_list_num == len(component_list): break
                
#                 component = component_list[component_list_num]
#                 name = component.component_type
#                 #print(name, component_list_num, component_dict_num)
#                 # For skipping Power and Fourier if already specified (save some iterations)
#                 #if component.param_values.get("skip", 0) == 1 and component.component_type in ("power", "fourier"):
#                     #continue
                    
#                 try:
#                     component_idx_start = input_keys.index(component.start_dict)
#                     component_idx_end   = input_keys.index(component.end_dict)
                    
#                 except ValueError as ve:
#                     # TODO: Don't default to spiral implementation! This is a hotfix for when no disk/spiral/etc.
#                     # i.e. in base container class, allow for an n component fit
                    
#                     # Sky defaults to component 3, this is why it fails for bulge only fits
#                     # or outputs the wrong component number
#                     if name == "sky":
#                         self.sky = Sky(2)
#                         component_list = self.to_list()
#                         send_to_helper = {}
#                         continue
                        
#                     elif "is not in list" in str(ve):
#                         component.add_skip(skip_val = 1)
#                         #components_to_pop.append(component.component_type)
                        
#                         #Any bending modes/rotations/etc. go here
#                         component_list_num += 1
                        
#                         # Assume if it fails it fails for the whole component
#                         send_to_helper = {}
#                         continue

#                     # End will *always* (header excluded) be #_param                
#                     component_end = [k for k in input_keys if k.endswith(component.end_dict[2:])][0]
#                     component_num = component_end[0]

#                     if component_num.isnumeric():
#                         # For COMP_#
#                         if component.start_dict[:-2] == "COMP":
#                             component_idx_start = input_keys.index(f"{component.start_dict[:-2]}_{component_end[0]}")
#                         else:
#                             component_idx_start = input_keys.index(f"{component_num}_{component.start_dict[2:]}")
                            
#                         component_idx_end   = input_keys.index(component_end)

#                     else:
#                         print(f"Can't find start/end of {component.component_type} segment.")
#                         print(f"Check the filename or start/end_dict variables.")
#                         print(f"Filename: {filename}")
#                         print(f"Start/End: {component.start_dict}/{component.end_dict}")
#                         raise ValueError(ve)

#                 if component_idx_start <= idx <= component_idx_end:
#                     send_to_helper[key] = value
                
#                 if idx == component_idx_end:
#                     print(send_to_helper)
#                     component.from_file_helper(send_to_helper)
                    
#                     send_to_helper = {}
#                     component_list_num += 1
#                     component_dict_num += 1
                    
            return #components_to_pop
        
        def from_text(self, input_file_obj):
            
            input_file = [line.rstrip("\n") for line in input_file_obj.readlines()]
            
            component_list = self.to_list()
            component_list_num = 0
            #send_to_helper = []
            #store = False
            #component_exists = False
            
            component_idx_nums = [(i + 1, line[-1]) for i, line in enumerate(input_file) if line.startswith("# Component number")]
            
            for idx, (component_begin, component_num) in enumerate(component_idx_nums):
                # White space *before* component number or end of file
                component_end = [i for i,line in enumerate(input_file) if line.strip().startswith("="*10)][-1] - 1
                
                if idx + 1 < len(component_idx_nums):
                    component_end = component_idx_nums[idx + 1][0] - 2
                
                chunk = input_file[component_begin : component_end]
                
                rotation_func  = [i for i, k in enumerate(chunk) if k.startswith("R0")]
                # Assume at least F1
                fourier_modes  = [i for i, k in enumerate(chunk) if k.startswith("F1")]
                
                component_name = input_file[component_begin].replace(" ", "").lstrip("0)").split("#")[0]
                component = component_list[component_list_num]
                
                # Aligning ourselves correctly
                count = 0
                while component.component_type != component_name or count == 100:
                    component.add_skip(skip_val = 1)
                    component_list_num += 1
                    component = component_list[component_list_num]
                    count += 1
                    
                component_list_num += 1
                
                 # Also assume only a single rotation function per component... which I think is fair
                if not rotation_func:
                    component.from_file_helper(chunk)
                    
                else:
                    r_start = rotation_func[0]
                    component.from_file_helper(chunk[:r_start - 1])
                    
                    component = component_list[component_list_num] #[c for c in component_list if c.component_type == component_name][0]
                    component_list_num += 1
                    
                    if not fourier_modes:
                        component.from_file_helper(chunk[r_start:])
                        
                    else:
                        f_start = fourier_modes[0]
                        component.from_file_helper(chunk[r_start : f_start])

                        component = component_list[component_list_num] #[c for c in component_list if c.component_type == "fourier"][0]
                        component_list_num += 1
                        
                        component.from_file_helper(chunk[f_start:])

# DEPRECATED
#             for line in input_file:
#                 component = component_list[component_num]
                
#                 if line.strip().startswith(component.start_text):
#                     store = True
#                     #component_exists = True
                    
#                 if store:
#                     send_to_helper.append(line)
#                     #component_exists = True
                    
#                 if line.strip().startswith(component.end_text):
#                     #component_exists = True
#                     store = False
#                     component.from_file_helper(send_to_helper)
#                     send_to_helper = []
#                     component_num += 1
                    
#                     if component_num == len(component_list): break
                    
                    #component_exists = False
                    #continue
                    
                # Basically if the if conditions are never met (which they should be)
                # if not component_exists:
                #     component.add_skip(skip_val = 1)
                #     continue
                    
                #components_to_pop.append(component.component_type)
        
            return #components_to_pop
        
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
    
    def new_container():
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
        return container

    container = new_container()

    print()
    print(container.to_pandas())
    
    #print(container)
    container.path_to_feedme = pj(TEST_OUTPUT_DIR, f"{base_out}_FeedmeContainer.in")
    
    container.to_file()


# In[10]:


# Testing FeedmeContainer from_file
if __name__ == "__main__":
    
    container = new_container()

    example_feedme = pj(TEST_DATA_DIR, "test-out", "1237667911674691747", "1237667911674691747.in")
    example_fits   = pj(TEST_DATA_DIR, "test-out", "1237667911674691747", "1237667911674691747_galfit_out.fits")
    
    print("These are feedme -> output")
    print("ignoring filepaths for reg tests...\n")
    
    container.from_file(example_feedme)
    print(iff(str(container)))
    
    print("*"*80)
    print("*"*80)
    
    container.from_file(example_fits)
    print(iff(str(container)))


# In[11]:


# Testing FeedmeContainer from_file with just bulge
if __name__ == "__main__":
    
    container = new_container()

    # This galaxy does not use the Power or Fourier functions
    example_feedme = pj(TEST_DATA_DIR, "test-out", "1237668589728366770", "1237668589728366770_galfit.01")
    example_fits   = pj(TEST_DATA_DIR, "test-out", "1237668589728366770", "1237668589728366770_galfit_out.fits")
    
    print("These are feedme -> output")
    print("ignoring filepaths for reg tests...\n")
    
    container.from_file(example_feedme)
    print(iff(str(container)))
    
    print("*"*80)
    print("*"*80)
    
    container.from_file(example_fits)
    print(iff(str(container)))


# In[12]:


# Testing FeedmeContainer from_file with no arms
if __name__ == "__main__":
    
    container = new_container()
    
    # This galaxy does not use the Power or Fourier functions
    example_feedme = pj(TEST_DATA_DIR, "test-out", "1237667912741355660", "1237667912741355660.in")
    example_fits   = pj(TEST_DATA_DIR, "test-out", "1237667912741355660", "1237667912741355660_galfit_out.fits")
    
    print("These are feedme -> output")
    print("ignoring filepaths for reg tests...\n")
    
    container.from_file(example_feedme)
    print(iff(str(container)))
    
    print("*"*80)
    print("*"*80)
    
    container.from_file(example_fits)
    print(iff(str(container)))


# In[15]:


# Testing extraction into FeedmeContainer attributes
if __name__ == "__main__":
    container = new_container()
    example_feedme = FeedmeContainer(path_to_feedme = "somewhere/out_there", 
                                     header         = container.header, 
                                     bulge          = container.bulge, 
                                     disk           = container.disk, 
                                     arms           = container.arms, 
                                     fourier        = container.fourier, 
                                     sky            = container.sky
                                    )
    
    feedme_components = example_feedme.extract_components()
    #print(feedme_components.to_list())
    #print()
    _ = [print("Key:", k) for k in example_feedme.to_dict().keys()]
    print()
    print(iff(str(example_feedme)))


# In[16]:


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


# In[ ]:


if __name__ == "__main__":
    export_to_py("Containers", pj(_MODULE_DIR, "Classes", "Containers"))

