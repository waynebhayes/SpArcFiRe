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

#from Classes.Parameters import *
from Classes.Components import *
# Relies on relative file hierarchy
from Functions.helper_functions import *


# In[4]:


class ComponentContainer:
    def __init__(self, **kwargs):
        # TODO: When split from SpArcFiRe, default to just a single Sersic
        self.bulge   = kwargs.get("bulge", Sersic(1))
        self.disk    = kwargs.get("disk", Sersic(2))
        self.arms    = kwargs.get("arms", Power(2))
        self.fourier = kwargs.get("fourier", Fourier(2))
        self.sky     = kwargs.get("sky", Sky(3))
        
        # Dict
        # Mutability strikes again
        self._components = deepcopy(vars(self))        
        self._components.update(kwargs)
        
# ==========================================================================================================

    @property
    def components(self):
        return self._components
    
    @components.setter
    def components(self, new_dict):
        self.check_component_types(new_dict)
        self._components = deepcopy(new_dict)
        
# ==========================================================================================================

    def check_component_types(self):
        for k, comp in self.components.items():
            assert isinstance(comp, GalfitComponent), f"The component fed into the ComponentContainer, {k}, is not a valid type."
        
# ==========================================================================================================
    
    def to_tuple(self):
        #self.header,
        # return (self.bulge,
        #         self.disk,
        #         self.arms,
        #         self.fourier,
        #         self.sky
        #        )
        return tuple(self.components.values())
    
    def to_list(self):
        #self.header,
        # return [self.bulge,
        #         self.disk,
        #         self.arms,
        #         self.fourier,
        #         self.sky
        #        ]
        return list(self.components.values())
    
# ==========================================================================================================
    
    def to_pandas(self):
        
        return pd.concat(
            [
                comp.to_pandas().reset_index() 
                for comp in self.to_list()
                if comp.component_type != "header"
            ], 
            axis = 1).drop(columns = ["index"])
    
    
    #TODO
    def from_pandas(self, input_df):
        pass
    
# ==========================================================================================================

    # def update_components(self, **kwargs):
    #     for key, value in kwargs.items():
    #         setattr(self, key, value)
       
    # This is more for the daughter classes
    #def extract_components(self):
    #    return self.components

# ==========================================================================================================

    def __str__(self):
        # Skipping the output of power and fourier if 'skipped'
        # i.e. don't exist in this current implementation
        out_str = "\n".join(
            str(comp) for comp in self.to_list()
             # if comp.param_values.get("skip", 0) != 1 or
             #    comp.component_type not in ("power", "fourier")
            )
        return out_str
    
    def __repr__(self):
        # Skipping the output of power and fourier if 'skipped'
        # i.e. don't exist in this current implementation
        out_str = "\n".join(
            repr(comp) for comp in self.to_list()
             # if comp.param_values.get("skip", 0) != 1 or
             #    comp.component_type not in ("power", "fourier")
            )
        return out_str


# In[5]:


class FeedmeContainer(ComponentContainer):
    def __init__(self, **kwargs):
        
        # Use vars to our advantage and declare this first
        # If it is then specified in kwargs, it will overwrite this one
        self.header = GalfitHeader()
        path_to_feedme = kwargs.pop("path_to_feedme","")
        
        # The path to the feedme that *generated* the components
        ComponentContainer.__init__(self, **kwargs)
        
        self.path_to_feedme = path_to_feedme

# ==========================================================================================================

    # def to_dict(self):
    #     return self.components.update({"header" : self.header})
    
    # def to_list(self):
    #     return [self.header] + self.components
    
# ==========================================================================================================
    
#     def __str__(self):
#         out_str = f"{str(self.header)}\n" + \
#                     "\n".join(
#                         str(comp) for comp in ComponentContainer.to_list(self)
#                         # if comp.param_values.get("skip", 0) != 1 or
#                         #    comp.component_type not in ("power", "fourier")
#                     )
#         return out_str
    
#     def __repr__(self):
#         out_str = f"{str(self.header)}\n" + \
#                     "\n".join(
#                         repr(comp) for comp in ComponentContainer.to_list(self)
#                         # if comp.param_values.get("skip", 0) != 1 or
#                         #    comp.component_type not in ("power", "fourier")
#                     )
#         return out_str
    
# ==========================================================================================================
    
    # def to_pandas(self):
    #     return pd.concat([comp.to_pandas().reset_index() 
    #                       for comp in ComponentContainer.to_list(self)
    #                       if not isinstance(comp, GalfitHeader)
    #                      ], axis = 1).drop(columns = ["index"])

# ==========================================================================================================

    def to_file(self, *args, filename = ""):
        if not filename:
            filename = self.path_to_feedme
            
        if args:
            self.header.to_file(filename, *args)
        else:
            self.header.to_file(filename, *self.to_list())
            
# ==========================================================================================================

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
            
            # Header does not fit the COMP paradigm used below so we explicitly account for it here
            # Reinclude magzpt since it the end is exclusive
            header_keys = input_keys[input_keys.index("INITFILE"):input_keys.index("MAGZPT")] + \
                          ["MAGZPT"]
                
            header_dict = {k:v for k,v in input_dict.items() if k in header_keys}
            
            component_list[0].from_file_helper_dict(header_dict)
            
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
                    # For when default components are used
                    # TODO: when generalized, take this out
                    if "skip" in component.parameters.keys():
                        component.skip.value = 1 #add_skip(skip_val = 1)
                        
                    component_list_num += 1
                    component           = component_list[component_list_num]
                    count              += 1
                
                component_list_num += 1
                
                # Also assume only a single rotation function per component... which I think is fair
                if not rotation_func:
                    send_to_helper = dict(input_list[i_begin : i_end])
                    component.from_file_helper_dict(send_to_helper)
                    
                else:
                    r_start = rotation_func[0]
                    send_to_helper = dict(chunk[:r_start])
                    component.from_file_helper_dict(send_to_helper)
                    
                    component = component_list[component_list_num] #[c for c in component_list if c.component_type == component_name][0]
                    component_list_num += 1
                    
                    if not fourier_modes:
                        send_to_helper = dict(chunk[r_start:])
                        component.from_file_helper_dict(send_to_helper)
                        
                    else:
                        f_start = fourier_modes[0]
                        send_to_helper = dict(chunk[r_start : f_start])
                        component.from_file_helper_dict(send_to_helper)

                        component = component_list[component_list_num] #[c for c in component_list if c.component_type == "fourier"][0]
                        component_list_num += 1
                        
                        send_to_helper = dict(chunk[f_start:])
                        component.from_file_helper_dict(send_to_helper)
        
        def from_text(self, input_file_obj):
            
            input_file = [line.rstrip("\n") for line in input_file_obj.readlines()]
            
            # Header is guaranteed to come first***
            component_list = self.to_list()
            
            header_begin_end_idx = [i + 1 
                                    for i, line in enumerate(input_file) 
                                    if line.startswith("# IMAGE and GALFIT CONTROL PARAMETERS")
                                    or line.startswith("P)")
                                   ]
            
            component_list[0].from_file_helper_list(
                input_file[header_begin_end_idx[0]:header_begin_end_idx[1]]
            )
            
            component_list_num = 1
            #send_to_helper = []
            #store = False
            #component_exists = False
            
            component_idx_nums = [(i + 1, line[-1]) 
                                  for i, line in enumerate(input_file) 
                                  if line.startswith("# Component number")
                                  #or line.startswith("# IMAGE and GALFIT CONTROL PARAMETERS")
                                 ]
            
            for idx, (component_begin, component_num) in enumerate(component_idx_nums):
                # White space *before* component number or end of file
                component_end = [i for i,line in enumerate(input_file) 
                                 if line.strip().startswith("="*10)
                                 #or line.strip().startswith("# INITIAL FITTING PARAMETERS")
                                 ][-1] - 1
                
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
                    component_list_num += 1
                    component = component_list[component_list_num]
                    count += 1
                    
                component_list_num += 1
                
                 # Also assume only a single rotation function per component... which I think is fair
                if not rotation_func:
                    component.from_file_helper_list(chunk)
                    
                else:
                    r_start = rotation_func[0]
                    component.from_file_helper_list(chunk[:r_start - 1])
                    
                    component = component_list[component_list_num] #[c for c in component_list if c.component_type == component_name][0]
                    component_list_num += 1
                    
                    if not fourier_modes:
                        component.from_file_helper_list(chunk[r_start:])
                        
                    else:
                        f_start = fourier_modes[0]
                        component.from_file_helper_list(chunk[r_start : f_start])

                        component = component_list[component_list_num] #[c for c in component_list if c.component_type == "fourier"][0]
                        component_list_num += 1
                        
                        component.from_file_helper_list(chunk[f_start:])

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
        
        store_text = kwargs.pop("store_text", False)
        
        FeedmeContainer.__init__(self, **kwargs)
        
        galfit_out_text    = galfit_out_obj.stdout        
        galfit_err_text    = galfit_out_obj.stderr
        galfit_return_code = galfit_out_obj.returncode
        
        # Default to this so it doesn't break if no text is fed in
        self.success = False
        
        # 0 is success
        if not galfit_return_code:
            self.success = True
            
            # Pop previous galfit out if stored on the assumption that we don't need it anymore
            # even if we decide not to store the next iteration
            # This way we can leave the call explicit in update_components and not have to worry
            # about using an old fit
            kwargs.pop("galfit_out_text", None)
            
        elif galfit_return_code < 1:
            # per subprocess documentation
            # A negative value -N indicates that the child was terminated by signal N (POSIX only).
            print(f"GALFIT was terminated by signal {galfit_return_code}")
            print(f"{galfit_err_text}")
        
        # else:
        #     print(f"GALFIT failed!")

        # For reading from galfit stdout to update classes
        def update_components(self, galfit_out_text, **kwargs) -> None: #, bulge, disk, arms, fourier, sky):
            
            #if not kwargs.get("galfit_out_text"):
            #    raise("Cannot update components, no output text provided.")
                
            last_it = galfit_out_text.split("Iteration")[-1]

            s_count = 0
            by_line = last_it.splitlines()
            for line in by_line:
                if line.strip().startswith("sersic"):
                    if kwargs.get("sersic_order"):
                        comp = eval("self." + kwargs.get("sersic_order")[s_count])
                        s_count += 1
                        
                    elif s_count == 0:
                        comp = self.bulge
                        s_count += 1 
                        
                    elif s_count == 1:
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
            update_components(self, galfit_out_text, **kwargs)
            
        if store_text:
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
    from RegTest.RegTest import *


# In[8]:


if __name__ == "__main__":
    # Testing basic functionality
    
    container = ComponentContainer()
    print(container)
    container_df = container.to_pandas()
    print()
    print(container_df)


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
    
    container = new_container()
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
    
    container = new_container()
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


# In[13]:


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
    
    #feedme_components = example_feedme.components
    #print(feedme_components.to_list())
    #print()
    _ = [print("Key:", k) for k in example_feedme.components.keys()]
    print()
    print(iff(str(example_feedme)))


# In[14]:


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
    _ = [print(str(comp)) for comp in good_output.to_list()]
    
    output_filename = pj(TEST_OUTPUT_DIR, f"{base_out}_OutputContainer.in")
    good_output.to_file(filename = output_filename)
    #good_output.header.to_file(output_filename, good_output.bulge, good_output.disk, good_output.arms, good_output.fourier, good_output.sky)


# In[15]:


if __name__ == "__main__":
    export_to_py("Containers", pj(_MODULE_DIR, "Classes", "Containers"))

