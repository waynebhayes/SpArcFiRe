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
import warnings

import pandas as pd
import numpy as np


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

from Functions.helper_functions import *
from Classes.Parameters import *


# In[20]:


# TODO: MAKE ITERABLE via generator
# def reverse(data):
    # for index in range(len(data)-1, -1, -1):
    #     yield data[index]
class GalfitComponent:
    def __init__(self,
                 parameters       = {}, #load_default_parameters(),
                 component_type   = "",
                 component_name   = "",
                 component_number = 0,
                 param_prefix     = " ",
                 **kwargs
                ):
        
        if parameters:
            # This should ONLY occur here since the other classes will
            # give a default to use in case parameters is a dict of 
            # keys and values
            assert all(isinstance(parameter, BaseParameter) for parameter in parameters.values()),\
                   f"Parameters in dictionary fed into {type(self).__name__} must be a parameter class.\n" \
                   f"We recommend feeding those in as key word arguments to the instantiation instead."

            for k,v in kwargs:
                if k in parameters:
                    parameters[k].value = v
                
        self.parameters       = parameters
        self.component_type   = component_type
        # If you decide to name your component *besides and/or including* the
        # variable itself
        self.component_name   = component_name
        self.component_number = component_number
        self.param_prefix     = param_prefix       
        
        # For reading from file
        self.start_dict = f"COMP_{self.component_number}"
        self.end_dict   = f"COMP_{self.component_number + 1}"
        
        #self.start_dict_value = self.component_type
        
        # This should work but clarity is probably best here
        # self.start_text = str(ComponentType(
        #     self.component_type, 
        #     parameter_number = f"{self.param_prefix}0", 
        #     component_number = component_number
        # ))
        self.start_text = f"# Component number: {self.component_number}"
        self.end_text   = f"{self.param_prefix.strip()}10"    
        
        self._section_sep = "="*80
        
# ==========================================================================================================

    @staticmethod
    def component_get_set(component_dict = None, component_type = ""):
        
        exec_dict = {}
        
        # TODO: Set up logging
        if not component_dict:
            exec_dict = load_default_parameters()
            
            warnings.warn(f"Component type: {component_type} not properly specified " \
                      f"upon initializing {GalfitComponent.__name__}.\n" \
                      f"Will initalize all component types currently set in 'load_default_parameters()' " \
                      f"from 'Parameters' module." \
                      "\n".join(exec_dict.keys())
                     )
        else:
            # Give it an empty key so that the following loop works
            exec_dict[""] = component_dict
        
        full_str = ""
            
        return "\n".join(
            [
                generate_get_set(
                    # inner dict looks like this: 
                    # "position" : parameters["position"]
                    {k : f"parameters[\"{k}\"]" 
                     for k in e_dict.keys()
                    }
                )
                for e_dict in exec_dict.values()
            ]
        )        
        
# ==========================================================================================================

    # def update_values(self, **kwargs):
    #     for key, value in kwargs.items():
    #         setattr(self, key, value)
    
    # Function must return a dict otherwise I can't use 
    # my fancy decorator
    def update_parameter_dict_with_dict(func):
        
        def func_wrapper(self, input_dict): #*args, **kwargs):
            input_dict = func(self, input_dict)
            
            assert isinstance(input_dict, dict), \
            f"update_parameter_dict_with_dict decorator from GalfitComponent " \
            f"improperly used for a {type(self).__name__} component. " \
            f"No dictionary was given."

            for pname, pval in input_dict.items():
                if pname in self.parameters:
                    self.parameters[pname].value = pval
                    
        return func_wrapper
    
    # This is by far the more dangerous way
    def update_parameter_dict_with_list(func):
        def func_wrapper(self, input_list): #*args, **kwargs):
            input_list = func(self, input_list)
            
            assert isinstance(input_list, list), \
            f"update_parameter_dict_with_list decorator from GalfitComponent " \
            f"improperly used for a {type(self).__name__} component. " \
            f"No list was given."

            # Must check for position to see if we need to add or subtract a list element
            # since position is multivalued in the log line
            parameter_names = self.parameters.keys()
            if "position" in parameter_names:
                input_list[1] = (input_list[0], input_list[1])
                input_list[0] = self.component_type
            else:
                input_list    = [self.component_type] + input_list
                
            if "skip" in parameter_names:
                input_list.append("0")
            
            assert len(input_list) == len(self.parameters), \
            f"update_parameter_dict_with_list decorator from GalfitComponent " \
            f"improperly used for a {type(self).__name__} component. " \
            f"List is not the same length as the dictionary of parameters."
            
            for pname, pval in zip(parameter_names, input_list):
                self.parameters[pname].value = pval
                    
        return func_wrapper
          
# ==========================================================================================================

    # TODO: Add comparison function via subtract(?)
    def __sub__(self):
        pass
    
# ==========================================================================================================
  
    def __str__(self):
        return f"# Component number: {self.component_number}\n" + \
                "\n".join([str(v) for v in self.parameters.values()]) + \
                "\n"

# ==========================================================================================================
    
    @update_parameter_dict_with_list
    def update_from_log(self, in_line:str):
        # Used to update from stdout i.e. what we see in fit.log
        # Requires outside function to loop through log file
        
        # NOTE: These necessarily round to two digits because that's
        # all Galfit outputs to stdout
        #print("Did this get properly overwritten from base in GalfitComponent?")
        
        # Extra empty slot to account for component name in parameter dictionary
        return [
            i.strip("[*,]() ") for i in in_line.split() 
            if any(map(str.isdigit, i))
        ]
        
# ==========================================================================================================

    # TODO: Update to work with series 
    def from_pandas(self, input_df):
        param_names  = [n.split(f"_{self.component_type}")[0] for n in input_df.columns]
        param_values = input_df.iloc[0].values.astype(float).round(4)
        new_param_dict = dict(zip(param_names, param_values))
        
        pos = "position"
        if pos in self.parameters.keys():
            new_param_dict[pos] = (new_param_dict[f"{pos}_x"], new_param_dict[f"{pos}_y"])
            new_param_dict.pop(f"{pos}_x")
            new_param_dict.pop(f"{pos}_y")
        
        # No graceful way to do this...
        # TODO: Can this be used for bending modes as well?
        if self.component_type in ("fourier"):
            f_modes = set([pn.split("_")[0] for pn in param_names])
            a   = "amplitude"
            pha = "phase_angle"
            
            for mode in f_modes:
                new_param_dict[mode] = (new_param_dict[f"{mode}_{a}"], new_param_dict[f"{mode}_{pha}"])
                new_param_dict.pop(f"{mode}_{a}")
                new_param_dict.pop(f"{mode}_{pha}")
        
        for k in self.parameters.keys():
            if k.startswith("_"):
                continue
                
            self.parameters[k].value = new_param_dict[k]

# ==========================================================================================================

    def to_pandas(self):
        name = f"{self.component_type}_{self.component_number}"
        parameter_dict = deepcopy(self.parameters)
                
        for pname, pval in self.parameters.items():
            if pname.startswith("_"):
                parameter_dict.pop(pname)
                continue
            
            # Split multivalued parameters like position
            # briefly convert to NumParameter to coincide with others
            if isinstance(pval, MultiParameter):
                old_keys = pval.value._asdict()
                parameter_dict.pop(pname)
                parameter_dict.update({f"{pname}_{k}" : NumParameter(v) for k, v in old_keys.items()})
        
        parameter_dict = {f"{k}_{name}" : v.value for k, v in parameter_dict.items()}
        
        all_data = pd.DataFrame(
            parameter_dict, 
            index = [name], 
            dtype = np.float32
        )
        
        # Move skip to end for reasons
        skip_col = f"skip_{name}"
        if skip_col in all_data.columns:
            all_data.insert(len(all_data.columns) - 1, skip_col, all_data.pop(skip_col))
                
        return all_data
# ==========================================================================================================

    def from_file_helper(*args, **kwargs):
        # This function is a placeholder
        # All components will have this, it'll make it easier and less
        # confusing to debug/parse the inputs
        raise Exception("Wrong file helper function called! It's either from_file_helper[_list | _dict].")
        
# ==========================================================================================================

    def update_parameters_file_helper(self, file_dict):

        for k, v in file_dict.items():
            if k.startswith("_"):
                continue
            
            if issubclass(type(self.parameters[k]), MultiParameter):
                value     = v[:2]
                fix_value = v[2:]
                
            else:
                value     = v[0]
                fix_value = v[1]
            
            self.parameters[k].value = value
            self.parameters[k].fix   = fix_value
            
# ==========================================================================================================    

    # File dict is only for fits
    def from_file_helper_dict(self, file_in):
        # Feed in just the chunk from the main 'from_file' caller
        # This requires determining param_begin/end in that caller
        # This can handle component type but it is unnecessary
        
        assert self.component_type, f"Component type must be specified to read from file."
        
        # to be abundantly safe
        file_in = deepcopy(file_in)
        
        p_names = [
            k for k in load_default_parameters()[self.component_type].keys()
            if not k.startswith("_") and k != "position"
        ]
        
        # File dict is only for fits          
        file_in   = {k : v for k,v in file_in.items() if v != self.component_type}
        list_vals = list(file_in.values())

        #file_in   = {k:v for k,v in file_in.items() if not k.endswith("_YC") and not k.endswith("_XC")}

        file_dict = {}
        count = 0
        for k, v in file_in.items():
            if not k.endswith("_YC") and not k.endswith("_XC"):
                file_dict[p_names[count]] = v.strip("[]").split()[0], 0 if ("[" in v) and ("]" in v) else 1
                
                count += 1

        if "position" in load_default_parameters()[self.component_type].keys():
            # Assume position is always the first and second index after the component
            file_dict["position"] = [i.strip("[]") for i in list_vals[:2]] + \
                                    [0 if ("[" in i) and ("]" in i) else 1 for i in list_vals[:2]]
        
        self.update_parameters_file_helper(file_dict)
            
# ==========================================================================================================
# These used to be unified but Fourier threw a wrench in all that 
# ==========================================================================================================

    def from_file_helper_list(self, file_in):
        # Feed in just the chunk from the main 'from_file' caller
        # This requires determining param_begin/end in that caller
        # This can handle component type but it is unnecessary
        
        assert self.component_type, f"Component type must be specified to read from file."
        
        # to be abundantly safe
        file_in = deepcopy(file_in)
        
        # Excludes 0
        #p_numbers = list(self.param_numbers.keys())[1:]
        # Invert keys and values for the dict comp a few lines down
        p_numbers = {
            str(v.parameter_number) : k 
            for k, v in load_default_parameters()[self.component_type].items() 
            if not k.startswith("_")
        }
                
        file_list = file_in
        file_dict = {
            line[:line.index(")")].strip(f" {self.param_prefix}") : line[line.index(")") + 1 : line.index("#")].strip()
            for line in file_list
            if line.strip()[0] not in ("#") #, "6", "7", "8")
        } #, "Z")}
        
        # join split split gets rid of all extra spacing except for one between everything
        file_dict = {
            p_numbers[num] : " ".join(v.split()).split() 
            for num, v in file_dict.items()           
            # Sometimes GALFIT outputs empty lines so make sure the line is valid
            if str(num) in p_numbers
        }
        
        # This is for setting the fix value later
        if "skip" in file_dict:
            file_dict["skip"].append(None)
                                    
        self.update_parameters_file_helper(file_dict)

# ==========================================================================================================
    
    def from_file(self, filename):
        # This function handles grabbing and storing the values from galfit files (input and output???)
        # It's written to generally handle both and stores everything in the respective component objects

        def from_fits(self, filename = "", image_num = 2):
            try: 
                # Grabbing the filename
                #input_filename = glob_name(galaxy_path, '', filename) 
                input_file = fits.open(filename)

            except FileNotFoundError:
                print(f"Can't open to read the file, {filename}. Check name/permissions/directory.")
                return None

            except OSError as ose:
                print(f"Something went wrong! {ose}")
                return None
        
            input_in = dict(input_file[image_num].header)
            keys = list(input_in.keys())
            
            # Power and Fourier now require component number (for the component they modify) at instantiation
            # Should not affect output so here try/except is for anything else
            try:
                feed_in  = {key : value for idx, (key, value) in enumerate(input_in.items()) 
                            if keys.index(self.start_dict) < idx <= keys.index(self.end_dict)}
                
            except ValueError as ve:
                # Trying to recover...
                # End will *always* (header excluded) be #_param                
                component_end = [k for k in keys if k.endswith(self.end_dict[2:])][0]
                             
                if component_end[0].isnumeric():
                    component_start = f"{component_end[0]}_{self.start_dict[2:]}"
                    
                else:
                    print(f"Can't find start/end of {self.component_type} segment.")
                    print(f"Check the filename or start/end_dict variables.")
                    print(f"Filename: {filename}")
                    print(f"Start/End: {self.start_dict}/{self.end_dict}")
                    raise ValueError(ve)
                    
                # Mix check with value (component type) and end key because that's usually known
                # Comp type should be handled now that we include the beginning
                feed_in  = {key : value for idx, (key, value) in enumerate(input_in.items()) 
                            if keys.index(component_start) <= idx <= keys.index(component_end)}   
            
            input_file.close()
            return feed_in
        
        def from_text(self, filename = ""):
            try: 
                # Grabbing the filename
                #input_filename = glob_name(galaxy_path, '', filename) 
                input_file = open(filename,'r')

            except FileNotFoundError:
                print(f"Can't open to read the file, {filename}. Check name/permissions/directory.")
                return None

            except OSError as ose:
                print(f"Something went wrong! {ose}")
                return None

            store   = False
            feed_in = []
            for line in input_file:
                if line.strip().startswith(self.start_text):
                    store = True
                    
                if store:
                    feed_in.append(line)
                    
                if line.strip().startswith(self.end_text):
                    store = False
                    
            input_file.close()
        
            return feed_in
        
        ext = os.path.splitext(filename)[1]
        if ext == ".fits":
            feed_in = from_fits(self, filename)
            self.from_file_helper_dict(feed_in)
        else:
            feed_in = from_text(self, filename)
            self.from_file_helper_list(feed_in)
        
        #raise Exception("Something went wrong importing from text!")
        #update_parameters_file_helper(self, file_dict)
        
# ==========================================================================================================

    def to_file(self, filename, *args):
        # For skipped power and fourier
        if self.parameters.get("skip", 0) == 1 and self.component_type in ("power", "fourier"):
            return None
        
        try:
            with open(filename, "w") as f:
                f.write("\n")
                f.write(str(self))
                f.write("\n")

                # *args for writing in additional classes at the same time (save I/O)
                comp_names = [c.component_type for c in args]
                with_fourier = "fourier" in comp_names

                # Arbitrary #
                fourier_index = 1000
                if with_fourier:
                    fourier_index = comp_names.index("fourier")

                for i, component in enumerate(args):
                    # For skipped power and fourier
                    if component.parameters.get("skip",0) == 1 and component.component_type in ("power", "fourier"):
                        continue
                        
                    f.write(str(component))
                    if i != fourier_index - 1:
                        f.write("\n")

                f.write("="*80 + "\n")
                
        except FileNotFoundError:
            print(f"Can't open to write the file, {filename}. Check permissions/directory.")
            
        except OSError as ose:
            print(f"Something went wrong! {ose}")
            


# In[21]:


class Sersic(GalfitComponent):
    def __init__(self, component_number, **kwargs):
        
        # SersicParameters               
        GalfitComponent.__init__(self, 
                                 load_default_sersic_parameters(component_number = component_number), 
                                 component_type = "sersic", 
                                 component_number = component_number,
                                 **kwargs
                                )
        
        # For reading from file
        self.start_dict = f"COMP_{self.component_number}"
        self.end_dict   = f"{self.component_number}_PA"
        
        # text kept at defaults
        #self.start_text = f"# Component number: {self.component_number}"
        #self.end_text   = f"{self.param_prefix.strip()}10"
    
    # Maybe it's silly to do it this way but in the future, it should be easier
    # to implement new components and it should be safer
    exec(GalfitComponent.component_get_set(load_default_sersic_parameters()))
    


# In[6]:


class Power(GalfitComponent):
    def __init__(self, component_number, **kwargs):
        #self.component_type = "power"
        #power_parameters = load_default_power_parameters(component_number = component_number)
        
        GalfitComponent.__init__(self, 
                                 load_default_power_parameters(component_number = component_number),
                                 component_type = "power", 
                                 component_number = component_number,
                                 param_prefix = "R",
                                 **kwargs
                                )
        
        # For reading from file
        # 2_ may not always be the case but that's why I have a try except in there ;)
        self.start_dict = f"{self.component_number}_ROTF"
        self.end_dict   = f"{self.component_number}_SPA"
        
        self.start_text = f"{self.param_prefix}0) power"
        # end kept at defult
        #self.end_text   = f"{self.param_prefix.strip()}10"
    
    
    
    # dict for get set looks like this: 
    # "position" : parameters["position"]
    exec(GalfitComponent.component_get_set(load_default_power_parameters()))  
    
    # Since this one does not have a component number, it gets special treatment
    def __str__(self):
        return "\n".join(GalfitComponent.__str__(self).split("\n")[1:])
    


# In[7]:


class Fourier(GalfitComponent):
    # kwargs is a placeholder
    def __init__(self, component_number, n = {1 : (0.05, 45), 3 : (0.05, 25)}, **kwargs):
        parameters = load_default_fourier_parameters(component_number = component_number)
        if n:
            for fnum, (amplitude, phase_angle) in n.items():
                parameters[f"F{fnum}"] = FourierMode(
                    mode = str(fnum),
                    amplitude = amplitude,
                    phase_angle = phase_angle,
                    component_number = component_number
                )
         
        GalfitComponent.__init__(self, 
                                 parameters,
                                 component_type = "fourier",
                                 param_prefix = "F",
                                 component_number = component_number,
                                 **kwargs
                                )
        
        self.sort_parameters()
        
        # normal rules don't apply here
        # Still use inheritance for the other functions
        
        # TODO: FIND SOME WAY TO UPDATE THIS WHEN OBJECT IS UPDATED
        # preferably without copying and pasting things
        # TODO: These do not update via update_param_values...
        self._amplitudes   = [mode.amplitude for pname, mode in self.parameters.items() if pname != "skip"]
        self._phase_angles = [mode.phase_angle for pname, mode in self.parameters.items() if pname != "skip"]
        
        p_numbers = list(self.parameters.keys())
        # For reading from file
        self.start_dict = f"{self.component_number}_F{p_numbers[0]}"
        self.end_dict   = f"{self.component_number}_F{p_numbers[-1]}PA"
        
        self.start_text = f"F{p_numbers[0]}"
        self.end_text   = f"{self.param_prefix}{p_numbers[-1]}"
        
    exec(GalfitComponent.component_get_set(load_default_power_parameters()))
    
    exec(
        generate_get_set(
            {
                "amplitudes" : "_amplitudes",
                "phase_angles" : "_phase_angles"
            }
        )
    )
    
# ==========================================================================================================

    def __str__(self):
        return "\n".join(GalfitComponent.__str__(self).split("\n")[1:])
    
# ==========================================================================================================

    # To keep things in proper order
    def sort_parameters(self):    
        self.parameters = dict(sorted(self.parameters.items()))
        
# ==========================================================================================================
    
    def include_fn(self, n:dict):
        for fnum, values in n.items():
            self.parameters[f"{self.param_prefix}{str(fnum)}"] = FourierMode(
                mode = str(fnum),
                amplitude = values[0],
                phase_angle = values[1],
                component_number = self.component_number
            )
            
        self.sort_parameters()
        
# ==========================================================================================================
        
    def from_file_helper_dict(self, file_in):
               
        # to be abundantly safe
        file_in = deepcopy(file_in)
        
        # Excludes 0
        #p_numbers = list(self.param_numbers.keys())[1:]
        n_dict = {}
            
        for k,v in file_in.items():

            if "+/-" in v:
                v = v.split("+/-")[0]

            #ex: 1_F1 -> 1
            # 1_F3PA -> 3
            k_num = int(k.split("_")[1][1])
            
            if k_num not in n_dict:
                n_dict[k_num] = [float(v.strip("[* ] "))]
            else:
                n_dict[k_num] += [float(v.strip("[* ] "))]
                
        self.include_fn(n_dict)
        
        for k,v in file_in.items():
            # For param_fix
            if "[" in v and "]" in v:
                if "PA" in k:
                    self.parameters[f"{self.param_prefix}{k_num}"].fix_phase_angle = "0"
                else:
                    self.parameters[f"{self.param_prefix}{k_num}"].fix_amplitude   = "0"

            else:
                if "PA" in k:
                    self.parameters[f"{self.param_prefix}{k_num}"].fix_phase_angle = "1"
                else:
                    self.parameters[f"{self.param_prefix}{k_num}"].fix_amplitude   = "1"

# ==========================================================================================================

    def from_pandas(self, input_df):
        param_names  = [n.split(f"_{self.component_type}")[0] for n in input_df.columns]
        param_values = input_df.iloc[0].values.astype(float)
        new_param_dict = dict(zip(param_names, param_values))
        
        # pos = "position"
        # if pos in self.parameters.keys():
        #     new_param_dict[pos] = (new_param_dict[f"{pos}_x"], new_param_dict[f"{pos}_y"])
        #     new_param_dict.pop(f"{pos}_x")
        #     new_param_dict.pop(f"{pos}_y")
        
        # No graceful way to do this...
        f_modes = set([pn.split("_")[0] for pn in param_names])
        a   = "amplitude"
        pha = "phase_angle"
            
        for mode in f_modes:
            if mode == "skip":
                continue
                
            new_param_dict[mode] = (new_param_dict[f"{mode}_{a}"], new_param_dict[f"{mode}_{pha}"])
            new_param_dict.pop(f"{mode}_{a}")
            new_param_dict.pop(f"{mode}_{pha}")
        
        for k in self.parameters.keys():
            if k.startswith("_"):
                continue
                
            self.parameters[k].value = new_param_dict[k]

# ==========================================================================================================

    def update_from_log(self, in_line):
        # example
        # fourier : (1:  0.06,   -6.67)   (3:  0.05,    0.18)
        
        # rstrip avoids a hanging ) later
        params = in_line.lstrip("fourier : ").replace(" ", "").rstrip(")").split(")(")
        
        for i, pname in enumerate(self.parameters.keys()):
            if pname != "skip":
                self.parameters[pname].value = eval(f"({params[i].split(':')[1].replace('*', '')})")


# In[22]:


class Sky(GalfitComponent):
    def __init__(self, component_number, **kwargs):
        
        GalfitComponent.__init__(self, 
                                 load_default_sky_parameters(component_number = component_number),
                                 component_type = "sky", 
                                 component_number = component_number,
                                 **kwargs
                                )
        
        # For reading from file
        self.start_dict = f"COMP_{self.component_number}"
        self.end_dict   = f"{self.component_number}_DSDY"
        
        self.start_text = f"# Component number: {self.component_number}"
        self.end_text   = f"{self.param_prefix.strip()}3"
        
    # dict for get set looks like this: 
    # "position" : parameters["position"]
    exec(GalfitComponent.component_get_set(load_default_sky_parameters()))
    
# ==========================================================================================================

    @GalfitComponent.update_parameter_dict_with_list
    def update_from_log(self, in_line):
#         # example
#         #sky       : [ 63.00,  63.00]  1130.51  -4.92e-02  1.00e-02

        # Ignore position
        return [i.strip("[*,]") for i in in_line.split() 
                  if any(map(str.isdigit, i))
               ][2:]


# In[9]:


class GalfitHeader(GalfitComponent):
    
    def __init__(self, parameters = {}, galaxy_name = "", **kwargs):
        
        # normal rules don't apply here
        # Still use inheritance for the other functions
        if not parameters:
            parameters = load_default_header_parameters(galaxy_name = galaxy_name)
        
        # If not fully specified, will use galaxy_name as default so it's good to use
        # it as an argument even if I am specifying each individually
        GalfitComponent.__init__(self, 
                                 parameters = parameters,
                                 component_type = "header",
                                 param_prefix   = "",
                                 **kwargs
                                )
        
        # For reading from file
        self.start_dict = "INITFILE"
        self.end_dict   = "MAGZPT"
        
        self.start_text = f"A" # {self.input_image}"
        self.end_text   = f"P" #{self.optimize}"
        
        # No newlines added so the strings can be added to directly
        self.input_menu_file   = f"#  Input menu file: {kwargs.get('input_menu_file', '')}.in"
        self.extra_header_info = f"#  {kwargs.get('extra_header_info', '')}"      
        
        # Don't mess with this tabbing
        self.post_header = """
# INITIAL FITTING PARAMETERS
#
#   For component type, the allowed functions are: 
#       sersic, expdisk, edgedisk, devauc, king, nuker, psf, 
#       gaussian, moffat, ferrer, and sky. 
#  
#   Hidden parameters will only appear when they're specified:
#       Bn (n=integer, Bending Modes).
#       C0 (diskyness/boxyness), 
#       Fn (n=integer, Azimuthal Fourier Modes).
#       R0-R10 (coordinate rotation, for creating spiral structures).
#       To, Ti, T0-T10 (truncation function).
# 
# ------------------------------------------------------------------------------
#   par)    par value(s)    fit toggle(s)    # parameter description 
# ------------------------------------------------------------------------------
"""

    # dict for get set looks like this: 
    # "position" : parameters["position"]
    exec(GalfitComponent.component_get_set(load_default_header_parameters()))
    
# ==========================================================================================================

    def __str__(self):
        return self.input_menu_file   + "\n\n" + \
               self.extra_header_info + "\n\n" + \
               self._section_sep      + "\n"   + \
               "# IMAGE and GALFIT CONTROL PARAMETERS\n" + \
               "\n".join(GalfitComponent.__str__(self).split("\n")[1:]) + \
               self.post_header
    
# ==========================================================================================================

    def update_parameters_file_helper(self, file_dict):
        
        for k,v in file_dict.items():
            if k.startswith("_"):
                continue
            
            value     = v
            if not issubclass(type(self.parameters[k]), MultiParameter):
                value = v[0]
            
            self.parameters[k].value = value
            
# ==========================================================================================================

    # This one is unlike the others so does not call
    # self.update_parameters_file_helper(file_dict)            
    @GalfitComponent.update_parameter_dict_with_dict
    def from_file_helper_dict(self, file_in, **kwargs):
        # Feed in just the chunk from the main 'from_file' caller
        # This requires determining param_begin/end in that caller
        
        # to be abundantly safe
        file_in = deepcopy(file_in)
        
        file_dict = kwargs
        # What can actually be gleaned from file
        file_dict["sigma_image"]         = file_in["SIGMA"]
        file_dict["psf"]                 = file_in["PSF"]
        file_dict["constraint_file"]     = file_in["CONSTRNT"]
        file_dict["region_to_fit"]       = tuple([int(v) for v in eval(file_in["FITSECT"].replace(":", ","))])
        file_dict["convolution_box"]     = tuple([int(v) for v in eval(f"({file_in['CONVBOX']})")])
        file_dict["mag_photo_zeropoint"] = float(file_in["MAGZPT"])
        
        return file_dict


# In[10]:


if __name__ == "__main__":
    from RegTest.RegTest import *


# In[11]:


# Unit Test for GalfitComponent
if __name__ == "__main__":
    component = GalfitComponent()
    print("Testing default values of base class GalfitComponent...")
    for k,v in component.__dict__.items():
        print(k,v)


# In[12]:


if __name__ == "__main__":
    bogus_list = """A) /home/portmanm/run6_1000_galfit_two_fit/sparcfire-in/1237667783385612464.fits      # Input data image (FITS file)
B) /home/portmanm/run6_1000_galfit_two_fit/sparcfire-tmp/galfits/1237667783385612464_out.fits      # Output data image block
C) none                # Sigma image name (made from data if blank or "none")
D) none                # Input PSF image and (optional) diffusion kernel
E) 1                   # PSF fine sampling factor relative to data
F) /home/portmanm/run6_1000_galfit_two_fit/sparcfire-tmp/galfit_masks/1237667783385612464_star-rm.fits      # Bad pixel mask (FITS image or ASCII coord list)
G) none                # File with parameter constraints (ASCII file)
H) 43   111  43   111  # Image region to fit (xmin xmax ymin ymax)
I) 50     50           # Size of the convolution box (x y)
J) 24.800              # Magnitude photometric zeropoint
K) 0.396  0.396        # Plate scale (dx dy)   [arcsec per pixel]
O) regular             # Display type (regular, curses, both)
P) 0                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps""".split("\n")

    bogus_dict = eval("""{'INITFILE': '/home/portmanm/testing_python_control/sparcfire-out/1237667429560025',
 'DATAIN': '/home/portmanm/testing_python_control/sparcfire-in/12376674295600251',
 'SIGMA': 'none',
 'PSF': 'none',
 'CONSTRNT': 'none',
 'MASK': '/home/portmanm/testing_python_control/sparcfire-tmp/galfit_masks/123',
 'FITSECT': '[36:98,37:99]',
 'CONVBOX': '51, 50',
 'MAGZPT': 24.8}""")

    header = GalfitHeader()
    header.from_file_helper_list(bogus_list)
    #header.update_param_values()
    print(header)
    header.to_file(f"{base_out}_header.txt")
    print()

    header.from_file_helper_dict(bogus_dict)
    #header.update_param_values()
    print(header)


# In[13]:


if __name__ == "__main__":
    bogus_list = """ 0) sersic                 #  Component type
 1) 76.7000  76.5000  0 0  #  Position x, y
 3) 12.9567     1          #  Integrated magnitude
 4) 18.5147     1          #  R_e (effective radius)   [pix]
 5) 0.6121      1          #  Sersic index n (de Vaucouleurs n=4)
 6) 0.0000      0          #     -----
 7) 0.0000      0          #     -----
 8) 0.0000      0          #     -----
 9) 0.3943      1          #  Axis ratio (b/a)
 10) -48.3372    1          #  Position angle (PA) [deg: Up=0, Left=90]
 Z) 0                      #  Skip this model in output image?  (yes=1, no=0)""".split("\n")

    bogus_dict = eval("""{'1_XC': '[67.3796]',
     '1_YC': '[67.7662]',
     '1_MAG': '13.1936 +/- 0.0257',
     '1_RE': '15.5266 +/- 0.1029',
     '1_N': '0.3433 +/- 0.0064',
     '1_AR': '0.6214 +/- 0.0039',
     '1_PA': '-19.1534 +/- 0.5867'}""")

    bulge = Sersic(1)
    bulge.from_file_helper_list(bogus_list)
    
    print(bulge)
    bulge.to_file(f"{base_out}_Sersic.txt")

    bulge.from_file_helper_dict(bogus_dict)
    print(bulge)
    
    bulge_df = bulge.to_pandas()
    print(bulge_df)
    print()
    bulge_df.iloc[0,0] = 111
    bulge_df.iloc[0,1] = 112
    bulge_df.iloc[0,2] = 113
    bulge.from_pandas(bulge_df)
    print(bulge)
    
    log_line = "sersic    : (  [62.90],  [62.90])  14.11*     13.75    0.30    0.63    60.82"
    
    bulge.update_from_log(log_line)
    print(bulge)


# In[14]:


# #Will use this once I diff check everything
# if __name__ == "__main__":
#     bogus_list = """ 0) sersic                 #  Component type
#  1) 76.7000  76.5000  0 0  #  Position x, y
#  3) 12.9567     1          #  Integrated magnitude
#  4) 18.5147     1          #  R_e (effective radius)   [pix]
#  5) 0.6121      1          #  Sersic index n (de Vaucouleurs n=4)
#  6) 0.0000      0          #     -----
#  7) 0.0000      0          #     -----
#  8) 0.0000      0          #     -----
#  9) 0.3943      1          #  Axis ratio (b/a)
#  10) -48.3372    1          #  Position angle (PA) [deg: Up=0, Left=90]
#  Z) 0                      #  Skip this model in output image?  (yes=1, no=0)""".split("\n")

#     bogus_dict = eval("""{'1_XC': '[67.3796]',
#      '1_YC': '[67.7662]',
#      '1_MAG': '13.1936 +/- 0.0257',
#      '1_RE': '15.5266 +/- 0.1029',
#      '1_N': '0.3433 +/- 0.0064',
#      '1_AR': '0.6214 +/- 0.0039',
#      '1_PA': '-19.1534 +/- 0.5867'}""")

#     bulge = Sersic(1)
#     print("Defaults--\n", bulge.parameters)
#     print()
    
#     bulge.magnitude.value = 5
#     print("Modifying magnitude directly--\n", bulge)
    
#     bulge.from_file_helper_list(bogus_list)
#     print("From file helper, list--\n", bulge)
    
#     print("Sending to file--\n", bulge)
#     bulge.to_file(f"{base_out}_Sersic.txt")

#     bulge = Sersic(1)
#     bulge.from_file_helper_dict(bogus_dict)
#     print("From file helper, dict--\n", bulge)
    
#     bulge_df = bulge.to_pandas()
#     print("To pandas--\n", bulge_df)
#     print()
    
#     bulge_df.iloc[0,0] = 111
#     bulge_df.iloc[0,1] = 112
#     bulge_df.iloc[0,2] = 113
#     bulge.from_pandas(bulge_df)
#     print("From pandas, modified parameters--\n", bulge)
    
#     log_line = "sersic    : (  [62.90],  [62.90])  14.11*     13.75    0.30    0.63    60.82"
#     bulge.update_from_log(log_line)
#     print("From log line--\n", bulge)


# In[15]:


if __name__ == "__main__":
    bogus_list = """ R0) power                  #  PA rotation func. (power, log, none)
R1) 0.0000      0          #  Spiral inner radius [pixels]
R2) 42.0200     0          #  Spiral outer radius [pixels]
R3) 595.0912    1          #  Cumul. rotation out to outer radius [degrees]
R4) -0.1961     1          #  Asymptotic spiral powerlaw
R9) 49.1328     1          #  Inclination to L.o.S. [degrees]
R10) 72.0972    1          #  Sky position angle""".split("\n")

    bogus_dict = eval("""{'2_ROTF': 'power',
 '2_RIN': '[0.0000]',
 '2_ROUT': '[22.0110]',
 '2_RANG': '79.0069 +/- 11.7225',
 '2_ALPHA': '-2.3697 +/- 0.0691',
 '2_INCL': '40.8043 +/- 2.7380',
 '2_SPA': '24.3010 +/- 4.5444'}""")

    arms = Power(2)
    arms.from_file_helper_list(bogus_list)
    #arms.update_param_values()
    arms.to_file(f"{base_out}_Power.txt")
    print(arms)

    arms.from_file_helper_dict(bogus_dict)
    #arms.update_param_values()
    print(arms)
    
    arms_df = arms.to_pandas()
    print(arms_df)
    print()
    arms_df.iloc[0,0] = 111
    arms_df.iloc[0,1] = 112
    arms_df.iloc[0,2] = 113
    arms.from_pandas(arms_df)
    print(arms)
    
    log_line = "power   :     [0.00]   23.51  219.64     -0.16*     ---  -44.95   -15.65"
    
    arms.update_from_log(log_line)
    print(arms)


# In[16]:


if __name__ == "__main__":
    bogus_list = """ F1) 0.2721   -56.9126 1 1  #  Azim. Fourier mode 1, amplitude, & phase angle
F3) -0.0690  -31.8175 1 1  #  Azim. Fourier mode 3, amplitude, & phase angle""".split("\n")

    bogus_dict = eval("""{'2_F1': '0.1449 +/- 0.0123',
 '2_F1PA': '44.3015 +/- 7.1154',
 '2_F3': '0.0979 +/- 0.0104',
 '2_F3PA': '-35.1366 +/- 4.4060'}""")

    fourier = Fourier(2)
    fourier.from_file_helper_list(bogus_list)
    #fourier.update_param_values()
    fourier.to_file(f"{base_out}_Fourier.txt")
    print(fourier)

    fourier.from_file_helper_dict(bogus_dict)
    #fourier.update_param_values()
    print(fourier)
    print()
    
    print(fourier.parameters)
    print()
    fourier_df = fourier.to_pandas()
    print(fourier_df)
    print()
    fourier_df.iloc[0,0] = 111
    fourier_df.iloc[0,1] = 112
    fourier_df.iloc[0,2] = 113
    fourier.from_pandas(fourier_df)
    print(fourier)
    print()
    
    log_line = "fourier : (1:  0.06,   -6.67)   (3:  0.05,    0.18)"
    
    fourier.update_from_log(log_line)
    print(fourier)


# In[17]:


# if __name__ == "__main__":
#     arms.add_skip(skip_val = 1)
#     fourier.add_skip(skip_val = 1)
    
#     # We should see an RZ, FZ which is nonsensical but when we write to file
#     # the power and fourier functions simply won't show

#     print(arms)
#     print(fourier)
    
#     bulge.to_file(f"{base_out}_PowerFourierSkip.txt", arms, fourier)


# In[18]:


if __name__ == "__main__":
    bogus_list = """  0) sky                    #  Component type
 1) 1112.1005    1          #  Sky background at center of fitting region [ADUs]
 2) 1.264e-02      1       #  dsky/dx (sky gradient in x)     [ADUs/pix]
 3) 1.813e-02      1       #  dsky/dy (sky gradient in y)     [ADUs/pix]
 Z) 0                      #  Skip this model in output image?  (yes=1, no=0)""".split("\n")

    sky = Sky(3)
    sky.from_file_helper_list(bogus_list)
    #sky.update_param_values()
    sky.to_file(f"{base_out}_Sky.txt")
    print(sky)
    
    bogus_dict = eval("""{'COMP_3': 'sky',
 '3_XC': '[67.0000]',
 '3_YC': '[68.0000]',
 '3_SKY': '1133.4166 +/- 0.1595',
 '3_DSDX': '0.0119 +/- 0.0048',
 '3_DSDY': '-0.0131 +/- 0.0047'}""")

    sky.from_file_helper_dict(bogus_dict)
    #sky.update_param_values()
    print(sky)
    
    sky_df = sky.to_pandas()
    print(sky_df)
    print()
    sky_df.iloc[0,0] = 111
    sky_df.iloc[0,1] = 112
    sky_df.iloc[0,2] = 113
    sky.from_pandas(sky_df)
    print(sky)
    
    log_line = "sky       : [ 63.00,  63.00]  1130.51  -4.92e-02  1.00e-02"
    
    sky.update_from_log(log_line)
    print(sky)


# In[19]:


if __name__ == "__main__":
    export_to_py("Components", pj(_MODULE_DIR, "Classes", "Components"))

