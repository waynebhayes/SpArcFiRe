#!/usr/bin/env python
# coding: utf-8

# In[969]:


import os
import sys
from os.path import join as pj
from os.path import exists
import subprocess
from copy import deepcopy
from IPython import get_ipython
from astropy.io import fits


# In[970]:


class GalfitComponent:
    def __init__(self,
                 component_type = "", 
                 component_number = 0, 
                 param_prefix = " "
                ):
        
        self.component_type   = component_type
        self.component_number = component_number
        self.param_prefix     = param_prefix
        
        key = "Component type"
        self.param_desc = {key : key}
        
        # tuple in value for printing
        self.param_numbers = {0 : key}
        self.param_values = {key : component_type} #(component_type, "")}
        self.param_fix = {key : ""}
        
        # For reading from file
        self.start_dict = f"COMP_{self.component_number}"
        self.end_dict   = f"COMP_{self.component_number + 1}"
        
        #self.start_dict_value = self.component_type
        
        self.start_text = f"# Component number: {self.component_number}"
        self.end_text   = f"{self.param_prefix.strip()}10"
        
# ==========================================================================================================

    def add_skip(self, skip_val = 0):
        key = "skip"
        self.param_numbers["Z"] = key
        self.param_values[key] = skip_val
        self.param_desc[key] = "Skip this model in output image?  (yes=1, no=0)"
        self.param_fix[key] = ""
        
# ==========================================================================================================

    def update_values(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)
        
# ==========================================================================================================

    def update_param_values(self):
        # Values are the only thing that will change via instance.param_name
        # We need to update the dictionary that holds all of them since I use
        # a deepcopy to save myself some lines of code/typing (cheeky)
        self.param_values.update({k:v for k,v in self.__dict__.items() if k in self.param_values})
        
# ==========================================================================================================

    def check_lengths(self):
        
        if not all(len(ec) == len(self.param_numbers) for ec in (self.param_numbers, self.param_values, self.param_desc, self.param_fix)):
            print("Length of array containing parameter numbers, values, and descriptions must be the same!")
            
            print("param_numbers:", end = " ")
            print(*self.param_numbers.keys(), sep = "\n")
            
            print("param_values:", end = " ")
            print(*self.param_values.keys(), sep = "\n")
            
            print("param_fix:", end = " ")
            print(*self.param_fix.keys(), sep = "\n")
            
            print("param_desc:", end = " ")
            print(*self.param_desc.keys(), sep = "\n")
            
            raise(AssertionError())
            
# ==========================================================================================================

    # TODO: Add comparison function via subtract
    def __sub__(self):
        pass
    
# ==========================================================================================================
  
    def __str__(self):
        output_str = ""
        num_type = ".4f"
        l_align = "<11"
        
        if self.component_type == "header":
            output_str += "\n\n".join((self.input_menu_file, self.extra_header_info, "="*80))
            output_str += "\n# IMAGE and GALFIT CONTROL PARAMETERS\n"
            num_type = "" #".0f"
            
        elif self.component_type not in ("power", "fourier", "header"):
            output_str = f"# Component number: {self.component_number}\n"
        
        self.check_lengths()
        
        for num, val, fix, desc in zip(self.param_numbers.keys(), self.param_values.values(), self.param_fix.values(), self.param_desc.values()):
            # Skip component type
            if isinstance(val, (int, float)):
                if isinstance(num, str):
                    line = f"{self.param_prefix}{num}) {val:{l_align}} {fix}"
                else:
                    line = f"{self.param_prefix}{num}) {val:{l_align}{num_type}} {fix}"
                
                # For 10) in sersic
                if len(f"{self.param_prefix}{num}") >= 3: 
                    line = line.lstrip()
                
            else:
                # position #, fourier, etc.
                if isinstance(val[0], (int, float)):
                    line = f"{self.param_prefix}{num}) {val[0]:<7{num_type}} {val[1]:{num_type}} {fix}"
                    
                # comp name
                else:
                    line = f"{self.param_prefix}{num}) {val}" #[0]}"
                
            output_str += f"{line:<23} # {desc}\n"
            
        if self.component_type == "header":
            output_str += self.post_header

        return output_str

# ==========================================================================================================
    
    def update_from_log(self, in_line:str):
        # This is a dummy function which is overwritten in all
        # sub-classes. Just leaving here for documentation
        # Used to update from stdout i.e. what we see in fit.log
        # Requires outside function to loop through log file
        
        # NOTE: These necessarily round to two digits because that's
        # all Galfit outputs to stdout
        print("Did this get properly overwritten?")
        
# ==========================================================================================================

    # TODO: Import parameters from pandas
    def from_pandas(self):
        pass

# ==========================================================================================================

    # TODO: Export parameters from pandas
    def to_pandas(self):
        pass
# ==========================================================================================================

    def from_file_helper(self):
        # This function is a placeholder
        # All components will have this, it'll make it easier and less
        # confusing to debug/parse the inputs
        pass

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
            
            # Really just for power and fourier function but technically this makes it more generalizable
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
                
                # print(ve)
                # for k,v in input_in.items():
                #     print(k,v)
                
                # 
                # 
                # # Usually of form #_param
                # component_end   = f"{component_start[0]}{self.end_dict[1:]}" 
                # print(component_start, component_end)               
            
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
        else:
            feed_in = from_text(self, filename)
        
        self.from_file_helper(feed_in)
        self.update_param_values()

# ==========================================================================================================

    def to_file(self, filename, *args):
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
                    f.write(str(component))
                    if i != fourier_index - 1:
                        f.write("\n")

                f.write("="*80 + "\n")
                
        except FileNotFoundError:
            print(f"Can't open to write the file, {filename}. Check permissions/directory.")
            
        except OSError as ose:
            print(f"Something went wrong! {ose}")
    
    
    # Either 
    # instance.param_name = #
    # instance.update_param_values()
    # or
    # instance.param_values["param_name"] = #
    # works
            
    # #https://stackoverflow.com/a/35282351
    # def __iter__(self):
    #     return self.__dict__.iteritems()


# In[971]:


class Sersic(GalfitComponent):
    def __init__(self, component_number, **kwargs):
        
        self.position         = kwargs.get("position", (0.0000, 0.0000))
        self.magnitude        = float(kwargs.get("magnitude", 13))
        self.effective_radius = float(kwargs.get("effective_radius", 10))
        self.sersic_index     = float(kwargs.get("sersic_index", 4))
        self.axis_ratio       = float(kwargs.get("axis_ratio", 0.6))
        self.position_angle   = float(kwargs.get("position_angle", 0))
        #self.skip = kwargs.get("skip", 0)
        
        param_dict = deepcopy(vars(self))
        
        GalfitComponent.__init__(self, component_type = "sersic", component_number = component_number)
        
        self.param_numbers.update(dict(
                                       zip([1,3,4,5,9,10], param_dict.keys())
                                      )
                                  )
        
        self.param_values.update(param_dict)
                
        self.param_fix.update({k : 1 for k in param_dict.keys()})
        self.param_fix["position"] = "0 0"
        
        pd = self.param_desc
        pd["position"] = "Position x, y"
        pd["magnitude"] = "Integrated magnitude"
        pd["effective_radius"] = "R_e (effective radius)   [pix]"
        pd["sersic_index"] = "Sersic index n (de Vaucouleurs n=4)"
        pd["axis_ratio"] = "Axis ratio (b/a)"
        pd["position_angle"] = "Position angle (PA) [deg: Up=0, Left=90]"
        
        # For reading from file
        self.start_dict = f"COMP_{self.component_number}"
        self.end_dict   = f"{self.component_number}_PA"
        
        # text kept at defaults
        #self.start_text = f"# Component number: {self.component_number}"
        #self.end_text   = f"{self.param_prefix.strip()}10"
        
    def from_file_helper(self, file_in):
        # Feed in just the chunk from the main 'from_file' caller
        # This requires determining param_begin/end in that caller
        # This can handle component type but it is unnecessary 
        
        # to be abundantly safe
        file_in = deepcopy(file_in)
        
        # Excludes 0
        p_numbers = list(self.param_numbers.keys())[1:]
        
        # File dict is only for fits
        if isinstance(file_in, dict):
            file_dict = {}
            file_in = {k:v for k,v in file_in.items() if v.lower() != "sersic"}
            
            x_key = [name for name in file_in.keys() if name.endswith("_XC")][0]
            x_pos = file_in.pop(x_key).strip("[]")
            
            for num, v in zip(p_numbers, file_in.values()): # file_dict.keys()
                if "+/-" in v:
                    v = v.split("+/-")[0]
                    
                file_dict[num] = [float(v.strip("[]"))]
                
                # For param_fix
                if "[" in v and "]" in v:
                    file_dict[num].append(0)

                    if num == 1:
                        file_dict[1] += [0, 0]
                        
                else:
                    file_dict[num].append(1)
            
            file_dict[1] = [float(x_pos)] + file_dict[1]
                        
        elif isinstance(file_in, list):
            file_list = file_in
            file_dict = {int(line[:line.index(")")].strip()) : line[line.index(")") + 1 : line.index("#")].strip()
                         for line in file_list if line.strip()[0] not in ("#", "Z")}
            
            file_dict = {k: " ".join(v.split()).split() for k,v in file_dict.items()}
            
            if file_list[-1].strip().startswith("Z"):
                self.add_skip()
                            
        else:
            print(f"Could not grab {self.component_type} components, no file list/dict specified.")
            return 
        
        self.position         = (float(file_dict[1][0]), float(file_dict[1][1]))
        position_fix          = f"{file_dict[1][2]} {file_dict[1][3]}"
        self.magnitude        = float(file_dict[3][0])
        self.effective_radius = float(file_dict[4][0])
        self.sersic_index     = float(file_dict[5][0])
        self.axis_ratio       = float(file_dict[9][0])
        self.position_angle   = float(file_dict[10][0])

        self.param_fix.update({self.param_numbers[k] : int(file_dict[k][1]) for k in p_numbers if k not in [1, "Z"]})
        self.param_fix["position"] = position_fix
        
        # This will happen in the main call
        #self.update_param_values()
        return
        
    def update_from_log(self, in_line):
        
        # example
        # sersic    : (  [62.90],  [62.90])  14.11     13.75    0.30    0.63    60.82
        # does not include position
        params = in_line.split("])")[1]
        params = " ".join(params.split())
        params = params.split()
        params = [i.strip("*") for i in params]
        self.magnitude        = float(params[0])
        self.effective_radius = float(params[1])
        self.sersic_index     = float(params[2])
        self.axis_ratio       = float(params[3])
        self.position_angle   = float(params[4])


# In[972]:


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
    bulge.from_file_helper(bogus_list)
    bulge.update_param_values()
    print(bulge)

    bulge.from_file_helper(bogus_dict)
    bulge.update_param_values()
    print(bulge)


# In[973]:


class Power(GalfitComponent):
    def __init__(self, **kwargs):
        
        #self.component_type = "power"
        self.inner_rad          = float(kwargs.get("inner_rad", 0))
        self.outer_rad          = float(kwargs.get("outer_rad", 20))
        self.cumul_rot          = float(kwargs.get("cumul_rot", 90))
        self.powerlaw           = float(kwargs.get("powerlaw", 1))
        self.inclination        = float(kwargs.get("inclination", 15))
        self.sky_position_angle = float(kwargs.get("sky_position_angle", 45))
        #self.skip = kwargs.get("skip", 0)
        
        param_dict = deepcopy(vars(self))
        
        GalfitComponent.__init__(self, component_type = "power", param_prefix = "R")
        
        self.param_numbers.update(dict(
                                       zip([1,2,3,4,9,10], param_dict.keys())
                                      )
                                  )
        
        self.param_values.update(param_dict)
        
        self.param_fix.update({k : 1 for k in param_dict.keys()})
        self.param_fix["inner_rad"] = 0
        self.param_fix["outer_rad"] = 0
        
        pd = self.param_desc
        pd["Component type"] = "PA rotation func. (power, log, none)"
        pd["inner_rad"] = "Spiral inner radius [pixels]"
        pd["outer_rad"] = "Spiral outer radius [pixels]"
        pd["cumul_rot"] = "Cumul. rotation out to outer radius [degrees]"
        pd["powerlaw"] = "Asymptotic spiral powerlaw"
        pd["inclination"] = "Inclination to L.o.S. [degrees]"
        pd["sky_position_angle"] = "Sky position angle"
        
        # For reading from file
        # 2_ may not always be the case but that's why I have a try except in there ;)
        self.start_dict = f"2_ROTF"
        self.end_dict   = f"2_SPA"
        
        self.start_text = f"{self.param_prefix}0) power"
        # end kept at defult
        #self.end_text   = f"{self.param_prefix.strip()}10"
        
    def from_file_helper(self, file_in):
        # Feed in just the chunk from the main 'from_file' caller
        # This requires determining param_begin/end in that caller
        # This can handle component type but it is unnecessary 
        
        # to be abundantly safe
        file_in = deepcopy(file_in)
        
        # Excludes 0
        p_numbers = list(self.param_numbers.keys())[1:]                 
            
        # File dict is only for fits
        if isinstance(file_in, dict):
            file_dict = {}
            file_in = {k:v for k,v in file_in.items() if v.lower() != "power"}
            
            for num, v in zip(p_numbers, file_in.values()): # file_dict.keys()
                if "+/-" in v:
                    v = v.split("+/-")[0]
                    
                file_dict[num] = [float(v.strip("[]"))]
                
                # For param_fix
                if "[" in v and "]" in v:
                    file_dict[num].append(0)

                else:
                    file_dict[num].append(1)
            
        elif isinstance(file_in, list):
            file_list = file_in
            file_dict = {int(line[:line.index(")")].strip(f"{self.param_prefix} ")) : line[line.index(")") + 1 : line.index("#")].strip()
                         for line in file_list if line.strip()[0] not in ("#", "Z")}
            
            file_dict = {k: " ".join(v.split()).split() for k,v in file_dict.items()}
            
            if file_list[-1].strip().startswith("Z"):
                self.add_skip()
                            
        else:
            print(f"Could not grab {self.component_type} components, no file list/dict specified.")
            return
                
        self.inner_rad          = float(file_dict[1][0])
        self.outer_rad          = float(file_dict[2][0])
        self.cumul_rot          = float(file_dict[3][0])
        self.powerlaw           = float(file_dict[4][0])
        self.inclination        = float(file_dict[9][0])
        self.sky_position_angle = float(file_dict[10][0])

        self.param_fix.update({self.param_numbers[k] : int(file_dict[k][1]) for k in p_numbers if k not in [1, "Z"]})
        
        # This will happen in the main call
        #self.update_param_values()
        return
        
    def update_from_log(self, in_line):
        # example
        # power   :     [0.00]   [23.51]  219.64     -0.16     ---  -44.95   -15.65
        # does not include inner/outer rad
        # Assumes those are the only ones fixed
        
        params = in_line.split("]")[2] 
        params = " ".join(params.split())
        params = params.split()
        params = [i.strip("*") for i in params]
        self.cumul_rot          = float(params[0])
        self.powerlaw           = float(params[1])
        # Some dashed line here...
        self.inclination        = float(params[3])
        self.sky_position_angle = float(params[4])


# In[974]:


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

    arms = Power()
    arms.from_file_helper(bogus_list)
    arms.update_param_values()
    print(arms)

    arms.from_file_helper(bogus_dict)
    arms.update_param_values()
    print(arms)


# In[975]:


class Fourier(GalfitComponent):
    def __init__(self, n = {1 : (0.05, 45), 3 : (0.05, 25)}):
        GalfitComponent.__init__(self, component_type = "fourier", param_prefix = "F")      
        # normal rules don't apply here
        # Still use inheritance for the other functions
        self.param_numbers = {}
        self.param_values = {}
        self.param_desc = {}
        self.param_fix = {}
        
        self.include_fn(n = n)
        
        p_numbers = list(self.param_numbers.keys())
        # For reading from file
        self.start_dict = f"2_F{p_numbers[0]}"
        self.end_dict   = f"2_F{p_numbers[-1]}PA"
        
        self.start_text = f"F{p_numbers[0]}"
        self.end_text   = f"{self.param_prefix}{p_numbers[-1]}"
        
    def include_fn(self, n:dict):
        for num, values in n.items():
            key = f"{self.param_prefix}{num}"
            self.param_numbers[num] = key
            self.param_values[key] = values
            self.param_desc[key] = f"Azim. Fourier mode {num}, amplitude, & phase angle"
            self.param_fix[key] = self.param_fix.get(key, "1 1")

    def from_file_helper(self, file_in):
        # Feed in just the chunk from the main 'from_file' caller
        # This requires determining param_begin/end in that caller
        
        # to be abundantly safe
        file_in = deepcopy(file_in)
        
        p_numbers = list(self.param_numbers.keys())
        
        # File dict is only for fits
        if isinstance(file_in, dict):
            n_dict = {}
            
            for k,v in file_in.items():
                
                if "+/-" in v:
                    v = v.split("+/-")[0]
                    
                #ex: 1_F1 -> 1
                # 1_F3PA -> 3
                k_num = int(k.split("_")[1][1])
                if k_num not in n_dict:
                    n_dict[k_num] = [float(v.strip("[] "))]
                else:
                    n_dict[k_num] += [float(v.strip("[] "))]
                
                # For param_fix
                if "[" in v and "]" in v:
                    self.param_fix[f"{self.param_prefix}{k_num}"] = "0 0"
                        
                else:
                    self.param_fix[f"{self.param_prefix}{k_num}"] = "1 1"
                
        elif isinstance(file_in, list):
            file_list = file_in
            file_dict = {int(line[:line.index(")")].strip(f"{self.param_prefix} ")) : line[line.index(")") + 1 : line.index("#")].strip()
                         for line in file_list if line.strip()[0] not in ("#", "Z")}
            
            n_dict = {}
            for k,v in file_dict.items():
                value = " ".join(v.split()).split()
                n_dict[k] = (float(value[0]), float(value[1]))
            
            if file_list[-1].strip().startswith("Z"):
                self.add_skip()
                            
        else:
            print(f"Could not grab {self.component_type} components, no file list/dict specified.")
            return
        
        self.include_fn(n_dict)
        
        # TODO: do I need these? Do I need to update param_fix? Fourier is different

        # self.param_fix.update({self.param_numbers[k] : int(file_dict[k][1]) for k in p_numbers if k not in [1, "Z"]})
        # self.param_fix["position"] = position_fix
        
        # This will happen in the main call
        #self.update_param_values()
        return
        
    def update_from_log(self, in_line):
        # example
        # fourier : (1:  0.06,   -6.67)   (3:  0.05,    0.18)
        
        # rstrip avoids a hanging ) later
        params = in_line.lstrip("fourier : ").replace(" ", "").rstrip(")").split(")(")
        
        self.param_values = {n: eval(f"({params[i].split(':')[1].replace('*', '')})")
                             for i, n in enumerate(self.param_values.keys())}


# In[976]:


if __name__ == "__main__":
    bogus_list = """ F1) 0.2721   -56.9126 1 1  #  Azim. Fourier mode 1, amplitude, & phase angle
F3) -0.0690  -31.8175 1 1  #  Azim. Fourier mode 3, amplitude, & phase angle""".split("\n")

    bogus_dict = eval("""{'2_F1': '0.1449 +/- 0.0123',
 '2_F1PA': '44.3015 +/- 7.1154',
 '2_F3': '0.0979 +/- 0.0104',
 '2_F3PA': '-35.1366 +/- 4.4060'}""")

    fourier = Fourier()
    fourier.from_file_helper(bogus_list)
    fourier.update_param_values()
    print(fourier)

    fourier.from_file_helper(bogus_dict)
    fourier.update_param_values()
    print(fourier)


# In[977]:


class Sky(GalfitComponent):
    def __init__(self, component_number, **kwargs):
        self.sky_background = float(kwargs.get("sky_background", 1000))
        self.dsky_dx        = float(kwargs.get("dsky_dx", 0))
        self.dsky_dy        = float(kwargs.get("dsky_dy", 0))
        
        param_dict = deepcopy(vars(self))
        
        GalfitComponent.__init__(self, component_type = "sky", component_number = component_number)
        
        self.param_numbers.update(dict(
                                       zip([1,2,3], param_dict.keys())
                                      )
                                  )
        
        self.param_values.update(param_dict)
        
        self.param_fix.update({k : 1 for k in param_dict.keys()})
        
        pd = self.param_desc
        pd["sky_background"] = "Sky background at center of fitting region [ADUs]"
        pd["dsky_dx"] = "dsky/dx (sky gradient in x)     [ADUs/pix]"
        pd["dsky_dy"] = "dsky/dy (sky gradient in y)     [ADUs/pix]"
        
        # For reading from file
        self.start_dict = f"COMP_{self.component_number}"
        self.end_dict   = f"{self.component_number}_DSDY"
        
        self.start_text = f"# Component number: {self.component_number}"
        self.end_text   = f"{self.param_prefix.strip()}3"
        
    def from_file_helper(self, file_in):
        # Feed in just the chunk from the main 'from_file' caller
        # This requires determining param_begin/end in that caller
        # This can handle component type but it is unnecessary 
        
        # to be abundantly safe
        file_in = deepcopy(file_in)
        
        # Excludes 0
        p_numbers = list(self.param_numbers.keys())[1:]                 
            
        # File dict is only for fits
        if isinstance(file_in, dict):
            file_dict = {}
            file_in = {k:v for k,v in file_in.items() if v.lower() != "sky" and k.split("_")[1] not in ["XC", "YC"]}
            
            for num, v in zip(p_numbers, file_in.values()): # file_dict.keys()
                if "+/-" in v:
                    v = v.split("+/-")[0]
                    
                file_dict[num] = [float(v.strip("[]"))]
                
                # For param_fix
                if "[" in v and "]" in v:
                    file_dict[num].append(0)

                else:
                    file_dict[num].append(1)
            
        elif isinstance(file_in, list):
            file_list = file_in
            file_dict = {int(line[:line.index(")")].strip(f"{self.param_prefix} ")) : line[line.index(")") + 1 : line.index("#")].strip()
                         for line in file_list if line.strip()[0] not in ("#", "Z")}
            
            file_dict = {k: " ".join(v.split()).split() for k,v in file_dict.items()}
            
            if file_list[-1].strip().startswith("Z"):
                self.add_skip()
                            
        else:
            print(f"Could not grab {self.component_type} components, no file list/dict specified.")
            return
                
        self.sky_background = float(file_dict[1][0])
        self.dsky_dx        = float(file_dict[2][0])
        self.dsky_dy        = float(file_dict[3][0])
        
        self.param_fix.update({self.param_numbers[k] : int(file_dict[k][1]) for k in p_numbers if k not in [1, "Z"]})
        
        # This will happen in the main call
        #self.update_param_values()
        return
        
    
    def update_from_log(self, in_line):
        # example
        #sky       : [ 63.00,  63.00]  1130.51  -4.92e-02  1.00e-02
        
        params = in_line.split("]")[1]
        params = " ".join(params.split())
        params = params.split()
        params = [i.strip("*") for i in params]
        self.sky_background = float(params[0])
        self.dsky_dx = float(params[1])
        self.dsky_dy = float(params[1])


# In[978]:


if __name__ == "__main__":
    bogus_list = """  0) sky                    #  Component type
 1) 1112.1005    1          #  Sky background at center of fitting region [ADUs]
 2) 1.264e-02      1       #  dsky/dx (sky gradient in x)     [ADUs/pix]
 3) 1.813e-02      1       #  dsky/dy (sky gradient in y)     [ADUs/pix]
 Z) 0                      #  Skip this model in output image?  (yes=1, no=0)""".split("\n")

    bogus_dict = eval("""{'COMP_3': 'sky',
 '3_XC': '[67.0000]',
 '3_YC': '[68.0000]',
 '3_SKY': '1133.4166 +/- 0.1595',
 '3_DSDX': '0.0119 +/- 0.0048',
 '3_DSDY': '-0.0131 +/- 0.0047'}""")

    sky = Sky(3)
    sky.from_file_helper(bogus_list)
    sky.update_param_values()
    print(sky)

    sky.from_file_helper(bogus_dict)
    sky.update_param_values()
    print(sky)


# In[979]:


class GalfitHeader(GalfitComponent):
    
    def __init__(self, galaxy_name = "", **kwargs):
        
        # If not fully specified, will use galaxy_name as default so it's good to use
        # it as an argument even if I am specifying each individually
        self.input_image     = kwargs.get("input_image", f"{galaxy_name}.fits")
        self.output_image    = kwargs.get("output_image", f"{galaxy_name}_galfit_out.fits")
        self.sigma_image     = kwargs.get("sigma_image", "none")
        self.psf             = kwargs.get("psf", "none") # May add gname to this
        self.fine_sampling   = kwargs.get("fine_sampling", 1)
        self.pixel_mask      = kwargs.get("pixel_mask", f"{galaxy_name}_star-rm.fits")
        self.constraints     = kwargs.get("constraints", "none")
        self.region_to_fit   = kwargs.get("region_to_fit", (0, 255, 0, 255))
        self.convolution_box = kwargs.get("convolution_box", (50, 50))
        self.mag_zeropoint   = kwargs.get("mag_zeropoint" , 24.800) # SDSS DR7
        self.plate_scale     = kwargs.get("plate_scale", (0.396, 0.396)) # SDSS DR7
        self.display_type    = kwargs.get("display_type", "regular")
        self.optimize        = kwargs.get("optimize", 0)
            
        header_dict = deepcopy(self.__dict__)
        # normal rules don't apply here
        # Still use inheritance for the other functions
        GalfitComponent.__init__(self, component_type = "header", param_prefix = "")

        self.param_numbers = dict(
                                  zip(["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "O", "P"], header_dict.keys())
                                  )
                                  
        self.param_values = deepcopy(header_dict)
        # For region to fit, otherwise the 'fix' is empty
        self.param_fix = {k: (f"  {v[2]}   {v[3]}" if k == "region_to_fit" else "") for k,v in deepcopy(header_dict).items()}
        
        self.param_desc = deepcopy(header_dict)
        pd = self.param_desc
        
        pd["input_image"]     = "Input data image (FITS file)"
        pd["output_image"]    = "Output data image block"
        pd["sigma_image"]     = "Sigma image name (made from data if blank or 'none')"
        pd["psf"]             = "Input PSF image and (optional) diffusion kernel"
        pd["fine_sampling"]   = "PSF fine sampling factor relative to data"
        pd["pixel_mask"]      = "Bad pixel mask (FITS image or ASCII coord list)"
        pd["constraints"]     = "File with parameter constraints (ASCII file)"
        pd["region_to_fit"]   = "Image region to fit (xmin xmax ymin ymax)"
        pd["convolution_box"] = "Size of the convolution box (x y)"
        pd["mag_zeropoint"]   = "Magnitude photometric zeropoint"
        pd["plate_scale"]     = "Plate scale (dx dy)   [arcsec per pixel]"
        pd["display_type"]    = "Display type (regular, curses, both)"
        pd["optimize"]        = "Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps"
        
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
        
    def from_file_helper(self, file_in, **kwargs):
        # Feed in just the chunk from the main 'from_file' caller
        # This requires determining param_begin/end in that caller
        
        # to be abundantly safe
        file_in = deepcopy(file_in)
        
        p_numbers = list(self.param_numbers.keys())
            
        # File dict is only for fits
        if isinstance(file_in, dict):
            file_dict = {}
            # From file
            file_dict["C"] = [file_in["SIGMA"]]
            file_dict["D"] = [file_in["PSF"]]
            file_dict["G"] = [file_in["CONSTRNT"]]
            file_dict["H"] = eval(file_in["FITSECT"].replace(":", ","))
            file_dict["I"] = eval(f"({file_in['CONVBOX']})")
            file_dict["J"] = [float(file_in["MAGZPT"])]
            
            # Everything else...
            file_dict["A"] = [kwargs.get("input_image", self.input_image)]
            file_dict["B"] = [kwargs.get("output_image", self.output_image)]
            file_dict["E"] = [kwargs.get("fine_sampling", self.fine_sampling)]
            file_dict["F"] = [kwargs.get("pixel_mask", self.pixel_mask)]
            file_dict["K"] = kwargs.get("plate_scale", self.plate_scale)
            file_dict["O"] = [kwargs.get("display_type", self.display_type)]
            file_dict["P"] = [kwargs.get("optimize", self.optimize)]
            
        elif isinstance(file_in, list):
            #print(file_in)
            file_list = file_in
            file_dict = {line[:line.index(")")].strip(f"{self.param_prefix} ") : line[line.index(")") + 1 : line.index("#")].strip()
                         for line in file_list if line.strip()[0] not in ("#")} #, "Z")}
            
            file_dict = {k: " ".join(v.split()).split() for k,v in file_dict.items()}
            #print(file_dict)
            # if file_list[-1].strip().startswith("Z"):
            #     self.add_skip()
                            
        else:
            print(f"Could not grab {self.component_type} components, no file list/dict specified.")
            return
        
        self.input_image     = file_dict["A"][0]
        self.output_image    = file_dict["B"][0]
        self.sigma_image     = file_dict["C"][0]
        self.psf             = file_dict["D"][0]
        self.fine_sampling   = int(file_dict["E"][0])
        self.pixel_mask      = file_dict["F"][0]
        self.constraints     = file_dict["G"][0]
        self.region_to_fit   = (int(file_dict["H"][0]), 
                                int(file_dict["H"][1]), 
                                int(file_dict["H"][2]), 
                                int(file_dict["H"][3]))
        self.param_fix["region_to_fit"] = f"{int(file_dict['H'][2])} {int(file_dict['H'][3])}"
        
        self.convolution_box = (int(file_dict["I"][0]), 
                                int(file_dict["I"][1]))
        self.mag_zeropoint   = float(file_dict["J"][0])
        self.plate_scale     = (float(file_dict["K"][0]),
                                float(file_dict["K"][0]))
        self.display_type    = file_dict["O"][0]
        self.optimize        = int(file_dict["P"][0])
        
        # This will happen in the main call
        #self.update_param_values()
        return


# In[980]:


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
 'CONVBOX': '50, 50',
 'MAGZPT': 24.8}""")

    header = GalfitHeader()
    header.from_file_helper(bogus_list)
    header.update_param_values()
    print(header)

    header.from_file_helper(bogus_dict)
    header.update_param_values()
    print(header)


# In[981]:


class ComponentContainer:
    def __init__(self, **kwargs):
        self.bulge   = kwargs.get("bulge", Sersic(1))
        self.disk    = kwargs.get("disk", Sersic(2))
        self.arms    = kwargs.get("arms", Power())
        self.fourier = kwargs.get("fourier", Fourier())
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
    
    def update_components(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)
       
    # This is more for the daughter classes
    def extract_components(self):
        return ComponentContainer(**vars(self))
    
    def __str__(self):
        out_str = "\n".join(str(comp) for comp in ComponentContainer.to_list(self))
        return out_str


# In[1017]:


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
        
    def to_file(self, *args):
        if args:
            self.header.to_file(self.path_to_feedme, *args)
        else:
            self.header.to_file(self.path_to_feedme, *ComponentContainer.to_list(self))
            
    def from_file(self, filename):
        # This function handles grabbing and storing the values from galfit files (input and output???)
        # It's written to generally handle both and stores everything in the respective component objects

        def from_fits(self, filename = filename):
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
        
            input_in = dict(input_file[2].header)
            input_keys = list(input_in.keys())
            
            component_list = self.to_list()
            component_num = 0
            send_to_helper = {}
            
            for idx, (key, value) in enumerate(input_in.items()):
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
        
        def from_text(self, filename = filename):
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
                    
            input_file.close()
        
            return
        
        ext = os.path.splitext(filename)[1]
        if ext == ".fits":
            from_fits(self)
        else:
            from_text(self)
            
        _ = [c.update_param_values() for c in self.to_list()]
        #_ = [c.update_param_values() for c in self.to_list()]        


# In[1018]:


class GalfitOutput(FeedmeContainer):
    def __init__(self, galfit_out_obj = subprocess.CompletedProcess("", 0), **kwargs):
        
        FeedmeContainer.__init__(self, **kwargs)
        
        galfit_out_text = galfit_out_obj.stdout        
        err_text = galfit_out_obj.stderr
        
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

            # Constrain to last 30 lines to save search time
            last_out_lines = galfit_out_text.split("\n")[-30:]
            if any(line.strip().startswith(success) for line in last_out_lines):
                self.success = True
                
            elif any(line.strip().startswith(failure) for line in last_out_lines):
                print(f"Galfit failed this run!")
                last_out_lines = '\n'.join(last_out_lines)
                print(f"{last_out_lines}")
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
                comp.update_param_values()

        if self.success:
            update_components(self, galfit_out_text)
            
    def __str__(self, galfit_out_text) -> str:
        return(galfit_out_text)


# In[1019]:


if __name__ == "__main__":
    # Testing basic functionality
    
    header = GalfitHeader(galaxy_name = "tester")
    bulge = Sersic(1, position = (25,25))
    disk  = Sersic(2, position = (25,25))
    arms  = Power()
    fourier = Fourier()
    sky   = Sky(3)
    
    container = FeedmeContainer(**{"header"  : header,
                                      "bulge"   : bulge,
                                      "disk"    : disk,
                                      "arms"    : arms,
                                      "fourier" : fourier,
                                      "sky"     : sky}
                                )
    
    #print(container)
    
    container.path_to_feedme = "tester.in"
    
    container.to_file() #, bulge, disk, arms, fourier, sky)


# In[1020]:


# Testing from_file
if __name__ == "__main__":
    
    header = GalfitHeader(galaxy_name = "fake_name")
    container.update_components(header = header)
    # bulge = Sersic(1)
    # disk  = Sersic(2)
    # arms  = Power()
    # fourier = Fourier()
    # sky   = Sky(3)
    
    real = "./sparcfire-out/1237667429560025179/1237667429560025179_galfit_out.fits"
    
    #for comp in (header, bulge, disk, arms, fourier, sky):
    #for comp in (fourier, sky):
        #comp.from_file("good_output.in")
        #print(comp.component_type)
        #comp.from_file(real)
        #print(comp)
    
    #container.from_file("good_output.in")
    #print(str(container))
    
    container.from_file(real)
    print(str(container))


# In[901]:


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
    


# In[707]:


if __name__ == "__main__":
    example_feedme = FeedmeContainer(path_to_feedme = "somewhere/out_there", 
                                     header         = header, 
                                     bulge          = bulge, 
                                     disk           = disk, 
                                     arms           = arms, 
                                     fourier        = fourier, 
                                     sky            = sky)
    feedme_components = example_feedme.extract_components()
    print(feedme_components.to_list())
    _ = [print("Key:", k) for k in example_feedme.to_dict().keys()]
    print(str(example_feedme))


# In[708]:


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
    
    bad_output  = GalfitOutput(dummy_obj)
    print(bad_output.bulge)
    print(bad_output.disk)
    print(bad_output.arms)
    print(bad_output.fourier)
    print(bad_output.sky)
    print(f"Did the bad example succeed? {bad_output.success}")
    
    print("="*80)
    
    print("And now checking the 'good' example (these should all be updated from the default values)\n")
    dummy_obj.stdout = good_example
    good_output = GalfitOutput(dummy_obj)
    
    print(good_output.bulge)
    print(good_output.disk)
    print(good_output.arms)
    print(good_output.fourier)
    print(good_output.sky)
    print(f"Did the good example succeed? {good_output.success}")
    
    print("="*80)
    print("Testing extraction into ComponentContainer...")
    _ = [print(str(comp)) for comp in good_output.extract_components().to_list()]
    
    good_output.header.to_file("good_output.in", good_output.bulge, good_output.disk, good_output.arms, good_output.fourier, good_output.sky)


# In[709]:


# For debugging purposes
def in_notebook():
    ip = get_ipython()
    
    if ip:
        return True
    else:
        return False


# In[710]:


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


# In[711]:


if __name__ == "__main__":
    export_to_py("galfit_objects")

