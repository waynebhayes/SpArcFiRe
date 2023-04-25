#!/usr/bin/env python
# coding: utf-8

# In[303]:


import os
import sys
from os.path import join as pj
from os.path import exists
import subprocess
from copy import deepcopy
from IPython import get_ipython


# In[304]:


class GalfitComponent:
    def __init__(self, 
                 component_name = "", 
                 component_number = 0, 
                 param_prefix = " "
                ):
        
        self.component_name = component_name
        self.component_number = component_number
        self.param_prefix = param_prefix
        
        key = "Component type"
        self.param_desc = {key : key}
        
        # tuple in value for printing
        self.param_numbers = {0 : key}
        self.param_values = {key : component_name} #(component_name, "")}
        self.param_fix = {key : ""}
        
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
        
        if self.component_name == "header":
            output_str += "\n\n".join((self.input_menu_file, self.extra_header_info, "="*80))
            output_str += "\n# IMAGE and GALFIT CONTROL PARAMETERS\n"
            num_type = "" #".0f"
            
        elif self.component_name not in ("power", "fourier", "header"):
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
                    
                # # For a bit of spacing in R10
                # This does not play nice with galfit *smacks forehead*
                # if len(line.split()[0]) >= 4: 
                #     line = line[:4] + line[5:]
                
            else:
                # position #, fourier, etc.
                if isinstance(val[0], (int, float)):
                    line = f"{self.param_prefix}{num}) {val[0]:<7{num_type}} {val[1]:{num_type}} {fix}"
                    
                # comp name
                else:
                    line = f"{self.param_prefix}{num}) {val}" #[0]}"
                
            output_str += f"{line:<23} # {desc}\n"
            
        if self.component_name == "header":
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
    
    def from_file(self, filename):
        # This function handles grabbing and storing the values from galfit files (input and output???)
        # It's written to generally handle both and stores everything in a dict of dicts
        # See parameter_names list for the corresponding labels
        # {Component: {Parameter : Value}}
    
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

        input_in = input_file.readlines()
        input_file.close()

        # Need component number for dual sersic
        try:
            c_num = self.component_number
        except:
            # For Power or Fourier
            c_num = 2

        component_start = False
        function_start  = False
        
        for line in input_in:
            # Empty
            if not line.strip():
                component_start = False
                function_start  = False
                continue
                
            c_str = f"# Component number: {c_num}"
            if line.startswith(c_str):
                component_start = True
                continue

            # This handy trick replaces all multi-spaces with a single space
            # Hence the weird syntax
            # Somewhere on stack overflow...
            # " ".join(line.split())
            try:
                split_line = " ".join(line.strip().split()).split()
            except ValueError:
                # Empty
                component_start = False
                function_start  = False
                continue
            
            line_num = split_line[0].strip(f"{self.param_prefix})")

            # Don't need the name since this is used independently by each component class
            if line_num == "0" and self.component_name == split_line[1].strip():
                # assert self.component_name == split_line[1].strip(), "The component name does not match the line " \ 
                #                                                      "being read in. Something went wrong!"
                function_start = True
                continue

            if self.component_name == "header" and "IMAGE and GALFIT CONTROL PARAMETERS" in line:
                component_start = True
                function_start  = True
                continue
            
            if not all((component_start, function_start)):
                continue

            # Special cases...
            if line_num == "1" and self.component_name == "sersic":
                line_num = int(line_num)
                value = (float(split_line[1]), float(split_line[2]))
                fix   = (int(split_line[3]), int(split_line[4]))

            elif line_num in ["1", "3"] and self.component_name == "fourier":
                line_num = int(line_num)
                value = (float(split_line[1]), float(split_line[2]))
                fix   = (int(split_line[3]), int(split_line[4]))
                
            elif self.component_name == "header":
                value = split_line[1]
                fix = ""
                if line_num in ("E", "J", "P"):
                    value = float(value)
                
                elif line_num == "H":
                    value = (int(split_line[1]), int(split_line[2]))
                    fix = (int(split_line[3]), int(split_line[4]))
                    
                elif line_num in ("I", "K"):
                    value = (float(split_line[1]), float(split_line[2]))

            else:
                line_num = int(line_num)
                value = float(split_line[1])
                fix   = int(split_line[2])

            key = self.param_numbers[line_num]
            self.param_values[key] = value
            self.param_fix[key] = fix

        self.update_values(**self.param_values)

# ==========================================================================================================

    def to_file(self, filename, *args):
        try:
            with open(filename, "w") as f:
                f.write("\n")
                f.write(str(self))
                f.write("\n")

                # *args for writing in additional classes at the same time (save I/O)
                comp_names = [c.component_name for c in args]
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


# In[305]:


class Sersic(GalfitComponent):
    def __init__(self, component_number, **kwargs):
        
        self.position = kwargs.get("position", (0.0000, 0.0000))
        self.magnitude = kwargs.get("magnitude", 13)
        self.effective_radius = kwargs.get("effective_radius", 10)
        self.sersic_index = kwargs.get("sersic_index", 4)
        self.axis_ratio = kwargs.get("axis_ratio", 0.6)
        self.position_angle = kwargs.get("position_angle", 0)
        #self.skip = kwargs.get("skip", 0)
        
        param_dict = deepcopy(vars(self))
        
        GalfitComponent.__init__(self, component_name = "sersic", component_number = component_number)
        
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


# In[306]:


class Power(GalfitComponent):
    def __init__(self, **kwargs):
        
        #self.component_type = "power"
        self.inner_rad = kwargs.get("inner_rad", 0)
        self.outer_rad = kwargs.get("outer_rad", 20)
        self.cumul_rot = kwargs.get("cumul_rot", 90)
        self.powerlaw = kwargs.get("powerlaw", 1)
        self.inclination = kwargs.get("inclination", 15)
        self.sky_position_angle = kwargs.get("sky_position_angle", 45)
        #self.skip = kwargs.get("skip", 0)
        
        param_dict = deepcopy(vars(self))
        
        GalfitComponent.__init__(self, component_name = "power", param_prefix = "R")
        
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


# In[307]:


class Fourier(GalfitComponent):
    def __init__(self, n = {1 : (0.05, 45), 3 : (0.05, 25)}):
        GalfitComponent.__init__(self, component_name = "fourier", param_prefix = "F")      
        # normal rules don't apply here
        # Still use inheritance for the other functions
        self.param_numbers = {}
        self.param_values = {}
        self.param_desc = {}
        self.param_fix = {}
        
        def include_fn(self, n:dict):

            for num, values in n.items():
                key = f"{self.param_prefix}{num}"
                self.param_numbers[num] = key
                self.param_values[key] = values
                self.param_desc[key] = f"Azim. Fourier mode {num}, amplitude, & phase angle"
                self.param_fix[key] = "1 1"
        
        include_fn(self, n = n)
        
    def update_from_log(self, in_line):
        # example
        # fourier : (1:  0.06,   -6.67)   (3:  0.05,    0.18)
        
        # rstrip avoids a hanging ) later
        params = in_line.lstrip("fourier : ").replace(" ", "").rstrip(")").split(")(")
        
        self.param_values = {n: eval(f"({params[i].split(':')[1].replace('*', '')})")
                             for i, n in enumerate(self.param_values.keys())}


# In[308]:


class Sky(GalfitComponent):
    def __init__(self, component_number, **kwargs):
        self.sky_background = kwargs.get("sky_background", 1000)
        self.dsky_dx = kwargs.get("dsky_dx", 0)
        self.dsky_dy = kwargs.get("dsky_dy", 0)
        
        param_dict = deepcopy(vars(self))
        
        GalfitComponent.__init__(self, component_name = "sky", component_number = component_number)
        
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


# In[309]:


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
        GalfitComponent.__init__(self, component_name = "header", param_prefix = "")

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


# In[310]:


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


# In[311]:


class FeedmeContainer(ComponentContainer):
    def __init__(self, **kwargs):
        ComponentContainer.__init__(self, **kwargs)
        # The path to the feedme that *generated* the components
        self.header  = kwargs.get("header", GalfitHeader())
        self.path_to_feedme = kwargs.get("path_to_feedme", "")
    
    def to_dict(self):
        return vars(self)
    
    def __str__(self):
        out_str = f"{str(self.header)}\n" + "\n".join(str(comp) for comp in ComponentContainer.to_list(self))
        return out_str
        
    def to_file(self, *args):
        if args:
            self.header.to_file(self.path_to_feedme, *args)
        else:
            self.header.to_file(self.path_to_feedme, *ComponentContainer.to_list(self))
        
    # TODO
    def from_file(self):
        pass


# In[312]:


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


# In[313]:


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
    
    print(container)
    
    container.to_file("tester.in") #, bulge, disk, arms, fourier, sky)


# In[316]:


# Testing from_file
if __name__ == "__main__":
    
    header = GalfitHeader(galaxy_name = "fake_name")
    container.update_components(header = header)
    # bulge = Sersic(1)
    # disk  = Sersic(2)
    # arms  = Power()
    # fourier = Fourier()
    # sky   = Sky(3)
    
    
    for comp in (header, bulge, disk, arms, fourier, sky):
        print(comp.component_name)
        comp.from_file("good_output.in")
    
    print(str(container))


# In[70]:


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
    


# In[56]:


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


# In[76]:


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


# In[15]:


# For debugging purposes
def in_notebook():
    ip = get_ipython()
    
    if ip:
        return True
    else:
        return False


# In[16]:


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


# In[ ]:


if __name__ == "__main__":
    export_to_py("galfit_objects")

