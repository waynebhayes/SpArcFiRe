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


# In[4]:


class BaseParameter():
    def __init__(self, value, **kwargs):
        
        self.value            = value
        self.fix_value        = kwargs.get("fix_value", 1)
        self.name             = kwargs.get("name", "")
        self.parameter_number = kwargs.get("parameter_number", "#")
        self.comment          = kwargs.get("comment", "")
        
        # Not sure if I'll need these but keeping them in for now
        self.component_name   = kwargs.get("component_name", "")
        self.component_number = kwargs.get("component_number", "")
        
# ==========================================================================================================

    # Formerly __str__
    # def __str__(self):
    #     return f"{self.value}"

    # Formerly __repr__
    def __str__(self):
        pre_comment = f"{self.parameter_number:>2}) {self.value:<11} {self.fix_value}"
        return f"{pre_comment:<23} # {self.comment}"
          
# ==========================================================================================================

    # def update_values(self, **kwargs):
    #     for key, value in kwargs.items():
    #         setattr(self, key, value)
        
# ==========================================================================================================

#     def from_pandas(self, input_df):
#         param_names  = [n.split(f"_{self.component_type}")[0] for n in input_df.columns]
#         param_values = input_df.iloc[0].values.astype(float)
#         new_param_dict = dict(zip(param_names, param_values))
        
#         pos = "position"
#         if pos in self.param_values:
#             new_param_dict[pos] = (new_param_dict[f"{pos}_x"], new_param_dict[f"{pos}_y"])
#             new_param_dict.pop(f"{pos}_x")
#             new_param_dict.pop(f"{pos}_y")
        
#         # No graceful way to do this...
#         # TODO: Can this be used for bending modes as well?
#         if self.component_type in ("fourier"):
#             f_modes = set([pn.split("_")[0] for pn in param_names])
#             a   = "amplitude"
#             pha = "phase_angle"
            
#             for mode in f_modes:
#                 new_param_dict[mode] = (new_param_dict[f"{mode}_{a}"], new_param_dict[f"{mode}_{pha}"])
#                 new_param_dict.pop(f"{mode}_{a}")
#                 new_param_dict.pop(f"{mode}_{pha}")
        
#         self.param_values.update(new_param_dict)
#         self.update_values()

# # ==========================================================================================================

#     def to_pandas(self):
#         name = f"{self.component_type}_{self.component_number}"
#         param_values = deepcopy(self.param_values)
        
        
#         if "Component type" in param_values:
#             param_values.pop("Component type")
            
#         multi_valued = [(k,v) for k,v in param_values.items() if isinstance(v, tuple)]
        
#         if multi_valued:
#             # For components which contain multi values, i.e. position
#             tuple_names = {"sersic"  : ("x","y"),
#                            "fourier" : ("amplitude", "phase_angle")}
            
#             # Usually only two but to keep things generic we can loop
#             for tup in multi_valued:
#                 tup_value = param_values.pop(tup[0])
                
#                 for i, val in enumerate(tup_value):
#                     param_values[f"{tup[0]}_{tuple_names[self.component_type][i]}"] = val
        
#         param_names  = [f"{i}_{name}" for i in param_values.keys()]
#         param_values = np.array(list(param_values.values()))
#         param_values = param_values.reshape(1, len(param_values))
        
#         # Redundancy in index name and column names is intentional!
#         all_data = pd.DataFrame(param_values, 
#                                 index = [name], 
#                                 columns = param_names,
#                                 dtype = np.float32)
        
#         return all_data

# ==========================================================================================================


# In[5]:


class HeaderParameter(BaseParameter, str):
    def __new__(cls, value = "", **kwargs):    
        return super(HeaderParameter, cls).__new__(cls, value)
    
    def __init__(self, value, **kwargs):
        BaseParameter.__init__(self, value, **kwargs)
    
        self.fix_value = ""


# In[6]:


class NumParameter(BaseParameter, float):
    def __new__(cls, value = 0, **kwargs):    
        return super(NumParameter, cls).__new__(cls, value)
    
    def __init__(self, value, **kwargs):
        
        BaseParameter.__init__(self, value, **kwargs)
        
        self.value = float(round(self.value, 4))
        
#         self.fix_value        = kwargs.get("fix_value", 0)
#         self.name             = kwargs.get("name", "")
#         self.parameter_number = kwargs.get("parameter_number", "#")
#         self.comment          = kwargs.get("comment", "")
        
#         # Not sure if I'll need these but keeping them in for now
#         self.component_name   = kwargs.get("component_name", "")
#         self.component_number = kwargs.get("component_number", "")
        
# ==========================================================================================================

    # Formerly __str__
    # def __str__(self):
    #     return f"{self.value}"

    # Formerly __repr__
    # def __str__(self):
    #     return f"{self.parameter_number:>2}) {self.value:<11.4f} {self.fix_value:<10d} #  {self.comment}"
    
# ==========================================================================================================


# In[7]:


class ComponentType(BaseParameter, str):
    def __new__(cls, name, **kwargs):    
        return super(ComponentType, cls).__new__(cls, name)
    
    def __init__(self, name, **kwargs):
        
        BaseParameter.__init__(self, name, **kwargs)
        self.fix_value = ""
        self.parameter_number = kwargs.get("parameter_number", 0)
        
        # self.value            = value
        self.comment          = "Component type"
        # self.component_number = kwargs.get("component_number", 0)
        
    # def __str__(self):
    #     empty_str = ""
    #     return f"{self.parameter_number:>2}) {self.name:<11} {empty_str:<10} #  {self.comment}"


# In[57]:


class MultiParameter(BaseParameter):
    def __init__(self, value = (0,0), **kwargs):
        
        BaseParameter.__init__(self, 0, **kwargs)
        
        # Generically use x and y
        if len(value) == 2:
            self.x = float(kwargs.get("x", value[0]))
            x_str = f"{self.x:.4f}"
            
            self.y = float(kwargs.get("y", value[1]))
            y_str = f"{self.y:.4f}"

            if self.x == int(self.x):
                x_str = f"{self.x:.0f}"
                
            if self.y == int(self.y):
                y_str = f"{self.y:.0f}"
                
            self.value_str  = f"{x_str} {y_str}"
            self.value = (self.x, self.y)
            
        elif len(value) == 4:
            self.xmin = int(kwargs.get("xmin", value[0]))
            self.xmax = int(kwargs.get("xmax", value[1]))
            self.ymin = int(kwargs.get("ymin", value[2]))
            self.ymax = int(kwargs.get("ymax", value[3]))

            self.value_str  = f"{self.xmin}   {self.xmax}   {self.ymin}   {self.ymax}"
            self.value = (self.xmin, self.xmax, self.ymin, self.ymax)
        
        self.fix_value        = ""
        self.fix_x            = ""
        self.fix_y            = ""
        
    # Override BaseParameter
    def __str__(self):
        pre_comment = f"{self.parameter_number:>2}) {self.value_str:<11} {self.fix_value}"
        return f"{pre_comment:<23} # {self.comment}"


# In[58]:


class Position(MultiParameter):
    def __init__(self, value = (0,0), **kwargs):
        
        MultiParameter.__init__(self, value, **kwargs)
    
        self.name             = "position"
        self.parameter_number = 1
        self.comment          = "Position x, y"
        self.fix_value        = kwargs.get("fix_value", "0 0")
        self.fix_x            = self.fix_value[0]
        self.fix_y            = self.fix_value[-1]


# In[93]:


class FourierMode(MultiParameter):
    def __init__(self, mode, amplitude = 0, phase_angle = 0, **kwargs):
        
        self.mode        = mode
        
        self.amplitude   = amplitude
        self.phase_angle = phase_angle
        
        if isinstance(amplitude, tuple):
            self.amplitude   = amplitude[0]
            self.phase_angle = amplitude[1]
            
        self.value_str  = f"{self.amplitude} {self.phase_angle}"
        self.value = (self.amplitude, self.phase_angle) #f"{value_str:<11}"
        MultiParameter.__init__(self, self.value, **kwargs)
        
        self.name             = "position"
        self.parameter_number = f"F{self.mode}"
        self.comment         = f"Azim. Fourier mode {self.mode}, amplitude, & phase angle"
        self.fix_value       = kwargs.get("fix_value", "1 1")
        self.fix_amplitude   = self.fix_value[0]
        self.fix_phase_angle = self.fix_value[-1]


# In[94]:


#class BendingModes(BaseParameter, float):
class BendingMode(NumParameter):
    def __init__(self, mode, amplitude, **kwargs):
        
        self.mode        = mode
        self.amplitude   = amplitude
        self.value       = amplitude
        
        NumParameter.__init__(self, self.value, **kwargs)
    
        self.name             = "bending mode"
        self.parameter_number = f"B{self.mode}"
        self.comment         = f"Bending mode {self.mode} amplitude"


# In[95]:


class ImageRegionToFit(HeaderParameter, MultiParameter):
    def __init__(self, value = (0, 256, 0, 256), **kwargs):
        
        # Header Parameter first so that we don't overwrite
        # multiparameter fix value
        HeaderParameter.__init__(self, value, **kwargs)
        MultiParameter.__init__(self, value, **kwargs)
        
        self.fix_value = ""
        self.parameter_number = "H"
        self.comment = "Image region to fit (xmin xmax ymin ymax)"

# Redundant because I may use different naming conventions elsewhere
class CropRegion(ImageRegionToFit):
    def __init__(self, *values, **kwargs):
        ImageRegionToFit.__init__(self, *values, **kwargs)


# In[96]:


if __name__ == "__main__":
    crop_region = ImageRegionToFit(
        (0, 100, 0, 100)
    )
    
    print(crop_region)
    
    crop_region = CropRegion(
        (45, 145, 45, 145)
    )
    
    print(crop_region)


# In[97]:


class ConvolutionBox(HeaderParameter, MultiParameter):
    def __init__(self, value = (50, 50), **kwargs):
        
        
        HeaderParameter.__init__(self, value, **kwargs)
        MultiParameter.__init__(self, value, **kwargs)
        
        self.fix_value = ""
        self.parameter_number = "I"
        self.comment = "Size of the convolution box (x y)"


# In[98]:


class PlateScale(HeaderParameter, MultiParameter):
    def __init__(self, value = (0.396, 0.396), **kwargs):
        HeaderParameter.__init__(self, value, **kwargs)
        MultiParameter.__init__(self, value, **kwargs)
        
        self.dx = self.x
        self.dy = self.y
        
        self.fix_value = ""
        self.parameter_number = "K"
        self.comment = "Plate scale (dx dy)   [arcsec per pixel]"


# In[99]:


if __name__ == "__main__":
    conv_box = ConvolutionBox(
        (100, 100)
    )
    
    print(conv_box)
    
    plate_scale = PlateScale(
        (0.396, 0.396)
    )
    
    print(plate_scale)


# In[100]:


# if __name__ == "__main__":
#     from RegTest.RegTest import *


# In[101]:


if __name__ == "__main__":
    bulge_line = ComponentType("sersic", component_number = 1)
    disk_line  = ComponentType("sersic", component_number = 2)
    arms_line  = ComponentType("power", parameter_number = "R0")
    sky_line   = ComponentType("sky", component_number = 3)

    print(bulge_line)
    print(disk_line)
    print(arms_line)
    print(sky_line)


# In[102]:


if __name__ == "__main__":
    position = Position(
        (100, 100),
        component_name = "Sersic",
        component_number = 1
    )
    
    print(position)
    
    position = Position(
        x = 101,
        y = 101,
        component_name = "Sersic",
        component_number = 1
    )
    
    print(position)


# In[103]:


if __name__ == "__main__":
    magnitude = NumParameter(
        16,
        name = "magnitude",
        parameter_number = 3,
        comment = "Integrated magnitude",
        component_name = "Sersic",
        component_number = 1
    )
    
    print(magnitude)
    print(magnitude + magnitude)


# In[107]:


if __name__ == "__main__":
    fourier1 = FourierMode(
        mode = 1,
        amplitude = 0.001, 
        phase_angle = 45,
        component_name = "Fourier"
    )
    
    print(fourier1)
    
    fourier3 = FourierMode(
        3,
        (0.002, 46),
        component_name = "Fourier"
    )
    
    print(fourier3)


# In[108]:


if __name__ == "__main__":
    bending2 = BendingMode(
        mode = 2,
        amplitude = 1.002
    )
    
    print(bending2)
    
    bending3 = BendingMode(
        mode = 3,
        amplitude = 1.003
    )
    
    print(bending3)


# In[109]:


# Parameters with defaults for Sersic profile
def load_default_sersic_parameters(component_number = None):
    sersic_line = ComponentType("sersic", component_number = component_number)
    
    position = Position(
        (100, 100),
        component_name = "Sersic",
        component_number = component_number
    )
    
    magnitude = NumParameter(
        16,
        name = "magnitude",
        parameter_number = 3,
        comment = "Integrated magnitude",
        component_name = "Sersic",
        component_number = component_number
    )
    
    effective_radius = NumParameter(
        10,
        name = "effective radius",
        parameter_number = 4,
        comment = "R_e (effective radius)   [pix]",
        component_name = "Sersic",
        component_number = component_number
    )
    
    sersic_index = NumParameter(
        1,
        name = "sersic index",
        parameter_number = 5,
        comment = "Sersic index n (de Vaucouleurs n=4)",
        component_name = "Sersic",
        component_number = component_number
    )
    
    axis_ratio = NumParameter(
        0.5,
        name = "axis ratio",
        parameter_number = 9,
        comment = "Axis ratio (b/a)",
        component_name = "Sersic",
        component_number = component_number
    )
    
    position_angle = NumParameter(
        90,
        name = "position_angle",
        parameter_number = 10,
        comment = "Position angle (PA) [deg: Up=0, Left=90]",
        component_name = "Sersic",
        component_number = component_number
    )
    
    skip = NumParameter(
        0,
        name = "skip",
        parameter_number = "Z",
        comment = "Skip this model in output image?  (yes=1, no=0)",
        component_name = "Sersic",
        component_number = component_number
    )
    
    loc = deepcopy(locals())
    loc.pop("component_number")
    
    return loc


# In[110]:


def load_default_power_parameters(component_number = None):
    
    param_prefix = "R"
    
    power_line = ComponentType("power", parameter_number = f"{param_prefix}0", component_number = component_number)
    
    inner_rad = NumParameter(
        0,
        name = "inner radius",
        parameter_number = f"{param_prefix}1",
        comment = "Spiral inner radius [pixels]",
        component_name = "Power",
        component_number = component_number
    )
    
    outer_rad = NumParameter(
        10,
        name = "outer radius",
        parameter_number = f"{param_prefix}2",
        comment = "Spiral outer radius [pixels]",
        component_name = "Power",
        component_number = component_number
    )
    
    cumul_rot = NumParameter(
        90,
        name = "cumulative rotation out",
        parameter_number = f"{param_prefix}3",
        comment = "Cumul. rotation out to outer radius [degrees]",
        component_name = "Power",
        component_number = component_number
    )
    
    powerlaw_index = NumParameter(
        0.5,
        name = "powerlaw index",
        parameter_number = f"{param_prefix}1",
        comment = "Asymptotic spiral powerlaw",
        component_name = "Power",
        component_number = component_number
    )
    
    inclination = NumParameter(
        0,
        name = "inclination",
        parameter_number = f"{param_prefix}9",
        comment = "Inclination to L.o.S. [degrees]",
        component_name = "Power",
        component_number = component_number
    )
    
    sky_position_angle = NumParameter(
        90,
        name = "sky_position_angle",
        parameter_number = f"{param_prefix}10",
        comment = "Sky position angle",
        component_name = "Power",
        component_number = component_number
    )
    
    # skip = Parameter(
    #     0,
    #     name = "skip",
    #     parameter_number = "Z",
    #     comment = "Skip this model in output image?  (yes=1, no=0)",
    #     component_name = "Power",
    #     component_number = component_number
    # )
    
    loc = deepcopy(locals())
    loc.pop("component_number")
    loc.pop("param_prefix")
    
    return loc


# In[111]:


def load_default_fourier_parameters(component_number = None):
    
    param_prefix = "F"
    
    fourier1 = FourierMode(
        mode = 1,
        amplitude = 0.05,
        phase_angle = 45,
        component_number = component_number
    )
    
    fourier3 = FourierMode(
        mode = 3,
        amplitude = 0.01,
        phase_angle = 25,
        component_number = component_number
    )
    
    skip = NumParameter(
        0,
        name = "skip",
        parameter_number = "Z",
        comment = "Skip this model in output image?  (yes=1, no=0)",
        component_name = "Fourier",
        component_number = component_number
    )
    
    loc = deepcopy(locals())
    loc.pop("component_number")
    loc.pop("param_prefix")
    
    return loc


# In[112]:


# Parameters with defaults for Sky profile
def load_default_sky_parameters(component_number = None):
    sky_line = ComponentType("sky", component_number = component_number)
    
    sky_background = NumParameter(
        1000,
        name = "sky background",
        parameter_number = 1,
        comment = "Sky background at center of fitting region [ADUs]",
        component_name = "Sky",
        component_number = component_number
    )
    
    dsky_dx = NumParameter(
        0,
        name = "dsky/dx",
        parameter_number = 2,
        comment = "dsky/dx (sky gradient in x)     [ADUs/pix]",
        component_name = "Sky",
        component_number = component_number
    )
    
    dsky_dy = NumParameter(
        0,
        name = "dsky/dy",
        parameter_number = 2,
        comment = "dsky/dy (sky gradient in y)     [ADUs/pix]",
        component_name = "Sky",
        component_number = component_number
    )
    
    loc = deepcopy(locals())
    loc.pop("component_number")
    
    return loc


# In[113]:


# Parameters with defaults for Sky profile
def load_default_header_parameters():
    
    input_image = HeaderParameter(
        "in.fits",
        name = "input image",
        parameter_number = "A",
        comment = "Input data image (FITS file)",
        component_name = "Header"
    )
    
    output_image = HeaderParameter(
        "out.fits",
        name = "output image",
        parameter_number = "B",
        comment = "Output data image block",
        component_name = "Header"
    )
    
    sigma_image = HeaderParameter(
        "none",
        name = "sigma image",
        parameter_number = "C",
        comment = "Sigma image name (made from data if blank or 'none')",
        component_name = "Header"
    )
    
    psf = HeaderParameter(
        "psf.fits",
        name = "PSF",
        parameter_number = "D",
        comment = "Input PSF image and (optional) diffusion kernel",
        component_name = "Header"
    )
    
    psf_fine_sampling = HeaderParameter(
        1,
        name = "PSF fine sampling factor",
        parameter_number = "E",
        comment = "PSF fine sampling factor relative to data",
        component_name = "Header"
    )
    
    pixel_mask = HeaderParameter(
        "none",
        name = "pixel mask",
        parameter_number = "F",
        comment = "Bad pixel mask (FITS image or ASCII coord list)",
        component_name = "Header"
    )
    
    constraint_file = HeaderParameter(
        "none",
        name = "constraint file",
        parameter_number = "G",
        comment = "File with parameter constraints (ASCII file)",
        component_name = "Header"
    )
    
    crop_region = ImageRegionToFit()
    
    conv_box    = ConvolutionBox()
    
    mag_photo_zeropoint = HeaderParameter(
        24.8,
        name = "Magnitude photometric zeropoint",
        parameter_number = "J",
        comment = "Magnitude photometric zeropoint",
        component_name = "Header"
    )
    
    plate_scale = PlateScale()
    
    display_type = HeaderParameter(
        "regular",
        name = "display type",
        parameter_number = "O",
        comment = "Display type (regular, curses, both)",
        component_name = "Header"
    )
    
    optimize = HeaderParameter(
        1,
        name = "optimize",
        parameter_number = "P",
        comment = "Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps",
        component_name = "Header"
    )
    
    loc = deepcopy(locals())
    
    return loc


# In[114]:


if __name__ == "__main__":
    _ = [print(v) for v in load_default_header_parameters().values()]
    print()
    _ = [print(v) for v in load_default_sersic_parameters().values()]
    print()
    _ = [print(v) for v in load_default_power_parameters().values()]
    print()
    _ = [print(v) for v in load_default_fourier_parameters().values()]
    print()
    _ = [print(v) for v in load_default_sky_parameters().values()]


# In[115]:


if __name__ == "__main__":
    _ = [print(type(v.value), v) for v in load_default_header_parameters().values()]
    print()
    _ = [print(type(v.value), v) for v in load_default_sersic_parameters().values()]
    print()
    _ = [print(type(v.value), v) for v in load_default_power_parameters().values()]
    print()
    _ = [print(type(v.value), v) for v in load_default_fourier_parameters().values()]
    print()
    _ = [print(type(v.value), v) for v in load_default_sky_parameters().values()]


# In[30]:


if __name__ == "__main__":
    export_to_py("Parameters", pj(_MODULE_DIR, "Classes", "Parameters"))

