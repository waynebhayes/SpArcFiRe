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

from collections import namedtuple

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
        
        # All underscored attributes are likely to change
        # And are thus relegated to @properties below for ease
        # of updating
        self._value            = value
        self._fix              = kwargs.get("fix", 1)
        
        self.name              = kwargs.get("name", "")
        self.parameter_number  = kwargs.get("parameter_number", "#")
        self.parameter_prefix  = kwargs.get("parameter_prefix", " ")
        self.comment           = kwargs.get("comment", "")
        
        # Not sure if I'll need these but keeping them in for now
        self.component_name    = kwargs.get("component_name", "")
        self.component_number  = kwargs.get("component_number", "")
          
# ==========================================================================================================

    exec(
        generate_get_set(
            {
                "value"     : "_value",
                "fix"       : "_fix"
            }
        )
    )
    
# ==========================================================================================================

    def __repr__(self):
        return repr(self.value)

    # Formerly __repr__
    def __str__(self):
        pre_fix = f"{self.parameter_prefix}{self.parameter_number}) {self.value:<10}"
        pre_comment = f"{pre_fix:<16}{self.fix}"
        return f"{pre_comment:<23} # {self.comment}"


# In[5]:


class ComponentType(BaseParameter, str):
    def __new__(cls, name, **kwargs):    
        return super(ComponentType, cls).__new__(cls, name)
    
    def __init__(self, name, **kwargs):
        
        BaseParameter.__init__(
            self, 
            name,
            fix              = "",
            parameter_number = 0,
            comment          = "Component type",
            **kwargs
        )
        
        #self.parameter_number = kwargs.get("parameter_number", 0)
        
        # self.value            = value
        #self.comment          = "Component type"
        # self.component_number = kwargs.get("component_number", 0)
        
    # def __str__(self):
    #     empty_str = ""
    #     return f"{self.parameter_number:>2}) {self.name:<11} {empty_str:<10} #  {self.comment}"


# In[6]:


class HeaderParameter(BaseParameter, str):
    def __new__(cls, value = "", **kwargs):    
        return super(HeaderParameter, cls).__new__(cls, value)
    
    def __init__(self, value, **kwargs):
        BaseParameter.__init__(
            self, 
            value, 
            fix              = "",
            parameter_prefix = "",
            **kwargs
        )


# In[7]:


class NumParameter(BaseParameter, float):
    def __new__(cls, value = 0, **kwargs):    
        return super(NumParameter, cls).__new__(cls, value)
    
    def __init__(self, value, **kwargs):
        
        BaseParameter.__init__(
            self, 
            float(round(value, 4)), 
            **kwargs
        )
        
        #self.value = float(round(self.value, 4))
        
#         self.fix        = kwargs.get("fix", 0)
#         self.name             = kwargs.get("name", "")
#         self.parameter_number = kwargs.get("parameter_number", "#")
#         self.comment          = kwargs.get("comment", "")
        
#         # Not sure if I'll need these but keeping them in for now
#         self.component_name   = kwargs.get("component_name", "")
#         self.component_number = kwargs.get("component_number", "")
        
# ==========================================================================================================

    @property
    def value(self):
        return self._value
    
    @value.setter
    def value(self, new_val):
        self._value = float(new_val)

    # Formerly __str__
    # def __str__(self):
    #     return f"{self.value}"

    # Formerly __repr__
    # def __str__(self):
    #     return f"{self.parameter_number:>2}) {self.value:<11.4f} {self.fix:<10d} #  {self.comment}"
    
# ==========================================================================================================


# In[8]:


class Skip(BaseParameter, int):
    def __new__(cls, value = 0, **kwargs):    
        return super(Skip, cls).__new__(cls, value)
    
    def __init__(self, value = 0, **kwargs):
        
        #self._value = value
        
        # In case it is accidentally given as a string such as 5.0
        try:
            value = int(value)
        except ValueError as ve:
            try:
                value = int(float(value))
            except ValueError as ve2:
                if "could not convert string to float" in ve2:
                    print("For Skip class:", ve2)
                    print("Resorting to default.")
        
        name             = "skip"
        parameter_number = f"Z"
        comment          = f"Skip this model in output image?  (yes=1, no=0)"
        fix              = ""
        
        BaseParameter.__init__(
            self, 
            value,
            name             = name,
            parameter_number = parameter_number,
            comment          = comment,
            fix              = fix
        )
        
    @property
    def value(self):
        return self._value
    
    @value.setter
    def value(self, new_val):
        try:
            self._value = int(new_val)
        except ValueError as ve:
            try:
                self._value = int(float(new_val))
            except ValueError as ve:
                raise Exception("Value given to Skip must be convertible to an int or a float.")
                
    @property
    def fix(self):
        return ""
    
    @fix.setter
    def fix(self, new_val):
        pass


# In[9]:


#class BendingModes(BaseParameter, float):
class BendingMode(NumParameter):
    def __init__(self, mode, amplitude, **kwargs):
        
        self._mode        = mode
        # self._amplitude   = amplitude
        
        NumParameter.__init__(
            self, 
            amplitude,
            name = "bending mode",
            parameter_number = f"{self.mode}",
            parameter_prefix = "B",
            comment = f"Bending mode {self.mode} amplitude",
            **kwargs)
    
# ==========================================================================================================
    
    exec(
        generate_get_set(
            {
                "amplitude" : "value",
                "mode" : "_mode"
            }
        )
    )
    
#     @property
#     def amplitude(self):
#         return self.value
    
#     @amplitude.setter
#     def amplitude(self, new_val):
#         self.value = new_val
        
#     @property
#     def mode(self):
#         return self._mode
    
#     @mode.setter
#     def mode(self, new_val):
#         self._mode = new_val


# In[10]:


# Use NamedTuples so we can access this like a dictionary
# while also using the indexing of a tuple where necessary

ntMultiParameter   = namedtuple("ntMultiParameter", "x y")
ntMultiParameterFix = namedtuple("ntMultiParameterFix", "fix_x fix_y")

class MultiParameter(BaseParameter):
    def __init__(self, value = (0,0), **kwargs):
        
        BaseParameter.__init__(self, 0, **kwargs)
        
        # Generically use x and y
        
        self._x = float(kwargs.get("x", value[0]))
        self._y = float(kwargs.get("y", value[1]))

#        self._value = ntMultiParameter(self._x, self._y)

        self._fix_x = kwargs.get("fix_x", "")
        self._fix_y = kwargs.get("fix_y", "")
    
    def check_type(self, iterable_in):
        assert isinstance(iterable_in, (list, tuple)), "The new value for this attribute must be a list/tuple!"
        
    def check_len(self, iterable_in, len_to_check = 2):
        assert len(iterable_in) == len_to_check, "The length of this iterable is not 2, MultiParameter cannot correctly update."
# ==========================================================================================================

    exec(
        generate_get_set(
            {
                "x" : "_x",
                "y" : "_y"
            }
        )
    )
        
# ==========================================================================================================

    @property
    def value(self):
        #self.check_len(self._value)
        return ntMultiParameter(self.x, self.y)
    
    @value.setter
    def value(self, new_tuple):
        self.check_type(new_tuple)
        self.check_len(new_tuple)
          
        self.x = float(new_tuple[0])
        self.y = float(new_tuple[1])
        
#        self._value = ntMultiParameter(self.x, self.y)

# ==========================================================================================================
        
    exec(
        generate_get_set(
            {
                "fix_x" : "_fix_x",
                "fix_y" : "_fix_y"
            }
        )
    )
    
    @property
    def fix(self):
        return ntMultiParameterFix(self.fix_x, self.fix_y)
    
    @fix.setter
    def fix(self, new_val):
        if isinstance(new_val, str):
            if len(new_val):
                self.fix_x = int(new_val[0])
                self.fix_y = int(new_val[-1])
            else:
                self.fix_x = ""
                self.fix_y = ""
            
        elif isinstance(new_val, (list, tuple)):
            if isinstance(new_val[0], (int, float)):
                self.fix_x = int(new_val[0])
                self.fix_y = int(new_val[1])
                
            else:
                self.fix_x = ""
                self.fix_y = ""
        else:
            raise Exception("The new value for this attribute must be a string, list, or tuple!")
        
# ==========================================================================================================

    def __repr__(self):
        return repr(self.value)

    # Override BaseParameter
    def __str__(self):

        x_str = f"{self.x:.4f}"
        y_str = f"{self.y:.4f}"

        if self.x == int(self.x):
            x_str = f"{self.x:.0f}"

        if self.y == int(self.y):
            y_str = f"{self.y:.0f}"

        value_str  = f"{x_str} {y_str}"
        
        # Check if 0 or 1
        # string likely indicates header parameter
        if isinstance(self.x, (int, float)):
            str_fix = "  ".join(str(i) for i in self.fix)
        else:
            str_fix = ""
            
        pre_comment = f"{self.parameter_prefix}{self.parameter_number}) {value_str:<11} {str_fix}"
        return f"{pre_comment:<23} # {self.comment}"


# In[11]:


#ntPosition    = namedtuple("ntPosition", "x y")

class Position(MultiParameter):
    def __init__(self, value = (0,0), **kwargs):
        
        MultiParameter.__init__(self, value, **kwargs)
    
        #self.value            = ntPosition(*self.value)
        
        self.name             = "position"
        self.parameter_number = 1
        self.comment          = "Position x, y"
        
        self.fix              = kwargs.get("fix", (0, 0))
        


# In[12]:


ntFourier = namedtuple("ntFourier", "amplitude phase_angle")

class FourierMode(MultiParameter):
    def __init__(self, mode, amplitude = 0, phase_angle = 0, **kwargs):
        
        self._mode      = mode
        
        if isinstance(amplitude, (list, tuple)):
            phase_angle = amplitude[1]
            amplitude   = amplitude[0]
            
        # Override these properties later
        #self.amplitude   = self._amplitude
        #self.phase_angle = self._phase_angle
        
        #self.value_str  = f"{self.amplitude} {self.phase_angle}"
        
        MultiParameter.__init__(
            self, 
            (amplitude, phase_angle),
            name = "fourier mode",
            parameter_number = f"{self.mode}",
            parameter_prefix = "F",
            comment = f"Azim. Fourier mode {self.mode}, amplitude, & phase angle",
            **kwargs
        )
        
        #self._value     = ntFourier(*self.value)
        
        # Use str instead of namedtuple for printing
        self.fix       = kwargs.get("fix", (1, 1))
        
        self.fix_x     = kwargs.get("fix_amplitude"  , self.fix[0])
        self.fix_y     = kwargs.get("fix_phase_angle", self.fix[1])
# ==========================================================================================================
        
    exec(
        generate_get_set(
            {
                "amplitude"   : "x",
                "phase_angle" : "y"
            }
        )
    )
    
# ==========================================================================================================   

    exec(
        generate_get_set(
            {
                "mode"            : "_mode",
                "fix_amplitude"   : "fix_x",
                "fix_phase_angle" : "fix_y"
            }
        )
    )

# ==========================================================================================================   

    @property
    def value(self):
        #self._value = ntFourier(self.amplitude, self.phase_angle)
        return ntFourier(self.amplitude, self.phase_angle)
    
    @value.setter
    def value(self, new_tuple):
        self.check_type(new_tuple)
        
        self.amplitude   = float(new_tuple[0])
        self.phase_angle = float(new_tuple[1])


# In[13]:


ntImageRegionToFit = namedtuple("ntImageRegionToFit", "x1 x2 y1 y2")

class ImageRegionToFit(MultiParameter, HeaderParameter):
    def __init__(self, value = (0, 256, 0, 256), **kwargs):
        
        # Declare MultiParameter even though we're pretty much
        # overriding everything but the explicit notion that it is 
        # multi valued comes in handy
        HeaderParameter.__init__(self, value, **kwargs)
        MultiParameter.__init__(self, value, **kwargs)
        
        #self.value = ntImageRegionToFit(*value)
        self.parameter_prefix = ""
        self.parameter_number = "H"
        self.comment          = "Image region to fit (xmin xmax ymin ymax)"
        
        self._x1 = int(kwargs.get("xmin", value[0]))
        self._x2 = int(kwargs.get("xmax", value[1]))
        self._y1 = int(kwargs.get("ymin", value[2]))
        self._y2 = int(kwargs.get("ymax", value[3]))
        
        # self._value = (self._x1, self._x2, self._y1, self._y2)
        
    # def check_type(self, iterable_in):
    #     assert isinstance(iterable_in, (list, tuple)), "The new value for this attribute must be a list/tuple!"
        
    # def check_len(self, iterable_in):
    #     assert len(iterable_in) == 4, "The length of this iterable is not 4, ImageRegionToFit cannot correctly update."

    exec(
        generate_get_set(
            {
                "x1" : "_x1",
                "x2" : "_x2",
                "y1" : "_y1",
                "y2" : "_y2"
            }
        )
    )
    
    exec(
        generate_get_set(
            {
                "xmin" : "_x1",
                "xmax" : "_x2",
                "ymin" : "_y1",
                "ymax" : "_y2"
            }
        )
    )
    
    @property
    def value(self):
        return ntImageRegionToFit(self.x1, self.x2, self.y1, self.y2)
    
    @value.setter
    def value(self, new_tuple):
        self.check_type(new_tuple)
        self.check_len(new_tuple, 4)
            
        self.x1     = int(new_tuple[0])
        self.x2     = int(new_tuple[1])
        self.y1     = int(new_tuple[2])
        self.y2     = int(new_tuple[3])
    
    def __str__(self):
        value_str = f"{self.x1:<5}{self.x2:<5}{self.y1:<5}{self.y2:<5}"
        
        pre_comment = f"{self.parameter_prefix}{self.parameter_number}) {value_str}" #{self.fix}"
        return f"{pre_comment:<23} # {self.comment}"
        
# Redundant because I may use different naming conventions elsewhere
class CropRegion(ImageRegionToFit):
    def __init__(self, *values, **kwargs):
        ImageRegionToFit.__init__(self, *values, **kwargs)


# In[14]:


# This is functionally the same as position but for clarity...
ntConvolutionBox = namedtuple("ntConvolutionBox", "x y")

class ConvolutionBox(HeaderParameter, MultiParameter):
    def __init__(self, value = (50, 50), **kwargs):
        
        
        HeaderParameter.__init__(self, value, **kwargs)
        MultiParameter.__init__(self, value, **kwargs)
        
        self.value            = ntConvolutionBox(*self.value)
        
        self.parameter_prefix = ""
        self.parameter_number = "I"
        self.comment          = "Size of the convolution box (x y)"
        
# ==========================================================================================================

    @property
    def value(self):
        return ntConvolutionBox(self.x, self.y)
    
    @value.setter
    def value(self, new_tuple):
        self.check_type(new_tuple)
        self.check_len(new_tuple)
            
        self.x = int(new_tuple[0])
        self.y = int(new_tuple[1])
        
# ==========================================================================================================

    def __str__(self):
        value_str = f"{self.x:<7}{self.y}"
        
        pre_comment = f"{self.parameter_prefix}{self.parameter_number}) {value_str}" #{self.fix}"
        return f"{pre_comment:<23} # {self.comment}"


# In[15]:


ntPlateScale = namedtuple("ntPlateScale", "dx dy")

class PlateScale(HeaderParameter, MultiParameter):
    def __init__(self, value = (0.396, 0.396), **kwargs):
        HeaderParameter.__init__(self, value, **kwargs)
        MultiParameter.__init__(self, value, **kwargs)
        
        #self._dx = self._x
        #self._dy = self._y
        #self.value = ntPlateScale(*self.value)
        self.parameter_prefix = ""
        self.parameter_number = "K"
        self.comment = "Plate scale (dx dy)   [arcsec per pixel]"
        
    exec(
        generate_get_set(
            {
                "dx" : "x",
                "dy" : "y"
            }
        )
    )
    
    @property
    def value(self):
        #self._value = ntPlateScale(self.dx, self.dy)
        return ntPlateScale(self.dx, self.dy)
    
    @value.setter
    def value(self, new_tuple):
        self.check_type(new_tuple)
            
        self.dx     = float(new_tuple[0])
        self.dy     = float(new_tuple[1])
        
# ==========================================================================================================

    def __str__(self):
        value_str = f"{self.x:<7}{self.y}"
        
        pre_comment = f"{self.parameter_prefix}{self.parameter_number}) {value_str}" #{self.fix}"
        return f"{pre_comment:<23} # {self.comment}"


# In[16]:


# if __name__ == "__main__":
#     from RegTest.RegTest import *


# In[17]:


if __name__ == "__main__":
    crop_region = ImageRegionToFit(
        (0, 100, 0, 100)
    )
    
    print(crop_region)
    print(repr(crop_region))
    print()
    
    crop_region = CropRegion(
        (45, 145, 45, 145)
    )
    
    print(crop_region)
    print(repr(crop_region))
    print()
    
    crop_region.xmin = 1
    crop_region.xmax = 155
    crop_region.ymin = 1
    crop_region.ymax = 155
    
    print(crop_region)
    print(repr(crop_region))


# In[18]:


if __name__ == "__main__":
    conv_box = ConvolutionBox(
        (100, 100)
    )
    
    print(conv_box)
    print(repr(conv_box))
    print()
    
    plate_scale = PlateScale(
        (0.396, 0.396)
    )
    
    print(plate_scale)
    print(repr(plate_scale))
    print()
    
    plate_scale.x = 0.4
    plate_scale.y = 0.4
    
    print(plate_scale)
    print(repr(plate_scale))
    print()


# In[19]:


if __name__ == "__main__":
    bulge_line = ComponentType("sersic", component_number = 1)
    disk_line  = ComponentType("sersic", component_number = 2)
    arms_line  = ComponentType("power" , parameter_prefix = "R")
    sky_line   = ComponentType("sky"   , component_number = 3)

    print(bulge_line)
    print(disk_line)
    print(arms_line)
    print(sky_line)


# In[20]:


if __name__ == "__main__":
    position = Position(
        (100, 100),
        component_name = "Sersic",
        component_number = 1
    )
    
    print(position)
    print(repr(position))
    print()
    
    position.x = 101
    position.y = 101
    print(position)
    print(repr(position))
    print()
    
    position.value = (102, 102)
    print(position)
    print(repr(position))


# In[21]:


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
    print("sum of magnitudes", magnitude + magnitude)
    
    magnitude.value = 5
    print(magnitude)


# In[22]:


if __name__ == "__main__":
    skip = Skip(
        component_name = "Sersic",
        component_number = 1
    )
    
    print(skip)
    skip.value = 1
    print(skip)


# In[23]:


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
    
    fourier1.amplitude   = 0.003
    fourier1.phase_angle = 47
    
    print(fourier1)


# In[ ]:


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
    
    bending3.mode = 4
    bending3.amplitude = 1.004
    
    print(bending3)


# In[ ]:


# Parameters with defaults for Sersic profile
def load_default_sersic_parameters(component_number = None):
    _sersic = ComponentType("sersic", component_number = component_number)
    
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
    
    skip = Skip(
        0,
        component_name = "Sersic",
        component_number = component_number
    )
    
    loc = deepcopy(locals())
    loc.pop("component_number")
    
    return loc


# In[ ]:


def load_default_power_parameters(component_number = None):
    
    param_prefix = "R"
    
    _power = ComponentType(
        "power", 
        parameter_prefix = param_prefix,
        component_number = component_number
    )
    
    inner_rad = NumParameter(
        0,
        name = "inner radius",
        parameter_number = "1",
        parameter_prefix = param_prefix,
        comment = "Spiral inner radius [pixels]",
        component_name = "Power",
        component_number = component_number
    )
    
    outer_rad = NumParameter(
        10,
        name = "outer radius",
        parameter_number = "2",
        parameter_prefix = param_prefix,
        comment = "Spiral outer radius [pixels]",
        component_name = "Power",
        component_number = component_number
    )
    
    cumul_rot = NumParameter(
        90,
        name = "cumulative rotation out",
        parameter_number = "3",
        parameter_prefix = param_prefix,
        comment = "Cumul. rotation out to outer radius [degrees]",
        component_name = "Power",
        component_number = component_number
    )
    
    powerlaw_index = NumParameter(
        0.5,
        name = "powerlaw index",
        parameter_number = "4",
        parameter_prefix = param_prefix,
        comment = "Asymptotic spiral powerlaw",
        component_name = "Power",
        component_number = component_number
    )
    
    inclination = NumParameter(
        0,
        name = "inclination",
        parameter_number = "9",
        parameter_prefix = param_prefix,
        comment = "Inclination to L.o.S. [degrees]",
        component_name = "Power",
        component_number = component_number
    )
    
    sky_position_angle = NumParameter(
        90,
        name = "sky_position_angle",
        parameter_number = "10",
        parameter_prefix = param_prefix,
        comment = "Sky position angle",
        component_name = "Power",
        component_number = component_number
    )
    
    # skip = Skip(
    #     0,
    #     component_name = "Power",
    #     component_number = component_number
    # )
    
    loc = deepcopy(locals())
    loc.pop("component_number")
    loc.pop("param_prefix")
    
    return loc


# In[ ]:


def load_default_fourier_parameters(component_number = None):
    
    F1 = FourierMode(
        mode = 1,
        amplitude = 0.05,
        phase_angle = 45,
        component_number = component_number
    )
    
    F3 = FourierMode(
        mode = 3,
        amplitude = 0.01,
        phase_angle = 25,
        component_number = component_number
    )
    
    skip = Skip(
        0,
        component_name = "Fourier",
        component_number = component_number
    )
    
    loc = deepcopy(locals())
    loc.pop("component_number")
    
    return loc


# In[ ]:


# Parameters with defaults for Sky profile
def load_default_sky_parameters(component_number = None):
    _sky = ComponentType("sky", component_number = component_number)
    
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
        parameter_number = 3,
        comment = "dsky/dy (sky gradient in y)     [ADUs/pix]",
        component_name = "Sky",
        component_number = component_number
    )
    
    skip = Skip(
        0,
        component_name = "Sky",
        component_number = component_number
    )
    
    loc = deepcopy(locals())
    loc.pop("component_number")
    
    return loc


# In[ ]:


# Parameters with defaults for Sky profile
def load_default_header_parameters(galaxy_name = ""):
    
    input_image = HeaderParameter(
        f"{galaxy_name}.fits",
        name = "input image",
        parameter_number = "A",
        comment = "Input data image (FITS file)",
        component_name = "Header"
    )
    
    output_image = HeaderParameter(
        f"{galaxy_name}_galfit_out.fits",
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
        f"{galaxy_name}_psf.fits",
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
        f"{galaxy_name}_star-rm.fits",
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
    
    region_to_fit = ImageRegionToFit()
    
    convolution_box    = ConvolutionBox()
    
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
        0,
        name = "optimize",
        parameter_number = "P",
        comment = "Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps",
        component_name = "Header"
    )
    
    loc = deepcopy(locals())
    loc.pop("galaxy_name")
    
    return loc


# In[ ]:


def load_default_parameters():
    return {
        "header"  : load_default_header_parameters(),
        "sersic"  : load_default_sersic_parameters(),
        "power"   : load_default_power_parameters(),
        "fourier" : load_default_fourier_parameters(),
        "sky"     : load_default_sky_parameters()
        #"bending"
    }


# In[ ]:


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


# In[ ]:


if __name__ == "__main__":
    _ = [
        print(*[parameter for parameter in component_dict.values()], sep = "\n", end = "\n\n")
        for component_dict in load_default_parameters().values()
    ]


# In[ ]:


if __name__ == "__main__":
    export_to_py("Parameters", pj(_MODULE_DIR, "Classes", "Parameters"))

