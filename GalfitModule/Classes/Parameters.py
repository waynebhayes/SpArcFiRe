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
        self._fix_value        = kwargs.get("fix_value", 1)
        
        self.name              = kwargs.get("name", "")
        self.parameter_number  = kwargs.get("parameter_number", "#")
        self.comment           = kwargs.get("comment", "")
        
        # Not sure if I'll need these but keeping them in for now
        self.component_name    = kwargs.get("component_name", "")
        self.component_number  = kwargs.get("component_number", "")
        
# ==========================================================================================================

    # Formerly __str__
    # def __repr__(self):
    #     return f"{self.value}"

    # Formerly __repr__
    def __str__(self):
        pre_comment = f"{self.parameter_number:>2}) {self.value:<11} {self.fix_value}"
        return f"{pre_comment:<23} # {self.comment}"
          
# ==========================================================================================================

    exec(
        generate_get_set(
            {
                "value"     : "_value",
                "fix_value" : "_fix_value"
            }
        )
    )
    
#     @property
#     def fix_value(self):
#         return self._fix_value
    
#     @fix_value.setter
#     def fix_value(self, new_val):
#         self._fix_value = new_val


# In[5]:


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


# In[6]:


class HeaderParameter(BaseParameter, str):
    def __new__(cls, value = "", **kwargs):    
        return super(HeaderParameter, cls).__new__(cls, value)
    
    def __init__(self, value, **kwargs):
        BaseParameter.__init__(self, value, **kwargs)
    
        self.fix_value = ""


# In[33]:


class NumParameter(BaseParameter, float):
    def __new__(cls, value = 0, **kwargs):    
        return super(NumParameter, cls).__new__(cls, value)
    
    def __init__(self, value, **kwargs):
        
        BaseParameter.__init__(self, value, **kwargs)
        
        self.value = float(round(self._value, 4))
        
#         self.fix_value        = kwargs.get("fix_value", 0)
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
    #     return f"{self.parameter_number:>2}) {self.value:<11.4f} {self.fix_value:<10d} #  {self.comment}"
    
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
        fix_value        = ""
        
        BaseParameter.__init__(
            self, 
            value,
            name             = name,
            parameter_number = parameter_number,
            comment          = comment,
            fix_value        = fix_value
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
    def fix_value(self):
        return ""
    
    @fix_value.setter
    def fix_value(self, new_val):
        pass


# In[9]:


class MultiParameter(BaseParameter):
    def __init__(self, value = (0,0,0,0), **kwargs):
        
        BaseParameter.__init__(self, 0, **kwargs)
        
        # Generically use x and y
        if len(value) == 2:
            self._x1 = float(kwargs.get("x", value[0]))
            x_str = f"{self._x1:.4f}"
            
            self._y1 = float(kwargs.get("y", value[1]))
            y_str = f"{self._y1:.4f}"

            if self._x1 == int(self._x1):
                x_str = f"{self._x1:.0f}"
                
            if self._y1 == int(self._y1):
                y_str = f"{self._y1:.0f}"
                
            #self.x = self._x1
            #self.y = self._y1
            self._value = (self._x1, self._y1)
            #self._value_str  = f"{x_str} {y_str}"
            
        elif len(value) == 4:
            self._x1 = int(kwargs.get("xmin", value[0]))
            self._x2 = int(kwargs.get("xmax", value[1]))
            self._y1 = int(kwargs.get("ymin", value[2]))
            self._y2 = int(kwargs.get("ymax", value[3]))
            
            # Redundancy for ease of use
            # self.xmin   = self._x1
            # self.xmax   = self._x2
            # self.ymin   = self._y1
            # self.ymax   = self._y2
            
            self._value = (self._x1, self._x2, self._y1, self._y2)
            #self._value_str  = f"{self._x1}   {self._x2}   {self._y1}   {self._y2}"
        
        #self._fix_value        = ""
        self._fix_x            = ""
        self._fix_y            = ""
    
    def check_type(*iterable_in):
        assert isinstance(iterable_in, (list, tuple)), "The new value for this attribute must be a list/tuple!"
    def check_len(*iterable_in):
        assert len(iterable_in) in (2,4), "The length of this iterable is not 2 or 4, MultiParameter cannot correctly update."
# ==========================================================================================================

#     @property
#     def x(self):
#         return self._x1
    
#     @x.setter
#     def x(self, new_val):
#         self._x1 = new_val
        
#     @property
#     def y(self):
#         return self._y1
    
#     @y.setter
#     def y(self, new_val):
#         self._y1 = new_val

    exec(
        generate_get_set(
            {
                "x" : "_x1",
                "y" : "_y1"
            }
        )
    )

# ==========================================================================================================

#     @property
#     def x1(self):
#         return self._x1
    
#     @x1.setter
#     def x1(self, new_val):
#         self._x1 = new_val
        
#     @property
#     def x2(self):
#         return self._x2
    
#     @x2.setter
#     def x2(self, new_val):
#         self._x2 = new_val
        
#     @property
#     def y1(self):
#         return self._y1
    
#     @y1.setter
#     def y1(self, new_val):
#         self._y1 = new_val
        
#     @property
#     def y2(self):
#         return self._y2
    
#     @y2.setter
#     def y2(self, new_val):
#         self._y2 = new_val
        
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

# ==========================================================================================================

#     @property
#     def xmin(self):
#         return self._x1
    
#     @xmin.setter
#     def xmin(self, new_val):
#         self._x1 = new_val
        
#     @property
#     def xmax(self):
#         return self._x2
    
#     @xmax.setter
#     def xmax(self, new_val):
#         self._x2 = new_val
        
#     @property
#     def ymin(self):
#         return self._y1
    
#     @ymin.setter
#     def ymin(self, new_val):
#         self._y1 = new_val
        
#     @property
#     def ymax(self):
#         return self._y2
    
#     @ymax.setter
#     def ymax(self, new_val):
#         self._y2 = new_val
        
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
        
# ==========================================================================================================

    @property
    def value(self):
        self.check_len(self._value)
        
        if len(self._value) == 2:
            return (self.x1, self.y1)
        elif len(self._value) == 4:
            return (self.x1, self.y1, self.x2, self.y2)
    
    @value.setter
    def value(self, new_tuple):
        self.check_type(new_tuple)
        self.check_len(new_tuple)
          
        new_tuple = (float(i) for i in new_tuple)
        if len(self._value) == 2:
            self._value = (self.x1, self.y1)
        elif len(self._value) == 4:
            self._value = (self.x1, self.y1, self.x2, self.y2)
        
#     @property
#     def value_str(self):
#         return self._value_str
    
#     @value_str.setter
#     def value_str(self, new_str):
#         self._value_str = new_str

# ==========================================================================================================

#     @property
#     def fix_x(self):
#         return self._fix_x
    
#     @fix_x.setter
#     def fix_x(self, new_str):            
#         self._fix_x = new_str
        
#     @property
#     def fix_y(self):
#         return self._fix_y
    
#     @fix_y.setter
#     def fix_y(self, new_str):            
#         self._fix_y = new_str
        
    exec(
        generate_get_set(
            {
                "fix_x" : "_fix_x",
                "fix_y" : "_fix_y"
            }
        )
    )
    
    @property
    def fix_value(self):
        return self._fix_value
    
    @fix_value.setter
    def fix_value(self, new_val):
        if isinstance(new_val, str):
            self._fix_value = new_val
            
        elif isinstance(new_val, (list, tuple)):
            self.fix_x = new_val[0]
            self.fix_y = new_val[1]
            
            self._fix_value = f"{new_val[0]} {new_val[1]}"
        else:
            raise Exception("The new value for this attribute must be a string, list, or tuple!")
        
# ==========================================================================================================

    def __repr__(self):
        return repr(self.value)

    # Override BaseParameter
    def __str__(self):
        if len(self.value) == 2:
            x_str = f"{self.x1:.4f}"
            y_str = f"{self.y1:.4f}"
            
            if self.x1 == int(self.x1):
                x_str = f"{self.x1:.0f}"
                
            if self.y1 == int(self.y1):
                y_str = f"{self.y1:.0f}"
                
            value_str  = f"{x_str} {y_str}"
            
        elif len(self.value) == 4:
            value_str  = f"{self.x1}   {self.x2}   {self.y1}   {self.y2}"
        
        else:
            raise Exception(f"Something went wrong! Value in multiparameter is {self.value}")
            
        pre_comment = f"{self.parameter_number:>2}) {value_str:<11} {self.fix_value}"
        return f"{pre_comment:<23} # {self.comment}"


# In[10]:


# Use NamedTuples so we can access this like a dictionary
# while also using the indexing of a tuple where necessary
ntPosition = namedtuple("ntPosition", "x y")

class Position(MultiParameter):
    def __init__(self, value = (0,0), **kwargs):
        
        MultiParameter.__init__(self, value, **kwargs)
    
        self.value             = ntPosition(*self.value)
        
        self.name              = "position"
        self.parameter_number  = 1
        self.comment           = "Position x, y"
        
        self.fix_value        = kwargs.get("fix_value", "0 0")
        self.fix_x            = kwargs.get("fix_x", self.fix_value[0])
        self.fix_y            = kwargs.get("fix_y", self.fix_value[-1])
        
# ==========================================================================================================

    @property
    def value(self):
        self._value = ntPosition(self.x, self.y)
        return self._value
    
    @value.setter
    def value(self, new_tuple):
        self.check_type(new_tuple)
        self.check_len(new_tuple)
        
        new_tuple = (float(i) for i in new_tuple)
        self._value = ntPosition(*new_tuple)
        self.x      = self._value[0]
        self.y      = self._value[1]
        
#     @property
#     def fix_x(self):
#         return self._fix_x
    
#     @fix_x.setter
#     def fix_x(self, new_val):
#         self._fix_x = new_val
    
#     @property
#     def fix_y(self):
#         return self._fix_y
    
#     @fix_y.setter
#     def fix_y(self, new_val):
#         self._fix_y = new_val


# In[11]:


ntFourier = namedtuple("ntFourier", "amplitude phase_angle")

class FourierMode(MultiParameter):
    def __init__(self, mode, amplitude = 0, phase_angle = 0, **kwargs):
        
        self._mode        = mode
        
#         self._x1 = amplitude        
#         self._y1 = phase_angle
        
        if isinstance(amplitude, (list, tuple)):
            phase_angle = amplitude[1]
            amplitude = amplitude[0]
            
        # Override these properties later
        #self.amplitude   = self._amplitude
        #self.phase_angle = self._phase_angle
        
        #self.value_str  = f"{self.amplitude} {self.phase_angle}"
        
        MultiParameter.__init__(self, (amplitude, phase_angle), **kwargs)
        
        self.value = ntFourier(*self.value)
        
        self.name             = "position"
        self.parameter_number = f"F{self.mode}"
        self.comment         = f"Azim. Fourier mode {self.mode}, amplitude, & phase angle"
        
        self.fix_value       = kwargs.get("fix_value", "1 1")
        
        self.fix_x = kwargs.get("fix_amplitude", self.fix_value[0])
        self.fix_y = kwargs.get("fix_phase_angle", self.fix_value[-1])

        
# ==========================================================================================================
#     @property
#     def amplitude(self):
#         return self.x1

#     @amplitude.setter
#     def amplitude(self, new_val):
#         self.x1 = new_val

#     @property
#     def phase_angle(self):
#         return self.y1

#     @phase_angle.setter
#     def phase_angle(self, new_val):
#         self.y1 = new_val
        
    exec(
        generate_get_set(
            {
                "amplitude" : "x1",
                "phase_angle" : "y1"
            }
        )
    )
    
# ==========================================================================================================   

    exec(
        generate_get_set(
            {
                "mode" : "_mode",
                "fix_amplitude" : "fix_x",
                "fix_phase_angle" : "fix_y"
            }
        )
    )
    
    @property
    def value(self):
        self._value = ntFourier(self.amplitude, self.phase_angle)
        return self._value
    
    @value.setter
    def value(self, new_tuple):
        self.check_type(new_tuple)
        
        new_tuple = (float(i) for i in new_tuple)
        self._value      = ntFourier(*new_tuple)
        self.amplitude   = self._value[0]
        self.phase_angle = self._value[1]
        
#     @property
#     def mode(self):
#         return self._mode
    
#     @mode.setter
#     def mode(self, new_val):
#         self._mode = new_val
        
#     @property
#     def fix_amplitude(self):
#         return self.fix_x
    
#     @fix_amplitude.setter
#     def fix_amplitude(self, new_val):
#         self.fix_x = new_val
        
#     @property
#     def fix_phase_angle(self):
#         return self.fix_y
    
#     @fix_phase_angle.setter
#     def fix_phase_angle(self, new_val):
#         self.fix_y = new_val


# In[12]:


#class BendingModes(BaseParameter, float):
class BendingMode(NumParameter):
    def __init__(self, mode, amplitude, **kwargs):
        
        self._mode        = mode
        # self._amplitude   = amplitude
        
        NumParameter.__init__(self, amplitude, **kwargs)
    
        self.name             = "bending mode"
        self.parameter_number = f"B{self.mode}"
        self.comment         = f"Bending mode {self.mode} amplitude"
        
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


# In[13]:


ntImageRegionToFit = namedtuple("ntImageRegionToFit", "x1 x2 y1 y2")

class ImageRegionToFit(HeaderParameter, MultiParameter):
    def __init__(self, value = (0, 256, 0, 256), **kwargs):
        
        # Header Parameter first so that we don't overwrite
        # multiparameter fix value
        HeaderParameter.__init__(self, value, **kwargs)
        MultiParameter.__init__(self, value, **kwargs)
        
        self.value = ntImageRegionToFit(*self.value)
        
        self.fix_value = ""
        self.parameter_number = "H"
        self.comment = "Image region to fit (xmin xmax ymin ymax)"
        
    @property
    def value(self):
        self._value = ntImageRegionToFit(self.x1, self.y1, self.x2, self.y2)
        return self._value
    
    @value.setter
    def value(self, new_tuple):
        self.check_type(new_tuple)
            
        new_tuple = (int(i) for i in new_tuple)
        self._value = ntImageRegionToFit(*new_tuple)
        self.x1     = self._value[0]
        self.y1     = self._value[1]
        self.x2     = self._value[2]
        self.y2     = self._value[3]
        
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
        
        self.value     = ntConvolutionBox(*self.value)
        
        self.fix_value = ""
        self.parameter_number = "I"
        self.comment = "Size of the convolution box (x y)"
        
# ==========================================================================================================

    @property
    def value(self):
        self._value = ntConvolutionBox(self.x, self.y)
        return self._value
    
    @value.setter
    def value(self, new_tuple):
        self.check_type(new_tuple)
        self.check_len(new_tuple)
            
        new_tuple = (float(i) for i in new_tuple)
        self._value = ntConvolutionBox(*new_tuple)
        self.x      = self._value[0]
        self.y      = self._value[1]


# In[15]:


ntPlateScale = namedtuple("ntPlateScale", "dx dy")

class PlateScale(HeaderParameter, MultiParameter):
    def __init__(self, value = (0.396, 0.396), **kwargs):
        HeaderParameter.__init__(self, value, **kwargs)
        MultiParameter.__init__(self, value, **kwargs)
        
        #self._dx = self._x
        #self._dy = self._y
        self.value = ntPlateScale(*self.value)
        
        self.fix_value = ""
        self.parameter_number = "K"
        self.comment = "Plate scale (dx dy)   [arcsec per pixel]"
        
    exec(
        generate_get_set(
            {
                "dx" : "x1",
                "dy" : "y1"
            }
        )
    )
    
    @property
    def value(self):
        self._value = ntPlateScale(self.dx, self.dy)
        return self._value
    
    @value.setter
    def value(self, new_tuple):
        self.check_type(new_tuple)
            
        self._value = ntPlateScale(*new_tuple)
        self.dx     = float(self._value[0])
        self.dy     = float(self._value[1])

#     @property
#     def dx(self):
#         return self.x1
    
#     @dx.setter
#     def dx(self, new_val):
#         self.x1 = new_val
    
#     @property
#     def dy(self):
#         return self.y1
    
#     @dy.setter
#     def dy(self, new_val):
#         self.y1 = new_val


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
    arms_line  = ComponentType("power", parameter_number = "R0")
    sky_line   = ComponentType("sky", component_number = 3)

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


# In[24]:


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


# In[25]:


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


# In[26]:


def load_default_power_parameters(component_number = None):
    
    param_prefix = "R"
    
    _power = ComponentType("power", parameter_number = f"{param_prefix}0", component_number = component_number)
    
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
    
    # skip = Skip(
    #     0,
    #     component_name = "Power",
    #     component_number = component_number
    # )
    
    loc = deepcopy(locals())
    loc.pop("component_number")
    loc.pop("param_prefix")
    
    return loc


# In[27]:


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
    
    skip = Skip(
        0,
        component_name = "Fourier",
        component_number = component_number
    )
    
    loc = deepcopy(locals())
    loc.pop("component_number")
    loc.pop("param_prefix")
    
    return loc


# In[28]:


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
        parameter_number = 2,
        comment = "dsky/dy (sky gradient in y)     [ADUs/pix]",
        component_name = "Sky",
        component_number = component_number
    )
    
    loc = deepcopy(locals())
    loc.pop("component_number")
    
    return loc


# In[29]:


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
        0,
        name = "optimize",
        parameter_number = "P",
        comment = "Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps",
        component_name = "Header"
    )
    
    loc = deepcopy(locals())
    
    return loc


# In[30]:


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


# In[31]:


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


# In[32]:


if __name__ == "__main__":
    export_to_py("Parameters", pj(_MODULE_DIR, "Classes", "Parameters"))

