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
import gc
import psutil

import numpy as np
import scipy.linalg as slg
from scipy.stats import norm, kstest
from skimage.draw import disk, ellipse
import matplotlib.pyplot as plt


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
#from Classes.Components import *
from Classes.Containers import *
from Functions.helper_functions import *


# In[4]:


class HDU:
    def __init__(self, 
                 name   = "observation",
                 header = {}, 
                 data   = None
                ):
        
        self._hdu_info = {
            "name"   : name,
            "header" : deepcopy(dict(header)),
            "data"   : deepcopy(np.array(data))
        }
        
# ==========================================================================================================

    @property
    def name(self):
        return self._hdu_info.get("name", "")
    
    @name.setter
    def name(self, new_name):
        self._hdu_info["name"] = new_name
        
    @property
    def header(self):
        return self._hdu_info.get("header", {})
    
    @header.setter
    def header(self, new_header):
        self._hdu_info["header"] = deepcopy(dict(new_header))
        
    @property
    def data(self):
        return self._hdu_info.get("data", "")
    
    @data.setter
    def data(self, new_data):
        self._hdu_info["data"] = deepcopy(np.array(new_data))
        
# ==========================================================================================================

    def __str__(self):
        header_str = ""
        for k,v in self.header.items():
            header_str += f"{k} = {v}\n"
            
        output_str = f"{self.name}\n{header_str}Img size: {np.shape(self.data)}"
        return output_str


# In[5]:


class FitsFile:
    def __init__(self,
                 filepath,
                 names       = ["observation"],
                 from_galfit = False,
                 wait        = False,
                 **kwargs
                ):
        
        self.filepath = filepath
        self.all_hdu  = {}
        
        # Use split over rstrip in case a waveband designation is given
        # (rstrip will remove any character that matches in the substring)
        # i.e. 12345678910_g would lose the _g for "_galfit_out.fits"
        # TODO: Replace rstrip with split in the rest of these scripts...
        self.gname    = kwargs.get("gname", os.path.basename(filepath).split("_galfit_out.fits")[0])
        
        assert os.path.splitext(filepath)[-1].lower() == ".fits", "File being passed into FitsHandler must be .fits!"
        
        try:
            file_in = fits.open(filepath)
            
        except FileNotFoundError:
            print(f"Can't open to read the file, {filepath}. Check name/permissions/directory.")
            raise(Exception()) 

        except OSError as ose:
            print(f"Something went wrong! {ose}")
            raise(Exception())
            
        # FITS starts the index at 0 but GALFIT outputs the observation image at 1
        # Also converting the header to a dict to save some trouble
        assert_str = f"Number of HDU names fed to object ({len(names)}) does not match number of HDUs in {filepath} ({len(file_in)})!"
        if from_galfit:
            assert len(names) + 1 == len(file_in), assert_str
            self.num_imgs = len(file_in) - 1
        else:
            assert len(names)     == len(file_in), assert_str
            self.num_imgs = len(file_in)
        
        for i, name in enumerate(names):
            
            index = i
            
            if from_galfit:
                index += 1

            header   = deepcopy(dict(file_in[index].header))
            data     = deepcopy(file_in[index].data)

            hdu = HDU(name = name, header = header, data = data)

            self.all_hdu[name] = hdu
            
        if self.num_imgs == 1:
            self.header = self.all_hdu[names[0]].header
            self.data   = self.all_hdu[names[0]].data
            
        self.file = file_in
        
        # Wait is for continuing to use the file in some other capacity
        # i.e. for outputfits below to grab more info
        if not wait:
            #file_in.close(verbose = True)
            self.close()
        
        #print("Did it close?", file_in.closed)
        # assert hdu_num == 4, "File being passed into FitsHandler has too few output HDUs."
# ==========================================================================================================

    def check_hdu_type(self, hdu):
        assert isinstance(hdu, HDU), "Input HDU is not an HDU class!"
        
# ==========================================================================================================
    @property
    def observation(self):
        return self.all_hdu.get("observation", None)
    
    @observation.setter
    def observation(self, new_hdu):
        self.check_hdu_type(new_hdu)
        self.all_hdu["observation"] = new_hdu

# ==========================================================================================================

    #def close(self):
    #    self.file.close(verbose = True)
    #    return

    def close(self):
        """Destructor for closing FITS files."""
        
        for ext in self.file:
            try:
                del ext.data
                del ext.header
            except AttributeError as ae:
                pass
            gc.collect()

        try:
            self.file.close(verbose = True)
        except Exception as ee:
            print(f'Failed to close FITS instance for {self.gname}: {ee}')
            print("May already be closed...")

# ==========================================================================================================

    def to_png(self, cleanup = True, **kwargs): #tmp_fits_path = "", tmp_png_path = "", out_filename = ""):
        
        gname          = kwargs.get("gname", self.gname)
        
        # TODO: BAD ASSUMPTION MOVING FORWARD
        #tmp_fits_path = kwargs.get("tmp_fits_path", self.filepath)
        fits_path      = kwargs.get("fits_path", self.filepath)
        # .../galfits -> galfit_png
        #tmp_png_dir   = os.path.split(tmp_fits_path)[0].rstrip("s") + "_png"
        tmp_png_dir    = kwargs.get("tmp_png_dir", "./")
        tmp_png_path   = pj(tmp_png_dir, gname)
        tmp_png_path   = kwargs.get("tmp_png_path", tmp_png_path)
        
        out_png_dir    = kwargs.get("out_png_dir", "./")
        
        capture_output = bool(kwargs.get("silent", False))
        
        fitspng_param  = "0.25,1" #1,150"
        
        # run_fitspng from helper_functions, string path to fitspng program
        fitspng_cmd1   = f"{run_fitspng} -fr \"{fitspng_param}\" -o {tmp_png_path}.png {fits_path}[1]"
        
        fitspng_cmd2   = f"{run_fitspng} -fr \"{fitspng_param}\" -o {tmp_png_path}_out.png {fits_path}[2]"
        
        fitspng_cmd3   = f"{run_fitspng} -fr \"{fitspng_param}\" -o {tmp_png_path}_residual.png {fits_path}[3]"
        
        cmds = [fitspng_cmd1, fitspng_cmd2, fitspng_cmd3]
        
        # sp is from helper_functions, subprocess.run call
        for cmd in cmds[:self.num_imgs]:
            # We must capture this call to check if the conversion worked
            fitspng_out = sp(cmd, capture_output = True)
            
            if "error" in fitspng_out.stderr:
                print("Skipping fitspng conversion... there is likely a library (libcfitsio) issue.")
                self.combined_png = ""
                return
        
        im1 = f"{tmp_png_path}.png"
        im2 = f"{tmp_png_path}_out.png"
        im3 = f"{tmp_png_path}_residual.png"
        
        combined = ""
        if self.num_imgs > 1:
            combined = "_combined"
        
        # Adding 'magick' to use the portable version in the GalfitModule
        montage_cmd = "magick montage " + \
                      " ".join(im_cmd for idx, im_cmd in enumerate([im1, im2, im3]) 
                               if idx <= self.num_imgs)
        
        tiling = f"1x{self.num_imgs}"
        if kwargs.get("horizontal", None):
            tiling = tiling = f"{self.num_imgs}x1"            
            
        # Combining the images using ImageMagick
        # If this is a single image, it'll also resize for me so that's why I leave it in
        montage_cmd += f" -tile {tiling} -geometry \"175x175+2+2\" \
                        {pj(out_png_dir, gname)}{combined}.png"
            
        _ = sp(montage_cmd, capture_output = capture_output)            
        
        if cleanup:
            _ = rm_files(im1, im2, im3)
        else:
            self.observation_png = im1
            self.model_png       = im2
            self.residual_png    = im3            
            
        self.combined_png    = f"{pj(out_png_dir, gname)}{combined}.png"

# ==========================================================================================================

    def __sub__(self, other):
        
        names = self.all_hdu.keys()
        
        assert_str1 = "Cannot subtract the data from these two FITS files, they do not contain the same number of HDUs!"
        assert len(self.all_hdu) == len(other.all_hdu), assert_str1
        
        assert_str2 = "Cannot subtract the data from these two FITS files, they do not have the same image dimensions!\n"
        for i, (a, b) in enumerate(zip(self.all_hdu.values(), other.all_hdu.values())):
            shape_a = np.shape(a.data)
            shape_b = np.shape(b.data)
            assert shape_a == shape_b, assert_str2 + f"At HDU {i}, the images have shapes {shape_a} & {shape_b}."
            
        # Python doesn't care if they're different lengths but
        # (for instance in the residual) we don't want to compare one to one
        result = {k : a.data - b.data for k, a, b in zip(names, self.all_hdu.values(), other.all_hdu.values())}
        
        return result

# ==========================================================================================================

    # # Use str to display feedme(?)
    # def __str__(self):
    #     pass
        
        #return output_str
        
# ==========================================================================================================

    # Use str to display feedme(?)
    def header_dict(self, name = ""):
        
        if name:
            output_dict = dict(self.all_hdu[name].header)
        else:
            output_dict = {name : dict(hdu.header) for name, hdu in self.all_hdu.items()}
            
        return output_dict
    
# ==========================================================================================================


# In[6]:


class OutputFits(FitsFile):

    def __init__(
        self, 
        filepath, 
        names = ["observation", "model", "residual"], 
        load_default = True, 
        **kwargs
    ):
            
        FitsFile.__init__(
            self, 
            filepath    = filepath,
            names       = names,
            wait        = True,
            from_galfit = True,
            **kwargs
        )
        
        # Dict is very redundant here but just for funsies
        # FITS header not Feedme header
        self.header = deepcopy(dict(self.model.header))
        
        # Can call the helper directly since we're just using the header dict
        #_header.from_file_helper_dict(self.header)
        self.feedme = FeedmeContainer(path_to_feedme = filepath, header = GalfitHeader(), load_default = load_default, **kwargs)
        self.feedme.from_file(self.header)
        
        self.data = deepcopy(self.model.data)
        
        self.bulge_mask = np.ones(np.shape(self.data))
        
        self.close()
        
# ==========================================================================================================

    @property
    def model(self):
        return self.all_hdu.get("model", None)
    
    @model.setter
    def model(self, new_hdu):
        self.check_hdu_type(new_hdu)
        self.all_hdu["model"] = new_hdu
        
    @property
    def residual(self):
        return self.all_hdu.get("residual", None)
    
    @residual.setter
    def residual(self, new_hdu):
        self.check_hdu_type(new_hdu)
        self.all_hdu["residual"] = new_hdu
        
# ==========================================================================================================
        
    def generate_bulge_mask(self, sparcfire_csv):
        
        # Thanks Azra!
        bulge_mask = np.ones(np.shape(self.model.data))
        
        try:
            info = pd.read_csv(sparcfire_csv, dtype = str).dropna()
        except FileNotFoundError as fe:
            print(fe)
            return bulge_mask
        
        if "rejected" in info[' fit_state']:
            print(f"SpArcFiRe fit_state 'rejected'. Cannot determine the bulge mask for {self.gname}.")
            return bulge_mask
            
        try:
            input_size = float(str(info[' iptSz'][0]).split()[0][1:])
        except Exception as e:
            print(f"There is an issue determining the bulge mask for {self.gname}.")
            return bulge_mask
        
        bulge_rad  = float(info[' bulgeMajAxsLen'][0])
        # In radians
        bulge_angle = float(info[' bulgeMajAxsAngle'][0])
        axis_ratio  = float(info[' bulgeAxisRatio'][0]) # Maj/minor
        if axis_ratio < 0.5:
            axis_ratio = 0.5
        
        # + 1 added for effect
        major_rad = int(bulge_rad * len(self.model.data[0]) // input_size) + 1
        minor_rad = int(major_rad/axis_ratio)

        crop_box = self.feedme.header.region_to_fit.value
        # To adjust for python indexing
        xbox_min, ybox_min = crop_box[0] - 1, crop_box[2] - 1
        
        # Shifting everything to origin
        center_x, center_y = np.array(self.feedme.bulge.position.value, dtype = int) - 1 -\
                             np.array((xbox_min, ybox_min), dtype = int)
        
        # center_x2, bulge_y2 = np.array(self.feedme.bulge.position, dtype = int) - \
        #                      np.array((xbox_min, ybox_min), dtype = int) + \
        #                      rad

        #xx, yy = disk((center_x, center_y), major_rad)
        try:
            xx, yy = ellipse(center_x, center_y, major_rad, minor_rad, rotation = bulge_angle, shape = np.shape(self.model.data))
        except Exception as e:
            print(e)
            print(self.gname)
            print(center_x, center_y, major_rad, minor_rad, rotation)
            return bulge_mask
        
        bulge_mask[xx, yy] = 0
        self.bulge_mask = bulge_mask
        
#         temp = self.model.data
#         plt.imshow(temp, origin = "lower")
#         plt.show()
        
#         temp[xx, yy] = np.min(self.model.data[np.nonzero(self.model.data)])
#         plt.imshow(temp, origin = "lower")
#         plt.show()

        #self.close()
        return bulge_mask
    
# ==========================================================================================================
        
    def generate_masked_residual(self, mask, use_bulge_mask = True, update_fits_header = True):

        small_number = 1e-8
        
        crop_box = self.feedme.header.region_to_fit.value
        # To adjust for python indexing
        # Also, reminder, non-inclusive of end
        xbox_min, xbox_max, ybox_min, ybox_max = crop_box[0] - 1, crop_box[1], crop_box[2] - 1, crop_box[3]
        # To invert the matrix since galfit keeps 0 valued areas
        crop_mask = 1
        if mask is not None and np.any(mask.data):
            cropped_mask = mask.data[xbox_min:xbox_max, ybox_min:ybox_max]
            mask_shape = np.shape(mask.data)
            # Fixing two common issues, hence two if statements
            if np.shape(cropped_mask) != np.shape(self.model.data):
                # The issue seems to pop up at the max border when either the min or the max runs up
                # against the original image size
                print("Shape mismatch between crop mask and model. Likely an indexing issue due to crop mask running into image bounds. Attempting to fix.")
                diff = np.array(np.shape(cropped_mask)) - np.array(np.shape(self.model.data))
                # One or the other... No graceful way to do this
                if xbox_min == 0:
                    xbox_max -= diff[0]
                else:
                    xbox_min -= diff[0]
                    
                if ybox_min == 0:
                    ybox_max -= diff[0]
                else:
                    ybox_min -= diff[0]
                
                cropped_mask = mask.data[xbox_min:xbox_max, ybox_min:ybox_max]
                
            # Now we're just trying to brute force it and hope for the best
            if np.shape(cropped_mask) != np.shape(self.model.data):
                print("Shape mismatch (again) between crop mask and model. Likely an indexing issue due to crop mask running into image bounds. Attempting to fix.")
                diff = np.abs(np.array(np.shape(cropped_mask)) - np.array(np.shape(self.model.data)))
                cropped_mask = np.pad(cropped_mask, ((diff[0],0), (diff[1],0)), 'constant')
                
            # Giving up and proceeding without
            if np.shape(cropped_mask) != np.shape(self.model.data):
                print("Shape mismatch. Proceeding without crop mask.")
                cropped_mask = 0
                
            crop_mask = 1 - cropped_mask
            
        if use_bulge_mask:
            feedme_dir, feedme_file = os.path.split(self.feedme.path_to_feedme)
            
            if exists(pj(feedme_dir, f"{self.gname}.csv")):
                crop_mask = self.generate_bulge_mask(pj(feedme_dir, f"{self.gname}.csv")) * crop_mask
            else:
                # REQUIRES GENERATE_BULGE_MASK TO BE RUN SEPARATE WITH CSV FILE SPECIFIED 
                try:
                    crop_mask = self.bulge_mask * crop_mask
                except AttributeError:
                    print(f"Could not generate bulge mask for {self.gname}. Check location of csv or run generate_bulge_mask with a specified csv file.")
                except ValueError:
                    print(f"Could not generate bulge mask for {self.gname}. There may be an issue with sparcfire output (broadcast issue).")
        
        try:
            # compare to gaussian with same mean, std via kstest
            # if p value high, not that different
            
            self.masked_residual  = (self.observation.data - self.model.data)*crop_mask
            exclude_masked_pixels = self.masked_residual[np.abs(self.masked_residual) > 0]
            mean                  = np.mean(exclude_masked_pixels)
            std                   = np.std(exclude_masked_pixels)
            gaussian              = norm.rvs(size = len(exclude_masked_pixels), loc = mean, scale = std, random_state = 0)
            self.kstest           = kstest(gaussian, exclude_masked_pixels.flatten())
            pvalue                = self.kstest.pvalue
            #statistic             = self.kstest.statistic
            # gaussian = norm.rvs(size = len(self.masked_residual)**2, loc = mean, scale = std, random_state = 0)
            # noised_masked_pixels = np.where(np.abs(self.masked_residual.flatten()) > 0, self.masked_residual.flatten(), gaussian)
            # self.kstest = kstest(gaussian, noised_masked_pixels)

            self.norm_observation           = slg.norm(crop_mask*self.observation.data)
            self.norm_model                 = slg.norm(crop_mask*self.model.data)
            self.norm_residual              = slg.norm(crop_mask*self.residual.data)
            self.masked_residual_normalized = self.masked_residual/min(self.norm_observation, self.norm_model)
            self.wayne_residual             = self.norm_residual/self.norm_observation
            self.wayne_quality              = pvalue/self.wayne_residual # bigger is better
            
#             obs_model = 1 - np.divide(
#                                 crop_mask*self.observation.data/self.norm_observation, 
#                                 crop_mask*self.model.data/self.norm_model + small_number
#                                      )

#             model_obs = 1 - np.divide( 
#                                 crop_mask*self.model.data/self.norm_model,
#                                 crop_mask*self.observation.data/self.norm_observation + small_number
#                                     )
            # Replace negative values with 1 - reciprocal
#            self.masked_residual_ratio = np.where(obs_model >= 0, obs_model, model_obs)
            # Masked residual normalized
            # I seem to use this acronym a lot
            self.nmr  = slg.norm(self.masked_residual_normalized)
#            self.nmrr = slg.norm(self.masked_residual_ratio)

            if update_fits_header:
                with fits.open(self.filepath, mode='update', output_verify='ignore') as hdul:
                    hdul[2].header["NMR"]   = (round(self.nmr, 8), "Norm of the masked residual")
                    hdul[2].header["W_NMR"] = (round(self.wayne_residual, 8), "Wayne's residual")
                    #hdul[2].header["W_Q"]   = (round(self.nmr, 8), "Wayne's quality measure")

                    # pvalue is sometimes none but round can't handle it
                    if isinstance(pvalue, float): # and isinstance(statistic, float):
                        hdul[2].header["KS_P"]    = (round(pvalue, 8), "p value of kstest vs noise")
                        #hdul[2].header["KS_STAT"] = (round(statistic, 8), "statistic value of kstest vs noise")
                    else:
                        hdul[2].header["KS_P"]    = (None, "p value of kstest vs noise")
                        #hdul[2].header["KS_STAT"] = (None, "statistic value of kstest vs noise")

        except ValueError:
            print(f"There may be a broadcast issue, observation, model, crop mask: ", end = "")
            print(f"{np.shape(self.observation.data)}, {np.shape(self.model.data)}, {np.shape(crop_mask)}")
            # print(np.shape(mask_fits_file.data))
            # print(np.shape(fits_file.data))
            # print(crop_box)
            #self.close()
            return None
        
        #self.close()
        
        return self.masked_residual_normalized
    
# ==========================================================================================================


# In[7]:


if __name__ == "__main__":
    from RegTest.RegTest import *


# In[8]:


# Testing from_file
if __name__ == "__main__":
    
    gname = "1237671124296532233"
    obs   = pj(SAMPLE_DIR, f"{gname}.fits")
    model = pj(SAMPLE_DIR, f"{gname}_galfit_out.fits")
    mask  = pj(SAMPLE_DIR, f"{gname}_star-rm.fits")
    
    test_obs   = FitsFile(obs)
    test_model = OutputFits(model)
    test_mask  = FitsFile(mask)
    
    print(test_obs.observation)
    print()
    print(test_model.feedme)
    print()
    print(test_model.model)
    
    # Purposefully do not fill in some of the header parameters
    # since those do not exist in the output FITS header
    # This is done to remind the user/programmer that the 
    # OutputFits object only serves to represent the header
    # nothing more, nothing less and so also reminds them to
    # use a different method to fill in the header.
    #print(test_model.feedme.header)
    
#     _header = GalfitHeader()
#     _header.from_file_helper(test_out.header)
    
#     crop_box = _header.region_to_fit
#     # To adjust for python indexing
#     box_min, box_max = crop_box[0] - 1, crop_box[1]
        
#     print(np.shape(test_in.data[box_min:box_max, box_min:box_max]))
    print("\nThese should all be the same .")
    print(np.shape(test_model.observation.data))
    print(np.shape(test_model.data))
    print(np.shape(test_model.residual.data))
    crop_box = test_model.feedme.header.region_to_fit.value
    # + 1 to account for python indexing
    crop_rad = crop_box[1] - crop_box[0] + 1
    print(f"({crop_rad}, {crop_rad})")
    print("Andddd pre crop")
    print(np.shape(test_obs.observation.data))


# In[9]:


# Unit test to check value of masked residual
if __name__ == "__main__":
    
    # Turn off updating FITS header to avoid modifying test data
    print("No bulge mask")
    _ = test_model.generate_masked_residual(test_mask, use_bulge_mask = False, update_fits_header = False)
    print(f"Norm of the observation: {test_model.norm_observation:.4f}")
    print(f"Norm of the model: {test_model.norm_model:.4f}")
    print(f"Norm of the residual: {test_model.norm_residual:.4f}")
    print(f"Norm of the masked residual: {test_model.nmr:.4f}")
    print(f"Wayne's residual: {test_model.wayne_residual:.4f}")
    #print(f"Norm of the masked residual ratio: {test_model.nmrr:.8f}")
    print(f"kstest p value: {test_model.kstest.pvalue:.4f}")
    #print(f"kstest statistic: {test_model.kstest.statistic:.4f}")
    
    print("\nNow with bulge mask")
    _ = test_model.generate_bulge_mask(pj(TEST_DATA_DIR, "test-out", gname, f"{gname}.csv"))
    _ = test_model.generate_masked_residual(test_mask, update_fits_header = False)
    print(f"Norm of the observation: {test_model.norm_observation:.4f}")
    print(f"Norm of the model: {test_model.norm_model:.4f}")
    print(f"Norm of the residual: {test_model.norm_residual:.4f}")
    print(f"Norm of the masked residual: {test_model.nmr:.4f}")
    print(f"Wayne's residual: {test_model.wayne_residual:.4f}")
    #print(f"Norm of the masked residual ratio: {test_model.nmrr:.8f}")
    print(f"kstest p value: {test_model.kstest.pvalue:.4f}")
    #print(f"kstest statistic: {test_model.kstest.statistic:.4f}")
    #print(np.min(test_model.observation.data))


# In[10]:


if __name__ == "__main__":
    model           = pj(SAMPLE_DIR, f"sample_model_galfit_out.fits")
    model_to_update = pj(TEST_OUTPUT_DIR, f"temp_galfit_out.fits")
    
    if exists(model_to_update):
        sp(f"rm -f {model_to_update}")
        
    _ = sp(f"cp {model} {model_to_update}")

    test_model = OutputFits(model_to_update)
    print("Checking FITS header update with NMR")
    
    print("Does the updated FITS file contain NMR and KStest keys?")
    keys_to_check = ("NMR", "KS_P", "W_NMR")
    
    # TODO: replace fits file with one without those header options
    print("Before... (expect False)", all(k in test_model.header for k in keys_to_check))
    
    _ = test_model.generate_masked_residual(test_mask)
    test_model = OutputFits(model_to_update)

    print("After...", all(k in test_model.header for k in keys_to_check))


# In[11]:


if __name__ == "__main__":
    print("Checking if all FITS files are closed...")
    print("Expect True:", not any("fits" in pof.path for pof in psutil.Process().open_files()))


# In[12]:


if __name__ == "__main__":
    export_to_py("FitsHandlers", pj(_MODULE_DIR, "Classes", "FitsHandlers"))

