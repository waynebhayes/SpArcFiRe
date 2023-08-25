#!/usr/bin/env python
# coding: utf-8

# **Author: Matthew Portman**
# 
# **Date (Commit date will likely be more accurate): 4/17/23**

# # MINI README
# 
# Please place this script into the directory which contains the folders for your input, temporary, and output files for SpArcFiRe typically denoted: sparcfire-in, sparcfire-tmp, and sparcfire-out. GALFIT will also be run from that same directory if you choose to use the control script included in the repo. I recommend running this alone *once* before running the big script to ensure everything is in its right place (https://www.youtube.com/watch?v=sKZN115n6MI). After you confirm that this works without any errors the first time around, feel free to run the control script from then on. 
# 
# Running from the overarching directory is a temporary measure which will be remedied upon completion of the entire control script and full integration with SpArcFiRe. 
# 
# TO RUN: `python3 sparc_to_galfit_feedme_gen.py`

# In[1]:


import numpy as np
from math import log
import glob
import csv
import subprocess
import random
import pandas as pd
from copy import deepcopy

from astropy.io import fits


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


import sys
import os
from os.path import join as pj
from os.path import exists

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
from Classes.Containers import FeedmeContainer
from Functions.helper_functions import *


# # To convert autocrop image to fits
# import subprocess
# 
# try: 
#     autocrop_img_filename = glob.glob('*_autoCrop.png')[0]
#     #convert_command = ["convert ", autocrop_img_filename,"input.fits"]
#     convert_command = "convert " + autocrop_img_filename + " input.fits"
# 
# except:
#     print("Cannot find ", autocrop_img_filename, " to convert to fits.")
#     print("Check Sparcfire output or directories. Cannot proceed. Quitting.")
#     raise SystemExit("Quitting.")
#     
# else:
#     try: 
#         subprocess.run(convert_command, shell=True)
#     
#     except:
#         print("Cannot convert to FITS:", autocrop_img_filename)
#         print("Check Sparcfire output or directories. Cannot proceed.")
#         raise SystemExit("Quitting.")

# In[ ]:


# # Grabbing filepath from command line
# def command_line(top_dir = os.getcwd(), **kwargs): # = True):
    
#     in_dir_out  = kwargs.get("in_dir", pj(top_dir, "sparcfire-in"))
#     tmp_dir_out = kwargs.get("tmp_dir", pj(top_dir, "sparcfire-tmp"))
#     out_dir_out = kwargs.get("out_dir", pj(top_dir, "sparcfire-out"))
        
#     #if not run_as_script:
#     #    return in_dir_out, tmp_dir_out, out_dir_out
    
#     try:
#         if len(sys.argv) != 4: # including name of python script
#             print(f"Using: \n{in_dir_out}\n{tmp_dir_out}\n{out_dir_out}\nto generate feedme.")

#         else:
#             in_dir_out = argv[1]
#             tmp_dir_out = argv[2]
#             out_dir_out = argv[3]
#             print(f"Using: \n{in_dir_out}\n{tmp_dir_out}\n{out_dir_out}\nto generate feedme.")
#     except:
#         pass

#     return absp(in_dir_out), absp(tmp_dir_out), absp(out_dir_out)


# In[ ]:


def get_galaxy_names_list(in_dir, tmp_dir, out_dir, galaxy_names = []):

    gnames_out = galaxy_names
#     try:
#         #filenames_read = find_files(in_dir, "*.fits", "f")
    
#     #except:
#         # Feel free to hardcode this 
#         #print("Please input the full path to the input directory for the SpArcFiRe run and include the trailing /")
#         #print("e.g. /home/usr/sparcfire-in/")
#         #in_dir = input()
    
#         #try:
#             #filenames_in = glob.glob(in_dir + "*.fits")

#         #except:
#             #print("See instructions above or check your directory path. Could not glob.")
#             #raise SystemExit("Exitting.")
    
#     except:
#         print("Please copy me into the directory which contains the folders for your")
#         print("input, temporary, and output files for SpArcFiRe denoted:")
#         print("sparcfire-in, sparcfire-tmp, and sparcfire-out.")
#         raise SystemExit("Exitting.")
        
#     else:
    if not gnames_out:
        print("No galaxy names fed in. Finding them to generate feedmes.")
    #     print("Something went wrong! No galaxies fed into feedme generator.")
    #     print("Please copy me into the directory which contains the folders for your")
    #     print("input, temporary, and output files for SpArcFiRe denoted:")
    #     print("sparcfire-in, sparcfire-tmp, and sparcfire-out.")
    #     raise SystemExit("Exitting.")
        gnames_out  = [os.path.basename(s).replace(".fits", "")
                       for s in find_files(in_dir, "*.fits", "f")
                       if exists(pj(out_dir, os.path.basename(s).replace(".fits", "")))
                      ]

    folders_out = [pj(out_dir, gname) for gname in gnames_out]
        
    return gnames_out, folders_out


# In[ ]:


def path_join(path='.', name='', file_ext=''):
    
    file_path = pj(path, name + file_ext)

    # deprecated
    # "./" + path + '/' + name + file_ext
    # print(file_path)
    # print(glob.glob(file_path))
    # file_name = glob.glob(file_path)[0]
    
    return file_path


# In[ ]:


def scale_var(x, scale = 1):
    # Scales
    if not x:
        return None
    
    return float(x)*scale


# In[ ]:


def galaxy_information(galaxy_name, galaxy_path):
     
    kwargs_out = {
        "bulge_maj_axs_len" : 2,
        "bulge_axis_ratio" : 0.5,
        "bulge_rot_angle" : 30,
        "crop_rad" : 30, # New!
        "center_pos_x" : 30,
        "center_pos_y" : 30,
        "disk_maj_axs_len" : 30,
        "disk_rot_angle" : 30,
        #"pos_angle_sersic" : 1,
        "pos_angle_power" : 30,
        "disk_axis_ratio" : 0.5,
        "avg_arc_length" : 30,
        "max_arc_length" : 20,
        "alpha" : 1,
        "est_arcs" : 2,
        "inclination" : 30,
        "bar_candidate" : 'FALSE',
        "bar_len"  : 5, 
        # spin_parity handled by if 
        "spin_parity" : '' #random.choice(['','-'])
                    }
           
    # Making this global since I'm now grabbing the necessary info from the csv
    # std for standardized image
    global scale_fact_std
    scale_fact_std = 2*kwargs_out["crop_rad"]/256
        
    chirality = 0
    chirality_2 = 0
    chirality_3 = 0
    pitch_angle = 20
    spin_parity = ""
    
    failure_modes = ["input rejected", 
                     "Subscript indices must either be real positive integers or logicals.",
                     "SIGMA must be a square  symmetric  positive definite matrix."
                    ]
    
    try:
        # Using this for the *big* runs because *somebody* (Wayne)
        # didn't validate and overwrite the original csv's to include CR (crop radius)
        csv_CR = path_join(galaxy_path, galaxy_name, '.csv+CR')
        csv_filename = path_join(galaxy_path, galaxy_name, '.csv')
        if exists(csv_CR):
            csv_filename = csv_CR
            
        csv_file = open(csv_filename, 'r')
        
    except:
        print("Can't open to read: ", csv_filename)
        print("Check Sparcfire output or directories. Proceeding with default values.")
        return kwargs_out

    
    reader = csv.DictReader(csv_file, skipinitialspace = True)
    for row in reader:
        # if sparcfire fails
        if row.get('fit_state', "") in failure_modes:
            return kwargs_out
        
        elif row.get('fit_state', "") != "OK":
            print(f"fit state is: {row.get('fit_state', '')}")
            return kwargs_out

        try:
            # Already accounted for, no need to scale
            center_pos_x = float(row.get('inputCenterR'))
            center_pos_y = float(row.get('inputCenterC'))

            crop_rad = float(row.get('cropRad'))

        except (ValueError, TypeError) as ve:
            print(f"SpArcFiRe likely failed on this galaxy, {ve}. Proceeding with default values...")
            return kwargs_out

        # Anything with converting to float/int should go in here since there are a couple issues besides
        # fit state which still result in an empty csv
        try:
            input_size    = int(row.get('iptSz', "[256 256]").strip().split()[0][1:])
            global scale_fact_ipt
            if input_size:
                scale_fact_ipt = 2*crop_rad/input_size
            else:
                scale_fact_ipt = 2*crop_rad/256

            bulge_maj_axs_len = scale_var(row.get('bulgeMajAxsLen'), scale_fact_ipt)
            # TODO: Take bulge axis ratio from SDSS? deVAB_r
            bulge_axis_ratio = float(row.get('bulgeAxisRatio'))

            # According to the dissertation this is given in terms of the input image not the standardized image
            disk_maj_axs_len = scale_var(row.get('diskMajAxsLen'), scale_fact_ipt) 

            # Angles don't need to scale
            # But sparcfire is flipped across the y-axis and uses mathematical convention
            # As opposed to astronomical (0 is vertical y axis, 90 is counter clockwise to x)
            bulge_rot_angle = 90 - np.degrees(float(row.get('bulgeMajAxsAngle')))
            
            # Power angle is mathematical convention, disk angle is astronomical
            pos_angle_power = np.degrees(float(row.get('diskMajAxsAngleRadians')))
            disk_rot_angle  = 90 - pos_angle_power
            
            disk_axis_ratio = float(row.get('diskAxisRatio'))

            # GRABBING ARC STATISTICS
            avg_arc_length = scale_var(row.get('avgArcLength'), scale_fact_std)
            #med_arc_length = row['medianArcLength'][:6]
            max_arc_length = scale_var(row.get('maxArcLength'), scale_fact_std)

            # Grabbing PA absolute average from dominant chirality
            pitch_angle = abs(float(row.get('pa_alenWtd_avg_domChiralityOnly')))

            # For estimating number of real arcs to influence fourier modes
            # This seems to be a much better way to do it
            est_arcs = int(row.get('rankAt75pct'))

        # Happens when row is empty... This is easier than messing with other conditionals for each scale factor/float conversion
        except ValueError as ve:
            return kwargs_out

        if disk_rot_angle < 0: 
            inclination = np.degrees(np.arccos(disk_axis_ratio))
            #bulge_rot_angle = 90 + bulge_angle #because sparcfire gets confused sometimes

            #pos_angle_power = -disk_rot_angle
            #pos_angle_sersic = 90 + disk_angle # Because sparcfire is checking negative direction too, 
                                                                      # assuming symmetry this should align us better
        else:
            inclination = -np.degrees(np.arccos(disk_axis_ratio)) 
            #bulge_rot_angle = 90 - bulge_angle       

            #pos_angle_power  = 180 - disk_rot_angle
            #pos_angle_sersic = 90 - disk_angle 

        #if disk_axis_ratio >= 0.5:
        #    pos_angle_power = pos_angle_power + bulge_angle # Big change, should this also apply to bulge?

        #pos_angle_power = scale_var(pos_angle_power)
        #pos_angle_sersic = pos_angle_sersic - 40 # 40 to adjust for coord rotation addition (see main gen)
        #inclination = inclination

        # Grabbing chirality
        chirality = row.get('chirality_maj', chirality)
        chirality_2 = row.get('chirality_alenWtd', chirality_2)
        chirality_3 = row.get('chirality_longestArc', chirality_3)

        #est_arcs = min(int(est_arcs[1]),int(est_arcs[-2]))

        # bar candidate: r_in = 0 according to Chien
        bar_candidate = row.get('bar_candidate_available')
        bar_len  = scale_var(row.get("bar_half_length_input_img"), scale_fact_ipt)

        #print(pitch_angle)

        # Based on linear regression from test set:
        #galaxy, alpha, pitch angle
        #6698 - 0.9319, 16.8303045
        #1248 - 2.0704, 28.10656079
        #9688 - 0.544, 12.11462932
        #9241 - 0.7194, 19.01479688
        #1827 - 1.2037, 21.55850102
        #4222 - 0.5625, 8.145239069
        #0761 - 1.1330, 19.53474969
        # anddd adjusted according to comparison results 
        alpha = max(0.07*pitch_angle - 0.3, 0.75)

    csv_file.close()
 
    if chirality == chirality_2 or chirality == chirality_3:
        if chirality == 'Z-wise':
            spin_parity = '-'
        elif chirality == 'S-wise':
            spin_parity = ''
        else:
            print(f"Couldn't use chirality.") #Proceeding with coin flip.")
    elif chirality_2 == chirality_3:
        if chirality_2 == 'Z-wise':
            spin_parity = '-'
        elif chirality_2 == 'S-wise':
            spin_parity = ''
        else:
            print(f"Couldn't use chirality.") # Proceeding with coin flip.")
    else:
        print("Something went wrong in choosing a chirality!") #Coin flip...")
        spin_parity = '' #random.choice(['','-'])
        
    loc = locals()
    kwargs_out.update({i: loc[i] for i in kwargs_out.keys() if i in loc})
    return kwargs_out


# In[ ]:


def arc_information(galaxy_name, galaxy_path, num_arms = 2, bulge_rad = 2):

    kwargs_out = {
        "inner_rad" : 0,
        "outer_rad" : 20,
        "cumul_rot" : 60
                 }
    
    try:
        #arcs_filename = glob.glob(galaxy_path + '/' + '*_arcs.csv')
        arcs_filename = path_join(galaxy_path, galaxy_name, '_arcs.csv')
        arcs_file = open(arcs_filename, 'r')
    except:
        print("Can't open to read: ", arcs_filename)
        print("Check Sparcfire output or directories. Proceeding with default values.")
        return kwargs_out

    else:
        reader = csv.DictReader(arcs_file)
        arcs_in = list(reader)
        
        # Nested dictionaries
        theta_sum = 0
        inner_rad = 0
        outer_rad = 0

        i = 0
        count = 0
        
        try:
            winding_dir = float(arcs_in[0]['pitch_angle']) > 0
        except IndexError as ie:
            # For when sparcfire fails
            # TODO: figure out where to put a function that checks if sparcfire failed so I only output
            # one message per galaxy instead of doing it for arcs and galaxy info
            return kwargs_out

        while i < num_arms:
            try:
                _ = arcs_in[i]['pitch_angle']
            except IndexError as ie:
                return kwargs_out
            
            if (float(arcs_in[i]['pitch_angle']) > 0) != winding_dir:
                i += 1
                continue

            # Rudimentary weighting scheme
            weight = (num_arms - i)/num_arms
            weight_2 = weight
            
            try:
                r_start = float(arcs_in[i]['r_start'])
                # Punish sparcfire's bulge-related misgivings
                if scale_var(r_start, scale_fact_std) <= bulge_rad:
                    weight_2 = 0.25
                    
                theta_sum += (np.degrees(float(arcs_in[i]['math_initial_theta'])) % 270) * weight_2
                theta_sum += (np.degrees(float(arcs_in[i]['relative_theta_end']))  % 270) * weight_2 #% 360

                inner_rad += r_start * weight_2
                outer_rad += float(arcs_in[i]['r_end']) * weight
            except ValueError as ve:
                break

            i += 1
            count += 1
        
        # How far the arm rotates in theta which approximates the point where sparcfire no longer traces the arm
        # and hence where the cumulative rotation out point is. Note, this already takes into account the relative
        # starting point of the arms which is the hard part... 
        
        # Averaging, tack on 180 to limit some of the craziness
        # Galfit tends to like the outer_rad smaller, hence larger divisor
        weight_div = 1/max(1, count + 1) #1/np.math.factorial(count)
        cumul_rot = abs(theta_sum)*weight_div # NOT A STRING

        inner_rad = inner_rad*weight_div # Averaging the inner distance to both arcs
        outer_rad = outer_rad*weight_div # Averaging outer distance
        
        arcs_file.close()
        
        try:
            inner_rad = scale_var(inner_rad, scale_fact_std)
            outer_rad = scale_var(outer_rad, scale_fact_std)
        except ValueError as ve:
            pass

    loc = locals()
    kwargs_out.update({i: loc[i] for i in kwargs_out.keys() if i in loc})
    return kwargs_out


# In[ ]:


# def csv_sdss_info(galaxy_names): # to grab petromag and also psf parameter things
    
#     gname_info = {}
    
#     # Defaults go here
#     run = 0
#     rerun = 0
#     cam = 0
#     field = 0
#     row = 0
#     col = 0
#     pmag = 16
    
#     try:
#         #star_dl_filename = glob_name('star_dl','star_dl','.csv')
#         star_dl_filename = path_join('star_dl','star_dl','.csv')
#         star_dl_file = open(star_dl_filename, 'r')
        
#     except:
#         print("Can't open to read star_dl.csv with star information (for PSF)")
#         print("Check Sparcfire output or directories. Proceeding with default values.")
    
#         for gname in galaxy_names:
#             gname_info[gname] = [run, rerun, cam, field, row, col, pmag]
            
#         return gname_info

#     else:
#         star_df = pd.read_csv(star_dl_file, index_col=0)
#         #csv_gal_name = star_df.columns[0]
        
#         for gname in galaxy_names:
#             try:
#                 temp = star_df.at[gname, 'run']

#             except:
#                 #print(star_df.loc[:, 'name'])
#                 print("Can't find the galaxy:", gname, "in our repository.")
#                 print("Proceeding with default values and no PSF.")
#                 run = 0
#                 rerun = 0
#                 camcol = 0
#                 field = 0
#                 rowc = 0
#                 colc = 0
#                 petromag = 16

#             else:
#                 galaxy_band = gname[-1]

#                 run = ' run: '     + str(star_df.at[gname, 'run'])
#                 rerun = ' rerun: ' + str(star_df.at[gname, 'rerun'])
#                 cam = ' camcol: '  + str(star_df.at[gname, 'camcol'])
#                 field = ' field: ' + str(star_df.at[gname, 'field'])
#                 row = ' row: '     + str(star_df.at[gname, 'rowC'])
#                 col = ' col: '     + str(star_df.at[gname, 'colC'])
#                 pmag = str(star_df.at[gname, 'petroMag_' + galaxy_band])
                
#             gname_info[gname] = [run, rerun, cam, field, row, col, pmag]
        
#     return gname_info


# In[ ]:


def write_starmask_ascii(starmask_filepath):
    with fits.open(starmask_filepath) as sm:
        mask_data = sm[0].data
        mask_indices = np.dstack(np.where(mask_data == 0))
        mask_indices = mask_indices.reshape(mask_indices.shape[1], 2).astype(str)

    gname = os.path.basename(starmask_filepath).split("_star-rm")[0]
    starmask_ascii_name = os.path.join(os.path.split(starmask_filepath)[0], f"{gname}_starmask.txt")

    # To avoid some file I/O
    if not mask_indices.size: return starmask_ascii_name

    with open(starmask_ascii_name, "w") as sma:
        _ = [sma.write(' '.join(i) + "\n") for i in mask_indices]

    return starmask_ascii_name

def write_to_feedme(path, list_in, feedme_name = "autogen_feedme_galfit.in"):
    
    file_path = os.path.join(".", path, feedme_name)
#    file_path = "/".join([".", path, ""]) + 'autogen_feedme_galfit.in'
    with open(file_path, 'w') as g:
        #for value in list_in:
        _ = [g.write(f"{value}\n") for value in list_in]
        
    return file_path
# In[ ]:


def write_to_feedmes(in_dir, tmp_dir, out_dir, **kwargs): # single_galaxy_name = "", **kwargs):
    
    # if top_dir:
    #     in_dir, tmp_dir, out_dir = command_line(top_dir)
    # else:
    #     in_dir, tmp_dir, out_dir = command_line()
    # cwd = os.getcwd()
    # in_dir = kwargs.get("in_dir", pj(cwd, "sparcfire-in"))
    # tmp_dir = kwargs.get("tmp_dir", pj(cwd, "sparcfire-tmp"))
    # out_dir = kwargs.get("out_dir", pj(cwd, "sparcfire-out"))
    
    galaxy_names, gfolders = get_galaxy_names_list(in_dir, tmp_dir, out_dir, galaxy_names = kwargs.get("galaxy_names", []))
    
    #psf_info = csv_sdss_info(galaxy_names)
    
    feedme_info_out = {}
    
    # Not worth it speed-wise
    #petromags = kwargs.get("petromags", 16*np.ones(len(gfolders)))
    #bulge_axis_ratios = kwargs.get("bulge_axis_ratios", [])
    
    for i, gfolder in enumerate(gfolders):
    
        # if single_galaxy_name:
        #     gname = single_galaxy_name
        #     #gfolder = pj(out_dir, gname)
        #     # Alternatively
        #     # galaxy = folders_out[galaxy_names.index(gname)]
        # else:
        #     gname = galaxy_names[count]
            
        gname = os.path.basename(gfolder)
        if not exists(gfolder):
            print(f"{gfolder} cannot be found. Continuing...")
            continue
            
        print(gname)
        
        if(os.path.basename(gfolder) != gname):
            print("uh oh naming went wrong in feedme generator.")
            print(f"basename of {gfolder} != {gname}")
            continue
        
        # From old implementation - 1/19/21
        # ************
        # x1crop, x2crop, y1crop, y2crop = autocrop_grab(gname, gfolder)
        # scale = (float(x2crop)-float(x1crop))/256 - from old implementation
        # ************
        
        galaxy_dict = galaxy_information(gname, gfolder)
        
        center_pos_x = float(galaxy_dict["center_pos_x"])
        center_pos_y = float(galaxy_dict["center_pos_y"])
        crop_rad = float(galaxy_dict["crop_rad"])
        
        x1crop = round(center_pos_x - 2*crop_rad)
        x2crop = round(center_pos_x + 2*crop_rad)        
        y1crop = round(center_pos_y - 2*crop_rad)
        y2crop = round(center_pos_y + 2*crop_rad)
    
        arc_dict = arc_information(gname, gfolder, num_arms = galaxy_dict["est_arcs"], bulge_rad = galaxy_dict["bulge_maj_axs_len"])
    
        if galaxy_dict["bar_candidate"].upper() == "FALSE": # According to Chien, if no bar, then r_in = 0 since r_in is more a mathematical construct relating to the bar
            in_rad = 0 #unsure if I'll keep this...
            #pass
        elif galaxy_dict["bar_candidate"].upper() == "TRUE":
            in_rad = min(galaxy_dict["bar_len"], arc_dict["inner_rad"])
        else:
            print(bar_candidate)
            print("bar_candidate is neither TRUE nor FALSE in the CSV. Check sparcfire output.")
            print("Defaulting the average inner distance to the arms.")
        
        tmp_dir_basename = os.path.basename(tmp_dir)
        header = GalfitHeader(input_menu_file = gname,
                              #extra_header_info = f"{run}{camcol}{field}; HDU: z{psf_row}{psf_col}",
                              galaxy_name = gname,
                              input_image = pj(in_dir, f"{gname}.fits"),
                              output_image = pj(tmp_dir, "galfits", f"{gname}_galfit_out.fits"),
                              # Unfortunately have to use a relative path here since GALFIT
                              # breaks when the filename is too long
                              #psf = pj("..", "..", tmp_dir_basename, "psf_files", f"{gname}_psf.fits"),
                              psf = f"{gname}_psf.fits",
                              pixel_mask = pj(tmp_dir, "galfit_masks", f"{gname}_star-rm.fits"),
                              region_to_fit = (x1crop, x2crop, y1crop, y2crop),
                              optimize = 0
                             )
        
        
        petromag = 16 #float(petromags[i])
        bulge_axis_ratio = float(galaxy_dict["bulge_axis_ratio"])
        # Take mag from SDSS and bulge_axis_ratio from SDSS
        with fits.open(pj(in_dir, f"{gname}.fits")) as gf:
            try:
                #SURVEY = SDSS-r  DR7
                color = gf[0].header["SURVEY"].split()[0][-1]
                petromag_str = f"pMag_{color}"
                devab_str = f"devAB_{color}"
                if petromag_str in gf[0].header and devab_str in gf[0].header:
                    petromag = gf[0].header[petromag_str]
                    bulge_axis_ratio = gf[0].header[devab_str]
            except KeyError:
                print(f"There is something wrong with the header for {gname}.")
        
        bulge = Sersic(component_number = 1, 
                       position = (center_pos_x, center_pos_y),
                       magnitude = float(petromag) - 1,
                       # Sometimes sparcfire messes this up
                       effective_radius = min(max(galaxy_dict["bulge_maj_axs_len"], 2), 0.2*crop_rad),
                       # According to other paper GALFIT usually doesn't have a problem with the index
                       sersic_index = 4,
                       axis_ratio = bulge_axis_ratio,
                       position_angle = galaxy_dict["bulge_rot_angle"]
                      )
        
        disk  = Sersic(component_number = 2, 
                       position = (center_pos_x, center_pos_y),
                       magnitude = float(petromag) - 3,
                       effective_radius = galaxy_dict["disk_maj_axs_len"],
                       # According to comparison tests, this usually ends up much higher than classical probably due to the spiral.
                       sersic_index = 1,
                       axis_ratio = galaxy_dict["disk_axis_ratio"],
                       position_angle = galaxy_dict["disk_rot_angle"]
                      )
            
        arms  = Power(component_number = 2,
                      inner_rad = in_rad, # Chosen based on where *detection* of arms usually start
                      outer_rad = arc_dict["outer_rad"],
                      cumul_rot = float(f"{galaxy_dict['spin_parity']}{arc_dict['cumul_rot']}"),
                      powerlaw = galaxy_dict["alpha"],
                      inclination = galaxy_dict["inclination"],
                      sky_position_angle = (galaxy_dict["pos_angle_power"] - galaxy_dict["disk_rot_angle"]) % 180 #90 
                     )
        
        # Take 90 pixels (in the 256x256 image) to be the cutoff for an arm
        # Use a simple cut off for now
        # Looking at the first two arms may be too unreliable
        #print("Max arc length in 256 img", scale_var(max_arc, 0.5*256/crop_rad))
        if scale_var(galaxy_dict["max_arc_length"], 1/scale_fact_std) < 75:
            print("Skipping Arms, max arc len is", galaxy_dict["max_arc_length"]/scale_fact_std)
            arms.add_skip(skip_val = 1)
            
        # May not do the fixing since we run the disk first
        # may however overwrite the output from the disk for this parameter
        # to help the arms in go_go_galfit
        #else:
            # Fixing this to 0.6 to give the arms the best chance to form
        #    disk.axis_ratio = 0.6
        #    disk.param_values["axis_ratio"] = 0.6

        fourier = Fourier(component_number = 2)
            
        sky   = Sky(component_number = 3)
        
        # Previously used for Fourier modes
        
#         f1 = 0.05
#         f3 = 0.02
#         f4 = 0.005
#         f5 = 0.001
#         if est_arcs <= 2:
#             f3 -= 0.015
#             f5 = 0
#             f4 = 0
#         elif est_arcs == 3:
#             f3 += 0.03
#             f4 += 0.005
#             f5 = 0
#         elif est_arcs >= 4:
#             f3 += 0.03
#             f4 += 0.015
#             f5 += 0.014
#         else:
#             print("Something went wrong with the estimation of arcs! Check chirality_votes_maj in csv. Proceeding")
        
        
        # count += 1
        
        container = FeedmeContainer(path_to_feedme = pj(gfolder, f"{gname}.in"),
                                    header         = header,
                                    bulge          = bulge, 
                                    disk           = disk,
                                    arms           = arms,
                                    fourier        = fourier,
                                    sky            = sky)
        
        if arms.param_values.get("skip", False):
            # By default includes the header
            container.to_file(bulge, disk, sky)
        else:
            container.to_file()
        
        feedme_info_out[gname] = container
        
    return feedme_info_out


# In[ ]:


if __name__ == "__main__":
    
    # FOR NOW (aka TODO) force >python 3.7 for f string and subprocess compatibility
    # out_str = """\t Python3.7 or greater required! Exitting without generating feedmes... 
    #             if feedmes have already been generated, galfit will run with those.\n"""
    # assert sys.version_info >= (3, 7), out_str
    
    
    # if in_notebook():
    #     cwd = cwd.replace("ics-home", username)
    #cwd = os.getcwd()
    if len(sys.argv) == 4:
        in_dir = sys.argv[1]
        tmp_dir = sys.argv[2]
        out_dir = sys.argv[3]
    else:
        cwd = os.getcwd()
        in_dir = pj(cwd, "sparcfire-in")
        tmp_dir = pj(cwd, "sparcfire-tmp")
        out_dir = pj(cwd, "sparcfire-out")
    
    print(f"Using: \n{in_dir}\n{tmp_dir}\n{out_dir}\nto generate feedme.")
    write_to_feedmes(
                     in_dir = in_dir, 
                     tmp_dir = tmp_dir, 
                     out_dir = out_dir
                    )


# In[14]:


if __name__ == "__main__":
    # in_notebook() is checked in the export function
    export_to_py("sparc_to_galfit_feedme_gen_debug", output_filename = pj(_MODULE_DIR, "sparc_to_galfit_feedme_gen"))

