#!/usr/bin/env python
# coding: utf-8

# **Author: Matthew Portman**
# 
# **Date (Github date will likely be more accurate): 4/17/23**

# # MINI README
# 
# Please place this script into the directory which contains the folders for your input, temporary, and output files for SpArcFiRe typically denoted: sparcfire-in, sparcfire-tmp, and sparcfire-out. GALFIT will also be run from that same directory if you choose to use the control script included in the repo. I recommend running this alone *once* before running the big script to ensure everything is in its right place (https://www.youtube.com/watch?v=sKZN115n6MI). After you confirm that this works without any errors the first time around, feel free to run the control script from then on. 
# 
# Running from the overarching directory is a temporary measure which will be remedied upon completion of the entire control script and full integration with SpArcFiRe. 
# 
# TO RUN: `python3 sparc_to_galfit_feedme_gen.py`
# 
# To run the control script: `bash control_script.sh`

# In[1]:


import numpy as np
from math import log
import glob
import csv
import subprocess
import random
import pandas as pd
import os
from os.path import join as pj

import sys
from astropy.io import fits

from galfit_objects import *
from copy import deepcopy


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

# In[2]:


# Grabbing filepath from command line
def command_line(top_dir = os.getcwd()): # = True):
    
    in_dir_out = os.path.join(top_dir, "sparcfire-in")
    tmp_dir_out = os.path.join(top_dir, "sparcfire-tmp")
    out_dir_out = os.path.join(top_dir, "sparcfire-out")
        
    #if not run_as_script:
    #    return in_dir_out, tmp_dir_out, out_dir_out
    
    try:
        if len(sys.argv) != 4: # including name of python script
            print(f"No path given. Defaulting to using: \n{in_dir_out}\n{tmp_dir_out}\n{out_dir_out}")

        else:
            in_dir_out = argv[1]
            tmp_dir_out = argv[2]
            out_dir_out = argv[3]
    except:
        pass

    return in_dir_out, tmp_dir_out, out_dir_out


# In[3]:


# Grabbing the file names
def get_galaxy_names_list(path_to_sparc_in):

    try:
        filenames_read = glob.glob(path_to_sparc_in + "/*.fits") # Previously hardcoded
    
    #except:
        # Feel free to hardcode this 
        #print("Please input the full path to the input directory for the SpArcFiRe run and include the trailing /")
        #print("e.g. /home/usr/sparcfire-in/")
        #path_to_sparc_in = input()
    
        #try:
            #filenames_in = glob.glob(path_to_sparc_in + "*.fits")

        #except:
            #print("See instructions above or check your directory path. Could not glob.")
            #raise SystemExit("Exitting.")
    
    except:
        print("Please copy me into the directory which contains the folders for your")
        print("input, temporary, and output files for SpArcFiRe denoted:")
        print("sparcfire-in, sparcfire-tmp, and sparcfire-out.")
        raise SystemExit("Exitting.")
        
    else:
        filenames_out = [os.path.splitext(s)[0] for s in filenames_read]
        galaxy_names_out = [os.path.basename(s) for s in filenames_out]
        filenames_out = [s.replace("-in", "-out") for s in filenames_out]
        
    return filenames_read, galaxy_names_out, filenames_out


# In[4]:


def path_join(path='.', name='', file_ext=''):
    
    file_path = os.path.join(path, name + file_ext)

    # deprecated
    # "./" + path + '/' + name + file_ext
    # print(file_path)
    # print(glob.glob(file_path))
    # file_name = glob.glob(file_path)[0]
    
    return file_path


# In[5]:


def scale_var(x, scale = 1):
    # Scales
    return float(x)*scale


# In[6]:


def galaxy_information(galaxy_name, galaxy_path):
   
    bulge_rad_out = 2
    bulge_axis_ratio_out = 0.5
    bulge_rot_angle_out = 1

    crop_rad_out = 30 # New!

    center_pos_x_out = 30
    center_pos_y_out = 30
    disk_maj_axs_len_out = 30
    pos_angle_sersic_out = 1
    pos_angle_power_out = 30
    axis_ratio_out = 0.5
    max_arc_length_out = 30
    chirality = 0
    chirality_2 = 0
    a_ratio = 1
    alpha_out = 1
    # spin_parity handled by if 
    est_arcs_out = 2
    inclination = 30
    bar_cand = 'FALSE'
    spin_parity = random.choice(['','-'])
    
    try:
        #tsv_filename = glob.glob("./" + galaxy_path + '/' + galaxy_name + '.tsv')[0]
        #csv_filename = glob_name(galaxy_path, galaxy_name, '.csv') # Changing to csv to be consistent
        csv_filename = path_join(galaxy_path, galaxy_name, '.csv') # Changing to csv to be consistent
        csv_file = open(csv_filename, 'r')
        
    except:
        print("Can't open to read: ", csv_filename)
        print("Check Sparcfire output or directories. Proceeding with default values.")

    else:
        reader = csv.DictReader(csv_file, skipinitialspace = True)
#        csv_in = list(reader)
        for row in reader:
            # if sparcfire fails
            if "input rejected" in row['fit_state']:
                break
                
            # Already accounted for, no need to scale
            center_pos_x_out = row['inputCenterR']
            center_pos_y_out = row['inputCenterC']
            
            crop_rad_out = row['cropRad']
            
            global scale_fact # Making this global since I'm now grabbing the necessary info from the csv
            scale_fact = 2*float(crop_rad_out)/256
            
            bulge_rad_out = scale_var(row['bulgeMajAxsLen'], scale_fact)
            bulge_axis_ratio_out = row['bulgeAxisRatio']
            
            # Scaled down, may scale down a bit further to better reflect the *half* radius
            disk_maj_axs_len_out = scale_var(row['diskMajAxsLen'], scale_fact) 
            
            # Angles don't need to scale
            bulge_angle = np.degrees(float(row['bulgeMajAxsAngle']))
            disk_angle = np.degrees(float(row['diskMajAxsAngleRadians']))
            
            a_ratio = float(row['diskAxisRatio'])
            
            if disk_angle < 0: 
                inclination = np.degrees(np.arccos(a_ratio))
                bulge_rot_angle_out = 90 + bulge_angle #because sparcfire gets confused sometimes
                
                pos_angle_power_out = -disk_angle
                pos_angle_sersic_out = 90 + disk_angle # Because sparcfire is checking negative direction too, 
                                                                          # assuming symmetry this should align us better
            else:
                inclination = -np.degrees(np.arccos(a_ratio)) 
                bulge_rot_angle_out = 90 - bulge_angle       
                
                pos_angle_power_out = 180 - disk_angle
                pos_angle_sersic_out = 90 - disk_angle 
                
            if a_ratio >= 0.5:
                pos_angle_power_out = pos_angle_power_out + bulge_angle # Big change, should this also apply to bulge?
            
            pos_angle_power_out = scale_var(pos_angle_power_out)
            pos_angle_sersic_out = pos_angle_sersic_out - 40 # 40 to adjust for coord rotation addition (see main gen)
            #inclination = inclination

            # GRABBING ARC STATISTICS
            avg_arc_length_out = scale_var(row['avgArcLength'], scale_fact)
            #med_arc_length = row['medianArcLength'][:6]
            max_arc_length_out = scale_var(row['maxArcLength'], scale_fact)
            
            # Grabbing chirality
            chirality = row['chirality_maj']
            chirality_2 = row['chirality_alenWtd']
            chirality_3 = row['chirality_longestArc']
            
            # For estimating number of real arcs to influence fourier modes
            # This seems to be a much better way to do it
            est_arcs_out = int(row['rankAt50pct'])
            #est_arcs_out = min(int(est_arcs_out[1]),int(est_arcs_out[-2]))
            
            # bar candidate: r_in = 0 according to Chien
            bar_cand = row['bar_candidate_available']
            
            # Grabbing PA absolute average from dominant chirality
            pitch_angle = abs(float(row['pa_alenWtd_avg_domChiralityOnly']))
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
            alpha_out = 0.07*pitch_angle - 0.3
            
        csv_file.close()
 
    if chirality == chirality_2 or chirality == chirality_3:
        if chirality == 'Z-wise':
            spin_parity = '-'
        elif chirality == 'S-wise':
            spin_parity = ''
        else:
            print(f"Couldn't use {chirality}. Proceeding with coin flip.")
    elif chirality_2 == chirality_3:
        if chirality_2 == 'Z-wise':
            spin_parity = '-'
        elif chirality_2 == 'S-wise':
            spin_parity = ''
        else:
            print(f"Couldn't use {chirality}. Proceeding with coin flip.")
    else:
        print("Something went wrong in choosing a chirality! Coin flip...")
        spin_parity = random.choice(['','-'])
        
    return (
        bulge_rad_out,
        bulge_axis_ratio_out,
        bulge_rot_angle_out,
        crop_rad_out,
        center_pos_x_out,
        center_pos_y_out,
        disk_maj_axs_len_out,
        pos_angle_sersic_out,
        pos_angle_power_out,
        a_ratio,
        max_arc_length_out,
        spin_parity,
        est_arcs_out,
        inclination,
        bar_cand,
        alpha_out
            )


# In[7]:


def arc_information(galaxy_name, galaxy_path, num_arms = 2):

    inner_rad = 0
    outer_rad = 20
    cumul_rot_out = 60
    
    try:
        #arcs_filename = glob.glob(galaxy_path + '/' + '*_arcs.csv')
        arcs_filename = path_join(galaxy_path, galaxy_name, '_arcs.csv')
        arcs_file = open(arcs_filename, 'r')
    except:
        print("Can't open to read: ", arcs_filename)
        print("Check Sparcfire output or directories. Proceeding with default values.")

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
            return inner_rad, outer_rad, cumul_rot_out #, alpha_out

        while count < num_arms:
            try:
                _ = arcs_in[i]['pitch_angle']
            except IndexError as ie:
                return inner_rad, outer_rad, cumul_rot_out #, alpha_out
            
            if (float(arcs_in[i]['pitch_angle']) > 0) != winding_dir:
                i += 1
                continue

            # Rudimentary weighting scheme
            weight = num_arms - count
            theta_sum += (np.degrees(float(arcs_in[i]['math_initial_theta'])) % 360) * weight
            theta_sum += np.degrees(float(arcs_in[i]['relative_theta_end']))*weight #% 360
            
            inner_rad += float(arcs_in[i]['r_start'])*weight
            outer_rad += float(arcs_in[i]['r_end'])*weight

            i += 1
            count += 1
        
        # How far the arm rotates in theta which approximates the point where sparcfire no longer traces the arm
        # and hence where the cumulative rotation out point is. Note, this already takes into account the relative
        # starting point of the arms which is the hard part... 
        
        # Averaging, tack on 180 to limit some of the craziness
        weight_div = 1/np.math.factorial(num_arms)
        cumul_rot_out = max((abs(theta_sum)*weight_div) % 180, 60.0) # NOT A STRING
 
        inner_rad = inner_rad*weight_div # Averaging the inner distance to both arcs
        outer_rad = outer_rad*weight_div # Averaging outer distance
        
        arcs_file.close()
                              
        inner_rad = scale_var(inner_rad, scale_fact)
        outer_rad = scale_var(outer_rad, scale_fact)

    return inner_rad, outer_rad, cumul_rot_out #, alpha_out


# In[8]:


def csv_sdss_info(galaxy_names): # to grab petromag and also psf parameter things
    
    gname_info = {}
    
    # Defaults go here
    run = 0
    rerun = 0
    cam = 0
    field = 0
    row = 0
    col = 0
    pmag = 16
    
    try:
        #star_dl_filename = glob_name('star_dl','star_dl','.csv')
        star_dl_filename = path_join('star_dl','star_dl','.csv')
        star_dl_file = open(star_dl_filename, 'r')
        
    except:
        print("Can't open to read star_dl.csv with star information (for PSF)")
        print("Check Sparcfire output or directories. Proceeding with default values.")
    
        for gname in galaxy_names:
            gname_info[gname] = [run, rerun, cam, field, row, col, pmag]
            
        return gname_info

    else:
        star_df = pd.read_csv(star_dl_file, index_col=0)
        #csv_gal_name = star_df.columns[0]
        
        for gname in galaxy_names:
            try:
                temp = star_df.at[gname, 'run']

            except:
                #print(star_df.loc[:, 'name'])
                print("Can't find the galaxy:", gname, "in our repository.")
                print("Proceeding with default values and no PSF.")
                run = 0
                rerun = 0
                camcol = 0
                field = 0
                rowc = 0
                colc = 0
                petromag = 16

            else:
                galaxy_band = gname[-1]

                run = ' run: '     + str(star_df.at[gname, 'run'])
                rerun = ' rerun: ' + str(star_df.at[gname, 'rerun'])
                cam = ' camcol: '  + str(star_df.at[gname, 'camcol'])
                field = ' field: ' + str(star_df.at[gname, 'field'])
                row = ' row: '     + str(star_df.at[gname, 'rowC'])
                col = ' col: '     + str(star_df.at[gname, 'colC'])
                pmag = str(star_df.at[gname, 'petroMag_' + galaxy_band])
                
            gname_info[gname] = [run, rerun, cam, field, row, col, pmag]
        
    return gname_info


# In[9]:


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
# In[17]:


def write_to_feedmes(top_dir = ""):
    
    if top_dir:
        in_dir, tmp_dir, out_dir = command_line(top_dir)
    else:
        in_dir, tmp_dir, out_dir = command_line()

    filenames_fits_in, galaxy_names, folders_out = get_galaxy_names_list(in_dir)
    
    psf_info = csv_sdss_info(galaxy_names)
    
    count = 0
    paths_to_feedme = {}
    
    for galaxy in folders_out:
    
        gname = galaxy_names[count]
        print(gname)
        
        if(os.path.basename(galaxy) != gname):
            print("uh oh naming went wrong")
            sys.exit()
        
        # From old implementation - 1/19/21
        # ************
        # x1crop, x2crop, y1crop, y2crop = autocrop_grab(gname, galaxy)
        # scale = (float(x2crop)-float(x1crop))/256 - from old implementation
        # ************
        
        bulge_rad, bulge_axis_ratio, pos_angle_bulge, \
            crop_rad, center_pos_x, center_pos_y, \
            disk_maj_axs_len, pos_angle_disk, pos_angle_power, \
            axis_ratio, max_arc, spin_dir, \
            est_arcs, inclination, bar_candidate, \
            alpha = galaxy_information(gname, galaxy)
        
        center_pos_x = float(center_pos_x)
        center_pos_y = float(center_pos_y)
        crop_rad = float(crop_rad)
        
        x1crop = round(center_pos_x - crop_rad)
        x2crop = round(center_pos_x + crop_rad)        
        y1crop = round(center_pos_y - crop_rad)
        y2crop = round(center_pos_y + crop_rad)
    
        in_rad, out_rad, cumul_rot = arc_information(gname, galaxy, num_arms = est_arcs)
    
        #cumul_rot = (cumul_rot + float(pos_angle_disk)) % 360 #+ float(pos_angle_power) 
    
        if bar_candidate.upper() == "FALSE": # According to Chien, if no bar, then r_in = 0 since r_in is more a mathematical construct relating to the bar
            in_rad = 0 #unsure if I'll keep this...
            #pass
        elif bar_candidate.upper() == "TRUE":
            pass
        else:
            print(bar_candidate)
            print("bar_candidate is neither TRUE nor FALSE in the TSV. Check sparcfire output.")
            print("Defaulting the average inner distance to the arms.")
        
        #To reconstruct the z PSF (i.e., the 5th HDU) at the position (row, col) = (500, 600) from run 1336, column 2, field 51 you’d say:
        #read_PSF psField-001336-2-0051.fit 5 500.0 600.0 foo.fit
        run, rerun, camcol, field, psf_row, psf_col, petromag = psf_info[gname]
        
        header = GalfitHeader(input_menu_file = gname,
                              extra_header_info = f"{run}{camcol}{field}; HDU: z{psf_row}{psf_col}",
                              galaxy_name = gname,
                              input_image = f"{filenames_fits_in[count]}",
                              output_image = pj(tmp_dir, "galfits", f"{gname}_out.fits"),
                              pixel_mask = pj(tmp_dir, "galfit_masks", f"{gname}_star-rm.fits"),
                              region_to_fit = (x1crop, x2crop, y1crop, y2crop),
                              optimize = 0
                             )      
        
        bulge = Sersic(component_number = 1, 
                       position = (center_pos_x, center_pos_y),
                       magnitude = float(petromag) - 1,
                       effective_radius = bulge_rad,
                       # According to other paper GALFIT usually doesn't have a problem with the index
                       sersic_index = 4,
                       axis_ratio = axis_ratio,
                       position_angle = pos_angle_bulge
                      )
        
        disk  = Sersic(component_number = 2, 
                       position = (center_pos_x, center_pos_y),
                       magnitude = float(petromag) - 3,
                       effective_radius = disk_maj_axs_len,
                       # According to comparison tests, this usually ends up much higher than classical probably due to the spiral.
                       sersic_index = 4,
                       # Fixing this to 0.6 to give the arms the best chance to form
                       axis_ratio = 0.6,
                       position_angle = pos_angle_disk
                      )       
        
        arms  = Power(inner_rad = in_rad, # Chosen based on where *detection* of arms usually start
                      outer_rad = out_rad,
                      cumul_rot = float(f"{spin_dir}{cumul_rot}"),
                      powerlaw = alpha,
                      inclination = inclination,
                      sky_position_angle = 90 # pos_angle_power
                     )
        
        fourier = Fourier()
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
        
        
        count += 1
        
        header.to_file(pj(galaxy, f"{gname}.in"), 
                       bulge, 
                       disk, 
                       arms, 
                       fourier, 
                       sky)
        
        paths_to_feedme[gname] = pj(galaxy, f"{gname}.in")
        
    return paths_to_feedme
        #write_to_feedme(galaxy, bulge_feedme, feedme_name = gname + "_bulge.in")
        #write_to_feedme(galaxy, disk_feedme, feedme_name = gname + "_disk.in")
        #paths_to_feedme.append(write_to_feedme(galaxy, formatted_feedme, feedme_name = gname + ".in")) # do I need paths_to_feedme? I used to use it for something...


# In[ ]:


if __name__ == "__main__":
    
    # FOR NOW (aka TODO) force >python 3.6 for f string compatibility
    out_str = """\t Python3.6 or greater required! Exitting without generating feedmes... 
                if feedmes have already been generated, galfit will run with those.\n"""
    assert sys.version_info >= (3, 6), out_str
    
    write_to_feedmes(top_dir = "/home/portmanm/run6_1000_galfit_two_fit")


# In[15]:


if __name__ == "__main__":
    export_to_py("notebook_feedme_gen", output_filename = "sparc_to_galfit_feedme_gen.py")

