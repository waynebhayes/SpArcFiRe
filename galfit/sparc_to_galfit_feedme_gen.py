#!/usr/bin/env python
# coding: utf-8

# **Author: Matthew Portman**
# 
# **Date (Github date will likely be more accurate): 10/22/20**

# # MINI README
# 
# Please place this script into the directory which contains the folders for your input, temporary, and output files for SpArcFiRe typically denoted: sparcfire-in, sparcfire-tmp, and sparcfire-out. GALFIT will also be run from that same directory if you choose to use the control script included in the repo. I recommend running this alone *once* before running the big script to ensure everything is in its right place. After you confirm that this works without any errors the first time around, feel free to run the control script from then on. 
# 
# Running from the overarching directory is a temporary measure which will be remedied upon completion of the entire control script and full integration with SpArcFiRe. 
# 
# TO RUN: `python sparc_to_galfit_feedme_gen.py`
# 
# To run the control script: `bash control_script.sh`

# In[ ]:


import numpy as np
from math import log
import glob
import csv
import subprocess
import random
import pandas as pd


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


# Grabbing the file names
def get_galaxy_names_list():

    try:
        filenames_read = glob.glob("sparcfire-in/*.fits") # Hardcoding is a temporary measure.
    
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
        print("input, temporary, and output files for SpArcFiRe typically denoted:")
        print("sparcfire-in, sparcfire-tmp, and sparcfire-out.")
        raise SystemExit("Exitting.")
        
    else:
        filenames_out = [s.split(".")[0] for s in filenames_read]
        galaxy_names_out = [s.split("/")[1] for s in filenames_out]
        filenames_out = [s.replace("in", "out") for s in filenames_out]
        
    return filenames_read, galaxy_names_out, filenames_out


# In[ ]:


def glob_name(path='', name='', desired_file=''):
    
    file_path = "./" + path + '/' + name + desired_file
    #print(file_path)
    file_name = glob.glob(file_path)[0]
    
    return file_name


# In[ ]:


def autocrop_grab(galaxy_name, galaxy_path):
    # for autocrop coordinates
    
    try: 
        cropcoord_filename = glob_name(galaxy_path, galaxy_name, '_crop_coord.txt') #glob.glob('*coord.txt')[0]
        cropcoord_file = open(cropcoord_filename,'r')

    except:
        print("Can't open to read: ", cropcoord_filename)
        print("Check Sparcfire output or directories. Proceeding without crop (good luck).")

        crop_x1 = 0
        crop_x2 = 256
        crop_y1 = 0
        crop_y2 = 256

    else:
        cropcoord_in = cropcoord_file.read()
        cropcoord_file.close()

        cropcoord_in = cropcoord_in.split('\n')

        crop_x1 = cropcoord_in[0]
        crop_x2 = cropcoord_in[1]
        crop_y1 = cropcoord_in[2]
        crop_y2 = cropcoord_in[3]
        
    return crop_x1, crop_x2, crop_y1, crop_y2


# In[ ]:


def scale_to_string(x, scale = 1):
    # Scales and chops
    return str(float(x)*scale)[:6]


# In[ ]:


def tsv_information(galaxy_name, galaxy_path, scale_fact):
   
    try:
        #tsv_filename = glob.glob("./" + galaxy_path + '/' + galaxy_name + '.tsv')[0]
        tsv_filename = glob_name(galaxy_path, galaxy_name, '.tsv')
        tsv_file = open(tsv_filename, 'r')
        
    except:
        print("Can't open to read: ", tsv_filename)
        print("Check Sparcfire output or directories. Proceeding with default values.")

        bulge_rad_out = '2'
        bulge_axis_ratio_out = '0.5'
        bulge_rot_angle_out = '1'

        center_pos_x_out = '30'
        center_pos_y_out = '30'
        disk_maj_axs_len_out = "30"
        pos_angle_sersic_out = "1"
        pos_angle_power_out = '30'
        axis_ratio_out = "0.5"
        max_arc_length_out = '30'
        chirality = 0
        chirality_2 = 0
        # spin_parity handled by if 
        est_arcs_out = 2
        inclination = '30'
        bar_cand = 'FALSE'

    else:
        reader = csv.DictReader(tsv_file, dialect='excel-tab')
        tsv_in = list(reader)
        for row in tsv_in:
            bulge_rad_out = scale_to_string(row['bulgeMajAxsLen'], scale_fact)
            bulge_axis_ratio_out = row['bulgeAxisRatio']         

            # Already accounted for, no need to scale
            center_pos_x_out = row['inputCenterR']
            center_pos_y_out = row['inputCenterC']
            
            # Scaled down, may scale down a bit further to better reflect the *half* radius
            disk_maj_axs_len_out = scale_to_string(row['diskMajAxsLen'], scale_fact) 
            
            # Angles don't need to scale
            bulge_angle = np.degrees(float(row['bulgeMajAxsAngle']))
            disk_angle = np.degrees(float(row['diskMajAxsAngleRadians']))
            
            a_ratio = float(row['diskAxisRatio']) # I use a_ratio later hence why it's not immediately converted
            axis_ratio_out = scale_to_string(a_ratio)
            
            if disk_angle < 0: 
                inclination = np.degrees(np.arccos(a_ratio))
                bulge_rot_angle_out = scale_to_string(90 + bulge_angle) #because sparcfire gets confused sometimes
                
                pos_angle_power_out = -disk_angle
                pos_angle_sersic_out = 90 + disk_angle # Because sparcfire is checking negative direction too, 
                                                                          # assuming symmetry this should align us better
            else:
                inclination = -np.degrees(np.arccos(a_ratio)) 
                bulge_rot_angle_out = scale_to_string(90 - bulge_angle)            
                
                pos_angle_power_out = 180 - disk_angle
                pos_angle_sersic_out = 90 - disk_angle 
                
            if a_ratio >= 0.5:
                pos_angle_power_out = pos_angle_power_out + bulge_angle # Big change, should this also apply to bulge?
            
            pos_angle_power_out = scale_to_string(pos_angle_power_out)
            pos_angle_sersic_out = scale_to_string(pos_angle_sersic_out)
            inclination = scale_to_string(inclination)

            # GRABBING ARC STATISTICS
            avg_arc_length_out = scale_to_string(row['avgArcLength'], scale_fact)
            #med_arc_length = row['medianArcLength'][:6]
            max_arc_length_out = scale_to_string(row['maxArcLength'], scale_fact)
            
            # Grabbing chirality
            chirality = row['chirality_maj']
            chirality_2 = row['chirality_alenWtd']
            
            # For estimating number of real arcs to influence fourier modes
            est_arcs_out = row['chirality_votes_maj']
            est_arcs_out = min(int(est_arcs_out[1]),int(est_arcs_out[-2]))
            
            # bar candidate: r_in = 0 according to Chien
            bar_cand = row['bar_candidate_available']
            
            # Grabbing PA absolute average from dominant chirality
            pitch_angle = abs(float(row['pa_alenWtd_avg_domChiralityOnly']))
            #print(pitch_angle)
            
            # Based on linear regression from test set: galaxy, alpha, pitch angle
            #6698 - 0.9319, 16.8303045
            #1248 - 2.0704, 28.10656079
            #9688 - 0.544, 12.11462932
            #9241 - 0.7194, 19.01479688
            #1827 - 1.2037, 21.55850102
            #4222 - 0.5625, 8.145239069
            #0761 - 1.1330, 19.53474969
            alpha_out = scale_to_string(0.07*pitch_angle - 0.3)
            
        tsv_file.close()
 
    if chirality == 'Zwise':
        spin_parity = '-'
    elif chirality == 'Swise':
        spin_parity = ''
    else:
        if chirality_2 == 'Zwise':
            spin_parity = '-'
        elif chirality_2 == 'Swise':
            spin_parity = ''
        else:
            print("Something went wrong in choosing a chirality! Coin flip...")
            spin_parity = random.choice(['','-'])
        
    return (
        bulge_rad_out,
        bulge_axis_ratio_out,
        bulge_rot_angle_out,
        center_pos_x_out,
        center_pos_y_out,
        disk_maj_axs_len_out,
        pos_angle_sersic_out,
        pos_angle_power_out,
        axis_ratio_out,
        max_arc_length_out,
        spin_parity,
        est_arcs_out,
        inclination,
        bar_cand,
        alpha_out
            )


# In[ ]:


def csv_information(galaxy_name, galaxy_path, scale_fact):

    try:
        #arcs_filename = glob.glob(galaxy_path + '/' + '*_arcs.csv')
        arcs_filename = glob_name(galaxy_path, galaxy_name, '_arcs.csv')
        arcs_file = open(arcs_filename, 'r')
    except:
        print("Can't open to read: ", arcs_filename)
        print("Check Sparcfire output or directories. Proceeding with default values.")
        
        inner_rad = '0'
        outer_rad = '20'
        cumul_rot_out = '60'

    else:
        reader = csv.DictReader(arcs_file)
        arcs_in = list(reader)
        arm_1_dict = arcs_in[0]
        
        try:
            arm_2_dict = arcs_in[1]
        
        except:
            print("Only found one arm. Proceeding.")
            arm_2_dict = arm_1_dict
        
        arc_init_theta = np.degrees(float(arm_1_dict['math_initial_theta'])) % 360
        arc_init_theta_2 = np.degrees(float(arm_2_dict['math_initial_theta'])) % 360
        
        arc_end_theta = np.degrees(float(arm_1_dict['relative_theta_end'])) % 360
        arc_end_theta_2 = np.degrees(float(arm_2_dict['relative_theta_end'])) % 360
        
        # How far the arm rotates in theta which approximates the point where sparcfire no longer traces the arm
        # and hence where the cumulative rotation out point is. Note, this already takes into account the relative
        # starting point of the arms which is the hard part... 
        # JK sparc is measuring as from its auto-crop and face-on transform. I can't use that!
        
        cumul_rot_out = 0.5*abs((arc_end_theta + arc_end_theta_2) - (arc_init_theta + arc_init_theta_2)) # NOT A STRING
        
        inner_rad = 0.5*(float(arm_1_dict['r_start']) + float(arm_2_dict['r_start'])) # Averaging the inner distance to both arcs            
        outer_rad = 0.5*(float(arm_1_dict['r_end']) + float(arm_2_dict['r_start'])) # Averaging outer distance
        
        #pitch_angle = 0.5*abs(float(arm_1_dict['pitch_angle']) + float(arm_2_dict['pitch_angle']))
        #print(pitch_angle)
        
        #alpha_out = scale_to_string(0.06*pitch_angle - 0.12)
        
    
        try:
            cumul_rot_out = abs((arc_init_theta % 360)/(log(inner_rad)/log(outer_rad)))
        except:
            print("Couldn't calculate theoretical theta_out. Proceeding.")
        
        arcs_file.close()
                              
        inner_rad = scale_to_string(inner_rad, scale_fact)
        outer_rad = scale_to_string(outer_rad, scale_fact)

    return inner_rad, outer_rad, cumul_rot_out #, alpha_out


# In[ ]:


def csv_sdss_info(galaxy_name): # to grab petromag and also psf parameter things
    
    try:
        #star_dl_filename = glob_name('star_dl','star_dl','.csv')
        star_dl_filename = glob_name('star_dl','star_dl','.csv')
        star_dl_file = open(star_dl_filename, 'r')
        
    except:
        print("Can't open to read csv with star information used for Galaxy ", galaxy_name)
        print("Check Sparcfire output or directories. Proceeding with default values.")

        # Defaults go here
        run_out = '0'
        rerun_out = '0'
        camcol_out = '0'
        field_out = '0'
        rowc_out = '0'
        colc_out = '0'
        petromag_out = '16'

    else:
        star_df = pd.read_csv(star_dl_file, index_col=0)
        #csv_gal_name = star_df.columns[0]
        try:
            temp = star_df.at[galaxy_name, 'run']
            
        except:
            #print(star_df.loc[:, 'name'])
            print("Can't find the galaxy: ", galaxy_name, " in our repository.")
            print("Proceeding with default values and no PSF.")
            run_out = '0'
            rerun_out = '0'
            camcol_out = '0'
            field_out = '0'
            rowc_out = '0'
            colc_out = '0'
            petromag_out = '16'
            
        else:
            galaxy_band = galaxy_name[-1]
            
            run_out = ' run: ' + scale_to_string(star_df.at[galaxy_name, 'run'])
            rerun_out = ' rerun: ' + scale_to_string(star_df.at[galaxy_name, 'rerun'])
            camcol_out = ' camcol: ' + scale_to_string(star_df.at[galaxy_name, 'camcol'])
            field_out = ' field: ' + scale_to_string(star_df.at[galaxy_name, 'field'])
            rowc_out = ' row: ' + scale_to_string(star_df.at[galaxy_name, 'rowC'])
            colc_out = ' col: ' + scale_to_string(star_df.at[galaxy_name, 'colC'])
            petromag_out = scale_to_string(star_df.at[galaxy_name, 'petroMag_' + galaxy_band])
        
        #csv_in = list(reader)
        #if galaxy_name in csv_in:
        #    gal_index = csv_in.index(galaxy_name)
        
        #    for row in csv_in:
        #        if galaxy_name in row:
        #            run_out = row['run']
        #            rerun_out = row['rerun']
        #            camcol_out = ['camcol']
        #            field_out = ['field']
        #            petromag = row['PetroMag_r']
        
    return run_out, rerun_out, camcol_out, field_out, rowc_out, colc_out, petromag_out


# In[ ]:


def write_to_feedme(path, list_in):
    
    file_path = "/".join([".", path, ""]) + 'autogen_feedme_galfit.in'
    galfit_in = open(file_path, 'w')

    for value in feedme_list:
        galfit_in.write(value + '\n')

    galfit_in.close()
    return file_path


# In[ ]:


if __name__ == "__main__":

    count = 0
    paths_to_feedme = []
    
    filenames_in, galaxy_names, filenames_out = get_galaxy_names_list()
    
    for galaxy in filenames_out:
    
        gname = galaxy_names[count]
        #print(gname)
        
        x1crop, x2crop, y1crop, y2crop = autocrop_grab(gname, galaxy)
    
        scale = (float(x2crop)-float(x1crop))/256
    
        bulge_rad, bulge_axis_ratio, pos_angle_bulge,             center_pos_x, center_pos_y,             disk_maj_axs_len, pos_angle_disk, pos_angle_power,             axis_ratio, max_arc, spin_dir,             est_arcs, inclination, bar_candidate,             alpha = tsv_information(gname, galaxy, scale)
    
        in_rad, out_rad, cumul_rot = csv_information(gname, galaxy, scale)
    
        run, rerun, camcol, field, psf_row, psf_col, petromag = csv_sdss_info(gname) 
    
        cumul_rot = scale_to_string((cumul_rot + float(pos_angle_disk)) % 360) #+ float(pos_angle_power) 
    
        if bar_candidate.upper() == "FALSE": # According to Chien, if no bar, then r_in = 0 since r_in is more a mathematical construct relating to the bar
            in_rad = "0" #unsure if I'll keep this...
            #pass
        elif bar_candidate.upper() == "TRUE":
            pass
        else:
            print(bar_candidate)
            print("bar_candidate is neither TRUE nor FALSE in the TSV. Check sparcfire output.")
            print("Defaulting the average inner distance to the arms.")
    
    
        # Initializing Feedme
        feedme_list = []
    
        #To reconstruct the z PSF (i.e., the 5th HDU) at the position (row, col) = (500, 600) from run 1336, column 2, field 51 you’d say:
        #read_PSF psField-001336-2-0051.fit 5 500.0 600.0 foo.fit

        feedme_list.append("#" + run + camcol + field + "; HDU: z" + psf_row + psf_col)
        feedme_list.append("")
        # Image and Galfit Control Param
        feedme_list.append("A) ./" + filenames_in[count]) 
        feedme_list.append("B) ./sparcfire-tmp/galfits/" + gname + "_out.fits")
        feedme_list.append("C) none")
        feedme_list.append("D) none")
        feedme_list.append("E) 1")
        feedme_list.append("F) ./sparcfire-tmp/galfit_masks/" + gname + "_star-rm.fits") #I hate to make this verbose but I have no other choice in this case...
        feedme_list.append("G) ./constraints.txt")
        feedme_list.append("H) " + x1crop + " " + x2crop + " " + y1crop + " " + y2crop)
        feedme_list.append("K) 10  10") 
        feedme_list.append("J) 24.8") # SDSS
        feedme_list.append("K) 0.396  0.396") # SDSS
        feedme_list.append("O) regular")
        feedme_list.append("P) 0")
        feedme_list.append("")
    
        # Sersic 1
        # Fixing as much as I can here... it's not exactly a priority.
        feedme_list.append("# Component number: 1")
        feedme_list.append("0) sersic")
        feedme_list.append("1) " + center_pos_x + " " + center_pos_y + "   0  0")
        feedme_list.append("3) " + str(float(petromag) + 3.0) + " 1") # Initial guess goes here
        feedme_list.append("4) " + bulge_rad + "0") 
        feedme_list.append("5) 4  1") # According to other paper GALFIT usually doesn't have a problem with the index
        feedme_list.append("6) 0  0")    
        feedme_list.append("7) 0  0")    
        feedme_list.append("8) 0  0")    
        # According to other papers, bulge (esp. in spiral galaxies) averages to about 2 so this is a good starting place
        # see https://ned.ipac.caltech.edu/level5/Sept11/Buta/Buta9.html
        feedme_list.append("9) " + axis_ratio + " 0")  
        feedme_list.append("10) " + pos_angle_bulge + " 0") 
        feedme_list.append("")
    
        # Sersic 2
        feedme_list.append("# Component number: 2")
        feedme_list.append("0) sersic")
        feedme_list.append("1) " + center_pos_x + " " + center_pos_y + "   0  0")
        feedme_list.append("3) " + petromag + "  1") # Initial guess goes here
        feedme_list.append("4) " + disk_maj_axs_len + " 1") # Use this for effective radius? Also 0 this one out? Will have to see how well it works in practice
        feedme_list.append("5) 1  1") # Classical disk follows Sersic n = 1 so good place to start (per Readme Exponential profile)
        feedme_list.append("6) 0  0")    
        feedme_list.append("7) 0  0")    
        feedme_list.append("8) 0  0")    
        feedme_list.append("9) " + axis_ratio + " 1") # Exactly what we're looking for... may have to adjust due to spiral arm degeneracy concerns (m=2)
        feedme_list.append("10) " + pos_angle_disk + " 1") #90  1") # fixing this to 'normal' 0 so that we can JUST rotate power function
        feedme_list.append("")
    
        # Power
        feedme_list.append("R0) power")
        feedme_list.append("R1) " + in_rad + "  0") # Chosen based on where *detection* of arms usually start
        feedme_list.append("R2) " + out_rad + "  0")
        feedme_list.append("R3) " + spin_dir + cumul_rot + "  1") # See calc above
        feedme_list.append("R4) " + alpha + "  1") # Another good thing to automate via Sparcfire 
        feedme_list.append("R9) " + inclination + " 1") # see if can't 0 this one out... 
        feedme_list.append("R10) 0  1")#+ pos_angle_power + " 1") # Always more to discover, looks like all the images are mirrored across the y axis.

        f1 = 0.05
        f3 = 0.02
        f4 = 0.005
        f5 = 0.001
        if est_arcs <= 2:
            f3 -= 0.015
            f5 = 0
            f4 = 0
        elif est_arcs == 3:
            f3 += 0.03
            f4 += 0.005
            f5 = 0
        elif est_arcs >= 4:
            f3 += 0.03
            f4 += 0.015
            f5 += 0.014
        else:
            print("Something went wrong with the estimation of arcs! Check chirality_votes_maj in csv. Proceeding")
        # ---- Fourier modes. May need to add more at some point (?)
        feedme_list.append("F1) " + str(f1) + " 45  1  1") # Need to experiment with amplitude and phase angle for better understanding of this
        feedme_list.append("F3) " + str(f3) + " 25  1  1")
        feedme_list.append("#F4) " + str(f4) + " 4  1  1")
        feedme_list.append("#F5) " + str(f5) + " 6  1  1")  
        feedme_list.append("")
    
        # Sky -- Necessary?
        feedme_list.append("# Component number: 3")
        feedme_list.append("0) sky")
        feedme_list.append("1) 1000  1")
        feedme_list.append("2) 0  1")
        feedme_list.append("3) 0  1")
                       
        count += 1
    
        paths_to_feedme.append(write_to_feedme(galaxy, feedme_list))

