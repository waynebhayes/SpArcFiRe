#!/usr/bin/env python
# coding: utf-8

# **Author: Matthew Portman**
# 
# **Date (Github date will likely be more accurate): 6/1/22**

# # MINI README
# 
# Please place this script into the directory which contains the folders for your input, temporary, and output files for SpArcFiRe typically denoted: sparcfire-in, sparcfire-tmp, and sparcfire-out. GALFIT will also be run from that same directory if you choose to use the control script included in the repo. I recommend running this alone *once* before running the big script to ensure everything is in its right place (https://www.youtube.com/watch?v=sKZN115n6MI). After you confirm that this works without any errors the first time around, feel free to run the control script from then on. 
# 
# Running from the overarching directory is a temporary measure which will be remedied upon completion of the entire control script and full integration with SpArcFiRe. 
# 
# TO RUN: `python sparc_to_galfit_feedme_gen.py`
# 
# To run the control script: `bash control_script.sh`

import numpy as np
from math import log
import glob
import csv
import subprocess
import random
import pandas as pd
import os
from sys import argv


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

# Grabbing filepath from command line
def command_line(run_as_script = True):
    
    in_dir_out = "sparcfire-in"
    tmp_dir_out = "sparcfire-tmp"
    out_dir_out = "sparcfire-out"
        
    if not run_as_script:
        return in_dir_out, tmp_dir_out, out_dir_out
    
    if len(argv) != 4: # including name of python script
        print(f"No path given. Defaulting to using {in_dir_out}\n{tmp_dir_out}\n{out_dir_out}")
        
    else:
        in_dir_out = argv[1]
        tmp_dir_out = argv[2]
        out_dir_out = argv[3]
        
    return in_dir_out, tmp_dir_out, out_dir_out


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
        filenames_out = [s.replace("in", "out") for s in filenames_out]
        
    return filenames_read, galaxy_names_out, filenames_out


def path_join(path='.', name='', file_ext=''):
    
    file_path = os.path.join(path, name + file_ext)

    # deprecated
    # "./" + path + '/' + name + file_ext
    # print(file_path)
    # print(glob.glob(file_path))
    # file_name = glob.glob(file_path)[0]
    
    return file_path

def scale_var(x, scale = 1):
    # Scales
    return float(x)*scale

def galaxy_information(galaxy_name, galaxy_path):
   
    spin_parity = random.choice(['','-'])
    try:
        #tsv_filename = glob.glob("./" + galaxy_path + '/' + galaxy_name + '.tsv')[0]
        #csv_filename = glob_name(galaxy_path, galaxy_name, '.csv') # Changing to csv to be consistent
        csv_filename = path_join(galaxy_path, galaxy_name, '.csv') # Changing to csv to be consistent
        csv_file = open(csv_filename, 'r')
        
    except:
        print("Can't open to read: ", csv_filename)
        print("Check Sparcfire output or directories. Proceeding with default values.")

        bulge_rad_out = '2'
        bulge_axis_ratio_out = '0.5'
        bulge_rot_angle_out = '1'
        
        crop_rad_out = '30' # New!
        
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
        reader = csv.DictReader(csv_file, skipinitialspace = True)
#        csv_in = list(reader)
        for row in reader:
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
        if chirality == 'Zwise':
            spin_parity = '-'
        elif chirality == 'Swise':
            spin_parity = ''
    elif chirality_2 == chirality_3:
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


def arc_information(galaxy_name, galaxy_path, num_arms = 2):

    try:
        #arcs_filename = glob.glob(galaxy_path + '/' + '*_arcs.csv')
        arcs_filename = path_join(galaxy_path, galaxy_name, '_arcs.csv')
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
        
        # Nested dictionaries
        theta_sum = 0
        inner_rad = 0
        outer_rad = 0

        i = 0
        count = 0
        winding_dir = float(arcs_in[0]['pitch_angle']) > 0

        while count < num_arms:
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


def csv_sdss_info(galaxy_names): # to grab petromag and also psf parameter things
    
    gname_info = {}
    
    try:
        #star_dl_filename = glob_name('star_dl','star_dl','.csv')
        star_dl_filename = path_join('star_dl','star_dl','.csv')
        star_dl_file = open(star_dl_filename, 'r')
        
    except:
        print("Can't open to read star_dl.csv with star information (for PSF)")
        print("Check Sparcfire output or directories. Proceeding with default values.")

        # Defaults go here
        run = '0'
        rerun = '0'
        cam = '0'
        field = '0'
        row = '0'
        col = '0'
        pmag = '16'
        
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
                run_out = '0'
                rerun_out = '0'
                camcol_out = '0'
                field_out = '0'
                rowc_out = '0'
                colc_out = '0'
                petromag_out = '16'

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


def write_to_feedme(path, list_in):
    
    file_path = os.path.join(".", path, "autogen_feedme_galfit.in")
#    file_path = "/".join([".", path, ""]) + 'autogen_feedme_galfit.in'
    with open(file_path, 'w') as g:
        #for value in list_in:
        _ = [g.write(f"{value}\n") for value in list_in]
        
    return file_path


def rebuild_template_dict(example_file_path):
    galfit_dict = {}
    sky = "#  Sky background at center of fitting region [ADUs]"
    add = ""
    with open(example_file_path, 'r') as e:
        for line in e.readlines():
            try:
                l_num = i.strip().split(")")
                if len(l_num) > 1:
                    comment = i[i.index("#"):]
                    if comment == sky:
                        add = "sky"
                    galfit_dict[add + l_num[0]] = comment
            except:
                pass
    
    # Taken from the template itself and copied/pasted 
    galfit_dict["fill"] = len("sersic                 ")-1 # -1 because I leave a space in format
    return galfit_dict


def quick_build_template():
    galfit_template_dict = eval(
        """{"# IMAGE and GALFIT CONTROL PARAMETERS" : "# IMAGE and GALFIT CONTROL PARAMETERS"
            "A" : "# Input data image (FITS file)"
            "B" : "# Output data image block"
            "C" : "# Sigma image name (made from data if blank or 'none') "
            "D" : "# Input PSF image and (optional) diffusion kernel"
            "E" : "# PSF fine sampling factor relative to data "
            "F" : "# Bad pixel mask (FITS image or ASCII coord list)"
            "G" : "# File with parameter constraints (ASCII file) "
            "H" : "# Image region to fit (xmin xmax ymin ymax)"
            "I" : "# Size of the convolution box (x y)"
            "J" : "# Magnitude photometric zeropoint "
            "K" : "# Plate scale (dx dy)   [arcsec per pixel]"
            "O" : "# Display type (regular, curses, both)"
            "P" : "# Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps"
            "# INITIAL FITTING PARAMETERS" : "# INITIAL FITTING PARAMETERS"
            "# Component number: 1" : "# Component number: 1"
            "0" : "#  Component type"
            "1" : "#  Position x, y"
            "3" : "#  Integrated magnitude "
            "4" : "#  R_e (effective radius)   [pix]"
            "5" : "#  Sersic index n (de Vaucouleurs n=4) "
            "6" : "#     ----- "
            "7" : "#     ----- "
            "8" : "#     ----- "
            "9" : "#  Axis ratio (b/a)  "
            "10" : "#  Position angle (PA) [deg: Up=0, Left=90]"
            "Z" : "#  Skip this model in output image?  (yes=1, no=0)"
            "R0" : "#  PA rotation func. (power, log, none)"
            "R1" : "#  Bar radius [pixels]"
            "R2" : "#  Spiral outer (i.e. asymptotic) radius [pixels]"
            "R3" : "#  Cumul. rotation out to outer radius [degrees]"
            "R4" : "#  Asymptotic spiral powerlaw "
            "R9" : "#  Inclination to L.o.S. [degrees]"
            "R10" : "#  Sky position angle"
            "F1" : "#  Azim. Fourier mode 1, amplitude, & phase angle"
            "F3" : "#  Azim. Fourier mode 3, amplitude, & phase angle"
            "F4" : "#  Azim. Fourier mode 4, amplitude, & phase angle"
            "F5" : "#  Azim. Fourier mode 5, amplitude, & phase angle"
            "sky1" : "#  Sky background at center of fitting region [ADUs]"
            "sky2" : "#  dsky/dx (sky gradient in x)     [ADUs/pix]"
            "sky3" : "#  dsky/dy (sky gradient in y)     [ADUs/pix]"
            "skyZ" : "#  Skip this model in output image?  (yes=1, no=0)"
            "fill" : 22}""".replace("\n",","))
    return galfit_template_dict


if __name__ == "__main__":

    count = 0
    paths_to_feedme = []
    
    in_dir, tmp_dir, out_dir = command_line() 
    filenames_in, galaxy_names, filenames_out = get_galaxy_names_list(in_dir)
    
    psf_info = csv_sdss_info(galaxy_names)
    
    for galaxy in filenames_out:
    
        gname = galaxy_names[count]
        
        # From old implementation - 1/19/21
        # ************
        # x1crop, x2crop, y1crop, y2crop = autocrop_grab(gname, galaxy)
        # scale = (float(x2crop)-float(x1crop))/256 - from old implementation
        # ************
        
        bulge_rad, bulge_axis_ratio, pos_angle_bulge,
        crop_rad, center_pos_x, center_pos_y,
        disk_maj_axs_len, pos_angle_disk, pos_angle_power,
        axis_ratio, max_arc, spin_dir,
        est_arcs, inclination, bar_candidate,
        alpha = galaxy_information(gname, galaxy)
        
        center_pos_x = float(center_pos_x)
        center_pos_y = float(center_pos_y)
        crop_rad = float(crop_rad)
        
        x1crop = center_pos_x - crop_rad
        x2crop = center_pos_x + crop_rad
        # It's a square so we can just use x1, x2
    
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
    
    
        # Initializing Feedme
        feedme_list = []
        
        # Initialize template dict
        gt = quick_build_template() # galfit template
        #gt = rebuild_template_dict("./m51.feedme")
    
        #To reconstruct the z PSF (i.e., the 5th HDU) at the position (row, col) = (500, 600) from run 1336, column 2, field 51 you’d say:
        #read_PSF psField-001336-2-0051.fit 5 500.0 600.0 foo.fit
        run, rerun, camcol, field, psf_row, psf_col, petromag = psf_info[gname]

        feedme_list.append(f"#{run}{camcol}{field}; HDU: z{psf_row}{psf_col}")
        feedme_list.append("")
        # Image and Galfit Control Param
        feedme_list.append(f"A) ./{filenames_in[count]}")
        feedme_list.append(f"B) {tmp_dir}/galfits/{gname}_out.fits")
        feedme_list.append(f"C) none")
        feedme_list.append(f"D) {tmp_dir}/psf_files/{gname}_psf.fits")
        feedme_list.append(f"E) 1")
        feedme_list.append(f"F) {tmp_dir}/galfit_masks/{gname}_star-rm.fits")
        feedme_list.append(f"G) none")  #./constraints.txt"
        feedme_list.append(f"H) {x1crop:.1f} {x2crop:.1f} {x1crop:.1f} {x2crop:.1f}") # Square!
        feedme_list.append(f"I) 50 50") # psf FWHM ~= 1, Chien recommends 40-80 times this value
        feedme_list.append(f"J) 24.8") # SDSS
        feedme_list.append(f"K) 0.396  0.396") # SDSS
        feedme_list.append(f"O) regular")
        feedme_list.append(f"P) 0")
        feedme_list.append("")
            
        # Sersic 1
        # Fixing as much as I can here... it's not exactly a priority.
        feedme_list.append(f"# Component number: 1")
        feedme_list.append(f"0) sersic")
        feedme_list.append(f"1) {center_pos_x:.1f} {center_pos_y:.1f} 0 0")
        feedme_list.append(f"3) {float(petromag) - 1:.2f} 1") # Initial guess goes here
        feedme_list.append(f"4) {bulge_rad:.2f} 1") 
        feedme_list.append(f"5) 4  1") # According to other paper GALFIT usually doesn't have a problem with the index
        feedme_list.append("6) 0  0")    
        feedme_list.append("7) 0  0")    
        feedme_list.append("8) 0  0")    
        # According to other papers, bulge (esp. in spiral galaxies) averages to about 2 so this is a good starting place
        # see https://ned.ipac.caltech.edu/level5/Sept11/Buta/Buta9.html
        feedme_list.append(f"9) {axis_ratio:.2f} 1")  
        feedme_list.append(f"10) {pos_angle_bulge:.2f} 1") 
        feedme_list.append("")
    
        # Sersic 2
        feedme_list.append("# Component number: 2")
        feedme_list.append(f"0) sersic")
        feedme_list.append(f"1) {center_pos_x:.1f} {center_pos_y:.1f} 0 0")
        feedme_list.append(f"3) {float(petromag) - 3:.2f} 1") 
        feedme_list.append(f"4) {disk_maj_axs_len:.2f} 1") # Use this for effective radius? Also 0 this one out? Will have to see how well it works in practice
        feedme_list.append(f"5) 4  1") # Classical disk follows Sersic n = 1 so good place to start (per Readme Exponential profile)
                                      # According to comparison tests, this usually ends up much higher probably due to the spiral.
        feedme_list.append("6) 0  0")    
        feedme_list.append("7) 0  0")    
        feedme_list.append("8) 0  0")    
        feedme_list.append(f"9) 0.6 1")  # Fixing this to 0.5 to give the arms the best chance to form
        #(f"9) {axis_ratio - 0.3} 1 {gt['9']}")
        feedme_list.append(f"10) {pos_angle_disk:.2f} 1") #90  1") # fixing this to 'normal' 0 so that we can JUST rotate power function
        feedme_list.append("")
    
        # Power
        feedme_list.append("R0) power")
        feedme_list.append(f"R1) {in_rad:.2f} 0") # Chosen based on where *detection* of arms usually start
        feedme_list.append(f"R2) {out_rad:.2f} 0")
        feedme_list.append(f"R3) {spin_dir}{cumul_rot:.2f} 1") # See calc above
        feedme_list.append(f"R4) {alpha:.2f} 1") # Another good thing to automate via Sparcfire 
        feedme_list.append(f"R9) {inclination:.2f} 1") # see if can't 0 this one out... 
        feedme_list.append(f"R10) 40  1")#+ pos_angle_power + " 1") # Always more to discover, looks like all the images are mirrored across the y axis.

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
        feedme_list.append(f"F1) {f1:.2f} 45  1  1") # Need to experiment with amplitude and phase angle for better understanding of this
        feedme_list.append(f"F3) {f3:.3f} 25  1  1")
        feedme_list.append(f"#F4) {f4:.3f} 4  1  1")
        feedme_list.append(f"#F5) {f5:.3f} 6  1  1")  
        feedme_list.append("")
    
        # Sky -- Necessary?
        feedme_list.append(f"# Component number: 3")
        feedme_list.append(f"0) sky")
        feedme_list.append(f"1) 1000  1")
        feedme_list.append(f"2) 0  1")
        feedme_list.append(f"3) 0  1")
                       
        count += 1
        
        formatted_feedme = []
        extra = ""
        for i in feedme_list:
            if i and not i.startswith("#"):
                str_split = i.split(")")
                component = extra + str_split[0]
                formatted_feedme.append(f"{i:<{gt['fill']}} {gt[component]}")
                
                # Sneakily do this at the end since 0) sky is just component name
                if "sky" in str_split[1]:
                    extra = "sky"
    
        paths_to_feedme.append(write_to_feedme(galaxy, formatted_feedme)) # do I need paths_to_feedme? I used to use it for something...

