#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Tabs I have open
# https://xgboost.readthedocs.io/en/stable/python/examples/multioutput_regression.html#sphx-glr-python-examples-multioutput-regression-py
# https://xgboost.readthedocs.io/en/stable/python/python_api.html#xgboost.XGBRegressor
# https://xgboost.readthedocs.io/en/stable/parameter.html
# https://xgboost.readthedocs.io/en/stable/prediction.html
# https://scikit-learn.org/stable/auto_examples/ensemble/plot_random_forest_regression_multioutput.html#sphx-glr-auto-examples-ensemble-plot-random-forest-regression-multioutput-py
# https://machinelearningmastery.com/xgboost-for-regression/
# https://machinelearningmastery.com/regression-metrics-for-machine-learning/


# In[2]:


#!pip install --upgrade pandas
#!pip install graphviz
#!pip install xgboost
#!pip install hyperopt
#!pip install kaleido


# In[2]:


from xgboost import XGBClassifier
import xgboost as xgb
from sklearn.model_selection import train_test_split, cross_val_score, RepeatedKFold
from sklearn.metrics import accuracy_score, SCORERS
from sklearn.metrics import mean_squared_error as MSE
#from sklearn.model_selection import GridSearchCV
from hyperopt import fmin, Trials, hp, tpe, STATUS_OK


# In[2]:


import numpy as np
from scipy.stats import kstest, anderson

from math import ceil
import itertools

from PIL import Image
import matplotlib.pyplot as plt

import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots

#import graphviz
import ssl
from glob import glob

#from annoy import AnnoyIndex
import random
import pandas as pd
from scipy import spatial

from shutil import copy2
import json
from copy import deepcopy

import pickle
import argparse

from joblib import Parallel, delayed


# In[5]:


# For debugging purposes
from IPython import get_ipython
def in_notebook():
    ip = get_ipython()
    
    if ip:
        return True
    else:
        return False


# In[6]:


import sys
import os
from os.path import join as pj

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
from Classes.Containers import *
from Classes.FitsHandlers import *
from Functions.helper_functions import *


# #!jupyter nbconvert --to script testing_xgboost.ipynb

# In[9]:


if __name__ == "__main__":
    
    # Force >python 3.10 for various compatabilities
    out_str = "\t Python3.10 or greater required! Exitting without generating feedmes..."
    assert sys.version_info >= (3, 10), out_str
    
    cwd = absp(os.getcwd()) # Doesn't work *in* notebook
    old_cwd = absp(cwd) # Strings are immutable
    
    username = os.environ["USER"]
    
    USAGE = f"""USAGE:

    python3 ./{sys.argv[0]} [OPTION] [[RUN-DIRECTORY] IN-DIRECTORY TMP-DIRECTORY OUT-DIRECTORY]
    
    OPTIONS => [-v | --verbose]

    This script is used to train XGBoost to feed better input to GALFIT via SpArcFiRe. 
    By default, it runs from the RUN (or current) directory and uses the
    '-in' '-tmp' and '-out' directories as specified or otherwise defaults to 
    'sparcfire-in', 'sparcfire-tmp', 'sparcfire-out'. 

    Please do not specify symlinks for the above, they discomfort the programmer.
    """
    
    parser = argparse.ArgumentParser(description = USAGE)
    
    parser.add_argument('-v', '--verbose',
                        dest     = 'verbose', 
                        action   = 'store_const',
                        const    = True,
                        default  = False,
                        help     = 'Verbose output for all bash commands in control script.'
                       )
    
    parser.add_argument(dest     = 'out_path',
                        action   = 'store',
                        type     = str,
                        help     = "[OUT-DIRECTORY] from SpArcFiRe. \
                                    SpArcFiRe directories should follow -in, -tmp, out or this probably won't work."
                       )
    
    if not in_notebook():
        args              = parser.parse_args() # Using vars(args) will call produce the args as a dict
        
        verbose           = args.verbose
        capture_output    = not args.verbose
        
        sparc_out_dir = args.out_path
            
    else:
        verbose = False
        capture_output = True
        
        cwd = cwd.replace("ics-home", username)
        sparc_out_dir = pj(_HOME_DIR, "run2_1000_galfit", "sparcfire-out") #pj(cwd, "sparcfire-out")
        
        sys.path.append(pj(_HOME_DIR, ".local", "bin"))
        
    # Making these absolute paths
    cwd     = absp(cwd)
    #in_dir  = absp(in_dir)
    #tmp_dir = absp(tmp_dir)
    sparc_out_dir = absp(sparc_out_dir)
    
    # Changing to specified working dir
    os.chdir(cwd)


# In[10]:


if __name__ == "__main__":
    png_tiled_dir  = pj(cwd, "labeled_png_out")
    good_labeled_dir  = pj(png_tiled_dir, "good")

    inputs_dir = pj(cwd, "galfit_inputs")
    good_inputs_dir = pj(inputs_dir, "good")
    not_good_inputs_dir = pj(inputs_dir, "not_good")

    outputs_dir = pj(cwd, "galfit_outputs")
    good_outputs_dir = pj(outputs_dir, "good")
    not_good_outputs_dir = pj(outputs_dir, "not_good")


# In[8]:


# Used to get organized
def generate_labels_from_imgs(all_galaxies_path, good_folder_path):
    
    # good galaxies are named ###_combined.png
    # inputs are named ###.in
    
    good_galaxy_names = list(map(lambda i: i.split("_")[0], os.listdir(good_folder_path)))
    all_galaxy_names = [os.path.basename(i) for i in os.listdir(all_galaxies_path) if os.path.basename(i).startswith("123")]
    
    label_dict = {i:(1 if i in good_galaxy_names else 0) for i in all_galaxy_names}
    
    return label_dict


# In[9]:


def organize_files(label_dict, sparcfire_out_dir, **kwargs):
    
    labels_top_dir  = kwargs.get("labels_top_dir", os.getcwd())
    good_name       = kwargs.get("good_name", "good")
    not_good_name   = kwargs.get("not_good_name", "not_good")
    label_dump_path = kwargs.get("label_dump_path", pj(labels_top_dir, "labeled_galaxies.json"))
    
    temp_dict = deepcopy(label_dict)
    for gname, label in temp_dict.items():
        if label:
            label_name = good_name
        else:
            label_name = not_good_name
        
        # We will extract feedme from output fits file
        # There should always be a generated output (because labeled) but just in case...
        # Copy that first so the except catches before trying to copy the input that way we don't have a mismatch
        try:
            #copy2(pj(sparcfire_out_dir, gname, f"{gname}_galfit_out.fits"), pj(labels_top_dir, "galfit_outputs", label_name))
            #copy2(pj(sparcfire_out_dir, gname, f"{gname}.in"), pj(labels_top_dir, "galfit_inputs", label_name))
            copy2(pj(sparcfire_out_dir, gname, f"{gname}_out.fits"), pj(labels_top_dir, "galfit_outputs", label_name, f"{gname}_galfit_out.fits"))
            copy2(pj(sparcfire_out_dir, gname, "autogen_feedme_galfit.in"), pj(labels_top_dir, "galfit_inputs", label_name, f"{gname}.in"))
        except FileNotFoundError:
            print(f"No output found for {gname}, continuing to copy...")
            label_dict.pop(gname)
            
    with open(label_dump_path, "w") as lg:
        json.dump(label_dict, lg)
    
    return label_dict


# In[10]:


if __name__ == "__main__":
    
    label_json = pj(cwd, 'labeled_galaxies.json')
    
    if exists(label_json):
        label_dict = json.load(open(label_json, 'r'))
    else:
        label_dict = generate_labels_from_imgs(sparc_out_dir, good_labeled_dir)
        label_dict = organize_files(label_dict, pj(_HOME_DIR, "run2_1000_galfit", "sparcfire-out"), label_dump_path = label_json)


# In[11]:


# class LabeledModel(OutputFits):
#     def __init__(self, 
#                  filepath = "",
#                  label = 0,
#                  **kwargs
#                 ):
        
#         OutputFits.__init__(self, filepath)
        
#         self.label = label        
#         self.df    = self.feedme.to_pandas(self)


# In[12]:


# # Non-parallelized
# def build_df(label_dict, label_dirs, file_suffix: str):
    
#     out_df = pd.DataFrame()
    
#     for gname, label in label_dict.items():

#         # 0 bad, 1 good
#         if label == 0:
#             input_dir = label_dirs[0]
#             #str_label = "not_good"
#         else:
#             input_dir = label_dirs[1]
#             #str_label = "good"

#         #gal_dict, param_names = #galfit_param_grab(pj(input_dir, gname + file_suffix))
# #         if not param_names: continue
        
# #         for i, param in enumerate(param_names):
# #             if param in param_names[:i]:
# #                 param += ""

#         #gname_df = flatten_to_pandas(gal_dict, param_names, gname)
#         path_to_output = pj(input_dir, f"{gname}{file_suffix}")
        
#         if file_suffix.endswith(".fits"):
#             feedme = OutputFits(path_to_output).feedme
#             gname_df = feedme.to_pandas()
            
#             region_to_fit = feedme.header.region_to_fit
#             gname_df["crop_rad"] = region_to_fit[1] - region_to_fit[0]
            
#         elif file_suffix.endswith(".in"):
#             feedme = FeedmeContainer(path_to_feedme = path_to_output)
#             feedme.from_file(path_to_output)
#             gname_df = feedme.to_pandas()
            
#             region_to_fit = feedme.header.region_to_fit
#             gname_df["crop_rad"] = region_to_fit[1] - region_to_fit[0]
            
#         gname_df["label"] = label #str_label
#         # TODO: Put this rename function in a to_pandas() function of OutputFits
#         # Make note in documentation: name is incoherent in everything but an outputfits
#         # since a feedme does not necessarily imply a galaxy
#         # That being said, make an optional gname parameter in feedme
#         gname_df.rename(index = {0:gname}, inplace = True)
#         out_df = pd.concat([out_df, gname_df])

#     # For the input galaxies, we have a lot of held values, these are uneccessary
#     # https://stackoverflow.com/a/39658662
#     nunique = out_df.nunique()
#     cols_to_drop = nunique[nunique == 1].index
#     out_df.drop(columns = cols_to_drop, inplace = True)
    
#     return out_df


# In[13]:


def build_df(gname, label, count, label_dirs, file_suffix: str):
    
    if not count % 100:
        print(gname, count)
        
    # 0 bad, 1 good
    if label == 0:
        input_dir = label_dirs[0]
        #str_label = "not_good"
    else:
        input_dir = label_dirs[1]
        #str_label = "good"

    #gal_dict, param_names = #galfit_param_grab(pj(input_dir, gname + file_suffix))
#         if not param_names: continue

#         for i, param in enumerate(param_names):
#             if param in param_names[:i]:
#                 param += ""

    #gname_df = flatten_to_pandas(gal_dict, param_names, gname)
    path_to_output = pj(input_dir, f"{gname}{file_suffix}")

    if file_suffix.endswith(".fits"):
        feedme = OutputFits(path_to_output).feedme
        gname_df = feedme.to_pandas()

        region_to_fit = feedme.header.region_to_fit
        gname_df["crop_rad"] = region_to_fit[1] - region_to_fit[0]

    elif file_suffix.endswith(".in"):
        feedme = FeedmeContainer(path_to_feedme = path_to_output)
        feedme.from_file(path_to_output)
        gname_df = feedme.to_pandas()

        region_to_fit = feedme.header.region_to_fit
        gname_df["crop_rad"] = region_to_fit[1] - region_to_fit[0]

    gname_df["label"] = label #str_label
    # TODO: Put this rename function in a to_pandas() function of OutputFits
    # Make note in documentation: name is incoherent in everything but an outputfits
    # since a feedme does not necessarily imply a galaxy
    # That being said, make an optional gname parameter in feedme
    gname_df.rename(index = {0:gname}, inplace = True)
    
    return gname_df


# In[14]:


def convert_angles(in_df, **kwargs):
    angles = ["position_angle_sersic_1",
              "position_angle_sersic_2",
              "cumul_rot_power_2",
              "inclination_power_2",
              "sky_position_angle_power_2"]
    
    angles = kwargs.get("angles", angles)
    
    out_df = in_df.copy()
    
    for col_name in angles:
        if col_name in out_df.columns:
            # Take advantage of symmetry across axes
            # Will have to be careful when outputting new template to retain correct direction
            # will likely pull this from original data
            out_df.loc[out_df[col_name] < 0, col_name] += 180
            out_df[col_name] *= np.pi/180
            
    return out_df


# In[15]:


def prepare_df(picklepath, label_dict, label_dirs, file_suffix):
    if exists(picklepath):
        df = pickle.load(open(picklepath, "rb"))
    else:
        if file_suffix.endswith(".fits"):
            df = Parallel(n_jobs = -2, timeout = 30)(
            delayed(build_df)(
                              gname,
                              label,
                              count,
                              label_dirs=label_dirs, 
                              file_suffix="_galfit_out.fits"
                             ) 
            for count, (gname, label) in enumerate(label_dict.items())
                                                    )
        elif file_suffix.endswith(".in"):
            df = [build_df(gname, label, count, label_dirs, file_suffix) 
                  for count, (gname, label) in enumerate(label_dict.items())]   
            
        else:
            raise Exception(f"file suffix {file_suffix} not valid.")
        
        df = pd.concat(df)
        
        # For the input galaxies, we have a lot of held values, these are uneccessary
        # https://stackoverflow.com/a/39658662
        nunique = df.nunique()
        cols_to_drop = nunique[nunique == 1].index
        df.drop(columns = cols_to_drop, inplace = True)
        
        df.to_pickle(picklepath)

    # Save *before* conversion for data posterity
    df = convert_angles(df)
    
    return df


# In[16]:


if __name__ == "__main__":
    print("Grabbing GALFIT input data...")
    galfit_in_df = prepare_df("galfit_in_df.pkl", label_dict, [not_good_inputs_dir, good_inputs_dir], ".in")
    print("Done!\n")
    
    print("Grabbing GALFIT output data, this may take awhile (if not already saved)...")
    galfit_out_df = prepare_df("galfit_out_df.pkl", label_dict, [not_good_outputs_dir, good_outputs_dir], "_galfit_out.fits")
    print("Done!\n")


# In[17]:


# if __name__ == "__main__":
#     in_pkl = pj(cwd, "galfit_in_df.pkl")
#     if exists(in_pkl):
#         galfit_in_df = pickle.load(open(in_pkl, "rb"))
#     else:   
#         galfit_in_df = [build_df(gname, label, count, [not_good_inputs_dir, good_inputs_dir], ".in") 
#                         for count, (gname, label) in enumerate(label_dict.items())]   
        
#         galfit_in_df = pd.concat(galfit_in_df)
        
#         # For the input galaxies, we have a lot of held values, these are uneccessary
#         # https://stackoverflow.com/a/39658662
#         nunique = galfit_in_df.nunique()
#         cols_to_drop = nunique[nunique == 1].index
#         galfit_in_df.drop(columns = cols_to_drop, inplace = True)
        
#         galfit_in_df.to_pickle(in_pkl)
        
#     galfit_in_df = convert_angles(galfit_in_df)


# In[18]:


# # TODO: Parallelize this since FITS take so much longer
# # I could alternatively output all the .in files from those FITS files separate and then build the DFs... a thought
# if __name__ == "__main__":
#     out_pkl = pj(cwd, "galfit_out_df.pkl")
#     if exists(out_pkl):
#         galfit_out_df = pickle.load(open(out_pkl, "rb"))
#     else:
#         galfit_out_df = Parallel(n_jobs = -2, timeout = 30)(
#             delayed(parallel_build_df)(
#                                        gname,
#                                        label,
#                                        count,
#                                        label_dirs=[not_good_outputs_dir, good_outputs_dir], 
#                                        file_suffix="_galfit_out.fits"
#                                       ) 
#             for count, (gname, label) in enumerate(label_dict.items())
#                                                                            )
#         galfit_out_df = pd.concat(galfit_out_df)
        
#         # For the input galaxies, we have a lot of held values, these are uneccessary
#         # https://stackoverflow.com/a/39658662
#         nunique = galfit_out_df.nunique()
#         cols_to_drop = nunique[nunique == 1].index
#         galfit_out_df.drop(columns = cols_to_drop, inplace = True)
        
#         galfit_out_df.to_pickle(out_pkl)
        
#     galfit_out_df = convert_angles(galfit_out_df)
#     #galfit_out_df['crop_rad'] = galfit_in_df['crop_rad']


# In[19]:


# TODO: Check deprecated then delete
# positions = ["position_x_sersic_1",
#              "position_y_sersic_1",
#              "position_x_sersic_2",
#              "position_y_sersic_2"]

# fourier = ["F1_amplitude_fourier_2",
#            "F1_phase_angle_fourier_2",
#            "F3_amplitude_fourier_2",
#            "F3_phase_angle_fourier_2"
#           ]

# file_prefixes = ["reduced_galfit_in", "reduced_galfit_out"]
# ignore_galfit_in  = ["inner_rad_power_2"]
# ignore_galfit_in.extend(positions)
# ignore_galfit_in.extend(fourier)

# ignore_galfit_out = ["crop_rad",
#                      "effective_radius_sersic_1", # Bulge radius goes crazzyyyyyy
#                      "inclination_power_2", # Trust sparc/galfit on this one
#                      "inner_rad_power_2",
#                      "outer_rad_power_2",
#                      "sky_background_sky_3", # Trust galfit on this one
#                      "dsky_dx_sky_3", # if we include them, does it get better? for future testing
#                      "dsky_dy_sky_3"
#                     ]
# ignore_galfit_out.extend(positions)
# ignore_galfit_out.extend(fourier)




# in_filter = []
# # filter_1 = "sersic_index_sersic_1"
# # filter_2 = "sersic_index_sersic_1"
# # filter_3 = "magnitude_sersic_1"
# # filter_4 = "effective_radius_sersic_1"
# # filter_5 = "effective_radius_sersic_2"




# In[20]:


def export_filter():
    positions = ["position_x_sersic_1",
                 "position_y_sersic_1",
                 "position_x_sersic_2",
                 "position_y_sersic_2"]

    fourier = ["F1_amplitude_fourier_2",
               "F1_phase_angle_fourier_2",
               "F3_amplitude_fourier_2",
               "F3_phase_angle_fourier_2"
              ]
    
    # Not all df may have columns I'd otherwise check
    # For instance, I hold the input constant for some variables (sersic index) so the input
    # df won't have those columns. So be careful here
    ignore_galfit_in  = ["inner_rad_power_2"]
    ignore_galfit_in.extend(positions)
    ignore_galfit_in.extend(fourier)
    
    ignore_galfit_out = ["crop_rad",
                         "effective_radius_sersic_1", # Bulge radius goes crazzyyyyyy
                         "inclination_power_2", # Trust sparc/galfit on this one
                         "inner_rad_power_2",
                         "outer_rad_power_2",
                         "F1_amplitude_fourier_2",
                         "F3_amplitude_fourier_2",
                         "F1_phase_angle_fourier_2",
                         "F3_phase_angle_fourier_2",
                         "sky_background_sky_3", # Trust galfit on this one
                         "dsky_dx_sky_3", # if we include them, does it get better? for future testing
                         "dsky_dy_sky_3"
                        ]
    ignore_galfit_out.extend(positions)
    ignore_galfit_out.extend(fourier)

    in_filter = []
    
    out_filter = [#"`Asymptotic spiral powerlaw disk` <= 1",
              f"`sersic_index_sersic_1` > 14",
              f"`sersic_index_sersic_1` < 0.05",
              f"`magnitude_sersic_1` > 26",
              f"`effective_radius_sersic_1` > `crop_rad`",
              f"`effective_radius_sersic_2` > `crop_rad`"
              #"`R_e (effective radius)   (pix) bulge` < `Crop Rad`"
              ]
    
    return ignore_galfit_in, ignore_galfit_out, in_filter, out_filter


# In[21]:


if __name__ == "__main__":
    # in_train/test should always have fewer exclusions and those exclusions must also be in the out_train/test
    # so that the algo isn't predicting things sight unseen
    sklearn = True

    ignore_galfit_in, ignore_galfit_out, in_filter, out_filter = export_filter()
    
    file_prefixes = ["reduced_galfit_in", "reduced_galfit_out"]
    
    ignore_galfit = (ignore_galfit_in, ignore_galfit_out)
    col_ignore = {fp:ig for fp, ig in zip(file_prefixes, ignore_galfit)}
    
    filter_strings = (in_filter, out_filter)
    gal_filter = {fp:fs for fp, fs in zip(file_prefixes, filter_strings)}


# In[22]:


# Per https://xgboost.readthedocs.io/en/stable/python/python_intro.html
# For training *regressor* on good data
def split_save_df_reg(*args, file_prefixes = [], col_ignore = {}, gal_filter = {}, sklearn = False):
    # col_ignore must be list of list of columns since we drop different things for different df
    # gal_filter must be given as: `column name` cond value, i.e. "`Asymptotic spiral powerlaw disk` <= 1"
    # Note, the conditions are what we *don't want*
    
    assert len(file_prefixes) == len(args), "File prefixes must be same length as # of dataframes being passed in! Try again"
    #assert col_ignore == len(args), "col_ignore must be same length as # of dataframes being passed in! Try again"
    
    return_dict = {}
    
    # Specifically looking for good galaxies since we're not classifying
    #gal_filter.append("label != 1") 
    #gal_filter = " or ".join(gal_filter)
    
    galaxies_to_drop = [in_df.query(
                                    " or ".join(gf + ["label != 1"])
                                    ).index 
                        for in_df, gf in zip(args, gal_filter.values()) if gf]
    
    # Unpacking to keep the list comp ;)
    galaxies_to_drop = list(itertools.chain.from_iterable(galaxies_to_drop))
    print(f"Keeping {len(args[0]) - len(galaxies_to_drop)} galaxies (out of {len(args[0].query('label == 1'))})")
    
    for in_df, file_prefix in zip(args, file_prefixes):
        # Filter and exclude
        
        exclude = ['label'] + col_ignore.get(file_prefix, [])
        
        exclude = list(set(exclude).intersection(set(in_df.columns)))
        in_df_good = in_df.drop(index = galaxies_to_drop, columns = exclude)
        
        label = pd.DataFrame(np.ones(len(in_df_good)), columns = ["label"], dtype = "int")
            
        # Pass random_state = 0 to guarantee input and output are lined up
        X_train, X_test, y_train, y_test = train_test_split(in_df_good, label, test_size=.3, random_state = 0)
    
#         dtrain = xgb.DMatrix(X_train, label = y_train)
#         dtest  = xgb.DMatrix(X_test,  label = y_test)
#         ddata  = xgb.DMatrix(in_df_good, label = label)

#         print(f"Saving dmatrices to file: {file_prefix}.train/test/data")
#         dtrain.save_binary(f'{file_prefix}.train')
#         dtest.save_binary(f'{file_prefix}.test')
#         ddata.save_binary(f'{file_prefix}.data')
        
        in_df_good["label"] = list(label.label)
        if sklearn:
            return_dict[file_prefix] = X_train, X_test, y_train, y_test, in_df_good
        else:
            return_dict[file_prefix] = dtrain, dtest, ddata, in_df_good, label

    return return_dict


# In[23]:


#TODO: FILTER OUT PREVIOUS GALFIT RESULTS FOR NON-PHYSICAL DATA, I.E. SERSIC INDEX TOO HIGH OR TOO LOW IN GOOD DATA
# Convert angles to radians -- done
# weight spiral law(?), look up how weighting works -- doesn't work
# Make better plots (seaborn?)
# Use physics based constraints to filter results
# Determine better loss method? Hmmmmm https://stats.stackexchange.com/a/445454
# For +/- things (like angle) find some way to use absolute value while retaining distribution, see if that makes things better -- done


# In[24]:


if __name__ == "__main__":
    split_save_out = split_save_df_reg(
                                       galfit_in_df, galfit_out_df,
                                       file_prefixes = file_prefixes,
                                       col_ignore    = col_ignore,
                                       gal_filter    = gal_filter,
                                       sklearn       = sklearn
                                      )

    # New...df is a combo of train/test, no different otherwise
    in_train,  in_test,  _, _, new_in_df    = split_save_out[file_prefixes[0]]
    out_train, out_test, _, _, new_out_df   = split_save_out[file_prefixes[1]]

    # Could train on bad galaxies too... but only for classification purposes
    
    assert len(galfit_out_df) - len(galfit_out_df.query(" or ".join(out_filter + ["label != 1"]))) == len(new_out_df)


# In[25]:


def load_data(file_prefix):
    dtrain = xgb.DMatrix(f'{file_prefix}.train')
    dtest  = xgb.DMatrix(f'{file_prefix}.test')
    ddata  = xgb.DMatrix(f'{file_prefix}.data')
    
    return dtrain, dtest, ddata


# In[26]:


# if __name__ == "__main__":
#     if not sklearn:
#         dtrain_in, dtest_in, ddata_in = load_data(file_prefixes[0])
#         dtrain_out, dtest_out, ddata_out = load_data(file_prefixes[1])


# In[27]:


def loss_methods(test_df, pred_df, method = "kstest"):
    method = method.lower()
    
    if method == "kstest":
        loss = np.average([kstest(test_df[col].values, pred_df[col].values).statistic
                                      for col in pred_df])
        method_str = "KStest Stat"


    elif method == "rmse":
        # squared = False gives RMSE
        loss = MSE(test_df, pred_df, squared = False)
        method_str = "RMSE"
        
    # Average of per galaxy rmse
    elif method == "galaxy_rmse":
        comparison_df = np.square(test_df - pred_df)
        loss = np.mean(np.sqrt(comparison_df.mean(axis=1)))
        method_str = "Per Galaxy RMSE"
        
    print(f"Average {method_str}: {loss:.3f}")

    return loss


# In[28]:


# Objective function
def objective(space): #, data, label, test_size = 0.3):
    clf = xgb.XGBRegressor(
                    learning_rate = space['learning_rate'],
                    n_estimators = space['n_estimators'],
                    #subsample = space['subsample'],
                    max_depth = int(space['max_depth']))
                    #gamma = space['min_split_loss'])
                    #reg_alpha = int(space['reg_alpha']),
                    #min_child_weight=int(space['min_child_weight']),
                    #colsample_bytree=int(space['colsample_bytree']))
    
    clf.set_params(#eval_metric="auc",
                   tree_method = "hist",
                   early_stopping_rounds=10,
                   objective='reg:squarederror'
                   #objective='reg:pseudohubererror'
                   #'eval_metric' : 'error' # Binary classification error rate
                   )
    
    clf.fit(in_train, 
            out_train, 
            eval_set = [(in_train, out_train)],
            verbose = False)

    pred = clf.predict(in_test)
    out_pred_df = pd.DataFrame(pred, columns = out_test.columns, index = in_test.index)

    # use weights to lower contribution from values with wider spread since these aren't all normalized
    # alternatively normalize everything to [0,1]
    #rmse = np.sqrt(MSE(out_test, pred))#, multioutput = "raw_values"))
    #print(f"RMSE: {rmse:.3f}")
    
#     acc = []
#     for col in out_test:
#         acc.append(kstest(out_test[col].values, out_pred_df[col].values).statistic)
        
#     print(np.average(acc))
        
    loss = loss_methods(out_test, out_pred_df, method = "kstest")
    
    return {'loss': loss, 'status': STATUS_OK }
    #accuracy = accuracy_score(y_test, pred>0.5)
    #print (f"SCORE: {accuracy:.3f}")
    #return {'loss': -accuracy, 'status': STATUS_OK }


# In[32]:


# Optimizing hyperparameters using hyperopt???
# Bayesian Optimization
if __name__ == "__main__":
    space = {'learning_rate': hp.uniform('learning_rate', 0.01, 0.2),
             'max_depth': hp.quniform('max_depth', 2, 18, 1),
             #'subsample' : hp.uniform('subsample', 0.5, 1),
             #'min_child_weight' : hp.quniform('min_child_weight', 0, 5, 1),
             #'max_delta_step' : hp.quniform('max_delta_step', 0, 10, 1),
             #'min_split_loss': hp.uniform ('min_split_loss', 1,9),
             #'reg_alpha' : hp.quniform('reg_alpha', 40, 180,1),
             #'reg_lambda' : hp.uniform('reg_lambda', 0, 1),
             #'colsample_bytree' : hp.uniform('colsample_bytree', 0.5,1),
             #'min_child_weight' : hp.quniform('min_child_weight', 0, 10, 1),
             'n_estimators' : hp.uniformint('n_estimators', 5, 100),
             'seed': 0
            }
    
    trials = Trials()

    data = new_in_df.drop(columns = "label")
    label = new_in_df["label"]

    best_hyperparams = fmin(fn = objective,
                            space = space,
                            algo = tpe.suggest,
                            max_evals = 100,
                            trials = trials
                           )
    
    best_hyperparams_int = {k:(int(v) if int(float(v)) == v else v) for k,v in best_hyperparams.items()}       
    print("The best hyperparameters are : ","\n")
    print(best_hyperparams_int)


# In[ ]:


# https://davetang.org/muse/2012/04/17/comparing-different-distributions/
# https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test#Two-sample_Kolmogorov.E2.80.93Smirnov_test
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.kstest.html
# or https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ks_2samp.html
# ?
# Also https://en.wikipedia.org/wiki/Anderson%E2%80%93Darling_test
# https://asaip.psu.edu/articles/beware-the-kolmogorov-smirnov-test/
# OK Darling compares against a known distribution so that won't work ha(!)
# The null hypothesis is that both samples come from the same distribution and is not rejected (p-value = 0.5361) since they do come from the exact same distribution.

# We can try making this the metric but according to this answer on stackexchange https://stats.stackexchange.com/a/511714
# we might be better served leaving it as MSE (especially if it's not differentiable) and use that for evaluation only
#scipy.stats.kstest(rvs, cdf, args=(), N=20, alternative='two-sided', method='auto')


# In[33]:


# Regresssion
if __name__ == "__main__":

    # If chosen by hand
    param = {"tree_method"   : "hist",
             "n_estimators"  : 120, 
             "max_depth"     : 6, 
             "learning_rate" : 0.15}
             #"objective"     : "binary:logistic"}

    reg = xgb.XGBRegressor(#**param
                           **best_hyperparams_int
                           )
                           
    reg.fit(in_train, 
            out_train, 
            eval_set = [(in_train, out_train)],
            verbose = False)

    reg_pred = reg.predict(in_test)
    out_pred_df = pd.DataFrame(reg_pred, 
                               columns = out_test.columns, 
                               index = in_test.index
                              )

    # squared = False gives RMSE
    rmse = MSE(out_test, 
               reg_pred, 
               multioutput = "raw_values", 
               squared = False
              )

    print("RMSE Eval")
    for col, score in zip(out_train.columns, rmse):
        if "cumul_rot" in col or "position_angle" in col.lower():
            score *= 180/np.pi
        print(f"{col}: {score:.4f}")

    print()
    print("KSTest eval")
    # Null hypothesis -- the samples are pulled from the same distribution
    # typical choice is <0.05, reject null hypothesis
    kstest_results = {col : kstest(out_test[col].values, out_pred_df[col].values)
                      for col in out_pred_df}

    insig_count = 0
    for col, result in kstest_results.items():
        print(col)
        print(result)
        if result.pvalue < 0.05:
            result_str = "Significant: Test and prediction distributions differ significantly."
        else:
            result_str = "Insignificant: Test and prediction distributions do not differ signficantly."
            insig_count += 1
        print(result_str)
        print()
    print(f"{insig_count}/{len(kstest_results)} prediction distributions match within p < 0.05")

    #plot_predt(out_train.to_numpy(), reg_pred, "multi")
    # Closer to 1 is better 
    # https://xgboost.readthedocs.io/en/stable/python/python_api.html#xgboost.XGBRegressor.score
    # I'm thinking R^2 isn't a good test here, just looking at the data
    # reg.score(in_test, out_test)

    # TODO: Can we average the difference of every parameter per galaxy and use that to evaluate the model on a per galaxy basis?
    # With some weighting I think


# In[34]:


def anderson_darling(in_df):
    # Null hypothesis -- the samples are pulled from the normal distribution
    # statistic is the 'result', critical_values are the values at which we can check each
    # significance level specified by significance_level
    # The statistic must be greater than a critical value to determine a significance value
    #
    # ex: statistic = 0.8
    # significance array (%) = [15. , 10. ,  5. ,  2.5,  1. ]
    # critical values        = [0.552, 0.629, 0.755, 0.88 , 1.047]
    # Then the result is significant at the 5%/0.05 confidence level
    
    anderson_results = {col : anderson(in_df[col].values, dist = 'norm')
                        for col in in_df}
    
    crit_sig = dict(zip(list(anderson_results.values())[0].critical_values, 
                        list(anderson_results.values())[0].significance_level))
    
    cv = list(crit_sig.keys())
    sl = list(crit_sig.values())
    
    print(f"Critical Values: {cv}")
    print(f"Significance Levels: {sl}")
    print()
    
    for col, result in anderson_results.items():
        rs = result.statistic
        # Checking highest crit val the significance is greater than
        cv_filter = list(filter(lambda x: rs > x, cv))
        if cv_filter:
            significance_value = crit_sig[cv_filter[-1]]
        else:
            significance_value = "N/A"
            
        print(f"{col}\n{rs}, satisfies {significance_value} confidence level\n", sep = "")
    
    return anderson_results


# In[35]:


# ************************************************************************
# Should we expect these subsets to be normal?
# ************************************************************************
if __name__ == "__main__":
    print("\nAnderson-Darling eval of prediction (vs Normal Distribution)")
    _ = anderson_darling(out_pred_df)


# In[36]:


# ************************************************************************
# Could be used for interesting statistics once I'm actually confident in the algorithm ;)
# ************************************************************************
if __name__ == "__main__":
    print("\nAnderson-Darling eval of test output (vs Normal Distribution)")
    _ = anderson_darling(out_test)


# In[37]:


def make_hist_plots(test_data, predicted_data, grid = False, bins = 30, adjust = 2.17, save = False):
    # Assume both test_data, predicted_data are pandas df
    # Number of bins arrived at empirically

    rows = len(test_data.columns)
    cols = 1
    
    if grid:
        rows = ceil(rows / 2)
        # Restricting to 2 for now since histograms are wide
        cols = 2
        
    fig = make_subplots(rows = rows, cols = cols, start_cell="top-left") 

    for i, col_name in enumerate(predicted_data.columns): 

        # plotly indexes at 1 
        if grid:
            row = ceil((i + 1)/rows) 
            col = i % cols + 1 
        else:
            row = i + 1
            col = 1    

        #print(row, col) 

        fig.add_trace(go.Histogram(x = test_data[col_name],
                                   nbinsx = bins,
                                   name = "Data",
                                   legendgroup = col_name,
                                   legendgrouptitle_text = col_name,
                                   marker_color = '#FFA15A'), #col_name + " test"), 
                      row = row, col = col) 
        
        fig.add_trace(go.Histogram(x = predicted_data[col_name],
                                   nbinsx = bins,
                                   name = "Predicted",
                                   legendgroup = col_name,
                                   legendgrouptitle_text = col_name,
                                   marker_color = 'cornflowerblue'), #col_name + " pred"), 
                      row = row, col = col)

        fig.update_layout(barmode='overlay')#,
                          #xaxis_title_text = col_name)
            
        # if save:
        #     filepath = pj(cwd, col_name.replace(" ", "_"))
        #     fig.write_image(f"{filepath}.svg")

    fig.update_traces(opacity = 0.6)
    height = rows*200
    width = 800
    # For legend placement
    # higher, less distance, lower, more distance
    adjust = adjust
    fig.update_layout(height = height, width = width,
                      legend_tracegroupgap = height / (adjust*len(in_test.columns)))
    
    if save:
        filepath = pj(os.getcwd, "hist_plots", "all_plots")
        fig.write_image(f"{filepath}.svg")
    else:
        fig.show()


# In[38]:


if __name__ == "__main__":
    make_hist_plots(out_test, out_pred_df, adjust=1.83, save = True)


# In[39]:


if __name__ == "__main__":
    reg.save_model(pj(cwd, "xgboost_model.json"))
    # reg.load_model("xgboost_model.json")


# In[11]:


if __name__ == "__main__":
    export_to_py("xgboost_train", pj(_MODULE_DIR, "XGBoost", "xgboost_train"))
    #export_to_py("xgboost_train", "xgboost_train")

