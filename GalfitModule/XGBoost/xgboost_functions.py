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

#!pip install xgboost
#!pip install hyperopt

from xgboost import XGBClassifier
import xgboost as xgb
from sklearn.model_selection import train_test_split, cross_val_score, RepeatedKFold
from sklearn.metrics import accuracy_score, SCORERS
from sklearn.metrics import mean_squared_error as MSE
#from sklearn.model_selection import GridSearchCV
from hyperopt import fmin, Trials, hp, tpe, STATUS_OK

import numpy as np
from math import ceil
from scipy.stats import kstest, anderson

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

# cwd = os.getcwd()
# home_dir = os.environ['HOME']


# In[6]:


# png_tiled_dir  = pj(cwd, "labeled_png_out")
# good_labeled_dir  = pj(png_tiled_dir, "good")

# inputs_dir = pj(cwd, "galfit_inputs")
# good_inputs_dir = pj(inputs_dir, "good")
# not_good_inputs_dir = pj(inputs_dir, "not_good")

# outputs_dir = pj(cwd, "galfit_outputs")
# good_outputs_dir = pj(outputs_dir, "good")
# not_good_outputs_dir = pj(outputs_dir, "not_good")

# sparc_out_dir = pj(home_dir, "run2_1000_galfit", "sparcfire-out")

# resize_obs_dir = os.path.join(png_resize_dir, "obs")
# resize_models_dir = os.path.join(png_resize_dir, "models")
# resize_residuals_dir = os.path.join(png_resize_dir, "residuals")

# features_files_dir = os.path.join(cwd, "image_feature_vectors")

# features_obs_dir = os.path.join(features_files_dir, "obs")
# features_models_dir = os.path.join(features_files_dir, "models")


# In[7]:


# Used to get organized
# Used to get organized
def generate_labels_from_imgs(all_galaxies_path, good_folder_path):
    
    # good galaxies are named ###_combined.png
    # inputs are named ###.in
    
    good_galaxy_names = list(map(lambda i: i.split("_")[0], os.listdir(good_folder_path)))
    all_galaxy_names = [os.path.basename(i) for i in os.listdir(all_galaxies_path) 
                        if os.path.basename(i).startswith("123")]
    
    label_dict = {i:(1 if i in good_galaxy_names else 0) for i in all_galaxy_names}
    
    return label_dict

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

# def organize_template_in_files(in_dict, top_dir, good_dir = "good", not_good_dir = "not_good"):
#     for gname, label in label_dict.items():
#         if label:
#             os.rename(pj(top_dir, gname + ".in"), pj(top_dir, good_dir, gname + ".in"))
#         else:
#             os.rename(pj(top_dir, gname + ".in"), pj(top_dir, not_good_dir, gname + ".in"))
            
#     return
# # In[8]:


# def organize_template_files(in_dict, sparcfire_out_dir, labels_top_dir = os.getcwd(), good_name = "good", not_good_name = "not_good"):
#     for gname, label in label_dict.items():
#         if label:
#             label_name = good_name
#         else:
#             label_name = not_good_name
        
#         # There will always be a generated feedme
#         # There should always be a generated output (because labeled) but just in case...
#         # Copy that first so the except catches before trying to copy the input that way we don't have a mismatch
#         try:
#             copy2(pj(sparcfire_out_dir, gname, "galfit.01"), pj(labels_top_dir, "galfit_outputs", label_name, gname + ".out"))
#             copy2(pj(sparcfire_out_dir, gname, "autogen_feedme_galfit.in"), pj(labels_top_dir, "galfit_inputs", label_name, gname + ".in"))
#         except FileNotFoundError:
#             print(f"No output found for {gname}, continuing to copy...")
#     return


# In[9]:


#label_dict = generate_labels_from_imgs(inputs_dir, good_labeled_dir)


# In[10]:


#organize_template_files(label_dict, pj(home_dir, "run2_1000_galfit", "sparcfire-out"))


# In[11]:


# Used if already organized
def generate_labels_from_paths(top_dir = os.getcwd(), good_dir = "good", not_good_dir = "not_good"):
    
    # good galaxies are named ###_combined.png
    # inputs are named ###.in
    
    good_fits = list(map(lambda i: i.split(".")[0], os.listdir(pj(top_dir, good_dir))))
    not_good_fits = list(map(lambda i: i.split(".")[0], os.listdir(pj(top_dir, not_good_dir))))
    all_fits = good_fits + not_good_fits
    
    # 1 good, 0 bad
    label_dict = {i:(1 if i in good_fits else 0) for i in all_fits}
    label_dict.pop("")
    
    return label_dict


# In[159]:


# # Borrowed and modified from in_out_comparison
# def galfit_param_grab(filename):
    
#     # This function handles grabbing and storing the values from galfit files (input and output???)
#     # It's written to generally handle both and stores everything in a dict of dicts
#     # See parameter_names list for the corresponding labels
#     # {Component: {Parameter : Value}}
    
#     try: 
#         # Grabbing the filename
#         #input_filename = glob_name(galaxy_path, '', filename) 
#         input_file = open(filename,'r')
    
#     except:
#         print("Can't open to read the " + filename + ". Continuing...")
    
#     else:
#         input_in = input_file.readlines()
#         input_file.close()
        
#         # Initialzing parameters.
#         input_dict = {}
#         #param_dict = {}
#         param_list = []
#         parameter_names = []
#         component_number = "0"
#         fourier_num = 1
#         #component_name = "bad programmer bad!"
        
#         # For convenient checking down the line.
#         n_range = [str(x) for x in range(3,11)] +                   ['R' + str(x) for x in range(1,11)] +                   ['F' + str(x) for x in range(1,6)]

#         # Looping through the galfit file.
#         for line in input_in[1:]:
#             line = line.replace("\n","")
            
#             # Grabbing crop coordinates for sanity checking
#             if line.startswith("H)"):
#                 values = line.split()
#                 input_dict["Crop Rad"] = [0.5*(float(values[2]) - float(values[1]))]
#                 parameter_names.append("Crop Rad")
            
#             # For those pesky blanks
#             if not line:
#                 # I don't use continue here I think because the blank lines indicate the end
#                 # of a component
#                 line = 'hi hi'
                
#             # Split for convenience and removing the ) for consistency among parameter numbers
#             values = line.split()
#             values[0] = values[0].replace(")", "")
                
#             # Grabbing the component number and innitializing the parameter dictionary.
#             # Essentially by doing it here, we guarantee that as we go through the lines
#             # we initialize a new param_dict to hold that component's parameters.
#             # This only triggers at every new component so there's no chance of triggering
#             # while grabbing the rest. 
#             if "Component" in line and "number:" in line:
#                 component_number = line[-1]
#                 #param_dict = {}
#                 param_list = []
#                 param_suffix = ""
                
#                 # To avoid repeating parameter names per sklearn/xgboost
#                 if component_number == "1":
#                     param_suffix = " bulge"
                    
#                 elif component_number == "2":
#                     param_suffix = " disk"
                    
#                 # For sky component
#                 elif component_number == "3": 
#                     n_range += ["1", "2"]
             
#             # Overwrites sersic unfortunately...
#             #elif values[0] == "0":
#                 #component_name = values[1]

#             # Storing the parameters themselves
#             elif values[0] in n_range:
                
#                 param_name = line.split("#")[-1].strip() + param_suffix
#                 param_name = param_name.replace("[", "(")
#                 param_name = param_name.replace("]", ")")

#                 # Accounting for the Fourier modes
#                 # Also we have to put everything in a list to appease pandas
#                 if line[0] == 'F':
#                     parameter_names.append(f"Fourier Amplitude {fourier_num}")
#                     parameter_names.append(f"Fourier Phase Angle {fourier_num}")
#                     fourier_num += 1
#                     #param_dict["Fourier Amplitude"] = [values[1]]
#                     #param_dict["Fourier Phase Angle"] = [values[2]]
#                     param_list.append(values[1])
#                     param_list.append(values[2])

#                 else:
#                     parameter_names.append(param_name)
#                     #param_dict[param_name] = [values[1]]
#                     # Convert to radians for testing loss function
#                     if "degree" in param_name.lower() or "angle" in param_name.lower():
#                         values[1] = str(float(values[1]) * np.pi/180)
                        
#                     param_list.append(values[1])

#             # There's a blank line between parameters so this is where we finalize the 
#             # parameter dictionary and place it in the main dict.
#             # Resetting component_number just in case and to store any extra faff which will then be popped.
#             else:
#                 #input_dict[component_number] = param_dict
#                 input_dict[component_number] = param_list
#                 component_number = "0"

#         # Popping extra faff
#         input_dict.pop("0")
#         #input_dict.pop("bad programmer bad!")
#         #input_dict.pop('3')
#         #print(input_dict)
#         return input_dict, parameter_names

#     return None, None

def flatten_to_pandas(in_dict, parameter_names, gname):
    
    param_values = np.array([value for sublist in in_dict.values() for value in sublist])
    param_values = param_values.reshape(1, len(param_values))
        
    all_data = pd.DataFrame(param_values, index = [gname], columns = parameter_names, dtype = np.float32)
    all_data = all_data.drop(columns=['----- bulge'])
    all_data = all_data.drop(columns=['----- disk'])
    
#     for component, sub_list in in_dict.items():
#         print(component)
#         print(sub_dict)
#         component_df = pd.DataFrame.from_dict(sub_dict, dtype = float)
        
#         if not all_data.empty:
#             #all_data.merge(component_df, how = "outer")
#         else:
#             all_data = component_df
        
    return all_data

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

def choose_df(top_dir, label_dirs, file_suffix):
    # top_dir = inputs_dir
    # label_dirs = [not_good_inputs_dir, good_inputs_dir] 
    # file_suffix = ".in", ".out"
    label_dict = generate_labels_from_paths(top_dir = top_dir)
    out_df = build_df(label_dict, label_dirs = label_dirs, file_suffix = file_suffix)
    return out_df

# Deprecated
# Per https://xgboost.readthedocs.io/en/stable/python/python_intro.html
def split_save_df(in_df, file_prefix = "xgboost", col_ignore = [], sklearn = False):
    exclude = ['label'] + col_ignore
    data = in_df.loc[:, in_df.columns.difference(exclude)]
    label = in_df['label']
    
    # Pass random_state = 0 to guarantee input and output are lined up
    X_train, X_test, y_train, y_test = train_test_split(data, label, test_size=.3, random_state = 0)
    
    dtrain = xgb.DMatrix(X_train, label = y_train)
    dtest  = xgb.DMatrix(X_test,  label = y_test)
    ddata  = xgb.DMatrix(data, label = label)
    
    dtrain.save_binary(f'{file_prefix}.train')
    dtest.save_binary(f'{file_prefix}.test')
    ddata.save_binary(f'{file_prefix}.data')
    
    if sklearn:
        return X_train, X_test, y_train, y_test, None
    
    return dtrain, dtest, ddata, data, label
# In[455]:


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


# In[18]:


def load_data(file_prefix):
    dtrain = xgb.DMatrix(f'{file_prefix}.train')
    dtest  = xgb.DMatrix(f'{file_prefix}.test')
    ddata  = xgb.DMatrix(f'{file_prefix}.data')
    
    return dtrain, dtest, ddata

# In[142]:


def make_hist_plots(test_data, predicted_data, grid = False, bins = 30):
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

    fig.update_traces(opacity = 0.6)
    height = rows*200
    width = 800
    # For legend placement
    # higher, less distance, lower, more distance
    adjust = 2.17
    fig.update_layout(height = height, width = width,
                      legend_tracegroupgap = height / (adjust*len(in_test.columns)))

    fig.show()
    
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
