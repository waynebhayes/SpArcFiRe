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


# In[163]:


#!pip install xgboost
#!pip install hyperopt


# In[240]:


from xgboost import XGBClassifier
import xgboost as xgb
from sklearn.model_selection import train_test_split, cross_val_score, RepeatedKFold
from sklearn.metrics import accuracy_score, SCORERS
from sklearn.metrics import mean_squared_error as MSE
from hyperopt import fmin, Trials, hp, tpe, STATUS_OK


# In[375]:


import numpy as np
from math import ceil
import itertools

from PIL import Image
import matplotlib.pyplot as plt

import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots

import graphviz
import ssl
from glob import glob

#from annoy import AnnoyIndex
import random
import pandas as pd
from scipy import spatial

import os
from os.path import join as pj
from shutil import copy2


# In[5]:


cwd = os.getcwd()
home_dir = os.environ['HOME']


# In[6]:


png_tiled_dir  = pj(cwd, "labeled_png_out")
good_labeled_dir  = pj(png_tiled_dir, "good")

inputs_dir = pj(cwd, "galfit_inputs")
good_inputs_dir = pj(inputs_dir, "good")
not_good_inputs_dir = pj(inputs_dir, "not_good")

outputs_dir = pj(cwd, "galfit_outputs")
good_outputs_dir = pj(outputs_dir, "good")
not_good_outputs_dir = pj(outputs_dir, "not_good")

sparc_out_dir = pj(home_dir, "run2_1000_galfit", "sparcfire-out")

# resize_obs_dir = os.path.join(png_resize_dir, "obs")
# resize_models_dir = os.path.join(png_resize_dir, "models")
# resize_residuals_dir = os.path.join(png_resize_dir, "residuals")

# features_files_dir = os.path.join(cwd, "image_feature_vectors")

# features_obs_dir = os.path.join(features_files_dir, "obs")
# features_models_dir = os.path.join(features_files_dir, "models")


# In[7]:


# Used to get organized
def generate_labels_from_imgs(all_inputs_path, good_folder_path):
    
    # good galaxies are named ###_combined.png
    # inputs are named ###.in
    
    good_galaxy_names = list(map(lambda i: i.split("_")[0], os.listdir(good_folder_path)))
    all_galaxy_names = list(map(lambda i: i.replace(".in", ""), os.listdir(all_inputs_path)))
    
    label_dict = {i:(1 if i in good_galaxy_names else 0) for i in all_galaxy_names}
    
    return label_dict

def organize_template_in_files(in_dict, top_dir, good_dir = "good", not_good_dir = "not_good"):
    for gname, label in label_dict.items():
        if label:
            os.rename(pj(top_dir, gname + ".in"), pj(top_dir, good_dir, gname + ".in"))
        else:
            os.rename(pj(top_dir, gname + ".in"), pj(top_dir, not_good_dir, gname + ".in"))
            
    return
# In[8]:


def organize_template_files(in_dict, sparcfire_out_dir, labels_top_dir = os.getcwd(), good_name = "good", not_good_name = "not_good"):
    for gname, label in label_dict.items():
        if label:
            label_name = good_name
        else:
            label_name = not_good_name
        
        # There will always be a generated feedme
        # There should always be a generated output (because labeled) but just in case...
        # Copy that first so the except catches before trying to copy the input that way we don't have a mismatch
        try:
            copy2(pj(sparcfire_out_dir, gname, "galfit.01"), pj(labels_top_dir, "galfit_outputs", label_name, gname + ".out"))
            copy2(pj(sparcfire_out_dir, gname, "autogen_feedme_galfit.in"), pj(labels_top_dir, "galfit_inputs", label_name, gname + ".in"))
        except FileNotFoundError:
            print(f"No output found for {gname}, continuing to copy...")
    return


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


# Borrowed and modified from in_out_comparison
def galfit_param_grab(filename):
    
    # This function handles grabbing and storing the values from galfit files (input and output???)
    # It's written to generally handle both and stores everything in a dict of dicts
    # See parameter_names list for the corresponding labels
    # {Component: {Parameter : Value}}
    
    try: 
        # Grabbing the filename
        #input_filename = glob_name(galaxy_path, '', filename) 
        input_file = open(filename,'r')
    
    except:
        print("Can't open to read the " + filename + ". Continuing...")
    
    else:
        input_in = input_file.readlines()
        input_file.close()
        
        # Initialzing parameters.
        input_dict = {}
        #param_dict = {}
        param_list = []
        parameter_names = []
        component_number = "0"
        fourier_num = 1
        #component_name = "bad programmer bad!"
        
        # For convenient checking down the line.
        n_range = [str(x) for x in range(3,11)] +                   ['R' + str(x) for x in range(1,11)] +                   ['F' + str(x) for x in range(1,6)]

        # Looping through the galfit file.
        for line in input_in[1:]:
            line = line.replace("\n","")
            
            # Grabbing crop coordinates for sanity checking
            if line.startswith("H)"):
                values = line.split()
                input_dict["Crop Rad"] = [0.5*(float(values[2]) - float(values[1]))]
                parameter_names.append("Crop Rad")
            
            # For those pesky blanks
            if not line:
                # I don't use continue here I think because the blank lines indicate the end
                # of a component
                line = 'hi hi'
                
            # Split for convenience and removing the ) for consistency among parameter numbers
            values = line.split()
            values[0] = values[0].replace(")", "")
                
            # Grabbing the component number and innitializing the parameter dictionary.
            # Essentially by doing it here, we guarantee that as we go through the lines
            # we initialize a new param_dict to hold that component's parameters.
            # This only triggers at every new component so there's no chance of triggering
            # while grabbing the rest. 
            if "Component" in line and "number:" in line:
                component_number = line[-1]
                #param_dict = {}
                param_list = []
                param_suffix = ""
                
                # To avoid repeating parameter names per sklearn/xgboost
                if component_number == "1":
                    param_suffix = " bulge"
                    
                elif component_number == "2":
                    param_suffix = " disk"
                    
                # For sky component
                elif component_number == "3": 
                    n_range += ["1", "2"]
             
            # Overwrites sersic unfortunately...
            #elif values[0] == "0":
                #component_name = values[1]

            # Storing the parameters themselves
            elif values[0] in n_range:
                
                param_name = line.split("#")[-1].strip() + param_suffix
                param_name = param_name.replace("[", "(")
                param_name = param_name.replace("]", ")")

                # Accounting for the Fourier modes
                # Also we have to put everything in a list to appease pandas
                if line[0] == 'F':
                    parameter_names.append(f"Fourier Amplitude {fourier_num}")
                    parameter_names.append(f"Fourier Phase Angle {fourier_num}")
                    fourier_num += 1
                    #param_dict["Fourier Amplitude"] = [values[1]]
                    #param_dict["Fourier Phase Angle"] = [values[2]]
                    param_list.append(values[1])
                    param_list.append(values[2])

                else:
                    parameter_names.append(param_name)
                    #param_dict[param_name] = [values[1]]
                    # Convert to radians for testing loss function
                    if "degree" in param_name.lower() or "angle" in param_name.lower():
                        values[1] = str(float(values[1]) * np.pi/180)
                        
                    param_list.append(values[1])

            # There's a blank line between parameters so this is where we finalize the 
            # parameter dictionary and place it in the main dict.
            # Resetting component_number just in case and to store any extra faff which will then be popped.
            else:
                #input_dict[component_number] = param_dict
                input_dict[component_number] = param_list
                component_number = "0"

        # Popping extra faff
        input_dict.pop("0")
        #input_dict.pop("bad programmer bad!")
        #input_dict.pop('3')
        #print(input_dict)
        return input_dict, parameter_names

    return None, None


# In[13]:


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


# In[399]:


def build_df(label_dict: dict, label_dirs, file_suffix: str):
    
    out_df = pd.DataFrame()

    for gname, label in label_dict.items():

        # 0 bad, 1 good
        if label == 0:
            input_dir = label_dirs[0]
            #str_label = "not_good"
        else:
            input_dir = label_dirs[1]
            #str_label = "good"

        gal_dict, param_names = galfit_param_grab(pj(input_dir, gname + file_suffix))
        if not param_names: continue
        
        for i, param in enumerate(param_names):
            if param in param_names[:i]:
                param += ""

        gname_df = flatten_to_pandas(gal_dict, param_names, gname)
        gname_df["label"] = label #str_label
        out_df = pd.concat([out_df, gname_df])

    # For the input galaxies, we have a lot of held values, these are uneccessary
    # https://stackoverflow.com/a/39658662
    nunique = out_df.nunique()
    cols_to_drop = nunique[nunique == 1].index
    out_df.drop(columns = cols_to_drop, inplace = True)
    
    return out_df


# In[15]:


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
    
    # Generator go fast brrrrr
    # _ = (
    # df.drop(index = galaxies_to_drop, inplace = True) for df in args
    # )
    
    for in_df, file_prefix in zip(args, file_prefixes):
        # Filter and exclude
        
        exclude = ['label'] + col_ignore.get(file_prefix, [])
        
        in_df_good = in_df.drop(index = galaxies_to_drop, columns = exclude)
        
        label = pd.DataFrame(np.ones(len(in_df_good)), columns = ["label"], dtype = "int")
            
        # Pass random_state = 0 to guarantee input and output are lined up
        X_train, X_test, y_train, y_test = train_test_split(in_df_good, label, test_size=.3, random_state = 0)
    
        dtrain = xgb.DMatrix(X_train, label = y_train)
        dtest  = xgb.DMatrix(X_test,  label = y_test)
        ddata  = xgb.DMatrix(in_df_good, label = label)

        print(f"Saving dmatrices to file: {file_prefix}.train/test/data")
        dtrain.save_binary(f'{file_prefix}.train')
        dtest.save_binary(f'{file_prefix}.test')
        ddata.save_binary(f'{file_prefix}.data')
        
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

def export_filter():
    ignore_galfit_in  = ["Bar radius (pixels) disk", "Fourier Amplitude 2"]
    ignore_galfit_out = ["Crop Rad",
                 "R_e (effective radius)   (pix) bulge", # Bulge radius goes crazzyyyyyy
                 "Inclination to L.o.S. (degrees) disk", # Trust sparc/galfit on this one
                 "Spiral inner radius (pixels) disk",
                 "Spiral outer radius (pixels) disk",
                 "Fourier Amplitude 1",
                 "Fourier Amplitude 2",
                 "Fourier Phase Angle 1",
                 "Fourier Phase Angle 2",
                 "Sky background at center of fitting region (ADUs)", # Trust galfit on this one
                 "dsky/dx (sky gradient in x)     (ADUs/pix)", # if we include them, does it get better? for future testing
                 "dsky/dy (sky gradient in y)     (ADUs/pix)"
                ]
    in_filter = []
    return ignore_galfit_in, ignore_galfit_out, in_filter