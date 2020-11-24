#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Author: Matthew Portman
# Date: 11/4/2020 - see GitHub for latest date

import glob
import csv


# In[2]:


# Grabbing the file names
def get_galaxy_names_list():
    
    # This function grabs the names of all input files in sparcfire-in and does some python string manipulation
    # to gab the names of the input files (with relative path), the numbers themselves, and the output folders
    # again with relative path

    try:
        filenames_read = glob.glob("sparcfire-in/*.fits") # Hardcoding is a temporary measure.
    
    except:
        print("Please copy me into the directory which contains the folders for your")
        print("input, temporary, and output files for SpArcFiRe denoted:")
        print("sparcfire-in, sparcfire-tmp, and sparcfire-out.")
        raise SystemExit("Exitting.")
        
    else:
        filenames_out = [s.split(".")[0] for s in filenames_read]
        galaxy_names_out = [s.split("/")[1] for s in filenames_out]
        filenames_out = [s.replace("in", "out") for s in filenames_out]
        
    return filenames_read, galaxy_names_out, filenames_out


# In[3]:


def glob_name(path='', name='', desired_file=''):
    
    # A helper function with glob to grab a desired file as input
    
    file_path = "./" + path + '/' + name + desired_file
    #print(file_path)
    file_name = glob.glob(file_path)[0]
    
    return file_name


# # Old version please ignore.
# 
# def input_grab(galaxy_name, galaxy_path):
#     # for auto-generated input
#     
#     try: 
#         input_filename = glob_name(galaxy_path, '', 'autogen_feedme_galfit.in') 
#         input_file = open(input_filename,'r')
#     
#     except:
#         raise SystemExit("Can't open to read the feedme. Exitting.")
#     
#     else:
#         input_in = input_file.read()
#         input_file.close()
#         
#         input_in = input_in.split("# Component number: ")[1:]
#         
#         input_dict = {}
#         count = 1
#         for component in input_in:
#             param_list = component.split("\n")
#             param_list = [x for x in param_list if x != '']
#             #print(param_list)
#             input_dict[count] = dict(x.split(')') for x in param_list[3:])
#             count += 1
# 
#     input_dict.pop(3)
#     return input_dict

# In[4]:


def input_grab(galaxy_path, filename):
    
    # This function handles grabbing and storing the values from galfit files (input and output)
    # It's written to generally handle both and stores everything in a nested dictionary
    # {Component: {Parmaeter : Value}}
    
    try: 
        # Grabbing the filename
        input_filename = glob_name(galaxy_path, '', filename) 
        input_file = open(input_filename,'r')
    
    except:
        raise SystemExit("Can't open to read the " + filename + ". Exitting.")
    
    else:
        input_in = input_file.readlines()
        input_file.close()
        
        # Initialzing parameters.
        input_dict = {} 
        param_dict = {} 
        component_number = 0
        
        # For convenient checking down the line.
        n_range = [str(x) for x in range(3,11)] +             ['R' + str(x) for x in range(1,11)] +             ['F' + str(x) for x in range(1,6)]

        # Looping through the galfit file.
        for line in input_in[1:]:
            line = line.replace("\n","")
            
            # For those pesky blanks
            if not line:
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
                param_dict = {}

            # Storing the parameters themselves
            elif values[0] in n_range:

                # Accounting for the Fourier modes
                if line[0] == 'F':
                    param_dict[values[0]] = values[1] + " " + values[2]

                else:
                    param_dict[values[0]] = values[1]

            # There's a blank line between parameters so this is where we finalize the 
            # parameter dictionary and place it in the main dict.
            # Resetting component_number just in case and to store any extra faff which will then be popped.
            else:
                input_dict[component_number] = param_dict
                component_number = 0

        # Popping extra faff and sky component
        input_dict.pop(0)
        input_dict.pop('3')
        #print(input_dict)
        return input_dict

    return None


# # Old version, please ignore.
# 
# def output_grab(galaxy_name, galaxy_path):
#     # for autocrop coordinates
#     
#     try: 
#         output_filename = glob_name(galaxy_path, '', 'galfit.01') 
#         output_file = open(output_filename,'r')
# 
#     except:
#         raise SystemExit("Can't open to read the GALFIT output. Exitting.")
#     
#     else:
#         output_in = output_file.read()
#         output_file.close()
#         
#         output_in = output_in.split("# Component number: ")[1:]
#         
#         output_dict = {}
#         count = 1
#         for component in output_in:
#             param_list = component.split("\n")
#             param_list = [x for x in param_list if x != '']
#             param_list = [x.split('#')[0].strip() for x in param_list[3:]]
#             output_dict[count] = dict(x.split(')') for x in param_list[:-1])
#             count += 1
#         
#     output_dict.pop(3)
#     return output_dict

# In[ ]:


def read_write_compare():
    
    # This function is the final bulk of the code which uses the previous functions and
    # writes the comparisons to two files, a text per galaxy in each galaxy folder in sparcfire-out 
    # and a file which contains just the difference in the output and input values (output - input)
    # placed in all_galfit_out.
    
    paths_to_feedme = []
    
    # Grabbing pertinent file directory info and galaxy names.
    filenames_in, galaxy_names, filenames_out = get_galaxy_names_list()
    
    compare_list = []
    
    # Looping through galaxy folders.
    for galaxy_path in filenames_out:
        
        # Grabbing the galaxy name.
        gname = galaxy_path.split('/')[-1]
        
        # Taking in the parameters from the input and output galfit files.
        feedme_dict = input_grab(galaxy_path, 'autogen_feedme_galfit.in')
        galfit_dict = input_grab(galaxy_path, 'galfit.01')

        # Initializing
        bulge_param = [] 
        disk_param = []
        power_param = []
        fourier_list = []
        
        bulge_diff = [gname]
        disk_diff = []
        power_diff = []
        
        # Convenient list for checking down the line
        sersic_params = ['Mag', 'R_e', 'n', '', '', '', 'axis ratio', 'PA']        
        p_params = ['R_in', 'R_out', 'Theta_out', 'Alpha', 'Inclination', 'Sky PA']
        
        # Looping through parameters
        for i in range(3, 11):
            
            # Taking the bulge parameters, input, output, and difference
            bulge_param.append(['Bulge ' + str(i) + ' - ' + sersic_params[i - 3],                                 feedme_dict['1'][str(i)], galfit_dict['1'][str(i)],                                 float(galfit_dict['1'][str(i)]) - float(feedme_dict['1'][str(i)])])
            
            # Putting the difference values in a separate list for those which aren't empty in sersic_params
            if sersic_params[i - 3]:
                bulge_diff.append(float(galfit_dict['1'][str(i)]) - float(feedme_dict['1'][str(i)]))
            
            # Disk parameters
            disk_param.append(['Disk ' + str(i) + ' - ' + sersic_params[i - 3],                                feedme_dict['2'][str(i)], galfit_dict['2'][str(i)],                                float(galfit_dict['2'][str(i)]) - float(feedme_dict['2'][str(i)])])
            
            # Difference in disk parameters
            if sersic_params[i - 3]:
                disk_diff.append(float(galfit_dict['2'][str(i)]) - float(feedme_dict['2'][str(i)]))

            # Some naming convention things
            if i < 7:
                k = 'R' + str(i - 2)
            else:
                k = 'R' + str(i + 2)

            # Power parameters
            if i < 9:
                power_param.append(['Power ' + k  + ' - ' + p_params[i - 3], feedme_dict['2'][k], galfit_dict['2'][k],                                     float(galfit_dict['2'][k]) - float(feedme_dict['2'][k])])
                
                power_diff.append(float(galfit_dict['2'][k]) - float(feedme_dict['2'][k]))
                
            # Fourier Parameters
            f1 = 'F' + str(i - 2)
            try:
                fourier_list.append(['Power ' + f1  + ' - Amp & Phase Angle', feedme_dict['2'][f1], galfit_dict['2'][f1], ''])
            except:
                pass
        
        # Updating comparison list and power parameter list with final values
        compare_list.append(bulge_diff + disk_diff + power_diff)
        power_param += fourier_list
        
        i_filename = './' + galaxy_path + '/' + 'galfit_io_compare.txt'
        
        # Writing to the individual galaxy comparison file
        with open(i_filename,'w') as i_file:
            
            i_file.write('Bulge Parameters:\n')
            i_file.writelines([', '.join(map(str, line))+'\n' for line in bulge_param])
            
            i_file.write('\nDisk Parameters:\n')
            i_file.writelines([', '.join(map(str, line))+'\n' for line in disk_param])
            
            i_file.write('\nPower Parameters:\n')
            i_file.writelines([', '.join(map(str, line))+'\n' for line in power_param])
    
    # Header for csv file and field names
    sersic_params = ['Mag', 'R_e', 'n', 'axis ratio', 'PA'] # Because python can't remove more than one instance easily
    field_headers = ['name'] + sersic_params + sersic_params + p_params 
    
    all_filename = './sparcfire-out/all_galfit_out/comparison_params.csv'
    
    # Writing csv to file
    with open(all_filename, 'w') as all_file:    
        csvwriter = csv.writer(all_file)  
        
        # writing the field headers && writing the data rows 
        csvwriter.writerow(field_headers)  
        csvwriter.writerows(compare_list)
        
    return None


# In[ ]:


if __name__ == "__main__":
    
    # Running all of the above. 
    
    read_write_compare()