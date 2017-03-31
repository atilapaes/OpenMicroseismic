#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 14:09:53 2017

@author: Atila Paes

Project:    OpenMicroseismic
Pack:       Potential-Event Detection
Method:     Energy-stack
File:       Parameters
"""

import multiprocessing, os

##################################################################
#%% Dataset I/O:
dataset_folder="Dataset_name"  #Folder contaning the MS files
save_plot=True                 #Save stack plots for future consult (Must create a folder with name dataset_folder+'_plots')

#Verifing if the result folder exists. In negative case, create it
if not os.path.exists(dataset_folder+"_plots"):
    os.makedirs(dataset_folder+"_plots")

### Data injection mode
# 0 - Default - need only the dataset_folder. Concatenate files with part of the previous one.
# 1 - Input a time list. The file is a csv with time in first column in the format YYYYmmDDhhMMSS
# 2 - Calculate file by file in a folder
data_injection_mode=2

##################################################################
#%% Verbose levels
# 0 - Everything. The code will output everything it is doing. Best for debug.
# 1 - Entering and leaving main processes. Avg performance with tracking of main activities.
# 2 - Just the highlights and finished % . Recommended after extensive testing.
verbose_level = 1

##################################################################
#%% Number of cores to be used
use_maximum_cores = True
core_to_be_used=1 #Define number of cores to use other than maximum

if use_maximum_cores == True:
    number_of_cores = multiprocessing.cpu_count()
else: 
    number_of_cores=core_to_be_used 

##################################################################
#%% Data processing parameters
 
#Number of samples in the moving avg- 200 is our golden value. In lower cases, the signal zig zag to fast
moving_average_samples=200

##################################################################
#%% Detecting peaks
#Peak identification parameters
threshold_stds=0    #3 Standart deviation is not a interesting parameter
threshold_means=2   # Treshold based in a multiple of the mean
mpd=500             # Minimum distance (in samples) between peaks. Due QC, 100 is minimum reasonable number.

# Peaks quality control: 
peak_width_minimum=500  # Minimum peak width
peak_width_snr=500  # Minimum peak width 
peak_snr_minimum=2  # Minimum peak SNR -> NOT USED IN THE CODE

# Data torrent parameters 
last_seconds=6  # Last_seconds = LTA +0.5 sec (Coef. time trust)+ 2 seconds (peak max width)
time_torrent_coef=2 #Time in seconds to use before the file start and before the file ends

##################################################################