#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 14:09:53 2017
Modified on Thu Feb 16

@author: Atila Paes

Project:    OpenMicroseismic
Pack:       Potential-Event Detection
Method:     Energy-stack
File:       Parameters
Version:    V2
"""

import multiprocessing

##################################################################
#%% Name of the dataset folder
dataset_folder="HornRiver-C88-S9F"
#dataset_folder="Dataset"


##################################################################
#%% Verbose levels
# 0 - Everything. The code will output everything it is doing. Best for debug.
# 1 - Entering and leaving main processes. Avg performance with tracking of main activities.
# 2 - No versose at all. Recommended after extensive testing.
verbose_level = 1

##################################################################
#%% Number of cores to be used
use_maximum_cores = False
core_to_be_used=1 #Define number of cores to use other than maximum

if use_maximum_cores == True:
    number_of_cores = multiprocessing.cpu_count()
else: 
    number_of_cores=core_to_be_used 

#%%
 
#Energy stack processing parameters
moving_average_samples=100

#%% Detecting peaks

threshold_stds=2 #3 is too restrictive cause it may not detect low-SNR event beside a high-SNR 


# Detect peaks that are at least separated by minimum peak distance (in  number of data). Due QC, 100 is reasonable number.
mpd=100 

# Quality control: 
peak_width_minimum=500  # Minimum peak width 
peak_snr_minimum=4      # Minimum peak SNR