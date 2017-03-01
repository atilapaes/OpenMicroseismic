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

import multiprocessing

##################################################################
#%% Name of the dataset folder
#dataset_folder="HornRiver-C88-S9F"
dataset_folder="Dataset2"

##################################################################
#%% Verbose levels
# 0 - Everything. The code will output everything it is doing. Best for debug.
# 1 - Entering and leaving main processes. Avg performance with tracking of main activities.
# 2 - No versose at all. Recommended after extensive testing.
verbose_level = 0

##################################################################
#%% Number of cores to be used
use_maximum_cores = False
core_to_be_used=1 #Define number of cores to use other than maximum

if use_maximum_cores == True:
    number_of_cores = multiprocessing.cpu_count()
else: 
    number_of_cores=core_to_be_used 

##################################################################
#%% Data processing parameters
 
#Energy stack processing parameters
moving_average_samples=200

##################################################################
#%% Detecting peaks
#Peak identification parameters
threshold_stds=0 #3 is too restrictive cause it may not detect low-SNR event beside a high-SNR 
threshold_means=2
mpd=500 # Minimum peak distance in samples. Due QC, 100 is reasonable number.

# Peaks quality control: 
peak_width_minimum=500  # Minimum peak width 
peak_snr_minimum=2      # Minimum peak SNR

# Data torrent parameters 
last_seconds=6  # Last_seconds = LTA +0.5 sec (Coef. time trust)+ 2 seconds (peak max width)
time_torrent_coef=2 #Time in seconds to use before the file start and before the file ends

##################################################################