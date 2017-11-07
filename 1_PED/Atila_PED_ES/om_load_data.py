#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 00:34:41 2017

@author: atilapaes

The main function (load_ms_data) is used to load MS data whitch is listed into a catalog.
The secondary function (split_3c_gph_signal) splits certain geophone's data from the whole stream. The info about the channels is imported from a catalog.
"""

#%% Auxiliar functions
def load_catalog(file_name):
    import pandas
    catalog = pandas.read_csv(file_name)
    return(catalog)
#%%
def load_ms_files(file1_name,file2_name, folder):
    """
    The file1_name and file2_name and folders are strings   
    """
    import obspy
    #loading first file
    ms_data_file1=obspy.read(folder+file1_name)

    #in case a second file is present
    if file2_name!='nan':
        ms_data_file2=obspy.read(folder+file2_name)       
        for channel in range(len(ms_data_file1)):
            ms_data_file1[channel] += ms_data_file2[channel]  
    return(ms_data_file1)
#%%
def slice_and_filter(ms_data,start_sec,end_sec):
    # Slice data
    ms_data=ms_data.slice(ms_data[0].stats.starttime+start_sec,ms_data[0].stats.starttime+end_sec)

    # Filtering
    ms_data.filter('bandpass', freqmin=30, freqmax=0.5*ms_data[0].stats.sampling_rate, corners=4, zerophase=True)
    return(ms_data)

#%%
def slice_and_filter_v2(ms_data,start_time,end_time):
    # Slice data
    ms_data=ms_data.slice(start_time,end_time)

    # Filtering
    ms_data.filter('bandpass', freqmin=30, freqmax=0.5*ms_data[0].stats.sampling_rate, corners=4, zerophase=True)
    return(ms_data)

#%% Main function for loading the data
"""
def load_ms_data(catalog_index,slice_sec_before,slice_sec_after,gph_set):

    Load and pre process a single or double MS files

    import obspy
    
    # Making the catalogs
    catalog_data=load_catalog(file_name_catalog_data)
    ms_data=''
    
    if catalog_data.valid_data[catalog_index]==True:
        # Event time    
        event_time=obspy.core.utcdatetime.UTCDateTime('20'+str(catalog_data.date[catalog_index])+str(catalog_data.time[catalog_index]))
        
        # Checking the gph_set to load correct files
        if gph_set=='3c':
            # Load raw data and merge
            ms_data=load_ms_files(file1_name=str(catalog_data.file_3c_identified[catalog_index]),file2_name=str(catalog_data.file_3c_auxiliar[catalog_index]), folder=folder_3c_files )    
            ms_data=slice_and_filter(ms_data,start_sec=event_time.second-slice_sec_before,end_sec=event_time.second+slice_sec_after)
        elif gph_set=='ch1':
            # Load raw data and merge
            ms_data=load_ms_files(file1_name=str(catalog_data.file_ch1_identified[catalog_index]),file2_name=str(catalog_data.file_ch1_auxiliar[catalog_index]), folder=folder_ch1_files )    
            ms_data=slice_and_filter(ms_data,start_sec=event_time.second-slice_sec_before,end_sec=event_time.second+slice_sec_after)
    else:
        print('====> THE INFORMED INDEX ('+str(catalog_index) +') CANNOT BE LOADED DUE MISSING FILES <=====')


    print('Data loaded. Date time: ', catalog_data.date_time_str[catalog_index])
    return(ms_data)
"""
#%% Secondary function for spliting the 
"""
def split_3c_gph_signal(ms_data,gph_number, gph_set):
    #
    NEEDS MAINTANACE 
    3C Case: Extract the 3C components. The new stream has the order z,h1,h2
    1C Case: Extract the 1C components.
    #
    import obspy
    catalog_gph=load_catalog(file_name_catalog_gph)    
    if gph_set=='3c':
        gph_signal=obspy.core.stream.Stream (traces=[ms_data[catalog_gph.gph_z[gph_number]], ms_data[catalog_gph.gph_h1[gph_number]],ms_data[catalog_gph.gph_h2[gph_number]]])
    elif gph_set=='ch1':
        gph_signal=obspy.core.stream.Stream (traces=[ms_data[gph_number]])
    return(gph_signal)
"""
#%% Split 3c data
def extract_3c_gph(ms_data,gph_number):
    import obspy, pandas    
    catalog_gph = pandas.read_csv('catalog_gph.csv')
    gph_signal=obspy.core.stream.Stream (traces=[ms_data[catalog_gph.gph_z[gph_number]], ms_data[catalog_gph.gph_h1[gph_number]],ms_data[catalog_gph.gph_h2[gph_number]]])  
    return(gph_signal)
