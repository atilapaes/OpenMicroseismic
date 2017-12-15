#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 18:34:48 2017

@author: atilapaes

This file is used to create the IO files for Open Microseismic project
"""


#%% Create files of filenames and time

dataset_name='sample' # Used to name of the CSV output file
files_path='Files/'   # Folder contning the Microseismic files
extension='dat'      # File extension without '.'. Ex: dat, SEGY, seg2...

def io_file_datetime(dataset_name,files_path,extension):
    """
    This function generates the file contaning the file names and its start times


    # For later importing the info in the output files
    #sample=pandas.read_csv('sample.csv')
    #sample['date_time'] = list(map(lambda x: obspy.core.utcdatetime.UTCDateTime(x), sample.date_time))

    """
    import obspy, pandas
    from os import listdir
    from tqdm import tqdm # Progress bar
    from os.path import isfile, join
    
    print('\n','---> Listing files...')
    ms_files = [f for f in listdir(files_path) if isfile(join(files_path, f)) and f.endswith("."+extension)]

    print('---> Reading files...')
    starttime=[]
    
    for index in tqdm(range(len(ms_files)),desc='Getting starttime'):
        ms_data=obspy.read(files_path+ms_files[index])
        starttime.append(ms_data[0].stats.starttime)
    
    print('\n','---> Exporting results...')
    file_time_df=pandas.DataFrame({'file_name':ms_files,'date_time':starttime})
    file_time_df=file_time_df.sort_values('date_time')
    file_time_df.to_csv(dataset_name+'.csv', sep=',', line_terminator='\n', index=False)
    print('Done!')
    
io_file_datetime(dataset_name,files_path,extension)

#%%