#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 18:34:48 2017

@author: atilapaes

This file is used to create the IO files for Open Microseismic project
"""


#%% Create files of filenames and time


# For later importing the info in the output files
#sample=pandas.read_csv('sample.csv')
#sample['date_time'] = list(map(lambda x: obspy.core.utcdatetime.UTCDateTime(x), sample.date_time))

dataset_name='sample'
files_path='Files/'
extension='dat'

def output_file_datetime(dataset_name,files_path,extension):
    """
    This function generates the file contaning the file names and its start times
    """
    import obspy, pandas
    from os import listdir
    from os.path import isfile, join
    ms_files = [f for f in listdir(files_path) if isfile(join(files_path, f)) and f.endswith("."+extension)]

    starttime=[]
    for index in range(len(ms_files)):
        ms_data=obspy.read(files_path+ms_files[index])
        starttime.append(ms_data[0].stats.starttime)

    file_time_df=pandas.DataFrame({'file_name':ms_files,'date_time':starttime})
    file_time_df=file_time_df.sort_values('date_time')
    file_time_df.to_csv(dataset_name+'.csv', sep=',', line_terminator='\n', index=False)

output_file_datetime(dataset_name,files_path,extension)

#%%

