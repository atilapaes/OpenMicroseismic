#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 04:01:07 2017

@author: atilapaes

Quality control filter for ATP CSV files
"""

import os, pandas, numpy, math

#%% Parameters
folder_csv_files='test_folder/'
file_list=[]

#%% List CSV files
for file in os.listdir(folder_csv_files):
    if file.endswith(".csv"):
        file_list.append(file)

#%%
for file_name in file_list:
    print(file_name)
    atp_data=pandas.read_csv(folder_csv_files+file_name)
    delta_t=[]
    delta_t_filtered=[]
    
    for gph_index in range(len(atp_data)):
        if (not(math.isnan(atp_data.pp[gph_index])) and not(math.isnan(atp_data.sp1[gph_index])) and not(math.isnan(atp_data.sp2[gph_index]))):
            delta_t.append(atp_data.sp1[gph_index]-atp_data.pp[gph_index])
        else:
            delta_t.append(None)
    
    delta_t_mean=numpy.mean(list(filter(lambda x: x!=None, delta_t)))
    delta_t_std=numpy.std(list(filter(lambda x: x!=None, delta_t)))
    
    #delta_t_mean=numpy.mean(delta_t)
    #delta_t_std=numpy.std(delta_t)
    
    for gph_index in range(len(atp_data)):
        if (delta_t[gph_index] != None):
            if ((delta_t[gph_index] < (delta_t_mean+delta_t_std)) or (delta_t[gph_index] > (delta_t_mean - delta_t_std))):
                delta_t_filtered.append(delta_t[gph_index])
            else:
                delta_t_filtered.append(None)
        else:
            delta_t_filtered.append(None)
            
            
    atp_data['delta_t']=delta_t
    atp_data['delta_t_filtered']=delta_t_filtered
    atp_data.to_csv(folder_csv_files+file_name, sep=',', line_terminator='\n', index=False)