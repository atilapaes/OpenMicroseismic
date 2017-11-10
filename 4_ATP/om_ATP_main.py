#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 03:58:22 2017

@author: atilapaes

This is the main script of Arrival-Time Picking (ATP) of the Open Microseismic project

Future activities:

1) Document the input files
"""
# Folder of shared modules
import sys
sys.path.insert(0,'')

#%%
import pandas, obspy,numpy
import om_load_data, om_ATP_functions

#%% Folders, files and output file
dataset_folder=''

atp_folder=''
catalog_gph=pandas.read_csv('catalog_gph.csv')

catalog_file_PE=''
catalog_PE=pandas.read_csv(catalog_file_PE)
#%%

# Loop over all Potential Events (PE)
for index_event in range(len(catalog_PE)):
    print('Processing ',("%.2f" % round((index_event/len(catalog_PE))*100,3)), '% complete')
    print('file'+str(catalog_PE.file[index_event]))
    
    # Load MS data 
    ms_data_3c=om_load_data.load_ms_files(file1_name=catalog_PE.file[index_event],file2_name=str(catalog_PE.file_aux[index_event]), folder=dataset_folder)
    dt=ms_data_3c[0].stats.delta

    list_pp=[]
    list_sp1=[]
    list_sp2=[]
    list_ps1pp=[]
    list_ps2pp=[]

    # Run method over each geophone        
    for gph_index in range(len(catalog_gph)):
        print('Processing gph ',gph_index)
        gph_data_3c=obspy.core.stream.Stream(traces=[ms_data_3c[catalog_gph.gph_z[gph_index]], ms_data_3c[catalog_gph.gph_h1[gph_index]],ms_data_3c[catalog_gph.gph_h2[gph_index]]])

        # Slice seismic signal using second info into the PE catalog
        gph_data_3c=gph_data_3c.slice(ms_data_3c[0].stats.starttime+catalog_PE.second[index_event]-2,ms_data_3c[0].stats.starttime+catalog_PE.second[index_event]+3)
                                                 
        # Filtering data
        gph_data_3c.filter('bandpass', freqmin=30, freqmax= 500, corners=4, zerophase=True) 
        
        # Torationg to ray-centered coordinates
        gph_data_3c=om_ATP_functions.rotate_ray_center(gph_data=gph_data_3c,print_results=True)
        
        # Calculating waves onset          
        [time_p,time_s1,time_s2]=om_ATP_functions.atp_kurt_diff(gph_data_3c=gph_data_3c,gph_number=gph_index, show=False, event_name='event '+str(index_event))
      
        list_pp.append(time_p)
        list_sp1.append(time_s1)
        list_sp2.append(time_s2)
        
        # Prevent error due mispicking
        if time_p !=None and time_s1 !=None and time_s2 !=None:
            list_ps1pp.append(time_s1-time_p)
            list_ps2pp.append(time_s2-time_p)
        else:
            list_ps1pp.append(None)
            list_ps2pp.append(None)
    

    # Export results. One CSV for each file
    atp=pandas.DataFrame({'gph':numpy.arange(0,len(catalog_gph)),'pp':list_pp,'sp1':list_sp1,'sp2':list_sp2})    
    atp.to_csv(atp_folder+str(catalog_PE.time_str[index_event])+'.csv', sep=',', line_terminator='\n', index=False)

#%%