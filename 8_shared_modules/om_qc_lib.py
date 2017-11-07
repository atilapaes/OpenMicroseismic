#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 14:03:42 2017

@author: atilapaes
Functions for event identification quality control

"""

import sys
sys.path.insert(0,'/Users/atilapaes/Documents/Python_projects/toc2me/shared_lib/')


def filter_dt(tp,ts, show_plot, lower_delta,min_len):
    """
    This function filter the ts - tp vector by removing outliers
    tp:         P- picking array
    ts:         S- picking array
    show_plot:  show the of the filtered delta_t and basic statistics
    lower_delta: minimum value of delta to be considered as valid picking
    
    """
    import numpy
    import matplotlib.pyplot as plt
    ####################################
    def array_none_removal(array_dt):
        output=[]
        for index in range(len(array_dt)):
            if array_dt[index]!=None:
                output.append(array_dt[index])
        return(output)
    
    ####################################
    # Calculate delta t based in P and S1 arrival
    delta_t=[]
    for gph_index in range(len(tp)):
        # verify is the pp and sp1 are NOT nan
        if tp[gph_index]!=None and ts[gph_index]!=None:
        #if (not math.isnan(tp[gph_index])) and (not math.isnan(ts[gph_index])):
            delta_t.append(ts[gph_index]-tp[gph_index])
        else:
            delta_t.append(None)

    ####################################
    #%% Filter delta_t
    
    # Remove Nones from delta t to calculate its stats
    dt_numbers=array_none_removal(delta_t)
    dt_numbers=list(filter(lambda x:x>0.25,  dt_numbers)) # remove values lower than 0.2
         
    # Build the new delta_t array with outliers removed
    dt_filtered=[]
    
    for dt_index in range(len(delta_t)):
        if delta_t[dt_index]!= None:
            if delta_t[dt_index] <= (numpy.median(dt_numbers)+numpy.std(dt_numbers)) and delta_t[dt_index] >= (numpy.median(dt_numbers)-numpy.std(dt_numbers)) and delta_t[dt_index]>=lower_delta:
                dt_filtered.append(delta_t[dt_index])
            else:
                dt_filtered.append(None)
        else:
            dt_filtered.append(None)
        
    dt_filtered_numbers=array_none_removal(dt_filtered)
    #print('Len filtered ',len(dt_filtered_numbers))
    
    ####################################
    if show_plot==True:# and len(dt_filtered_numbers)>=min_len:    
        plt.figure(1, figsize=(8,1))
        plt.plot(dt_filtered,'.k')
        plt.axhline(y=numpy.median(dt_numbers),c='g', linestyle='--')#, lw=0.5)
        plt.axhline(y=numpy.median(dt_numbers)+numpy.std(dt_numbers),c='r', linestyle=':')#, lw=0.5)
        plt.axhline(y=numpy.median(dt_numbers)-numpy.std(dt_numbers),c='r', linestyle=':')#,lw=0.5)
        plt.ylim(0.5,1.5)
        plt.xlabel('Gph id')
        plt.ylabel('ts -tp (s)')
        plt.title('ATP ts-tp')    
        plt.show()
        plt.close()
    ####################################    
    return(dt_filtered,dt_filtered_numbers)