#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 18 16:46:20 2017

@author: atilapaes
"""

#%%
import numpy
from scipy import stats

#%%
def cft_energy_stack_v2(stream):
    """    
    This function calculates the energy stack of the signal for N-dimensions streams
    
    """
    # Calculating Characteristic Function, SQUARE OF VELOCITY for each CHANNEL
    import numpy
    output=numpy.zeros(len(stream[0]))   
    for index in range(len(stream)):       
        output=numpy.add(output,numpy.square(stream[index].data))
    return(output)

#%%
def energy_stack_selected_gph(stream,valid_gph):
    """    
    This function calculates the energy-stack based in a list of valid gph
    
    """
    # Calculating Characteristic Function, SQUARE OF VELOCITY for each CHANNEL
    import numpy
    output=numpy.zeros(len(stream[0]))   
    
    for index in range(int(len(stream)/3)):       
        if valid_gph[index]==True:
            #print('============',valid_gph[index])
            output=numpy.add(output,numpy.square(stream[index].data))
            output=numpy.add(output,numpy.square(stream[index+69].data))
            output=numpy.add(output,numpy.square(stream[index+2*69].data))
    return(output)
#%%
def energy_stack_selected_gph_v2(stream,gph_catalog):
    """    
    This function calculates the energy-stack based in a list of valid gph and gph channels
    
    """
    # Calculating Characteristic Function, SQUARE OF VELOCITY for each CHANNEL
    import numpy
    output=numpy.zeros(len(stream[0]))   
    
    for index in range(int(len(stream)/3)):       
        if gph_catalog.valid_gph[index]==True:
            #print('============',valid_gph[index])
            output=numpy.add(output,numpy.square(stream[gph_catalog.gph_z[index]].data))
            output=numpy.add(output,numpy.square(stream[gph_catalog.gph_h1[index]].data))
            output=numpy.add(output,numpy.square(stream[gph_catalog.gph_h2[index]].data))
    return(output)

#%%
def cft_modified_energy(ms_data):    
    """    
    This function calculates the energy stack of the signal for 1C or 3C streams
    
    Parameters:
    -----------
    ms_data:        1c or 3C stream of data
    mavg_samples:   number of samples to use in the moving average    
    """
    # Calculating Characteristic Function, SQUARE OF VELOCITY for each CHANNEL
    if len(ms_data)==1:
        ms_data[0].data=numpy.square(ms_data[0].data)
    elif len(ms_data)==3:        
        ms_data[0].data=numpy.add(numpy.power(numpy.absolute(ms_data[0]),3),numpy.add(numpy.power(numpy.absolute(ms_data[1]),3),numpy.power(numpy.absolute(ms_data[2]),3)))

    return(ms_data[0].data)

#%% CALCULATE Kurtosis function
def cft_kurtosis_old(signal, kurt_samples):    
    #==Initialization of Variables
    signal_kurt=numpy.zeros(len(signal))
    
    # Calculation of kurtosis
    for index in range(kurt_samples,len(signal)):
        signal_kurt[index]=stats.kurtosis(signal[index:index+kurt_samples], axis=0, fisher=True, bias=True)
    signal_kurt=signal_kurt/signal_kurt.max()
    return(signal_kurt)
#%%
def cft_kurtosis_center(signal, kurt_samples):    
    #==Initialization of Variables
    signal_kurt=numpy.zeros(len(signal))
    
    # Calculation of kurtosis
    for index in range(int(kurt_samples/2),len(signal)-int(kurt_samples/2)):
        signal_kurt[index]=stats.kurtosis(signal[index-int(kurt_samples/2):index+int(kurt_samples/2)], axis=0, fisher=True, bias=True)
    signal_kurt=signal_kurt/signal_kurt.max()
    return(signal_kurt)

#%% CALCULATE Kurtosis function
def cft_kurtosis(signal, kurt_samples):    
    #==Initialization of Variables
    signal_kurt=numpy.zeros(len(signal))
    
    # Calculation of kurtosis
    for index in range(kurt_samples-1,len(signal)):
        signal_kurt[index]=stats.kurtosis(signal[index-kurt_samples+1:index], axis=0, fisher=True, bias=True)
    return(signal_kurt)
#%% Auxiliar functions
def aux_moving_avg(signal,samples):
    """
    This function calculates the moving average of a provided 1-C signal (array).
    
    signal_ma = numpy.zeros((len(signal),))
    for index in range(len(signal)):
         signal_ma[index] = numpy.sum(signal[index:(index+samples)])
    return (signal_ma/samples)
    
    """
    signal_ma = numpy.zeros((len(signal),))
    
    #Regular signal
    for index in range(samples//2, len(signal)-samples//2):
        signal_ma[index] = (numpy.sum(signal[index-samples//2:(index+samples//2)]))/samples
         
    #Borders of the signal
    
    signal_ma[0:samples//2] = 0
    signal_ma[len(signal)-samples//2:len(signal)] = 0
    signal_ma=signal_ma/signal_ma.max()
    
    return (signal_ma)   

#%%
def aux_moving_avg_arriving(signal,mavg_samples):
    """
    This function calculates the moving average of a provided 1-C signal (array).    
    """
    signal_ma = numpy.zeros((len(signal),))
    
    #Regular signal
    for index in range(mavg_samples, len(signal)-mavg_samples):
        signal_ma[index] = (numpy.sum(signal[index:(index+mavg_samples)]))/mavg_samples

    signal_mean=numpy.mean(signal_ma[mavg_samples:len(signal)-mavg_samples])
    
    # Signal borders
    for index in range(0,mavg_samples):
        signal_ma[index] = signal_mean

    for index in range(len(signal)-mavg_samples, len(signal)):
        signal_ma[index] = signal_mean
       
    return (signal_ma/signal_ma.max())   

#%% new function Aug 23, 2017
def diff_kurt(stream,kurt_samples):
    import numpy
    output = numpy.zeros((len(stream),len(stream[0]))) #utput with same dimention of input data

    for channel in range(3):    
        stream[channel].data=stream[channel].data/stream[channel].data.max()
        kurt=cft_kurtosis(stream[channel].data, kurt_samples)
        kurt=kurt/kurt.max()

        kurt_diff=numpy.diff(kurt)
        kurt_diff=numpy.append([0],kurt_diff) # appending [0] in the left
        kurt_diff=kurt_diff/kurt_diff.max()
        
        for index in range(len(kurt_diff)):                
            if kurt_diff[index] >0:
                output[channel][index]=kurt_diff[index]
                
    return(output)

#%%