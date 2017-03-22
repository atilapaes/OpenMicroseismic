#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  4 02:42:34 2017

@author: AtilaPaes

Suppose a file with sequence ABC ABC ABC [...] where ABC is any commutation of H1, H2 and Z

Future activities:
------------------
1) Include option for 1-C data.
2) Exclude dead of noisy channels or geophones. (Should I make them zero?)
3) Deal with null file of file slice

"""

import numpy, matplotlib, math

from scipy import stats


#%% PLOT 3C data        
def plot_3c(data, time, title='Test',define_time_range=False,time_range=[1,2]):        
    """"
    Future Activities: save figure {Dpi, name}, Name the channels if different than std
    """
    for geophone in range(int(len(data)/3)): #Range of number of geophones
        matplotlib.pyplot.plot(time,data[geophone*3+0] + (geophone+1),'r') #First channel  in Red
        matplotlib.pyplot.plot(time,data[geophone*3+1] + (geophone+1),'g') #Second channel in Green
        matplotlib.pyplot.plot(time,data[geophone*3+2] + (geophone+1),'b') #Third channel in in Blue
            
    matplotlib.pyplot.xlabel('time (s)',fontweight='bold',fontsize=12)
    matplotlib.pyplot.ylabel('geophones',fontweight='bold',fontsize=12)
    matplotlib.pyplot.title(title,fontweight='bold',fontsize=14)
    if define_time_range==True:
        matplotlib.pyplot.axis(time_range + [0, int(len(data)/3)+1])
    else:
        matplotlib.pyplot.axis([time[0], time[-1], 0, int(len(data)/3)+1])      
    matplotlib.pyplot.gca().invert_yaxis() #Invert Y Axis
    matplotlib.pyplot.show()        


#%% PLOT 1C data        
def plot_1c(data, time, title='Test',define_time_range=False,time_range=[1,2]):
    """"
    Future Activities: save figure {Dpi, name}, Name the channels if different than std
    """
    for geophone in range(len(data)): #Range of number of geophones
        matplotlib.pyplot.plot(time,data[geophone] + (geophone+1),'b') 
            
    matplotlib.pyplot.xlabel('Time (s)',fontweight='bold',fontsize=12)
    matplotlib.pyplot.ylabel('Geophones',fontweight='bold',fontsize=12)
    matplotlib.pyplot.title(title,fontweight='bold',fontsize=14)
    if define_time_range==True:
        matplotlib.pyplot.axis(time_range + [0, len(data)+1])
    else:
        matplotlib.pyplot.axis([time[0], time[-1], 0, len(data)+1])      
    matplotlib.pyplot.gca().invert_yaxis() #Invert Y Axis 
    matplotlib.pyplot.show()        

#%% PLOT stack
def plot_stack(data, time, title='Test',define_time_range=False,time_range=[1,2], save=False):        
    """"
    Future Activities: save figure {Dpi, name}, Name the channels if different than std
    """
    matplotlib.pyplot.figure()
    matplotlib.pyplot.plot(time,data,'b') 
            
    matplotlib.pyplot.xlabel('Time (s)',fontweight='bold',fontsize=12)
    matplotlib.pyplot.ylabel('stack',fontweight='bold',fontsize=12)
    matplotlib.pyplot.title(title,fontweight='bold',fontsize=14)
    if define_time_range==True:
        matplotlib.pyplot.axis(time_range + [0, numpy.max(data)])
    else:
        matplotlib.pyplot.axis([time[0], time[-1], 0, numpy.max(data)])      
    
    if save == True:
        matplotlib.pyplot.savefig(title+'.png',format='png', dpi=100)   
    matplotlib.pyplot.show()


#%% sta lta
def stalta(data,sta,lta):
   
    sta_data=numpy.zeros((len(data),len(data[0])))
    lta_data=numpy.zeros((len(data),len(data[0])))    
    stalta_data=numpy.zeros((len(data),len(data[0]))) 
    
    for geophone in range(len(data)):
        
        #== lta Calculation an normalization
        for lta_index in range(lta,len(data[0])):
            lta_data[geophone][lta_index]=numpy.mean(data[geophone][lta_index-lta:lta_index])
        
        #== STA Calculation and normalization
        for sta_index in range(lta-sta,len(data[0])):
        #for sta_index in range(lta-sta,len(data[0])-sta):
            sta_data[geophone][sta_index]=numpy.mean(data[geophone][sta_index:sta_index+sta])
        
        # Defining stalta (Traditional numpy.divide may returns Nan and Inf). 
        #Defining all conditions to calculate ratio. Otherwise, attribure value Zero
        for index in range(len(data[0])):
            if lta_data[geophone][index] != 0 and not(math.isnan(lta_data[geophone][index])) and sta_data[geophone][index] != 0 and not(math.isnan(sta_data[geophone][index])):
                stalta_data[geophone][index]=sta_data[geophone][index]/lta_data[geophone][index]
            else:
                stalta_data[geophone][index]=0

        #Normalization
        stalta_data[geophone]=stalta_data[geophone]/numpy.nanmax(stalta_data[geophone]) #nanmax IGNORES NAN IN THE ARRAY
        
    return(stalta_data)
            
#%% CALCULATE CHARACTERISTICS FUNCTIONS
def calculate_function(ms_data, function_kind=1, stack=False,
                       print_plot=False, print_log = False,
                       calculate_stalta=False, sta=1000, lta=5000, 
                       save_plot=False,plot_name='Test.png',
                       define_time_range=False, time_range=[1,1],
                       sample_kurt=50):
    """
    Future activities: implement save for other kinds of plot.
    
    This is the core function of the module. It assumes the signal in each channel as velocity.
    It calculates the main characteristic functions of the microseismic signal.
    
    Parameters:
    -----------
    ms_data:         data to be used as input.
    function_kind:  Define the kind of function the user wants:
                
                    1 - Filetered Waveform.
                    2 - V^2 in each channel
                    3 - V^2 in each geophone (Vx^2+Vy^2+Vz^2) -> MOST USED FUNCTION
                    4 - Modified Energy Ratio (by Wong 2009) |Vx|^3 + |Vy|^3 + |Vz|^3 -> Highlight peaks from background
                    5 - Kurtosis of {V^2 in each geophone}. It performs neither fast or accurate for event detection.
                        sample_kurt must be even.
    stack:          Calculate the stack of signal.

    print_plot:     Plots the Characteristic function of given signal.
    print_log:      Used to follow progress and troubleshooting.
    
    calculate_stalta:   Allow the stalta calculation.
    sta:                Number of samples in sta.
    lta:                Number of samples in lta.
        
    define_time_range:  Allow definition of the plot time-range.
    time_range:         List of two elements with initial and final time.
    
    sample_kurt:     Number of samples to calcuate Kurtosis. Use at least 32 to get suitable random sample.
    
    """
#=============================================================================    
    #%% INPUTS

    #==Loading file    
    if print_log == True:
        print('Loading SEGY...')
    
    #ms_data=obspy.read(FileName)
    #ms_data = obspy.read(FileName)
    #Future option: send a slice of the file contaning just the data
#=============================================================================        
    #%% FILTERING DATA

    #==Defining time array
    time=numpy.arange(0,ms_data[0].stats.npts*ms_data[0].stats.delta,ms_data[0].stats.delta)
    
    #==Signal treatment
    if print_log == True:
        print('Loading Signal treatment...')
    
    #==Filtering
    ms_data.filter('bandpass', freqmin=30, freqmax=0.5*ms_data[0].stats.sampling_rate, corners=4, zerophase=True)
    
    #==Remove of 60 Hz
    #ms_data.filter('bandstop', freqmin=59.5, freqmax=60.5, corners=4, zerophase=True)
    
    #==Normalize signal
    ms_data.normalize()
    #Multiply by 0.5 to fit in the plot
    if function_kind==1:  #Necessary just if the function is waveform, once it goes to negative and subsequent Geophones occupy the same space on graph
        for channel in range(int(len(ms_data))):
            ms_data[channel].data=ms_data[channel].data*0.5
#=============================================================================    
    #%% FUNCTION KIND =1 NORMALIZED FILTERED WAVEFORM    
           
    #==Function Filtered Waveform
    if function_kind == 1: #Verify if Function "1" (Filtered Waveform) was requested
        
        #==Calculate signal square
        ms_data_fw=numpy.copy(ms_data)
        
        #==PLOTTING 
        if print_log == True:
            print('Plotting...')
        if print_plot == True:
            plot_3c(data=ms_data_fw, time=time, title='Filtered Waveform',define_time_range=define_time_range,time_range=time_range)
        
        #==Error report about invalid options in waveform
        if calculate_stalta==True:
            print("Sorry, it's the Waveform. stalta wasn't calculated or plotted.")
        if stack==True:
            print("Sorry, it's the Waveform. Stack wasn't calculated or plotted.")
        
        
        #== Returning waveform
        if print_log == True:
            print('Returning data...')
        return(ms_data_fw, time)
#=============================================================================
    #%% FUNCTION KIND = 2, SQUARE OF VELOCITY for each CHANNEL
           
    #==Function Filtered Waveform
    if function_kind == 2: #Verify if Function "2" (channel V^2) was requested
        
        if print_log == True:
            print('Calculating Function V^2 for each channel...')
        #==Calculate Function V channel square
        ms_data_fvc2=numpy.copy(ms_data)
        
        for channel in range(len(ms_data)):
            ms_data_fvc2[channel]=numpy.square(ms_data[channel])

        #== Calculate stalta
        if calculate_stalta==True:
            ms_data_fvc2=stalta(data=ms_data_fvc2,sta=sta,lta=lta)
        
        #== stack data?
        if stack==False:
            #==PLOTTING         
            if print_plot == True:
                plot_3c(data=ms_data_fvc2, time=time, title='V^2 in each channel',define_time_range=define_time_range,time_range=time_range)
        else:
            stacked=numpy.zeros(len(ms_data_fvc2[0]))
            for channel in range(len(ms_data)):
                stacked = numpy.add(stacked,ms_data_fvc2[channel])
            ms_data_fvc2=stacked
            
            #==PLOTTING stack        
            if print_plot == True:
                plot_stack(data=stacked, time=time, title='stack of V^2 in each channel',define_time_range=define_time_range,time_range=time_range)
            
        #== Returning Function
        if print_log == True:
            print('Returning data...')
        return(ms_data_fvc2, time)
        
    #%% FUNCTION KIND = 3, SQUARE OF VELOCITY for each Geophone  = Vx^2 + Vy^2 + Vz^3
           
    #==Function Filtered Waveform
    if function_kind == 3: #Verify if Function "2" (Geophne V^2) was requested
        
        if print_log == True:
            print('Calculating Function V^2 for each Geophone...')
        #==Calculate Function V channel square
        ms_data_fvg2=numpy.zeros((int(len(ms_data)/3),len(ms_data[0])))
        
        for geophone in range(int(len(ms_data)/3)): # Signal of Each geophone is defined by sum of squares of each channel
            ms_data_fvg2[geophone]=numpy.add(numpy.square(ms_data[geophone*3+0]),numpy.add(numpy.square(ms_data[geophone*3+1]),numpy.square(ms_data[geophone*3+2])))
            
        #== Calculate stalta
        if calculate_stalta==True:
            ms_data_fvg2=stalta(data=ms_data_fvg2,sta=sta,lta=lta)
            
         #== stack data?
        if stack==False:
            #==PLOTTING         
            if print_plot == True:
                plot_1c(data=ms_data_fvg2, time=time, title='V^2 in each Geophone',define_time_range=define_time_range,time_range=time_range)  
        
        else:
            stacked=numpy.zeros(len(ms_data_fvg2[0]))
            for channel in range(len(ms_data_fvg2)):
                stacked = numpy.add(stacked,ms_data_fvg2[channel])
            ms_data_fvg2=stacked
            
            #==PLOTTING stack        
            if print_plot == True:
                plot_stack(data=ms_data_fvg2, time=time, title='stack of V^2 in each channel',define_time_range=define_time_range,time_range=time_range)   
            
        #== Returning Function
        if print_log == True:
            print('Returning data...')
        return(ms_data_fvg2, time)


        
    #%% FUNCTION KIND = 4, Cube OF ABSOLUTE VELOCITY for each Geophone  = |Vx|^3 + |Vy|^3 + |Vz|^3
           
    #==Function Filtered Waveform
    if function_kind == 4: #Verify if Function "4" (Geophne |V|^3) was requested
        
        if print_log == True:
            print('Calculating Function |V|^3 for each Geophone...')
        #==Calculate Function V channel square
        ms_data_fvg3=numpy.zeros((int(len(ms_data)/3),len(ms_data[0])))
        
        for geophone in range(int(len(ms_data)/3)): # Signal of Each geophone is defined by sum of squares of each channel
            ms_data_fvg3[geophone]=numpy.add(numpy.power(numpy.absolute(ms_data[geophone*3+0]),3),numpy.add(numpy.power(numpy.absolute(ms_data[geophone*3+1]),3),numpy.power(numpy.absolute(ms_data[geophone*3+2]),3)))
            
        #== Calculate stalta
        if calculate_stalta==True:
            ms_data_fvg3=stalta(data=ms_data_fvg3,sta=sta,lta=lta)
            
         #== stack data?
        if stack==False:
            #==PLOTTING         
            if print_plot == True:
                plot_1c(data=ms_data_fvg3, time=time, title='|V|^3 in each Geophone',define_time_range=define_time_range,time_range=time_range)  
        
        else:
            stacked=numpy.zeros(len(ms_data_fvg3[0]))
            for channel in range(len(ms_data_fvg3)):
                stacked = numpy.add(stacked,ms_data_fvg3[channel])
            ms_data_fvg3=stacked
            
            #==PLOTTING stack        
            if print_plot == True:
                plot_stack(data=ms_data_fvg3, time=time, title='stack of |V|^3 in each channel',define_time_range=define_time_range,time_range=time_range)   
            
        #== Returning Function
        if print_log == True:
            print('Returning data...')
        return(ms_data_fvg3, time)
################################################################################

#%% FUNCTION KIND = 5, KURTOSIS OF THE SQUARE VELOCITY OF EACH GEOPHONE -> KURT(V^2)
    
    #==Function of Kurtosis of V^2 through  geophones
    if function_kind == 5:     
        if print_log == True:
            print('Calculating Function Kurtosis for each Geophones...')
        #==Calculate Function V^2 in each geophone, then, V^3

        #==Initialization of Variables
        ms_data_kurt=numpy.zeros((int(len(ms_data)/3),len(ms_data[0]))) #Initializing the array for Kurtosis
        #[ms_data_fvg2,time]=Charac_Functions.calculate_function(ms_data,function_kind=3) #Can't calculte this simple until modify input mehotd in Charac_Functions 
                
        #==Calculate Function V channel square
        ms_data_fvg2=numpy.zeros((int(len(ms_data)/3),len(ms_data[0])))
        
        for geophone in range(int(len(ms_data)/3)): # Signal of Each geophone is defined by sum of squares of each channel
            ms_data_fvg2[geophone]=numpy.add(numpy.square(ms_data[geophone*3+0]),numpy.add(numpy.square(ms_data[geophone*3+1]),numpy.square(ms_data[geophone*3+2])))
             
        for geophone in range(int(len(ms_data)/3)): # Signal of Each geophone is defined by sum of squares of each channel
            for index in range(int(sample_kurt/2),len(ms_data[0])-int(sample_kurt/2)):
                ms_data_kurt[geophone][index]=stats.kurtosis(ms_data_fvg2[geophone][index-int(sample_kurt/2):index+int(sample_kurt/2)], axis=0, fisher=True, bias=True)
            max_kurt=ms_data_kurt[geophone].max()
            ms_data_kurt[geophone]=ms_data_kurt[geophone]/max_kurt
        
       #== Calculate stalta
        if calculate_stalta==True:
            ms_data_kurt=stalta(data=ms_data_kurt,sta=sta,lta=lta)
            
         #== stack data?
        if stack==False:
            #==PLOTTING         
            if print_plot == True:
                plot_1c(data=ms_data_kurt, time=time, title='Kurtosis in each Geophone',define_time_range=define_time_range,time_range=time_range)  
        
        else:
            stacked=numpy.zeros(len(ms_data_kurt[0]))
            for channel in range(len(ms_data_kurt)):
                stacked = numpy.add(stacked,ms_data_kurt[channel])
            ms_data_kurt=stacked
            
            #==PLOTTING stack        
            if print_plot == True:
                plot_stack(data=ms_data_kurt, time=time, title='stack of Kurtosis in each Geophone',define_time_range=define_time_range,time_range=time_range)   
                    
        #== Returning Function
        if print_log == True:
            print('Returning data...')
        return(ms_data_kurt, time)
################################################################################
#%%
#REVIEWED
def list_extensions(folder):
    """
    This function list all the (segy, sgy, dat) files into the speciefied folder.
    
    Parameters
    ----------
    folder : Name of the folder that contains the files to be listed.
    
    Returns
    -------
    list_of_files : List of files with the extensions listed in parameters file.

    Examples
    --------
    """
    import om_ped_es_parameters, os
    if om_ped_es_parameters.verbose_level <= 1:
        print('Listing MS files...')

    ### Entering in the folder contaning the MS Files
    os.chdir(folder)
    
    ## List all the files in the directory
    list_of_files=os.listdir()
    
    ## Initialing the list of SEGYs
    list_of_ms_files=[]
    
    ## Loop to identify the SEGYs and sgy 
    for index in range(len(list_of_files)):
        if list_of_files[index].endswith("."+"segy") or list_of_files[index].endswith("."+"sgy") or list_of_files[index].endswith("."+"dat"):
            list_of_ms_files.append(list_of_files[index])
    
    ## Debug: verify if the list is empty
    if len(list_of_ms_files)==0:
        print('The list of files is empty. The current compatible formats are: segy, sgy and dat')
            
    ## Leaving the folder contaning the MS Files
    os.chdir('..')
    if om_ped_es_parameters.verbose_level <= 1:
        print('Listing MS files -> Done!')
    return (list_of_ms_files)
################################################################################

#%% Math functions to manipulate data #########################################
""" REVIEWED"""   
def moving_avg(signal,samples):
    """
    This function calculates the moving average of a provided 1-C signal (array).
    
    Parameters
    ----------
    folder : Name of the folder that contains the files to be listed.
    
    Returns
    -------
    list_of_files : List of files with the extensions listed in parameters file.

    Examples
    --------
    
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
    signal_avg = numpy.average(signal_ma[samples//2:len(signal)-samples//2])
    signal_ma[0:samples//2] = signal_avg
    signal_ma[len(signal)-samples//2:len(signal)] = signal_avg 
        
    return (signal_ma)   
################################################################################

#%%% ##################################################################################################
def resume_array(original_array, first_sort, second_sort):
    import om_ped_es_parameters
    """
    This function get an array with repeated values in first_sort column and return 
    one element of same first_sort column with maximum value in second_sort column
    
    Ex: 
         [(16, 2), (16, 2), (16, 3), (17, 2), (18, 4), (18, 5), (18, 1)]
 
        1) Identifing the different single values of first_sort
            [16, 17, 18]


        2) Scan first array and look for elements with first_sort value
            Make the blocked array with blocks of same first_sort
            
            [(16, 3), (16, 2), (16, 2)]
            [(17, 2)]
            [(18, 5), (18, 4), (18, 1)]
           
            
        3) Sort the blocked array and get the ones with maximum second_sort parameters
            [(16, 3), (17, 2), (18, 5)]
            
    """
    #Identifing the different single values of first_sort
    block_values=[]
    block_values.append(original_array[0][first_sort])
    for count in range(1,(len(original_array)-1)):
        if original_array[count -1][first_sort] != original_array[count][first_sort]:
            block_values.append(original_array[count][first_sort])        
            
    resumed_array=[]
    for count_b in range(len(block_values)):
        blocked_array=[]
        # Scan first array and look for elements with first_sort value
        for count_a in range(len(original_array)):
            #Make the blocked array with blocks of same first_sort
            if original_array[count_a][first_sort] == block_values[count_b]: 
                blocked_array.append(original_array[count_a][:])
        # Sort the blocked array
    
        blocked_array=sorted(blocked_array, key=lambda blocked_array:blocked_array[:][second_sort], reverse=True)
        #print(blocked_array)
        resumed_array.append(blocked_array[0][:])
        
    if om_ped_es_parameters.verbose_level <= 1: #For debugging
            print('Identified peaks after QC: ',resumed_array) 
    
    return(resumed_array)          
    
##################################################################################################

#%% Peak evaluation
#Reviewed
def peak_evaluation(signal, peaks_positions, time_torrent_index):
    import numpy, om_ped_es_parameters, om_general_signal_processing
    """
    #Getting the left and right coordinates where the peak borders cross the limit condition
    signal:           1-C Vector with peaks.
    peaks_positions:  Vector contaning the index of maximum amplitude of peaks in the signal   
                        (Argument of "signal" that results in peak maximum)
    
    #Output: peaks properties of the signal
    
    Inner variables:
    threshold:  Signal level where the function should cross to trigger the width border
                Ex: numpy.mean(V) --> it really makes sense in this context
    """
    threshold=numpy.mean(signal)
    peaks_properties=numpy.zeros((len(peaks_positions),4))
    
    
    # Eliminationg the peaks too close from the borders
    if om_ped_es_parameters.verbose_level == 0:
        print('Peaks before border cutting: ',peaks_positions)
    
    peaks_positions_sliced=[]
    for index in range(len(peaks_positions)):
        if (peaks_positions[index] >= time_torrent_index) and (peaks_positions[index] <= (len(signal)-time_torrent_index)):
            peaks_positions_sliced.append(peaks_positions[index])
    
    if om_ped_es_parameters.verbose_level == 0:
        print('Peaks after border cutting: ',peaks_positions_sliced)


    peaks_properties_qc1=[]
    for peak_counter in range(len(peaks_positions_sliced)):
        
        # Bait for Peaks Quality-Control. It garantees all peaks identified have all parameters calculated correctly.
        right_idx=left_idx=0
        
        #Evaluation the peak position and width
        for index in range(peaks_positions_sliced[peak_counter],len(signal)):
            if signal[index]<threshold:
                right_idx=index
                break
        for index in range(peaks_positions_sliced[peak_counter],0,-1):
            if signal[index]<threshold:
                left_idx=index
                break        

        #Registering the properties of each peak
        peaks_properties[peak_counter][0]=left_idx                           #Arrival index
        peaks_properties[peak_counter][1]=(signal[peaks_positions_sliced[peak_counter]]/numpy.mean(signal))  #SNR
        #print('Peak level', signal[peaks_positions_sliced[peak_counter]])
        #print('Peak mean', numpy.mean(signal))
        peaks_properties[peak_counter][2]=right_idx-left_idx                 #Width of Peak
        peaks_properties[peak_counter][3]=right_idx                         #Right index
        
        
        # Quality control 1, baits-based: Use peaks with proper calculated info
        if (peaks_properties[peak_counter][0] !=0) or (peaks_properties[peak_counter][3] !=0):
            peaks_properties_qc1.append([peaks_properties[peak_counter][0], peaks_properties[peak_counter][1], peaks_properties[peak_counter][2],peaks_properties[peak_counter][3]])
    #################################
    
    # Quality control output
    if om_ped_es_parameters.verbose_level == 0:
            print("Peak properties after QC1. Left index, SNR, right index and width: ",peaks_properties)
    
            
    #At this point, all info is calculated. Now we have to make sure to exclude repetead peaks
    # Ideal situation: peak with "^"shape. Regular situation:  peak with "M" shape.
    #The method used is to compare the arrival index. If they are the same, keep the peak with maximum SNR
    
    peaks_properties_qc2=[] #Variable to be returned
    
    if len(peaks_properties_qc1)!=0: # Verifing if any peak last after quality control
        # Lets verify if any of the arraival times are equal:
            identified_peaks=om_general_signal_processing.resume_array(original_array=peaks_properties_qc1, first_sort=0, second_sort=1)
            
            #>>>>> peaks properties after resume array
            if om_ped_es_parameters.verbose_level == 0:
                print('Peaks properties after resume array',identified_peaks)
            
            # Second Quality control for the events based in peaks width and SNR
            for index in range(len(identified_peaks)):
                if identified_peaks[index][2] >= om_ped_es_parameters.peak_width_minimum and peaks_properties[index][1]>=om_ped_es_parameters.peak_snr_minimum :
                    peaks_properties_qc2.append(peaks_properties[index][:])
                
    if om_ped_es_parameters.verbose_level <= 1:
        print("Returned values from peak evaluation ",peaks_properties_qc2)               
    
    return(peaks_properties_qc2) 
################################################################################

#%%
#REVIEWED
def detect_peaks(x, mph=None, mpd=1, threshold=0, edge='rising',
                 kpsh=False, valley=False, show=False, ax=None):

    import numpy

    """Detect peaks in data based on their amplitude and other features.
    
    
    Detect peaks in data based on their amplitude and other features.
    http://nbviewer.jupyter.org/github/demotu/BMC/blob/master/notebooks/DetectPeaks.ipynb

    Atila's example
    u=detect_peaks.detect_peaks(data1, mph=numpy.mean(data1)+1.5*numpy.std(data1), mpd=100, show=True)
    
    
    Parameters
    ----------
    x : 1D array_like
        data.
    mph : {None, number}, optional (default = None)
        detect peaks that are greater than minimum peak height.
    mpd : positive integer, optional (default = 1)
        detect peaks that are at least separated by minimum peak distance (in
        number of data).
    threshold : positive number, optional (default = 0)
        detect peaks (valleys) that are greater (smaller) than `threshold`
        in relation to their immediate neighbors.
    edge : {None, 'rising', 'falling', 'both'}, optional (default = 'rising')
        for a flat peak, keep only the rising edge ('rising'), only the
        falling edge ('falling'), both edges ('both'), or don't detect a
        flat peak (None).
    kpsh : bool, optional (default = False)
        keep peaks with same height even if they are closer than `mpd`.
    valley : bool, optional (default = False)
        if True (1), detect valleys (local minima) instead of peaks.
    show : bool, optional (default = False)
        if True (1), plot data in matplotlib figure.
    ax : a matplotlib.axes.Axes instance, optional (default = None).

    Returns
    -------
    ind : 1D array_like
        indeces of the peaks in `x`.

    Notes
    -----
    The detection of valleys instead of peaks is performed internally by simply
    negating the data: `ind_valleys = detect_peaks(-x)`
    
    The function can handle NaN's 

    See this IPython Notebook [1]_.

    References
    ----------
    .. [1] http://nbviewer.ipython.org/github/demotu/BMC/blob/master/notebooks/DetectPeaks.ipynb

    """

    x = numpy.atleast_1d(x).astype('float64')
    if x.size < 3:
        return numpy.array([], dtype=int)
    if valley:
        x = -x
    # find indices of all peaks
    dx = x[1:] - x[:-1]
    # handle NaN's
    indnan = numpy.where(numpy.isnan(x))[0]
    if indnan.size:
        x[indnan] = numpy.inf
        dx[numpy.where(numpy.isnan(dx))[0]] = numpy.inf
    ine, ire, ife = numpy.array([[], [], []], dtype=int)
    if not edge:
        ine = numpy.where((numpy.hstack((dx, 0)) < 0) & (numpy.hstack((0, dx)) > 0))[0]
    else:
        if edge.lower() in ['rising', 'both']:
            ire = numpy.where((numpy.hstack((dx, 0)) <= 0) & (numpy.hstack((0, dx)) > 0))[0]
        if edge.lower() in ['falling', 'both']:
            ife = numpy.where((numpy.hstack((dx, 0)) < 0) & (numpy.hstack((0, dx)) >= 0))[0]
    ind = numpy.unique(numpy.hstack((ine, ire, ife)))
    # handle NaN's
    if ind.size and indnan.size:
        # NaN's and values close to NaN's cannot be peaks
        ind = ind[numpy.in1d(ind, numpy.unique(numpy.hstack((indnan, indnan-1, indnan+1))), invert=True)]
    # first and last values of x cannot be peaks
    if ind.size and ind[0] == 0:
        ind = ind[1:]
    if ind.size and ind[-1] == x.size-1:
        ind = ind[:-1]
    # remove peaks < minimum peak height
    if ind.size and mph is not None:
        ind = ind[x[ind] >= mph]
    # remove peaks - neighbors < threshold
    if ind.size and threshold > 0:
        dx = numpy.min(numpy.vstack([x[ind]-x[ind-1], x[ind]-x[ind+1]]), axis=0)
        ind = numpy.delete(ind, numpy.where(dx < threshold)[0])
    # detect small peaks closer than minimum peak distance
    if ind.size and mpd > 1:
        ind = ind[numpy.argsort(x[ind])][::-1]  # sort ind by peak height
        idel = numpy.zeros(ind.size, dtype=bool)
        for i in range(ind.size):
            if not idel[i]:
                # keep peaks with the same height if kpsh is True
                idel = idel | (ind >= ind[i] - mpd) & (ind <= ind[i] + mpd) \
                    & (x[ind[i]] > x[ind] if kpsh else True)
                idel[i] = 0  # Keep current peak
        # remove the small peaks and sort back the indices by their occurrence
        ind = numpy.sort(ind[~idel])

    if show:
        if indnan.size:
            x[indnan] = numpy.nan
        if valley:
            x = -x
        _plot(x, mph, mpd, threshold, edge, valley, ax, ind)

    return ind


def _plot(x, mph, mpd, threshold, edge, valley, ax, ind):
    """Plot results of the detect_peaks function, see its help."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print('matplotlib is not available.')
    else:
        if ax is None:
            _, ax = plt.subplots(1, 1, figsize=(8, 4))

        ax.plot(x, 'b', lw=1)
        if ind.size:
            label = 'valley' if valley else 'peak'
            label = label + 's' if ind.size > 1 else label
            ax.plot(ind, x[ind], '+', mfc=None, mec='r', mew=2, ms=8,
                    label='%d %s' % (ind.size, label))
            ax.legend(loc='best', framealpha=.5, numpoints=1)
        ax.set_xlim(-.02*x.size, x.size*1.02-1)
        ymin, ymax = x[numpy.isfinite(x)].min(), x[numpy.isfinite(x)].max()
        yrange = ymax - ymin if ymax > ymin else 1
        ax.set_ylim(ymin - 0.1*yrange, ymax + 0.1*yrange)
        ax.set_xlabel('data #', fontsize=14)
        ax.set_ylabel('Amplitude', fontsize=14)
        mode = 'Valley detection' if valley else 'Peak detection'
        ax.set_title("%s (mph=%s, mpd=%d, threshold=%s, edge='%s')"
                     % (mode, str(mph), mpd, str(threshold), edge))
        # plt.grid()
        plt.show()

 #%%