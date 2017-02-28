#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  4 02:42:34 2017

@author: AtilaPaes

Suppose a file with sequence ABC ABC ABC [...] where ABC is any commutation of H1, H2 and Z

Future activities:
------------------
1) Calculate functions and send data to plot in case of 1-C data.
2) Exclude channels or geophones. Make them zero.


"""

import numpy, matplotlib, math

from scipy import stats


#%% PLOT 3C DATA        
def Plot3C(Data, Time, Title='Test',Define_Time_Range=False,Time_Range=[1,2]):        
    """"
    Future Activities: save figure {Dpi, name}, Name the Channels if different than std
    """
    for geophone in range(int(len(Data)/3)): #Range of number of geophones
        matplotlib.pyplot.plot(Time,Data[geophone*3+0] + (geophone+1),'r') #First channel  in Red
        matplotlib.pyplot.plot(Time,Data[geophone*3+1] + (geophone+1),'g') #Second channel in Green
        matplotlib.pyplot.plot(Time,Data[geophone*3+2] + (geophone+1),'b') #Third channel in in Blue
            
    matplotlib.pyplot.xlabel('Time (s)',fontweight='bold',fontsize=12)
    matplotlib.pyplot.ylabel('Geophones',fontweight='bold',fontsize=12)
    matplotlib.pyplot.title(Title,fontweight='bold',fontsize=14)
    if Define_Time_Range==True:
        matplotlib.pyplot.axis(Time_Range + [0, int(len(Data)/3)+1])
    else:
        matplotlib.pyplot.axis([Time[0], Time[-1], 0, int(len(Data)/3)+1])      
    matplotlib.pyplot.gca().invert_yaxis() #Invert Y Axis
    matplotlib.pyplot.show()        


#%% PLOT 1C DATA        
def Plot1C(Data, Time, Title='Test',Define_Time_Range=False,Time_Range=[1,2]):        
    """"
    Future Activities: save figure {Dpi, name}, Name the Channels if different than std
    """
    for geophone in range(len(Data)): #Range of number of geophones
        matplotlib.pyplot.plot(Time,Data[geophone] + (geophone+1),'b') 
            
    matplotlib.pyplot.xlabel('Time (s)',fontweight='bold',fontsize=12)
    matplotlib.pyplot.ylabel('Geophones',fontweight='bold',fontsize=12)
    matplotlib.pyplot.title(Title,fontweight='bold',fontsize=14)
    if Define_Time_Range==True:
        matplotlib.pyplot.axis(Time_Range + [0, len(Data)+1])
    else:
        matplotlib.pyplot.axis([Time[0], Time[-1], 0, len(Data)+1])      
    matplotlib.pyplot.gca().invert_yaxis() #Invert Y Axis
    matplotlib.pyplot.show()        

#%% PLOT stack
def Plotstack(Data, Time, Title='Test',Define_Time_Range=False,Time_Range=[1,2]):        
    """"
    Future Activities: save figure {Dpi, name}, Name the Channels if different than std
    """
    matplotlib.pyplot.plot(Time,Data,'b') 
            
    matplotlib.pyplot.xlabel('Time (s)',fontweight='bold',fontsize=12)
    matplotlib.pyplot.ylabel('stack',fontweight='bold',fontsize=12)
    matplotlib.pyplot.title(Title,fontweight='bold',fontsize=14)
    if Define_Time_Range==True:
        matplotlib.pyplot.axis(Time_Range + [0, numpy.max(Data)])
    else:
        matplotlib.pyplot.axis([Time[0], Time[-1], 0, numpy.max(Data)])      
    matplotlib.pyplot.show()        


#%% STA LTA
def STALTA(Data,STA,LTA):
   
    STAData=numpy.zeros((len(Data),len(Data[0])))
    LTAData=numpy.zeros((len(Data),len(Data[0])))    
    STALTAData=numpy.zeros((len(Data),len(Data[0]))) 
    
    for Geophone in range(len(Data)):
        
        #== LTA Calculation an normalization
        for LTAIndex in range(LTA,len(Data[0])):
            LTAData[Geophone][LTAIndex]=numpy.mean(Data[Geophone][LTAIndex-LTA:LTAIndex])
        
        #== STA Calculation and normalization
        for STAIndex in range(LTA-STA,len(Data[0])):
        #for STAIndex in range(LTA-STA,len(Data[0])-STA):
            STAData[Geophone][STAIndex]=numpy.mean(Data[Geophone][STAIndex:STAIndex+STA])
        
        # Defining STALTA (Traditional numpy.divide may returns Nan and Inf). 
        #Defining all conditions to calculate ratio. Otherwise, attribure value Zero
        for index in range(len(Data[0])):
            if LTAData[Geophone][index] != 0 and not(math.isnan(LTAData[Geophone][index])) and STAData[Geophone][index] != 0 and not(math.isnan(STAData[Geophone][index])):
                STALTAData[Geophone][index]=STAData[Geophone][index]/LTAData[Geophone][index]
            else:
                STALTAData[Geophone][index]=0

        #Normalization
        STALTAData[Geophone]=STALTAData[Geophone]/numpy.nanmax(STALTAData[Geophone]) #nanmax IGNORES NAN IN THE ARRAY
        
    return(STALTAData)
            
#%% CALCULATE CHARACTERISTICS FUNCTIONS
def calculate_function(ms_data, function_kind=1, stack=False,
                       print_plot=False, print_log = False,
                       Calculate_STALTA=False, STA=1000, LTA=5000, 
                       Define_Time_Range=False, Time_Range=[1,1],
                       SampleKurt=50):
    """
    This is the core function of the module. It assumes the signal in each channel as velocity.
    It calculates the main characteristic functions of the microseismic signal.
    
    Parameters:
    -----------
    ms_data:         Data to be used as input.
    function_kind:  Define the kind of function the user wants:
                
                    1 - Filetered Waveform.
                    2 - V^2 in each channel
                    3 - V^2 in each geophone (Vx^2+Vy^2+Vz^2) -> MOST USED FUNCTION
                    4 - Modified Energy Ratio (by Wong 2009) |Vx|^3 + |Vy|^3 + |Vz|^3 -> Highlight peaks from background
                    5 - Kurtosis of {V^2 in each geophone}. It performs neither fast or accurate for event detection.
                        SampleKurt must be even.
    stack:          Calculate the stack of signal.

    print_plot:     Plots the Characteristic function of given signal.
    print_log:      Used to follow progress and troubleshooting.
    
    Calculate_STALTA:   Allow the STALTA calculation.
    STA:                Number of samples in STA.
    LTA:                Number of samples in LTA.
        
    Define_Time_Range:  Allow definition of the plot time-range.
    Time_Range:         List of two elements with initial and final time.
    
    SampleKurt:     Number of samples to calcuate Kurtosis. Use at least 32 to get suitable random sample.
    
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
    Time=numpy.arange(0,ms_data[0].stats.npts*ms_data[0].stats.delta,ms_data[0].stats.delta)
    
    #==Signal treatment
    if print_log == True:
        print('Loading Signal treatment...')
    
    #==Filtering
    ms_data.filter('bandpass', freqmin=30, freqmax=0.5*ms_data[0].stats.sampling_rate, corners=4, zerophase=True)
    
    #==Remove of 60 Hz
    ms_data.filter('bandstop', freqmin=59.5, freqmax=60.5, corners=4, zerophase=True)
    
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
        ms_dataFW=numpy.copy(ms_data)
        
        #==PLOTTING 
        if print_log == True:
            print('Plotting...')
        if print_plot == True:
            Plot3C(Data=ms_dataFW, Time=Time, Title='Filtered Waveform',Define_Time_Range=Define_Time_Range,Time_Range=Time_Range)
        
        #==Error report about invalid options in waveform
        if Calculate_STALTA==True:
            print("Sorry, it's the Waveform. STALTA wasn't calculated or plotted.")
        if stack==True:
            print("Sorry, it's the Waveform. Stack wasn't calculated or plotted.")
        
        
        #== Returning waveform
        if print_log == True:
            print('Returning data...')
        return(ms_dataFW, Time)
#=============================================================================
    #%% FUNCTION KIND = 2, SQUARE OF VELOCITY for each CHANNEL
           
    #==Function Filtered Waveform
    if function_kind == 2: #Verify if Function "2" (Channel V^2) was requested
        
        if print_log == True:
            print('Calculating Function V^2 for each channel...')
        #==Calculate Function V Channel square
        ms_dataFVC2=numpy.copy(ms_data)
        
        for Channel in range(len(ms_data)):
            ms_dataFVC2[Channel]=numpy.square(ms_data[Channel])

        #== Calculate STALTA
        if Calculate_STALTA==True:
            ms_dataFVC2=STALTA(Data=ms_dataFVC2,STA=STA,LTA=LTA)
        
        #== stack data?
        if stack==False:
            #==PLOTTING         
            if print_plot == True:
                Plot3C(Data=ms_dataFVC2, Time=Time, Title='V^2 in each channel',Define_Time_Range=Define_Time_Range,Time_Range=Time_Range)
        else:
            stacked=numpy.zeros(len(ms_dataFVC2[0]))
            for Channel in range(len(ms_data)):
                stacked = numpy.add(stacked,ms_dataFVC2[Channel])
            ms_dataFVC2=stacked
            
            #==PLOTTING stack        
            if print_plot == True:
                Plotstack(Data=stacked, Time=Time, Title='stack of V^2 in each channel',Define_Time_Range=Define_Time_Range,Time_Range=Time_Range)
            
        #== Returning Function
        if print_log == True:
            print('Returning data...')
        return(ms_dataFVC2, Time)
        
    #%% FUNCTION KIND = 3, SQUARE OF VELOCITY for each Geophone  = Vx^2 + Vy^2 + Vz^3
           
    #==Function Filtered Waveform
    if function_kind == 3: #Verify if Function "2" (Geophne V^2) was requested
        
        if print_log == True:
            print('Calculating Function V^2 for each Geophone...')
        #==Calculate Function V Channel square
        ms_dataFVG2=numpy.zeros((int(len(ms_data)/3),len(ms_data[0])))
        
        for Geophone in range(int(len(ms_data)/3)): # Signal of Each geophone is defined by sum of squares of each channel
            ms_dataFVG2[Geophone]=numpy.add(numpy.square(ms_data[Geophone*3+0]),numpy.add(numpy.square(ms_data[Geophone*3+1]),numpy.square(ms_data[Geophone*3+2])))
            
        #== Calculate STALTA
        if Calculate_STALTA==True:
            ms_dataFVG2=STALTA(Data=ms_dataFVG2,STA=STA,LTA=LTA)
            
         #== stack data?
        if stack==False:
            #==PLOTTING         
            if print_plot == True:
                Plot1C(Data=ms_dataFVG2, Time=Time, Title='V^2 in each Geophone',Define_Time_Range=Define_Time_Range,Time_Range=Time_Range)  
        
        else:
            stacked=numpy.zeros(len(ms_dataFVG2[0]))
            for Channel in range(len(ms_dataFVG2)):
                stacked = numpy.add(stacked,ms_dataFVG2[Channel])
            ms_dataFVG2=stacked
            
            #==PLOTTING stack        
            if print_plot == True:
                Plotstack(Data=ms_dataFVG2, Time=Time, Title='stack of V^2 in each channel',Define_Time_Range=Define_Time_Range,Time_Range=Time_Range)   
            
        #== Returning Function
        if print_log == True:
            print('Returning data...')
        return(ms_dataFVG2, Time)


        
    #%% FUNCTION KIND = 4, Cube OF ABSOLUTE VELOCITY for each Geophone  = |Vx|^3 + |Vy|^3 + |Vz|^3
           
    #==Function Filtered Waveform
    if function_kind == 4: #Verify if Function "4" (Geophne |V|^3) was requested
        
        if print_log == True:
            print('Calculating Function |V|^3 for each Geophone...')
        #==Calculate Function V Channel square
        ms_dataFVG3=numpy.zeros((int(len(ms_data)/3),len(ms_data[0])))
        
        for Geophone in range(int(len(ms_data)/3)): # Signal of Each geophone is defined by sum of squares of each channel
            ms_dataFVG3[Geophone]=numpy.add(numpy.power(numpy.absolute(ms_data[Geophone*3+0]),3),numpy.add(numpy.power(numpy.absolute(ms_data[Geophone*3+1]),3),numpy.power(numpy.absolute(ms_data[Geophone*3+2]),3)))
            
        #== Calculate STALTA
        if Calculate_STALTA==True:
            ms_dataFVG3=STALTA(Data=ms_dataFVG3,STA=STA,LTA=LTA)
            
         #== stack data?
        if stack==False:
            #==PLOTTING         
            if print_plot == True:
                Plot1C(Data=ms_dataFVG3, Time=Time, Title='|V|^3 in each Geophone',Define_Time_Range=Define_Time_Range,Time_Range=Time_Range)  
        
        else:
            stacked=numpy.zeros(len(ms_dataFVG3[0]))
            for Channel in range(len(ms_dataFVG3)):
                stacked = numpy.add(stacked,ms_dataFVG3[Channel])
            ms_dataFVG3=stacked
            
            #==PLOTTING stack        
            if print_plot == True:
                Plotstack(Data=ms_dataFVG3, Time=Time, Title='stack of |V|^3 in each channel',Define_Time_Range=Define_Time_Range,Time_Range=Time_Range)   
            
        #== Returning Function
        if print_log == True:
            print('Returning data...')
        return(ms_dataFVG3, Time)
################################################################################

#%% FUNCTION KIND = 5, KURTOSIS OF THE SQUARE VELOCITY OF EACH GEOPHONE -> KURT(V^2)
    
    #==Function of Kurtosis of V^2 through  geophones
    if function_kind == 5:     
        if print_log == True:
            print('Calculating Function Kurtosis for each Geophones...')
        #==Calculate Function V^2 in each geophone, then, V^3

        #==Initialization of Variables
        ms_dataKurt=numpy.zeros((int(len(ms_data)/3),len(ms_data[0]))) #Initializing the array for Kurtosis
        #[ms_dataFVG2,Time]=Charac_Functions.calculate_function(ms_data,function_kind=3) #Can't calculte this simple until modify input mehotd in Charac_Functions 
                
        #==Calculate Function V Channel square
        ms_dataFVG2=numpy.zeros((int(len(ms_data)/3),len(ms_data[0])))
        
        for Geophone in range(int(len(ms_data)/3)): # Signal of Each geophone is defined by sum of squares of each channel
            ms_dataFVG2[Geophone]=numpy.add(numpy.square(ms_data[Geophone*3+0]),numpy.add(numpy.square(ms_data[Geophone*3+1]),numpy.square(ms_data[Geophone*3+2])))
             
        for Geophone in range(int(len(ms_data)/3)): # Signal of Each geophone is defined by sum of squares of each channel
            for index in range(int(SampleKurt/2),len(ms_data[0])-int(SampleKurt/2)):
                ms_dataKurt[Geophone][index]=stats.kurtosis(ms_dataFVG2[Geophone][index-int(SampleKurt/2):index+int(SampleKurt/2)], axis=0, fisher=True, bias=True)
            MaxKurt=ms_dataKurt[Geophone].max()
            ms_dataKurt[Geophone]=ms_dataKurt[Geophone]/MaxKurt
        
       #== Calculate STALTA
        if Calculate_STALTA==True:
            ms_dataKurt=STALTA(Data=ms_dataKurt,STA=STA,LTA=LTA)
            
         #== stack data?
        if stack==False:
            #==PLOTTING         
            if print_plot == True:
                Plot1C(Data=ms_dataKurt, Time=Time, Title='Kurtosis in each Geophone',Define_Time_Range=Define_Time_Range,Time_Range=Time_Range)  
        
        else:
            stacked=numpy.zeros(len(ms_dataKurt[0]))
            for Channel in range(len(ms_dataKurt)):
                stacked = numpy.add(stacked,ms_dataKurt[Channel])
            ms_dataKurt=stacked
            
            #==PLOTTING stack        
            if print_plot == True:
                Plotstack(Data=ms_dataKurt, Time=Time, Title='stack of Kurtosis in each Geophone',Define_Time_Range=Define_Time_Range,Time_Range=Time_Range)   
                    
        #== Returning Function
        if print_log == True:
            print('Returning data...')
        return(ms_dataKurt, Time)
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
    import om_ped_es_parameters_v2, os
    if om_ped_es_parameters_v2.verbose_level <= 1:
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
    if om_ped_es_parameters_v2.verbose_level <= 1:
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
    """
    signal_ma = numpy.zeros((len(signal),))
    for index in range(len(signal)):
         signal_ma[index] = numpy.sum(signal[index:(index+samples)])
    return (signal_ma/samples)
################################################################################

#%%% ##################################################################################################
def resume_array(original_array, first_sort, second_sort):
    import om_ped_es_parameters_v2
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
        
    if om_ped_es_parameters_v2.verbose_level <= 1: #For debugging
            print('Identified peaks after QC: ',resumed_array) 
    
    return(resumed_array)          
    
##################################################################################################

#%% Peak evaluation
#Reviewed
def peak_evaluation(signal, peaks_positions, time_torrent_index):
    import numpy, om_ped_es_parameters_v2, om_general_signal_processing
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
    
    
    # Eliminationg the peaks to close from the borders
    peaks_positions_sliced=[]
    for index in range(len(peaks_positions)):
        if (peaks_positions[index] >= time_torrent_index) and (peaks_positions[index] <= (len(signal)-time_torrent_index)):
            peaks_positions_sliced.append(peaks_positions[index])

    
    for peak_counter in range(len(peaks_positions_sliced)):
        
        #Evaluation the peak position and width
        for index in range(peaks_positions_sliced[peak_counter],len(signal)):
            if signal[index]<threshold:
                right_idx=index
                break
        for index in range(peaks_positions_sliced[peak_counter],0,-1):
            if signal[index]<threshold:
                left_idx=index
                break        

        #Registering the properties of each peak -> Properties other than Position index and SNR can be used later
        peaks_properties[peak_counter][0]=left_idx                           #Arrival index
        peaks_properties[peak_counter][1]=(signal[peaks_positions_sliced[peak_counter]]/numpy.mean(signal))-1  #SNR
        peaks_properties[peak_counter][2]=right_idx-left_idx                 #Width of Peak
        
    
    #At this point, all info is calculated. Now we have to make sure to exclude repetead peaks
    # Ideal situation: peak with "^"shape. Regular situation:  peak with "M" shape.
    #The method used is to compare the arrival index. If they are the same, keep the peak with maximum SNR
    
    #lets verify if any of the arraival times are equal:
    identified_peaks=om_general_signal_processing.resume_array(original_array=peaks_properties, first_sort=0, second_sort=1)

    #Some Quality control for the events based in peaks width and SNR
    identified_peaks_qc=[]
    for index in range(len(identified_peaks)):
        if identified_peaks[index][2] >= om_ped_es_parameters_v2.peak_width_minimum and peaks_properties[index][1]>=om_ped_es_parameters_v2.peak_snr_minimum :
            identified_peaks_qc.append(peaks_properties[index][:])
      
    return(identified_peaks_qc) 
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
        ax.set_xlabel('Data #', fontsize=14)
        ax.set_ylabel('Amplitude', fontsize=14)
        mode = 'Valley detection' if valley else 'Peak detection'
        ax.set_title("%s (mph=%s, mpd=%d, threshold=%s, edge='%s')"
                     % (mode, str(mph), mpd, str(threshold), edge))
        # plt.grid()
        plt.show()

 #%%