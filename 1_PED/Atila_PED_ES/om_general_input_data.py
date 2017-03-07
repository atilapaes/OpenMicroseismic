#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 07:35:25 2017

@author: AtilaPaes

This module is used to import Microseismic data in several different ways.

"""

import obspy, numpy, datetime
#%% #Retuns content of a file
def input_file(dataset_folder,file_name):
    ms_data=obspy.read(dataset_folder+"/"+file_name)
    return(ms_data)

#================================================================================================================
#%% Returns the content of two concatenated files or a slice of the first + the second
def concatenate_files(file1,file2,partitial_file1=True,last_seconds=5):
    ms_data_file1=obspy.read(file1)  
    ms_data_file2=obspy.read(file2)  

    if partitial_file1==True:
        ms_data=ms_data_file2    
        for channel in range(len(ms_data_file1)):
            ms_data[channel] +=ms_data_file1[channel].slice(ms_data_file1[channel].stats.endtime-last_seconds,ms_data_file1[channel].stats.endtime)
    else:    
        ms_data=ms_data_file1
        for channel in range(len(ms_data_file1)):
            ms_data[channel] +=ms_data_file2[channel]   
    
    return(ms_data)
    
#================================================================================================================
#%% Inform a time and extract it from dataset
def input_time(input_time,dataset_name, load_before=2,load_after=5, print_log=False):
    import om_general_signal_processing
    """
    FUTURE IMPROVEMENTS:
    -------------------
    Make a mixed argument between print log, the verbose level and approach using case excelption
    Condition at the moment: delta_time <= file size
    -------------------
    
    This module imports ms_data from a specific time in a dataset. 
    By default, the slice starts 2 seconds before and ends 5 seconds after the specified time.
    Due the short time range, this script doesn't include slices bigger than 2 files.
    
    Parameters
    ----------
    input_time:     A specific time (%Y%m%d%H%M%S).
    dataset_name:    The name of folder file contaning the Dataset list. Result from "ls -1 >> Dataset.txt".
    delta_time:     The total monitoring time of each file.
    load_before:    Time in seconds to load before the input_time
    load_after:     Time in seconds to load after the input_time
    
    print_log:  {True or False}, optional (default=False)
                Print the log of processing.

    Returns a slice of the seismogram from input_time - load_before to input_time +load_after
    -------------------------------------------------------------------------------------------
    
    """
        
    #%% Variable initialization
    if print_log==True:
        print("Variable initialization...")
    
    #== List_of_files=numpy.genfromtxt(str(dataset_name+".txt"),dtype=str) -> Old version
    list_of_files=om_general_signal_processing.list_extensions(dataset_name)    
    
    input_time=datetime.datetime.strptime(input_time, '%Y%m%d%H%M%S')
    load_before=datetime.timedelta(seconds=load_before)
    load_after=datetime.timedelta(seconds=load_after)
    time_of_files=[]    
    identified_files=[]


    #== Convert the list of file names in list of time
    for ii in range(len(list_of_files)):
        time_of_files.append(datetime.datetime.strptime(list_of_files[ii].split(".")[0], '%Y%m%d%H%M%S'))

    # Consulting the time duration of each file
    sample_data=obspy.read(str(dataset_name+"/"+list_of_files[0]))
    delta_time=int(sample_data[0].stats.npts*sample_data[0].stats.delta)
        
    ### First test - Identify if the provided test is into the Dataset list
    log=[]
    identification = False
    
    #== Verify if the input_time is in the time range =========================
    if input_time <= time_of_files[-1] and input_time >= time_of_files[0]: #COndition to be in dataset range
        log.append('INPUT TIME IS IN DATASET RANGE')
        
        for file_number in range(len(list_of_files)): #Identifying the file
            if input_time >= time_of_files[file_number] and input_time <= time_of_files[file_number]+datetime.timedelta(seconds=delta_time):
                identification = True
                identified_file_number=file_number
                #identified_start_time=time_of_files[file_number]
                log.append("INPUT TIME WAS IDENTIFIED")
                
                #Identified file list: 0-Identified file, 1-Previous file, 2-Next file
                #identified_files=[list_of_files[file_number],list_of_files[file_number-1],list_of_files[file_number+1]]                
                
        if identification == False:
            log.append("INPUT TIME IS IN A GAP OF MONITORING")

    else:
        log.append('INPUT TIME IS NOT IN DATASET RANGE')

    if print_log == True:
        print("Log of file identification: ",log)    
    #==========================================================================    

    #%% Loading continuous dataset and slicing the data of interest 
    #== NOTE: OBSPY ONLY TAKES TIME ARGUMENTS IN FORMAT OF ITS OWN UTCDateTime
    if  identification == True: # In case the time is identified into the dataset, let's to slice it from the dataset and return to the workflow
        
        ms_data=obspy.read(str(dataset_name+"/"+identified_files[identified_file_number]))   
        
        # Case 1: Too close from the file beginning:
        if obspy.core.utcdatetime.UTCDateTime(input_time - load_before) < time_of_files[identified_file_number]:
    
            if identified_file_number == 0:
                if print_log==True:
                    print("Case 2: Slice too close to first file beginning.")
                
                # Slice dataset
                ms_data=ms_data.slice(time_of_files[identified_file_number],time_of_files[identified_file_number]+delta_time)
            
            else: # Merging two subsequent files
                if print_log==True:
                    print("Case 2: Slice too close to file beginning. Merging with previous file.")

                #Loading auxiliar file
                ms_data_aux=obspy.read(str(dataset_name+"/"+identified_files[identified_file_number+1]))

                # Merging files - OBSPY automatically reorder time position of arrays  
                for channel in range(len(ms_data)): 
                    ms_data[channel] += ms_data_aux[channel]

                # Slicing 
                ms_data=ms_data.slice(obspy.core.utcdatetime.UTCDateTime(input_time-load_before),obspy.core.utcdatetime.UTCDateTime(input_time+load_after))


        # Case 2: Too close from file end
        elif (obspy.core.utcdatetime.UTCDateTime(input_time + load_after)) > (time_of_files[identified_file_number]+ datetime.timedelta(seconds=delta_time)):
            
            if identified_file_number ==len(time_of_files):
                if print_log==True:
                    print("Case 2: Slice too close to first file beginning.")
                
                # Slice dataset
                ms_data=ms_data.slice(time_of_files[identified_file_number]+delta_time,time_of_files[identified_file_number])
            
            else: # Merging two subsequent files
                if print_log == True:
                    print("Case 3: Slice too close to file end. Merging with next file.")

                # Loading auxiliar file
                ms_data_aux=obspy.read(str(dataset_name+"/"+identified_files[identified_file_number-1]))
                
                # Merging files - OBSPY automatically reorder time position of arrays  
                for channel in range(len(ms_data)): 
                    ms_data[channel] += ms_data_aux[channel]
                
                # Slicing 
                ms_data=ms_data.slice(obspy.core.utcdatetime.UTCDateTime(input_time-load_before),obspy.core.utcdatetime.UTCDateTime(input_time+load_after))
        
        # Case 3: Data into a single file
        elif (obspy.core.utcdatetime.UTCDateTime(input_time - load_before) >= time_of_files[identified_file_number]) and (obspy.core.utcdatetime.UTCDateTime(input_time + load_after) <= time_of_files[identified_file_number]+ datetime.timedelta(seconds=delta_time)):
            if print_log==True:
                print("Case 3: Slice into single file.")

            # Slicing
            ms_data=ms_data.slice(obspy.core.utcdatetime.UTCDateTime(input_time-load_before),obspy.core.utcdatetime.UTCDateTime(input_time+load_after))

    else:   # Case 4: Aler the user that the time is not present 
        print('Error loading file: Input time is not into in the dataset.')

    return(ms_data)


#================================================================================================================
#%% Inform a time and extract it from dataset
def input_time_old(input_time,dataset_name,delta_time, load_before=2,load_after=5, print_log=False):

    """
    
    Future improvements: exclude delta_time and the txt file from arguments
    
    
    This module imports ms_data from a specific time in a dataset. 
    By default, the slice starts 2 seconds before and ends 5 seconds after the specified time.
    Due the short time range, this script doesn't include slices bigger than 2 files.
    
    Parameters
    ----------
    input_time:     A specific time (%Y%m%d%H%M%S).
    dataset_name:    The name of TXT file contaning the Dataset list. Result from "ls -1 >> Dataset.txt".
    delta_time:     The total monitoring time of each file.
    load_before:    Time in seconds to load before the input_time
    load_after:     Time in seconds to load after the input_time
    
    print_log:  {True or False}, optional (default=False)
                Print the log of processing.

    Returns a slice of the seismogram from input_time - load_before to input_time +load_after
    -------------------------------------------------------------------------------------------
    
    FUTURE IMPROVEMENTS:
    -------------------
    Calculate delta_time instead asking the use
    """
        
    #%%Variable initialization
    if print_log==True:
        print("Variable initialization...")
    list_of_files=numpy.genfromtxt(str(dataset_name+".txt"),dtype=str) 
    input_time=datetime.datetime.strptime(input_time, '%Y%m%d%H%M%S')
    load_before=datetime.timedelta(seconds=load_before)
    load_after=datetime.timedelta(seconds=load_after)
    time_of_files=[]    
    identified_files=[]
    
    #== Convert txt file in list of time
    for ii in range(len(list_of_files)):
        time_of_files.append(datetime.datetime.strptime(list_of_files[ii].split(".")[0], '%Y%m%d%H%M%S'))

    
    #==First test - Identify if the provided test is into the Dataset list
    log=[]
    identification = False
    
    #Verify if the input_time is in the time range
    if input_time <= time_of_files[-1] and input_time >= time_of_files[0]: 
        log.append('INPUT TIME IS IN DATASET RANGE')
        
        for file_number in range(len(list_of_files)): #Identifying the file
            if input_time >= time_of_files[file_number] and input_time <= time_of_files[file_number]+datetime.timedelta(seconds=delta_time):
                identification = True
                
                #Identified file list: 0-Identified file, 1-Previous file, 2-Next file
                identified_files=[list_of_files[file_number],list_of_files[file_number-1],list_of_files[file_number+1]]
                
        if identification == False:
            log.append("INPUT TIME IS IN A GAP OF MONITORING")

    else:
        log.append('INPUT TIME IS NOT IN DATASET RANGE')

    if print_log == True:
        print("Log of file identification: ",log)    
        
    #%% Loading file contaning the time    
    ms_data=obspy.read(str(dataset_name+"/"+identified_files[0]))
    
    #NOTE: OBSPY ONLY TAKES TIME ARGUMENTS IN FORMAT OF ITS OWN UTCDateTime
    #Case 1: Slice into single file:
    if obspy.core.utcdatetime.UTCDateTime(input_time - load_before) >= ms_data[0].stats.starttime and obspy.core.utcdatetime.UTCDateTime(input_time + load_after) <= ms_data[0].stats.endtime:
        if print_log==True:
            print("Case 1: Slice into single file.")
        #ms_data=ms_data.slice(obspy.core.utcdatetime.UTCDateTime(input_time-load_before),obspy.core.utcdatetime.UTCDateTime(input_time+load_after))
        
    else: #Cases 2 and 3: Time too close of file's border
        
        #Case 2 - Slice to close to file beginning
        if obspy.core.utcdatetime.UTCDateTime(input_time - load_before) < ms_data[0].stats.starttime and obspy.core.utcdatetime.UTCDateTime(input_time + load_after) < ms_data[0].stats.endtime:           
            if print_log==True:
                print("Case 2: Slice too close to file beginning. Merging with previous file.")
            ms_data_aux=obspy.read(str(dataset_name+"/"+identified_files[1])) #Loading auxiliar file
        
        #Case 3 - Slice too close to file end
        if obspy.core.utcdatetime.UTCDateTime(input_time - load_before) > ms_data[0].stats.starttime and obspy.core.utcdatetime.UTCDateTime(input_time + load_after) > ms_data[0].stats.endtime:
            if print_log==True:
                print("Case 3: Slice too close to file end. Merging with next file.")
            ms_data_aux=obspy.read(str(dataset_name+"/"+identified_files[2])) #Loading auxiliar file
            print("--->>> ms_data_aux",str(dataset_name+"/"+identified_files[2]))
            print(ms_data_aux)
        
        # Merging files - OBSPY automatically reorder time position of arrays  
        for channel in range(len(ms_data)): 
            ms_data[channel] += ms_data_aux[channel]
        
        if print_log==True:
            print("Slice after merge:")
            print(ms_data)
    
    #Slicing the requested sample
    ms_data=ms_data.slice(obspy.core.utcdatetime.UTCDateTime(input_time-load_before),obspy.core.utcdatetime.UTCDateTime(input_time+load_after))

    return(ms_data)

#================================================================================================================
#%%    Make a Ring buffer from the first to the last file specified.
def make_ring_buffer(first_file,last_file, dataset_name):
    """
    first_file:      First File to be loaded in the ring buffer.
    last_file:       Last File to be loaded in the ring buffer
    dataset_name:    Name of the dataset. Same name of folder and TXT file
    
    """
    #%%
    list_of_files=numpy.genfromtxt(str(dataset_name+".txt"),dtype=str) 
    print("Length of dataset", len(list_of_files))
    ListIndex=numpy.zeros(2)
    
    #Identifing the First and last file in the ringbuffer
    for file_number in range(len(list_of_files)):
        if list_of_files[file_number]==first_file:
            ListIndex[0]=file_number
        if list_of_files[file_number]==last_file:
            ListIndex[1]=file_number
    
    ms_data=obspy.read(str(dataset_name+"/"+list_of_files[int(ListIndex[0])])) #Loading auxiliar file    
    
    for file_number in range(int(ListIndex[0])+1,int(ListIndex[1])+1):
        print(list_of_files[file_number])
        ms_data_aux=obspy.read(str(dataset_name+"/"+list_of_files[file_number])) #Loading auxiliar file  
        
        # Merging files - OBSPY automatically reorder time position of arrays  
        for channel in range(len(ms_data)): 
            ms_data[channel] += ms_data_aux[channel]    
    
    return(ms_data)

#================================================================================================================
#%% Returns the content of two concatenated files or a slice of the first + the second

def input_torrent(dataset_folder, list_of_files, file_number, last_seconds):
    """
    For a list o file names, this function concatenate files with previous one, except for the first file.
    It is mainly used as data torrent for STA/LTA calculations.
    """
    
    if file_number == 0:
        ms_data=obspy.read(dataset_folder+"/"+list_of_files[0])
    
    else:
        ms_data_file1=obspy.read(dataset_folder+"/"+list_of_files[file_number-1])
        ms_data_file2=obspy.read(dataset_folder+"/"+list_of_files[file_number])

        ms_data=ms_data_file2    
        for channel in range(len(ms_data_file1)):
            ms_data[channel] +=ms_data_file1[channel].slice(ms_data_file1[channel].stats.endtime-last_seconds,ms_data_file1[channel].stats.endtime)
    return(ms_data)
    
#================================================================================================================
#%% Make a list of files to be used in input_torrent
def split_list_of_files(folder_name, number_of_cores):
    """
    List MS Data files in a given dataset. 
    Divide in N sub-lists.
    Include the last element in sub-list X as first element in list X+1.
    
    Original list: [a b c] [d e f] [g h i] [j l m]
    New list:      [a b c] [c d e f] [f g h i] [i j l m]
    """
    import numpy, om_general_signal_processing    
    
    ### List the MS files in the folder
    list_of_files=om_general_signal_processing.list_extensions(folder_name)
    
    ### Spliting the array in the number of cores to be used- Includes case of N cores =1
    file_list_for_cores=numpy.array_split(list_of_files,number_of_cores)
    
    if number_of_cores != 1:    
        # Insert the last file of list N-1 in the beginning of list N
        for list_index in range(1,number_of_cores):
            file_list_for_cores[list_index]=numpy.insert(file_list_for_cores[list_index],0,file_list_for_cores[list_index-1][-1])
            
    return(file_list_for_cores)
    
#================================================================================================================

#%% Make a list of files to be used in input_torrent
def split_list_of_times(folder_name, number_of_cores):
    """
    List MS Data files in a given dataset. 
    Divide in N sub-lists.
    Include the last element in sub-list X as first element in list X+1.
    
    Original list: [a b c d e f g h i j l m]
    New list:      [a b c] [d e f] [g h i] [j l m]
    """
    import numpy, csv
    file_name=folder_name+'.csv'
    
    # Read text file 
    exampleFile = open(file_name)
    exampleReader = csv.reader(exampleFile)
    exampleData = list(exampleReader)
    
    # Use the first column of the file as list of time
    list_of_times=[]
    for index in range(len(exampleData)):
        list_of_times.append(exampleData[index][0])
        
    # Spliting the array in the number of cores to be used- Includes case of N cores =1
    file_list_for_cores=numpy.array_split(list_of_times,number_of_cores)
                
    return(file_list_for_cores)
