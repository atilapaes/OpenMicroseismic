#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 07:35:25 2017

@author: AtilaPaes

This module is used to import Microseismic data in several different ways.

"""

import obspy, numpy, datetime
#%% #Retuns content of a file
""" Reviewed"""
def input_file(dataset_folder,file_name):
    ms_data=obspy.read(dataset_folder+"/"+file_name)
    return(ms_data)
#================================================================================================================

#%% Returns the content of two concatenated files or a slice of the first + the second
def InputTorrent(File1,File2,Partitial_File1=True,Last_Seconds=5):
    MSDataFile1=obspy.read(File1)  
    MSDataFile2=obspy.read(File2)  

    if Partitial_File1==True:
        MSData=MSDataFile2    
        for Channel in range(len(MSDataFile1)):
            MSData[Channel] +=MSDataFile1[Channel].slice(MSDataFile1[Channel].stats.endtime-Last_Seconds,MSDataFile1[Channel].stats.endtime)
    else:    
        MSData=MSDataFile1
        for Channel in range(len(MSDataFile1)):
            MSData[Channel] +=MSDataFile2[Channel]   
    
    return(MSData)

#================================================================================================================

#%% Inform a time and extract it from dataset
def input_time(Input_Time,DatasetName,Delta_Time, Load_Before=2,Load_After=5, Print_Log=False):

    """
    This module imports MSData from a specific time in a dataset. 
    By default, the slice starts 2 seconds before and ends 5 seconds after the specified time.
    Due the short time range, this script doesn't include slices bigger than 2 files.
    
    Parameters
    ----------
    Input_Time:     A specific time (%Y%m%d%H%M%S).
    DatasetName:    The name of TXT file contaning the Dataset list. Result from "ls -1 >> Dataset.txt".
    Delta_Time:     The total monitoring time of each file.
    Load_Before:    Time in seconds to load before the Input_Time
    Load_After:     Time in seconds to load after the Input_Time
    
    Print_Log:  {True or False}, optional (default=False)
                Print the Log of processing.

    Returns a slice of the seismogram from Input_Time - Load_Before to input_Time +Load_After
    -------------------------------------------------------------------------------------------
    
    FUTURE IMPROVEMENTS:
    -------------------
    Calculate Delta_Time instead asking the use
    """
        
    #%%Variable initialization
    if Print_Log==True:
        print("Variable initialization...")
    ListOfFiles=numpy.genfromtxt(str(DatasetName+".txt"),dtype=str) 
    Input_Time=datetime.datetime.strptime(Input_Time, '%Y%m%d%H%M%S')
    Load_Before=datetime.timedelta(seconds=Load_Before)
    Load_After=datetime.timedelta(seconds=Load_After)
    TimeOfFiles=[]    
    Identified_File=[]
    
    #== Convert txt file in list of time
    for ii in range(len(ListOfFiles)):
        TimeOfFiles.append(datetime.datetime.strptime(ListOfFiles[ii].split(".")[0], '%Y%m%d%H%M%S'))

    
    #==First test - Identify if the provided test is into the Dataset list
    Log=[]
    Indentification = False
    
    #Verify if the Input_Time is in the time range
    if Input_Time <= TimeOfFiles[-1] and Input_Time >= TimeOfFiles[0]: 
        Log.append('INPUT TIME IS IN DATASET RANGE')
        
        for FileNumber in range(len(ListOfFiles)): #Identifying the file
            if Input_Time >= TimeOfFiles[FileNumber] and Input_Time <= TimeOfFiles[FileNumber]+datetime.timedelta(seconds=Delta_Time):
                Indentification = True
                
                #Identified file list: 0-Identified file, 1-Previous file, 2-Next file
                Identified_File=[ListOfFiles[FileNumber],ListOfFiles[FileNumber-1],ListOfFiles[FileNumber+1]]
                
        if Indentification == False:
            Log.append("INPUT TIME IS IN A GAP OF MONITORING")

    else:
        Log.append('INPUT TIME IS NOT IN DATASET RANGE')

    if Print_Log == True:
        print("Log of file identification: ",Log)    
        
    #%% Loading file contaning the time    
    MSData=obspy.read(str(DatasetName+"/"+Identified_File[0]))
    
    #NOTE: OBSPY ONLY TAKES TIME ARGUMENTS IN FORMAT OF ITS OWN UTCDateTime
    #Case 1: Slice into single file:
    if obspy.core.utcdatetime.UTCDateTime(Input_Time - Load_Before) >= MSData[0].stats.starttime and obspy.core.utcdatetime.UTCDateTime(Input_Time + Load_After) <= MSData[0].stats.endtime:
        if Print_Log==True:
            print("Case 1: Slice into single file.")
        #MSData=MSData.slice(obspy.core.utcdatetime.UTCDateTime(Input_Time-Load_Before),obspy.core.utcdatetime.UTCDateTime(Input_Time+Load_After))
        
    else: #Cases 2 and 3: Time too close of file's border
        
        #Case 2 - Slice to close to file beginning
        if obspy.core.utcdatetime.UTCDateTime(Input_Time - Load_Before) < MSData[0].stats.starttime and obspy.core.utcdatetime.UTCDateTime(Input_Time + Load_After) < MSData[0].stats.endtime:           
            if Print_Log==True:
                print("Case 2: Slice too close to file beginning. Merging with previous file.")
            MSData_Aux=obspy.read(str(DatasetName+"/"+Identified_File[1])) #Loading auxiliar file
        
        #Case 3 - Slice too close to file end
        if obspy.core.utcdatetime.UTCDateTime(Input_Time - Load_Before) > MSData[0].stats.starttime and obspy.core.utcdatetime.UTCDateTime(Input_Time + Load_After) > MSData[0].stats.endtime:
            if Print_Log==True:
                print("Case 3: Slice too close to file end. Merging with next file.")
            MSData_Aux=obspy.read(str(DatasetName+"/"+Identified_File[2])) #Loading auxiliar file
            print("--->>> MSData_Aux",str(DatasetName+"/"+Identified_File[2]))
            print(MSData_Aux)
        
        # Merging files - OBSPY automatically reorder time position of arrays  
        for Channel in range(len(MSData)): 
            MSData[Channel] += MSData_Aux[Channel]
        
        if Print_Log==True:
            print("Slice after merge:")
            print(MSData)
    
    #Slicing the requested sample
    MSData=MSData.slice(obspy.core.utcdatetime.UTCDateTime(Input_Time-Load_Before),obspy.core.utcdatetime.UTCDateTime(Input_Time+Load_After))

    return(MSData)
#================================================================================================================

#%%    Make a Ring buffer from the first to the last file specified.
def MakeRingBuffer(FirstFile,LastFile, DatasetName):
    """
    FirstFile:      First File to be loaded in the ring buffer.
    LastFile:       Last File to be loaded in the ring buffer
    DatasetName:    Name of the dataset. Same name of folder and TXT file
    
    """
    #%%
    ListOfFiles=numpy.genfromtxt(str(DatasetName+".txt"),dtype=str) 
    print("Length of dataset", len(ListOfFiles))
    ListIndex=numpy.zeros(2)
    
    #Identifing the First and last file in the ringbuffer
    for FileNumber in range(len(ListOfFiles)):
        if ListOfFiles[FileNumber]==FirstFile:
            ListIndex[0]=FileNumber
        if ListOfFiles[FileNumber]==LastFile:
            ListIndex[1]=FileNumber
    
    MSData=obspy.read(str(DatasetName+"/"+ListOfFiles[int(ListIndex[0])])) #Loading auxiliar file    
    
    for FileNumber in range(int(ListIndex[0])+1,int(ListIndex[1])+1):
        print(ListOfFiles[FileNumber])
        MSData_Aux=obspy.read(str(DatasetName+"/"+ListOfFiles[FileNumber])) #Loading auxiliar file  
        
        # Merging files - OBSPY automatically reorder time position of arrays  
        for Channel in range(len(MSData)): 
            MSData[Channel] += MSData_Aux[Channel]    
    
    return(MSData)

