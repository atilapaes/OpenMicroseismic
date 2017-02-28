#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 11:18:23 2017

@author: Atila Paes

Project:    OpenMicroseismic
Pack:       Potential-Event Detection
Method:     Energy-stack
File:       Functions
"""


#%% Main function of the Energy Stack method
##################################################################################################    
def energy_stack(dataset_folder,list_of_files,core_counter,start_time):
    import numpy, time, datetime, om_ped_es_parameters_v2, om_general_input_data, om_general_signal_processing
    
    if om_ped_es_parameters_v2.verbose_level <= 2:
         print('Core number ',core_counter, '-> Initializing...')

    #Creation of identified_events List
    identified_events=[]
    peaks_properties=[]
    
    #==Loop over the file list=================================================
    for file_number in range(len(list_of_files)):
        #print('List of files ',list_of_files)
        
        ### Loading MS Data using torrent method ==============================
        ms_data=om_general_input_data.input_torrent(dataset_folder=dataset_folder, list_of_files=list_of_files, file_number=file_number, last_seconds=om_ped_es_parameters_v2.last_seconds)                

        #Software log
        if om_ped_es_parameters_v2.verbose_level <= 2 and core_counter==0 :
            print(round((file_number/len(list_of_files))*100,1),'% complete')    
        #======================================================================
        
        #### Processing center ================================================
        #Calculating the characteristic function 
        if om_ped_es_parameters_v2.verbose_level <= 1:
            print_log = print_plot = True
        else:
            print_log = print_plot = False

        [stacked,time_array]=om_general_signal_processing.calculate_function(ms_data=ms_data, print_log = print_log, function_kind=3, print_plot=print_plot, stack=True)
            
        #Curve smooth by Moving average
        if om_ped_es_parameters_v2.verbose_level <= 1:
            show=True
        else:
            show=False
        stacked_ma=om_general_signal_processing.moving_avg(stacked,om_ped_es_parameters_v2.moving_average_samples)

        #Detecting highest peaks-positioin into the signal
        peaks_positions=om_general_signal_processing.detect_peaks(x=stacked_ma, mph=om_ped_es_parameters_v2.threshold_means*numpy.mean(stacked_ma)+om_ped_es_parameters_v2.threshold_stds*numpy.std(stacked_ma), mpd=om_ped_es_parameters_v2.mpd, show=show)
        
        #Processing in case of at leat one identified peak
        if len(peaks_positions)!=0: 
            if om_ped_es_parameters_v2.verbose_level == 1:
                print('Peak position',peaks_positions) #For debugging
            
            #Calculationg the peaks properties and discarding the repeated events
            peaks_properties=om_general_signal_processing.peak_evaluation(signal=stacked_ma, peaks_positions=peaks_positions, time_torrent_index=int(om_ped_es_parameters_v2.time_torrent_coef/ms_data[0].stats.delta))
            #peaks_properties=om_general_signal_processing.peak_evaluation(signal=stacked_ma, peaks_positions=peaks_positions)
            if om_ped_es_parameters_v2.verbose_level == 0:
                print('Peaks properties',peaks_properties) #For debugging
            
            #Calculating time of occurence
            for index in range (len(peaks_properties)):
                # Output event time in format '%Y%m%d%H%M%S'
                 #Identifing the time in seconds in the time array. Truncked in seconds
                time_array[int(peaks_properties[index][0])]
                event_time=int((ms_data[0].stats.starttime+datetime.timedelta(seconds=int(time_array[int(peaks_properties[index][0])]))).strftime('%Y%m%d%H%M%S'))
                
                # Output the SNR with 2 decimal places
                snr=round(peaks_properties[index][1],2)
                
                identified_events.append([event_time,snr])
                if om_ped_es_parameters_v2.verbose_level <= 1:    
                    print("N.th event time and SNR", identified_events[-1])
        ### END OF PROCESSING CENTER ==========================================
    #==========================================================================

    ### Exporting results   ===================================================
    PE_Output='PE_Results_Core'+str(core_counter)+'_'+ dataset_folder+'.csv'
    numpy.savetxt(PE_Output, identified_events, delimiter=",",newline='\n')
    
    
    if om_ped_es_parameters_v2.verbose_level <= 1:
        print("Identified events. Time and SNR",identified_events)

    # Log
    #if om_ped_es_parameters_v2.verbose_level <= 1:
    print('Core', core_counter, '-> Processing Done!')
    #==========================================================================

    #=============================================
    ###Printing the processing time ###############################################
    #if core_counter+1==om_ped_es_parameters_v2.coreN: #Presenting the elapsed time after exporting the last Text file 
    elapsed_time = time.time() - start_time
    print('Elapsed time ', round(elapsed_time,3), 'sec',' in Core ',core_counter)

    return()

