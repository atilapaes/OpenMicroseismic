#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 02:34:08 2017

@author: atilapaes

Open Microseismic Project:
Main function for Potential Event detection

"""
#%% Include path to all modules into shared_lib folder

import pandas, numpy
import matplotlib.pyplot as plt
import om_load_data, om_char_function, om_ped_functions

#%% Paramethers

peak_distance=3         # Minimum peak distance in seconds
cut_index_sec =2.5      # Number of sec to be desconsidered at begining and end of file for processing peaks.
                        # Don'r worry, there are some data loading superposition to properly handle this.
peak_threshold_std=1    # Numbers of STD to add to the mean to set threshold
peak_width_min_sec=0.1  # Time in sec for min width of peak
peak_width_max_sec=5    # Time in sec for max width of peak

superposition_sec=6     # Number of seconds to superpose loading of adjacents files

# Tools to follow data processing
show_plot=False
print_plot=False
plot_title='day xx, hour yy' # Example day xx, hour yy

#%% File loading and folders definition

dataset_folder='/Volumes/toc2me/Fox_Creek_3D_2016_3C/Fox_Creek_3D_2016_3C/Fox_Creek_3D_2016_3C_SEG2/'

catalog_file='toc2me_3c_d305_test.csv'
file_catalog=pandas.read_csv('input_catalogs/'+catalog_file)
output_file='output_catalogs/output_'+catalog_file


gph_catalog=pandas.read_csv('catalog_gph.csv')
    
#%% Output file
output = open(output_file, 'w')
output.write('peak_time,peak_height,signal_mean,signal_std,width,peak_time_left,peak_time_right' + '\n')
output.close()

#%% Function that tranform a number in two digits
def int2two_digit(number):
    if number < 10:
    		output='0'+str(number)
    else:
    		output=str(number)
    return(output)
    
#%% Data loading and pre processing - including data counter
for index_file in range(len(file_catalog)):#(len(file_catalog)):
    print('Processing ',("%.2f" % round((index_file/len(file_catalog))*100,3)), '% complete')
    if index_file==len(file_catalog)-1:
        #load single file
        ms_data_3c=om_load_data.load_ms_files(file1_name=file_catalog.file_name[len(file_catalog)-1],file2_name='nan', folder=dataset_folder)
        ms_data_3c=om_load_data.slice_and_filter_v2(ms_data=ms_data_3c,start_time=ms_data_3c[0].stats.starttime,end_time=ms_data_3c[0].stats.endtime)    
    else:
        ms_data_3c=om_load_data.load_ms_files(file1_name=file_catalog.file_name[index_file],file2_name=file_catalog.file_name[index_file+1], folder=dataset_folder)      
  
        # File size inseconds
        file_size_sec=int(ms_data_3c[0].stats.npts*ms_data_3c[0].stats.delta/2)
        ms_data_3c=om_load_data.slice_and_filter_v2(ms_data=ms_data_3c,start_time=ms_data_3c[0].stats.starttime,end_time=ms_data_3c[0].stats.starttime+file_size_sec+superposition_sec)
        
#%% Characteristic function 
    es=om_char_function.energy_stack_selected_gph_v2(stream=ms_data_3c,gph_catalog=gph_catalog)
    es_mavg=om_char_function.aux_moving_avg(signal=es,samples=200)
    es_mavg=es_mavg/es_mavg.max()
    
    if print_plot==True:
        plt.figure(index_file)
        plt.plot(es_mavg)
        plt.title(str(index_file)+plot_title)
        plt.savefig(str(index_file),format='png')
        plt.close()    
    
    cf_mean=numpy.mean(es_mavg)
    cf_std=numpy.std(es_mavg) 
        
#%% C.F. Analysis
    # Analise of es_mavg and picking of index in max peaks
    peaks_positions=om_ped_functions.detect_peaks(x=es_mavg,mph=cf_mean+peak_threshold_std*cf_std,mpd=int(peak_distance/ms_data_3c[0].stats.delta), show=show_plot)        
        
    if len(peaks_positions)!=0: # Processing the case of at least one identified peak
        # print("Len peaks_position > 0") # For debug tests
        
        # This function cuts off edges and retunr just the properties of valid peaks
        peaks_properties=om_ped_functions.peak_evaluation(signal=es_mavg,peaks_positions=peaks_positions, cut_index=int(cut_index_sec/ms_data_3c[0].stats.delta),peak_threshold=cf_mean+peak_threshold_std*cf_std, peak_width_min=int(peak_width_min_sec/ms_data_3c[0].stats.delta), peak_width_max=int(peak_width_max_sec/ms_data_3c[0].stats.delta))
        
        for index_peak_prop in range(len(peaks_properties)):

            peak_time=ms_data_3c[0].stats.starttime+(peaks_properties[index_peak_prop][0])*ms_data_3c[0].stats.delta
            peak_time_str=str(peak_time.year)+(int2two_digit(peak_time.month))+(int2two_digit(peak_time.day))+(int2two_digit(peak_time.hour))+(int2two_digit(peak_time.minute))+(int2two_digit(peak_time.second))#+str(peak_time.microsecond)
                        
            peak_time_width=str('%.3f' %((peaks_properties[index_peak_prop][4])*ms_data_3c[0].stats.delta))

            peak_time_left=ms_data_3c[0].stats.starttime+((peaks_properties[index_peak_prop][5])*ms_data_3c[0].stats.delta)
            peak_time_left_str=str(peak_time_left.year)+(int2two_digit(peak_time_left.month))+(int2two_digit(peak_time_left.day))+(int2two_digit(peak_time_left.hour))+(int2two_digit(peak_time_left.minute))+(int2two_digit(peak_time_left.second))+(int2two_digit(peak_time_left.microsecond))
                      
            peak_time_right=ms_data_3c[0].stats.starttime+((peaks_properties[index_peak_prop][6])*ms_data_3c[0].stats.delta)
            peak_time_right_str=str(peak_time_right.year)+int2two_digit(peak_time_right.month)+int2two_digit(peak_time_right.day)+int2two_digit(peak_time_right.hour)+int2two_digit(peak_time_right.minute)+int2two_digit(peak_time_right.second)+int2two_digit(peak_time_right.microsecond)

            #%% Appending results in external file
            output=open(output_file, 'a')        
            output.write(peak_time_str +','+ str('%.5f' %(peaks_properties[index_peak_prop][1])) +','+ str('%.5f' %(peaks_properties[index_peak_prop][2])) +','+ str('%.5f' %(peaks_properties[index_peak_prop][3])) +','+ peak_time_width +','+ peak_time_left_str +','+  peak_time_right_str +'\n') #                     
            output.close()
    #%%
    else: 
        print('===> NO PEAKS IN THE SPECIFIED PARAMETERS <===')
print('===> END OF SCRIPT <===')        