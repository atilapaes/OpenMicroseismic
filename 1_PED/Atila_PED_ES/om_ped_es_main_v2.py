#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 16:15:53 2017

@author: Atila Paes

Project:    OpenMicroseismic
Pack:       Potential-Event Detection
Method:     Energy-stack
File:       Main
"""

import time, multiprocessing #Importing libraries from Python
import om_ped_es_functions_v2, om_ped_es_parameters_v2, om_general_input_data #Defined by user


#Initializing time of processing
start_time = time.time()

print('Initializing Potential-Event detection...')

### List files in folder and split in sub-list to send to cores
file_list_for_cores=om_general_input_data.split_list_of_files(folder_name=om_ped_es_parameters_v2.dataset_folder, number_of_cores=om_ped_es_parameters_v2.number_of_cores)


###############################################################################
### MULTICORE PROCESSING        
if __name__ == '__main__':
    jobs = []
    for core_counter in range(om_ped_es_parameters_v2.number_of_cores):
        #Defingn inputs and outputs of the task to the cores
        p = multiprocessing.Process(target=om_ped_es_functions_v2.energy_stack, args=(om_ped_es_parameters_v2.dataset_folder,file_list_for_cores[core_counter][:],core_counter, start_time))
        jobs.append(p)
        p.start()   
############################################################################### 
