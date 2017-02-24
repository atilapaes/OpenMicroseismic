#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 16:15:53 2017

Modified on Thu Feb 16

@author: Atila Paes

Project:    OpenMicroseismic
Pack:       Potential-Event Detection
Method:     Energy-stack
File:       Main
Version:    V2
"""

import time, multiprocessing, numpy #Importing libraries from Python
import om_ped_es_functions_v2, om_ped_es_parameters_v2,om_general_signal_processing #Defined by user


#Initializing time of processing
start_time = time.time()

print('Initializing Potential-Event detection...')

### List the MS files on the folder
list_of_files=om_general_signal_processing.list_extensions(om_ped_es_parameters_v2.dataset_folder)


### Spliting the array in the number of cores to be used
file_list_for_cores=numpy.array_split(list_of_files,om_ped_es_parameters_v2.number_of_cores)


###############################################################################
### MULTICORE PROCESSING        
if __name__ == '__main__':
    jobs = []
    for core_counter in range(om_ped_es_parameters_v2.number_of_cores):
        #Defingn inputs and outputs of the task to the cores
        p = multiprocessing.Process(target=om_ped_es_functions_v2.energy_stack, args=(om_ped_es_parameters_v2.dataset_folder,file_list_for_cores[core_counter][:],core_counter))
        jobs.append(p)
        p.start()   
############################################################################### 
