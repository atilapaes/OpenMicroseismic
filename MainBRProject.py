##DATA LOADING 

import os, datetime, time, shelve, obspy, numpy, FunctionsBRP


####### STEP 1 - DATA LOADING  #############
#"""
### Input data ########
wellnumber=3
injection=1
subrun=1

folder='Subrun1'
filestoload=10
#######################

### Creation of the filename
#fileheader="W"+ str(wellnumber)+"-"+"I"+str(injection)+"-"+"SR"+str(subrun) # Example W3-I1-SR1

 
###Loading SEGY Raw data and exporting to external file
#FunctionsBRP.Concatenation2(folder,filestoload,fileheader)
#"""


###########################################################################

####### STEP 2 - DATA TREATMENT  ###################
## Procedures in data treatment:

freq_min=50
freq_max=300


# Generate the list of files to load.
list_raw_concateneted=FunctionsBRP.listKindofFile(folder,'db')
# Excluding the extension .db from the list
for ii in range(len(list_raw_concateneted)):
	list_raw_concateneted[ii]=list_raw_concateneted[ii].split('.db')[0]


for ii in range(len(list_raw_concateneted)): 
	FunctionsBRP.FilterRawData(list_raw_concateneted[ii],folder,freq_min,freq_max)
	print('Data treatment on file ',ii+1,' of ', len(list_raw_concateneted))


# Into a loop, apply the filtering, name the channels and save an external file
## 1) Apply mean remove and filtering (Just Butterworth at first tests)





## 2) Naming the channels (cause the data in SEGYs from Petrobras has no complete headers)
