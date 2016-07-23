
##################################################################################################
def listKindofFile(folder,extension):
	"""Argunments: string (folder). This function list all the SEGY files into the speciefied folder
	"""

	print('Listing SEGYs -> Initializing...')
	### Function listSEGYs####
	### This function enters a directory, list the files and returns a list of SEGYs

	import os

	### Interaction with the user ########################
	
	folderofSEGYs = folder

	######################################################

	### Entering on the folder contaning the SEGYs
	os.chdir(folderofSEGYs)

	## List all the files in the directory
	listoffiles=os.listdir()

	## Initialing the list of SEGYs
	listofSEGYs=[]

	##Loop to identify the SEGYs and 
	for index in range(len(listoffiles)):
		if listoffiles[index].endswith("."+extension):
			listofSEGYs = listofSEGYs + [listoffiles[index]]
	
	### Leaving the folder contaning the SEGYs##################

	os.chdir('..')
	print('Listing SEGYs -> Done!')
	return(listofSEGYs)


##################################################################################################
##################################################################################################
def axisname(stream):
	import numpy

	totalchannels=len(stream) #total number of channels

	for index1 in range(int(totalchannels/3)):
		
		stream[(0+index1*3)].stats.channel='Z'
		stream[(0+index1*3)].stats.location='Geoph'+str(index1+1)
	
		stream[(1+index1*3)].stats.channel='N'
		stream[(1+index1*3)].stats.location='Geoph'+str(index1+1)

		stream[(2+index1*3)].stats.channel='E'
		stream[(2+index1*3)].stats.location='Geoph'+str(index1+1)

	return(stream)

##################################################################################################
##################################################################################################

def writefile(msdataRam,filename):
	"""This function will be used to write a certain file into the memory
	msdataRam -> the structure to be stored
	filename -> name of the file to store
	"""
	import shelve
	
	## Save the data in a external file
	##################################
	#Call shelve.open and pass it a filename to open/initialize
	database=shelve.open(filename) 
	# Store the variable msdata from .read as 'msdata' in the database struct
	database['msdataStored']=msdataRam
	#Close the file
	database.close() 


##################################################################################################
##################################################################################################

def Concatenation2(folder,filestoload,fileheader):
	""" folder -> name of folder to enter
		filestoload -> amount of files to load in each external file
	
		Output: file 1: segys from 0 to (files to load -1)
				file 2: segys from (files to load -1) to (2x files to load -1)
				last file: segys from (nx files to load -1) to last segy
	"""

	import os, datetime, time, shelve, obspy, numpy
	import FunctionsBRP

	## Input data  #################
	#folder="Subrun1"
	#filestoload=10
	################################


	#Listing the files
	
	listofSEGYs=FunctionsBRP.listKindofFile(folder,"segy")
	#listofSEGYs=FunctionsBRP.listSEGY(folder)
	numberoffiles=len(listofSEGYs)


	### Creation of the array to files to load
	#loadarray=numpy.arange(0,numberoffiles,filestoload)
	loadarray=numpy.arange(0,numberoffiles,filestoload-1)
	loadarray=numpy.append(loadarray,numberoffiles)
	#print(loadarray)
	numberoflooops=len(loadarray)-1

	#Entering the folder to load files
	os.chdir(folder)


	sample=obspy.read(listofSEGYs[0])
	nchannels=len(sample)
	print("======================================================")
	for ii in range(0,numberoflooops): # This loop runs the code over each loop of n files
		
		print("Range of files from ",loadarray[ii]," to ",loadarray[ii+1])
		print("################################")

		for jj in range(loadarray[ii],loadarray[ii+1]+1): # This loop loads the n files

			if jj == loadarray[ii]:
				print("Loading the first file ", jj)
				subrun=obspy.read(listofSEGYs[jj])
			else:
				if jj != loadarray[-1]: #load all files but the last on the list. First element is zero
					print("Loading file ", jj)
					filebuffer=obspy.read(listofSEGYs[jj])	
					for channel in range(nchannels): #This loop concatenate the channels of the filebuffer with the subrun
						subrun[channel]=subrun[channel] +filebuffer[channel]	
			#print("######")	
		
		print("----Writing the data on disk---File ", str(ii+1)," of ",str(numberoflooops))

		filename=str(fileheader)+"-"+"file-"+str(ii)+"-RAW"
		FunctionsBRP.writefile(subrun,filename)
		
		print("################################################################")		
	print("======================================================")


	os.chdir('..')

	####################################


##################################################################################################
##################################################################################################
## This is the test file for development of the data treatment on the raw concatenated data
## Input freq_min, freq_max, filename
## Output: another .db file filtered

def FilterRawData(filename,folder,freq_min,freq_max):
	import os, datetime, time, shelve, obspy
	import numpy as np
	import FunctionsBRP
	
	#Entering the folder to load files
	os.chdir(folder)

	shelffile=shelve.open(filename)
	msdata_raw=shelffile['msdataStored']
	shelffile.close()

	#### Collecting info about the file #####
	ntraces=len(msdata_raw) #Getting the number of traces on the stream

	msdata=msdata_raw.copy()

	# Removing the AVG
	for ii in range(ntraces):
		msdata[ii].data=msdata_raw[ii].data-np.mean(msdata_raw[ii].data)

	# Naming the axises
	FunctionsBRP.axisname(msdata)

	#### Filtering by BUTTERWORTH ####################################################
	# Band pass filter
	msdata.filter('bandpass', freqmin=freq_min, freqmax=freq_max, corners=4, zerophase=True)

	#Remove of 60 Hz
	msdata.filter('bandstop', freqmin=59.5, freqmax=60.5, corners=4, zerophase=True)

	#Remove of 180 Hz
	msdata.filter('bandstop', freqmin=179.5, freqmax=180.5, corners=4, zerophase=True)
	###########################################################################


	##### Folders navigation ##################################################
	## Leaving the folder of the raw files
	os.chdir('..')
	
	######Folder of the filtered files
	folderfiltered=folder+'-Filtered' # Name of the folder
	
	if not os.path.exists(folderfiltered):  #Verifying if the folder already exists
		os.makedirs(folderfiltered)
		print('Folder --> ', folderfiltered, ' <-- created')
	
	os.chdir(folderfiltered) # Entering the folder

	## 5) Save the filtered data in an external file 
	#Changing the end of the filename from 'RAW' to 'Filtered'
	filename_filtered=filename.split('RAW')[0]+'Filtered'

	FunctionsBRP.writefile(msdata,filename_filtered)

	os.chdir('..') # Leaving to the main folder
	
########  OBSOLETE SCRIPTS   ###################################################################
################################################################################################
################################################################################################
################################################################################################
def listSEGY(folder):
	"""Argunments: string (folder). This function list all the SEGY files into the speciefied folder
	"""

	print('Listing SEGYs -> Initializing...')
	### Function listSEGYs####
	### This function enters a directory, list the files and returns a list of SEGYs

	import os

	### Interaction with the user ########################
	
	folderofSEGYs = folder

	######################################################

	### Entering on the folder contaning the SEGYs
	os.chdir(folderofSEGYs)

	## List all the files in the directory
	listoffiles=os.listdir()

	## Initialing the list of SEGYs
	listofSEGYs=[]

	##Loop to identify the SEGYs and 
	for index in range(len(listoffiles)):
		if listoffiles[index].endswith(".segy"):
			listofSEGYs = listofSEGYs + [listoffiles[index]]
	
	### Leaving the folder contaning the SEGYs##################

	os.chdir('..')
	print('Listing SEGYs -> Done!')
	return(listofSEGYs)

##################################################################################################


##################################################################################################
def Concatenation(folder):
	"""Argument: two strings, first is folder name and second is the number of subrun.
	This function concatenate all the SEGYs into the folder and save a .db file of the final variable
	"""
	
	import os, datetime, time, shelve, obspy
	import FunctionsBRP
	print('Concatenation -> Initializing...')

	listofSEGYs=FunctionsBRP.listKindofFile(folder,"segy")
	#listofSEGYs=FunctionsBRP.listSEGY(folder)
	
	### Entering on the folder contaning the SEGYs
	os.chdir(folder)

	subrun=obspy.read(listofSEGYs[0])
	nchannels=len(subrun)


	for index1 in range(len(listofSEGYs)): #This loop loads the SEGY files from the second
		filebuffer=obspy.read(listofSEGYs[index1])

		for index2 in range(nchannels): #This loop concatenate the channels of the filebuffer with the subrun
			subrun[index2]=subrun[index2] +filebuffer[index2]

	### Leaving the folder contaning the SEGYs##################
	os.chdir('..')
	
	print('Concatenation -> Done!')	
	print('Writing external file -> Initializing...')
	
	## Save the SubRun in a external file at the same folder of the scripts
	################################################################################
	shelfFile=shelve.open(folder) #Call shelve.open and pass it a fileneme "***" 
	# A file named ***.db will be created
	shelfFile['folder']=subrun# Store the variable "subrun" as "subrun#" in the file
	shelfFile.close() #Close the file
	################################################################################
	print('Writing external file -> Done!')

##################################################################################################

	