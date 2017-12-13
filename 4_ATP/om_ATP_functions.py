#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 17:30:21 2017

@author: atilapaes

Functions for ATP module

"""
import sys
om_root_folder=''
sys.path.insert(0,om_root_folder+'8_shared_modules')


#%%
import numpy
import om_char_function

#%%
# Parameters
kurt_samples=50
es_peak_width=0.2 # in seconds
kurt_time_before_peak=0.5
kurt_time_after_peak=2

kurt_diff_window_left_shift=0.5 # in seconds


#%%
###############################################################################
# ATP SOFTARE VERSION 3
###############################################################################

#%%  Method for picking the Chacteristic function max in a specified narrow window
def kurt_picking(kurt_diff,es_h_arg,es_z_arg,start_second, dt):    
    # Version 3 improves the windows for searching max kurt_diff. 
    # The idea is define windows based in [max_energy_index - shift: max_energy_index]
    
    left_shift_samples=int(kurt_diff_window_left_shift/dt)
    
    # for debug
    #print('interval z',es_z_arg- left_shift_samples,es_z_arg)
    #print('interval h',es_h_arg- left_shift_samples,es_h_arg)
    
    if (es_z_arg- left_shift_samples) >= 0:
        
        arg_max_kurt_p=(numpy.argmax(kurt_diff[0][es_z_arg- left_shift_samples:es_z_arg]))+(es_z_arg- left_shift_samples)
        time_p=arg_max_kurt_p*dt + start_second        
        
        arg_max_kurt_s1=(numpy.argmax(kurt_diff[1][es_h_arg- left_shift_samples:es_h_arg]))+(es_h_arg- left_shift_samples)
        time_s1=arg_max_kurt_s1*dt + start_second
            
        arg_max_kurt_s2=(numpy.argmax(kurt_diff[2][es_h_arg- left_shift_samples:es_h_arg]))+(es_h_arg- left_shift_samples)
        time_s2=arg_max_kurt_s2*dt + start_second
    
    else:
        print('Picking error: Cannot locate max charac. function into the interval.')
        time_p = time_s1 = time_s2 = None
    
    return(time_p,time_s1,time_s2)

#%%
def atp_kurt_diff(gph_data_3c,gph_number, show, event_name):
    import numpy
    import matplotlib.pyplot as plt

    # Prevent method to run over null data
    if numpy.sum(gph_data_3c[0]) != float(0) and numpy.sum(gph_data_3c[0]) !=float(0) and numpy.sum(gph_data_3c[0]) !=float(0):
        dt=gph_data_3c[0].stats.delta        
        
        # Functions for  search-windows definition                

        # Energy stack for verical and horizontal channels
        es_z=numpy.square(gph_data_3c[0].data)
        es_h=numpy.add(numpy.square(gph_data_3c[1].data),numpy.square(gph_data_3c[2].data))
        
        # Initial arguments for kurt_diff search
        if (numpy.argmax(es_z) <= numpy.argmax(es_h)):
            es_h_arg=numpy.argmax(es_h)
            es_z_arg=numpy.argmax(es_z[0:es_h_arg-20]) # maybe introduce a safety coeff = -20
        else:
            print('Picking error: ATP for S-wave earlier than P-wave.')
            time_p = time_s1 = time_s2 = None
            return(time_p,time_s1,time_s2)
            
        # Characteristic function
        kurt_diff=om_char_function.diff_kurt(stream=gph_data_3c,kurt_samples=kurt_samples)
        
        start_second=gph_data_3c[0].stats.starttime.second
        [time_p,time_s1,time_s2]=kurt_picking(kurt_diff, es_h_arg,es_z_arg,start_second, dt)
        
        # Print plot
        if show==True:
            time=numpy.arange(gph_data_3c[0].stats.starttime.second, gph_data_3c[0].stats.endtime.second+0.5*dt, dt)
            
            plt.figure(gph_number, figsize=(10,5))
            plt.subplot(3,1,1)
            plt.title(event_name + '- Gph '+str(gph_number))
            plt.plot(time,gph_data_3c[0].data,'r',lw=0.5)
            plt.plot(time,kurt_diff[0]-1,'k', lw=0.5)
            plt.scatter(time_p,0,c='m',marker='*')
            plt.axvline(x=time_p,c='m',lw=0.5)
            plt.xlim(time_p-0.5,time_s1+0.5)
            plt.xticks([])
    
            plt.subplot(3,1,2)
            plt.plot(time,gph_data_3c[1].data,'g',lw=0.5)
            plt.plot(time,kurt_diff[1]-1,'k',lw=0.5)
            plt.axvline(x=time_s1,c='m',lw=0.5)
            plt.scatter(time_s1,0,c='m',marker='*')
            plt.xlim(time_p-0.5,time_s1+0.5)
            plt.xticks([])
    
            plt.subplot(3,1,3)
            plt.plot(time,gph_data_3c[2].data,'b',lw=0.5)
            plt.plot(time,kurt_diff[2]-1,'k',lw=0.5)
            plt.scatter(time_s2,0,c='m',marker='*')
            plt.axvline(x=time_s2,c='m',lw=0.5)
            plt.xlim(time_p-0.5,time_s1+0.5)
            plt.show()
            plt.close()
        
        return(time_p,time_s1,time_s2)
    
    else:
        print('C.F. skipping null data.')
        time_p = time_s1 = time_s2 = None
        return(time_p,time_s1,time_s2)    
#%%
def rotate_ray_center(gph_data,print_results=False,print_plot=False,gph_number='No Info'):
    """
    This function rotates the gph signal from ZH1H2 to LQT coordinates
    """
    import numpy
    import matplotlib.pyplot as plt
    from obspy.signal.polarization import particle_motion_odr
    from obspy.signal.rotate import rotate_zne_lqt
   
    rotated_gph_data=gph_data.copy()
    time=numpy.linspace(gph_data[0].stats.starttime.second, gph_data[0].stats.endtime.second,gph_data[0].stats.npts)
    
    # Calculating angles for non-null arrays
    if ((numpy.sum(gph_data[0].data) !=0) and (numpy.sum(gph_data[1].data)!=0) and (numpy.sum(gph_data[2].data)!=0)):
        odr_azimuth, odr_incidence, odr_azi_error, odr_inc_error=particle_motion_odr(gph_data)

        # Correcting angles
        inclination=90-odr_incidence
        if odr_azimuth <180:
            ba=odr_azimuth+180
        else:
            ba=odr_azimuth-180

        # Print results for debugging
        if print_results==True:
            print('Gph_number: '+str(gph_number)+'. BA: ',str(round(ba,2)),'. Inclination: ', str(round(inclination,2)))
 
        # Signal Rotation
        rotated_gph_data[0].data, rotated_gph_data[1].data, rotated_gph_data[2].data=rotate_zne_lqt(gph_data[0].data, gph_data[1].data, gph_data[2].data, ba=ba, inc=inclination)
    else:
        print('Null gph signal identified')
   
    if print_plot==True:
        plt.figure(1,figsize=(8,5))
        
        plt.subplot(3,2,1)
        plt.title('Original data- Gph '+str(gph_number))
        plt.ylabel('Z', weight='bold')
        plt.plot(time,gph_data[0].data,lw=0.5,c='r')
        
        plt.subplot(3,2,3)
        plt.ylabel('h1', weight='bold')
        plt.plot(time,gph_data[1].data,lw=0.5,c='g')
        
        plt.subplot(3,2,5)
        plt.ylabel('h2', weight='bold')
        plt.plot(time,gph_data[2].data,lw=0.5,c='b')
                
        plt.subplot(3,2,2)
        plt.title('Rotated. BA: '+str(round(ba,2))+'. Inc: '+str(round(inclination,2)))
        plt.ylabel('L', weight='bold')
        plt.plot(time,rotated_gph_data[0].data,c='k',lw=0.5)
        
        plt.subplot(3,2,4)
        plt.ylabel('Q', weight='bold')
        plt.plot(time,rotated_gph_data[1].data,lw=0.5,c='k')
        
        plt.subplot(3,2,6)
        plt.ylabel('T', weight='bold')
        plt.plot(time,rotated_gph_data[2].data,lw=0.5,c='k')
        plt.tight_layout()
        
        plt.show()
        plt.close()
        
    return(rotated_gph_data)
#%%