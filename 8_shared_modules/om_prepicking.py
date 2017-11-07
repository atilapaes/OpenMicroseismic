#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 19:17:32 2017
@author: atilapaes
"""

import numpy, math
import numpy as np
from scipy import signal
from math import pi, sin,cos
#%%
def moving_avg_arriving(signal,mavg_samples):
    """
    This function calculates the moving average of a provided 1-C signal (array).    
    """
    signal_ma = numpy.zeros((len(signal),))
    
    #Regular signal
    for index in range(mavg_samples, len(signal)-mavg_samples):
        signal_ma[index] = (numpy.sum(signal[index:(index+mavg_samples)]))/mavg_samples

    signal_mean=numpy.mean(signal_ma[mavg_samples:len(signal)-mavg_samples])
    
    for index in range(0,mavg_samples):
        signal_ma[index] = signal_mean

    for index in range(len(signal)-mavg_samples, len(signal)):
        signal_ma[index] = signal_mean
       
    return (signal_ma/signal_ma.max())   
################################################################################

#%% CALCULATE CHARACTERISTICS FUNCTIONS
def signal_energy_stack(ms_data, apply_mavg=True, mavg_samples=30):    
    """    
    This function calculates the energy stack of the signal for 1C or 3C streams
    
    Parameters:
    -----------
    ms_data:        1c or 3C stream of data
    mavg_samples:   number of samples to use in the moving average    
    """
    import numpy
    #%% Calculating Characteristic Function, SQUARE OF VELOCITY for each CHANNEL
    if len(ms_data)==1:
        ms_data[0].data=numpy.square(ms_data[0].data)
    elif len(ms_data)==3:        
        ms_data[0].data=numpy.add(numpy.square(ms_data[0].data),numpy.add(numpy.square(ms_data[1].data),numpy.square(ms_data[2].data)))
        
    #== Data smothing
    if apply_mavg==True:
        ms_data[0].data = moving_avg_arriving(signal=ms_data[0].data,mavg_samples=mavg_samples)
    
    ms_data[0].plot()
    return(ms_data[0].data)

#%%
def flinn(stream, noise_thres=0):
    """
    Computes the azimuth, incidence, rectilinearity and planarity after the
    eigenstructure decomposition method of [Flinn1965b]_.

    :param stream: ZNE sorted trace data
    :type stream: List of ZNE sorted numpy arrays
    :param noise_tresh: Variance of the noise sphere; data points are excluded
        when falling within the sphere of radius sqrt(noise_thres),
        default is set to 0.
    :type noise_thres: float
    :returns:  azimuth, incidence, rectilinearity, and planarity
    """
    import numpy as np
    mask = (stream[0][:] ** 2 + stream[1][:] ** 2 + stream[2][:] ** 2
            ) > noise_thres
    x = np.zeros((3, mask.sum()), dtype=np.float64)
    # East
    x[0, :] = stream[2][mask]
    # North
    x[1, :] = stream[1][mask]
    # Z
    x[2, :] = stream[0][mask]

    covmat = np.cov(x)
    eigvec, eigenval, v = np.linalg.svd(covmat)
    # Rectilinearity defined after Montalbetti & Kanasewich, 1970
    rect = 1.0 - np.sqrt(eigenval[1] / eigenval[0])
    # Planarity defined after [Jurkevics1988]_
    plan = 1.0 - (2.0 * eigenval[2] / (eigenval[1] + eigenval[0]))
    azimuth = math.degrees(math.atan2(eigvec[0][0], eigvec[1][0]))
    eve = np.sqrt(eigvec[0][0] ** 2 + eigvec[1][0] ** 2)
    incidence = math.degrees(math.atan2(eve, eigvec[2][0]))
    if azimuth < 0.0:
        azimuth = 360.0 + azimuth
    if incidence < 0.0:
        incidence += 180.0
    if incidence > 90.0:
        incidence = 180.0 - incidence
        if azimuth > 180.0:
            azimuth -= 180.0
        else:
            azimuth += 180.0
    if azimuth > 180.0:
        azimuth -= 180.0
    
    data_return_flinn=[azimuth, incidence, rect, plan]
    return(data_return_flinn)

#%%
def particle_motion_odr(stream, noise_thres=0):
    """
    Computes the orientation of the particle motion vector based on an
    orthogonal regression algorithm.

    :param stream: ZNE sorted trace data
    :type stream: :class:`~obspy.core.stream.Stream`
    :param noise_tres: variance of the noise sphere; data points are excluded
        when falling within the sphere of radius sqrt(noise_thres)
    :type noise_thres: float
    :returns: azimuth, incidence, error of azimuth, error of incidence
    """
    import scipy.odr
    import numpy as np
    z = []
    n = []
    e = []
    comp, npts = np.shape(stream)

    for i in range(0, npts):
        if (stream[0][i] ** 2 + stream[1][i] ** 2 + stream[2][i] ** 2) \
                > noise_thres:
            z.append(stream[0][i])
            n.append(stream[1][i])
            e.append(stream[2][i])

    def fit_func(beta, x):
        # XXX: Eventually this is correct: return beta[0] * x + beta[1]
        return beta[0] * x

    data = scipy.odr.Data(e, n)
    model = scipy.odr.Model(fit_func)
    odr = scipy.odr.ODR(data, model, beta0=[1.])
    out = odr.run()
    az_slope = out.beta[0]
    az_error = out.sd_beta[0]

    n = np.asarray(n)
    e = np.asarray(e)
    z = np.asarray(z)
    r = np.sqrt(n ** 2 + e ** 2)

    data = scipy.odr.Data(r, abs(z))
    model = scipy.odr.Model(fit_func)
    odr = scipy.odr.ODR(data, model, beta0=[1.0])
    out = odr.run()
    in_slope = out.beta[0]
    in_error = out.sd_beta[0]

    azimuth = math.atan2(1.0, az_slope)
    incidence = math.atan2(1.0, in_slope)

    az_error = 1.0 / ((1.0 ** 2 + az_slope ** 2) * azimuth) * az_error
    # az_error = math.degrees(az_error)
    in_error = 1.0 / ((1.0 ** 2 + in_slope ** 2) * incidence) * in_error
    # in_error = math.degrees(in_error)

    azimuth = math.degrees(azimuth)
    incidence = math.degrees(incidence)

    if azimuth < 0.0:
        azimuth = 360.0 + azimuth
    if incidence < 0.0:
        incidence += 180.0
    if incidence > 90.0:
        incidence = 180.0 - incidence
        if azimuth > 180.0:
            azimuth -= 180.0
        else:
            azimuth += 180.0
    if azimuth > 180.0:
        azimuth -= 180.0

    data_return_odr=[azimuth, incidence, az_error, in_error]
    return(data_return_odr)

#%%
def rotate_zne_lqt(z, n, e, ba, inc):
    """
    Rotates all components of a seismogram.

    The components will be rotated from ZNE (Z, North, East, left-handed) to
    LQT (e.g. ray coordinate system, right-handed). The rotation angles are
    given as the back-azimuth and inclination.

    The transformation consists of 3 steps::

        1. mirroring of E-component at ZN plain: ZNE -> ZNW
        2. negative rotation of coordinate system around Z-axis with angle ba:
           ZNW -> ZRT
        3. negative rotation of coordinate system around T-axis with angle inc:
           ZRT -> LQT

    :type z: :class:`~numpy.ndarray`
    :param z: Data of the Z component of the seismogram.
    :type n: :class:`~numpy.ndarray`
    :param n: Data of the North component of the seismogram.
    :type e: :class:`~numpy.ndarray`
    :param e: Data of the East component of the seismogram.
    :type ba: float
    :param ba: The back azimuth from station to source in degrees.
    :type inc: float
    :param inc: The inclination of the ray at the station in degrees.
    :return: L-, Q- and T-component of seismogram.
    """    
    if len(z) != len(n) or len(z) != len(e):
        raise TypeError("Z, North and East component have different length!?!")
    if ba < 0 or ba > 360:
        raise ValueError("Back Azimuth should be between 0 and 360 degrees!")
    if inc < 0 or inc > 360:
        raise ValueError("Inclination should be between 0 and 360 degrees!")
    ba *= 2 * pi / 360
    inc *= 2 * pi / 360
    l = z * cos(inc) - n * sin(inc) * cos(ba) - e * sin(inc) * sin(ba)
    q = z * sin(inc) + n * cos(inc) * cos(ba) + e * cos(inc) * sin(ba)
    t = n * sin(ba) - e * cos(ba)
    
    rotated=[l, q, t]
    return(rotated)

#%% Rotation Horizontal
def rotate_ne_rt(n, e, ba):
    """
    Rotates horizontal components of a seismogram.

    The North- and East-Component of a seismogram will be rotated in Radial
    and Transversal Component. The angle is given as the back-azimuth, that is
    defined as the angle measured between the vector pointing from the station
    to the source and the vector pointing from the station to the North.

    :type n: :class:`~numpy.ndarray`
    :param n: Data of the North component of the seismogram.
    :type e: :class:`~numpy.ndarray`
    :param e: Data of the East component of the seismogram.
    :type ba: float
    :param ba: The back azimuth from station to source in degrees.
    :return: Radial and Transversal component of seismogram.
    """
    if len(n) != len(e):
        raise TypeError("North and East component have different length.")
    if ba < 0 or ba > 360:
        raise ValueError("Back Azimuth should be between 0 and 360 degrees.")
    r = e * sin((ba + 180) * 2 * pi / 360) + n * cos((ba + 180) * 2 * pi / 360)
    t = e * cos((ba + 180) * 2 * pi / 360) - n * sin((ba + 180) * 2 * pi / 360)
    return [r, t]