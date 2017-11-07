#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 04:25:04 2017

@author: atilapaes

Open Microseismic Project:

Major functions used by Potential Event Detection by Energy Stack

"""
    #%%
def peak_evaluation(signal,peaks_positions, cut_index,peak_threshold, peak_width_min, peak_width_max):
    import numpy
    """
    This function calculates the ES properties of a list of peaks
    
    PARAMETERS:
    ==========
    
    signal              The es_mavg signal
    peaks_positions     The list of idenfified positions (in the array signal) of the peaks
    cut_index           Number of samples in the begining and end of signal array to cut off peaks from the list
    peak_threshold      The threshold used to identify peaks
    peak_width_min      Minimum number of samples of peak_width to validate a peak
    peak_width_limit    Maximum number of samples of peak_width to validate a peak    
    
    OUTPUT:
    =======    
    List contaning lists of event properties.
    
    """
    #%%
    def cut_borders(peaks_positions, index_start, index_end):
        """
        Cut peaks at the borders of Characteristic Function
        """
        output=[]
        for index in range(len(peaks_positions)):
            if (peaks_positions[index] >= index_start) and (peaks_positions[index] <= index_end):
                output.append(peaks_positions[index])
        return(output)
    
    #%%
    def peak_width(signal,signal_index):
        """
        Calculate both edges of the peak width
        """
        # Bait for Peaks Quality-Control. It garantees all peaks identified have all parameters calculated correctly.
        right_idx=left_idx=0
        
        #Evaluation the peak position and width
        for index_scan in range(peaks_positions[index_peak],len(signal)):
            if signal[index_scan]<= peak_threshold:
                right_idx=index_scan
                break
        for index_scan in range(peaks_positions[index_peak],0,-1):
            if signal[index_scan]<peak_threshold:
                left_idx=index_scan
                break
        if right_idx !=0 and left_idx!=0:
            return([left_idx,right_idx])
        else:
            return(None)
    #%%
    peaks_properties=[]
    
    # Eliminationg the peaks too close from the borders
    peaks_positions=cut_borders(peaks_positions, index_start=cut_index, index_end=len(signal)-cut_index)

    # Parameters calculation
    for index_peak in range(len(peaks_positions)):
        # Calculating peak width
        width = peak_width(signal,signal_index=peaks_positions[index_peak])
        
        # Quality control for width. (quality control for peak height was done in the detect_peak function)
        if width == None:
            print('===> Peak width NOT calculated correctly')
        else: 
            #print('===> Peak width calculated correctly')
            if (width[1]-width[0] >= peak_width_min) and (width[1]-width[0]) <= peak_width_max:
                #print('===> Peak width into parameter range')
                # Peak_position time, Peak height, signal mean, signal std , width, left_idx, right_idx
                peaks_properties.append([peaks_positions[index_peak],signal[peaks_positions[index_peak]],numpy.mean(signal),numpy.std(signal),width[1]-width[0],width[0], width[1]])
                #peaks_properties.append([peaks_positions[index_peak],str(round(signal[peaks_positions[index_peak]],4)),str(round(numpy.mean(signal),4)),str(round(numpy.std(signal),4)),str(round(width[1]-width[0],4)),str(round(width[0],4)), str(round(width[1],4))])
            #else:
                #print('===> Peak width OUT OF parameter range')

    return(peaks_properties)

################################################################################
#%%
#REVIEWED
def detect_peaks(x, mph=None, mpd=1, threshold=0, edge='rising',
                 kpsh=False, valley=False, show=False, ax=None):

    import numpy

    """Detect peaks in data based on their amplitude and other features.
    
    
    Detect peaks in data based on their amplitude and other features.
    http://nbviewer.jupyter.org/github/demotu/BMC/blob/master/notebooks/DetectPeaks.ipynb

    Atila's example
    u=detect_peaks.detect_peaks(data1, mph=numpy.mean(data1)+1.5*numpy.std(data1), mpd=100, show=True)
    
    
    Parameters
    ----------
    x : 1D array_like
        data.
    mph : {None, number}, optional (default = None)
        detect peaks that are greater than minimum peak height.
    mpd : positive integer, optional (default = 1)
        detect peaks that are at least separated by minimum peak distance (in
        number of data).
    threshold : positive number, optional (default = 0)
        detect peaks (valleys) that are greater (smaller) than `threshold`
        in relation to their immediate neighbors.
    edge : {None, 'rising', 'falling', 'both'}, optional (default = 'rising')
        for a flat peak, keep only the rising edge ('rising'), only the
        falling edge ('falling'), both edges ('both'), or don't detect a
        flat peak (None).
    kpsh : bool, optional (default = False)
        keep peaks with same height even if they are closer than `mpd`.
    valley : bool, optional (default = False)
        if True (1), detect valleys (local minima) instead of peaks.
    show : bool, optional (default = False)
        if True (1), plot data in matplotlib figure.
    ax : a matplotlib.axes.Axes instance, optional (default = None).

    Returns
    -------
    ind : 1D array_like
        indeces of the peaks in `x`.

    Notes
    -----
    The detection of valleys instead of peaks is performed internally by simply
    negating the data: `ind_valleys = detect_peaks(-x)`
    
    The function can handle NaN's 

    See this IPython Notebook [1]_.

    References
    ----------
    .. [1] http://nbviewer.ipython.org/github/demotu/BMC/blob/master/notebooks/DetectPeaks.ipynb

    """

    x = numpy.atleast_1d(x).astype('float64')
    if x.size < 3:
        return numpy.array([], dtype=int)
    if valley:
        x = -x
    # find indices of all peaks
    dx = x[1:] - x[:-1]
    # handle NaN's
    indnan = numpy.where(numpy.isnan(x))[0]
    if indnan.size:
        x[indnan] = numpy.inf
        dx[numpy.where(numpy.isnan(dx))[0]] = numpy.inf
    ine, ire, ife = numpy.array([[], [], []], dtype=int)
    if not edge:
        ine = numpy.where((numpy.hstack((dx, 0)) < 0) & (numpy.hstack((0, dx)) > 0))[0]
    else:
        if edge.lower() in ['rising', 'both']:
            ire = numpy.where((numpy.hstack((dx, 0)) <= 0) & (numpy.hstack((0, dx)) > 0))[0]
        if edge.lower() in ['falling', 'both']:
            ife = numpy.where((numpy.hstack((dx, 0)) < 0) & (numpy.hstack((0, dx)) >= 0))[0]
    ind = numpy.unique(numpy.hstack((ine, ire, ife)))
    # handle NaN's
    if ind.size and indnan.size:
        # NaN's and values close to NaN's cannot be peaks
        ind = ind[numpy.in1d(ind, numpy.unique(numpy.hstack((indnan, indnan-1, indnan+1))), invert=True)]
    # first and last values of x cannot be peaks
    if ind.size and ind[0] == 0:
        ind = ind[1:]
    if ind.size and ind[-1] == x.size-1:
        ind = ind[:-1]
    # remove peaks < minimum peak height
    if ind.size and mph is not None:
        ind = ind[x[ind] >= mph]
    # remove peaks - neighbors < threshold
    if ind.size and threshold > 0:
        dx = numpy.min(numpy.vstack([x[ind]-x[ind-1], x[ind]-x[ind+1]]), axis=0)
        ind = numpy.delete(ind, numpy.where(dx < threshold)[0])
    # detect small peaks closer than minimum peak distance
    if ind.size and mpd > 1:
        ind = ind[numpy.argsort(x[ind])][::-1]  # sort ind by peak height
        idel = numpy.zeros(ind.size, dtype=bool)
        for i in range(ind.size):
            if not idel[i]:
                # keep peaks with the same height if kpsh is True
                idel = idel | (ind >= ind[i] - mpd) & (ind <= ind[i] + mpd) \
                    & (x[ind[i]] > x[ind] if kpsh else True)
                idel[i] = 0  # Keep current peak
        # remove the small peaks and sort back the indices by their occurrence
        ind = numpy.sort(ind[~idel])

    if show:
        if indnan.size:
            x[indnan] = numpy.nan
        if valley:
            x = -x
        _plot(x, mph, mpd, threshold, edge, valley, ax, ind)

    return ind


def _plot(x, mph, mpd, threshold, edge, valley, ax, ind):
    """Plot results of the detect_peaks function, see its help."""
    import numpy
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print('matplotlib is not available.')
    else:
        if ax is None:
            _, ax = plt.subplots(1, 1, figsize=(8, 4))

        ax.plot(x, 'b', lw=1)
        if ind.size:
            label = 'valley' if valley else 'peak'
            label = label + 's' if ind.size > 1 else label
            ax.plot(ind, x[ind], '+', mfc=None, mec='r', mew=2, ms=8,
                    label='%d %s' % (ind.size, label))
            ax.axhline(mph,color='cyan')
            
            ax.legend(loc='best', framealpha=.5, numpoints=1)
        ax.set_xlim(-.02*x.size, x.size*1.02-1)
        ymin, ymax = x[numpy.isfinite(x)].min(), x[numpy.isfinite(x)].max()
        yrange = ymax - ymin if ymax > ymin else 1
        ax.set_ylim(ymin - 0.1*yrange, ymax + 0.1*yrange)
        ax.set_xlabel('data #', fontsize=14)
        ax.set_ylabel('Amplitude', fontsize=14)
        mode = 'Valley detection' if valley else 'Peak detection'
        ax.set_title("%s (mph=%s, mpd=%d, threshold=%s, edge='%s')"
                     % (mode, str(mph), mpd, str(threshold), edge))
        # plt.grid()
        plt.show()
#%%