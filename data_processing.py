# -*- coding: utf-8 -*-
#
# - Author: steve simmert
# - GitHub: https://github.com/stevosn
# - E-mail: steve.simmert@uni-tuebingen.de
# - Copyright: 2017#
################################################################################
# some data processing related stuff
################################################################################
from scipy.signal import convolve
from scipy import ones, sqrt

def block_avg(x, n, t=None, fs=1.0, method=None):
    """
    Calculate block averages of n data point of the given data in x.
    
    The function calulates non-overlapping block averages of subsequent n-tuples.
    If t is None, a time vector is generated according to the sampling frequency fs.
    
    Arguments
    ---------
    x : array
        Data array
    n : int
        block size, unmber of data point in the block
    t : None or array
        If None, a timevector is generated based on a sampling frequency fs.
        If t is given it should have the same length as x.
    fs : float
        Sampling frequency. Only relevant if t=None.
    method : None or 'convolve'
        For debugging. None let the function use 
        resize(array, (floor(len(x)/n), n)) and then calculates the means
        and stdev of the blocks of n data points. 'convolve' uses a convolution
        approach to calculate the block averages, only.
    
    Returns
    -------
    (t, x_ba, std) : tuple
    t : array
        Time points for d. Each time point corresponds to the center of the block.
    x_ba : array
        Block averaged data
    std : array
        Standard deviations of the blocks.
    """
    if method == 'convolve':
        d = convolve(x, ones(n), mode='valid')
        d_out = d[::n] / n
        sem = d_out * 0.0
    else:
        N = len(x)
        n = int(n)
        m = int(floor(N/n))
        d = resize(x.copy(), (m, n))
        d_out = d.mean(axis=1)
        sem = d.std(axis=1)
    if t is None:
        t = array([*range(len(data))]) / fs
    t_out = t[::n][:len(d_out)] + t[:n].mean()
    
    return (t_out, d_out, sem)
