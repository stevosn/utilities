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
from scipy import ones, sqrt, floor, resize, array

def block_avg(x, n, t=None, fs=1.0):
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
    N = len(x)
    n = int(n)
    m = int(floor(N/n))
    x_ = resize(x.copy(), (m, n))
    x_ba = x_.mean(axis=1)
    std = x_.std(axis=1)
    
    if t is None:
        t = array([*range(len(x))]) / fs

    t_out = t[::n][:len(x_ba)] + t[:n].mean()
    
    return (t_out, x_ba, std)
