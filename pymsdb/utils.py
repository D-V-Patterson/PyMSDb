#!/usr/bin/env python
""" utils.py
Copyright (C) 2020 Dale V. Patterson (dale.v.patterson@gmail.com)

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

Defines utility functions to be used with histograms, pivots and other data
returned from a MSDB object
"""

#__name__ = 'utils'
__license__ = 'GPLv3'
__version__ = '0.0.2'
__date__ = 'April 2021'
__author__ = 'Dale Patterson'
__maintainer__ = 'Dale Patterson'
__email__ = 'dale.v.patterson@gmail.com'
__status__ = 'Development'

import numpy as np

# helper function
def is_int(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

def percent_change(old,new): return (new-old) / old * 100

def lobf(ts):
    """
    calculates line of best fit from points in ts
    :param ts: list of tuples t = (x,y)
    :return: the lobf
    """
    xs = [x for x,_ in ts]
    ys = [y for _,y in ts]
    [m,b] = np.polyfit(xs,ys,1)
    return [(x,y) for x,y in zip(np.array(xs),np.array(xs)*m+b)]

def ema(xs,p=5):
    """
    exponential moving average
    :param xs: list of values
    :param p: period
    :return: the ema as numpy array
    """
    # create the ema array, set initial value and calculate smoothing variable
    es = np.zeros_like(xs,dtype=np.float32)
    k = 2/(p+1)

    # iterate over xs and calcuate ema
    for i,x in enumerate(xs):
        if i == 0: es[i] = sum(xs[:p])/p
        else: es[i] = k*x + (1-k)*es[i-1]

    # return it
    return es

def ewma(xs,l=0.5):
    """
    exponential weighted moving average
    :param xs: list of values
    :param l: lambda, the smoothing weight
    :return: the ewma
    """
    es = np.zeros_like(xs,dtype=np.float32)
    for i,x in enumerate(xs):
        if i == 0: es[i] = np.average(xs)
        else: es[i] = l*x + (1-l)*es[i-1]
    return es

def sma(xs,p=5):
    # TODO returns list that is shorter than xs based on the value for p
    return np.convolve(xs,np.ones(p) / p,mode='valid')
