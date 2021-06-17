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

# degree/radian/minute of angle(moa)
def d2m(d): return d*60.0             # degrees to MOA
def d2r(d): return d*np.pi/180.0      # degrees to radians
def r2d(r): return r*180.0/np.pi      # radians to degrees
def r2m(r): return r*60.0*180.0/np.pi # radians to MOA
def m2d(m): return m/60.0             # MOA to degrees
def m2r(m): return m/60.0*np.pi/180.0 # MOA to radians

# hack trigonometry fcts to return 0. rather than very small number for certain
# degrees i.e. sin of 180 degrees

# all trig fcts expect values in radians but have built-in convert-from parameter
# where convert-from is one of {'r':radians i.e. do nothing,'d'=degrees.'m'=moa}
# TODO:
#  is there a better way of setting elements to 0
def cos(x,m='d'):
    if m == 'd': x = d2r(x)
    elif m == 'm': x = m2r(x)
    x = np.cos(x)
    z = np.isclose(x,0.)
    if isinstance(z,np.ndarray):
        if np.any(z):
            for i,zi in enumerate(z):
                if zi: x[i] = 0.
    elif z: x = 0.
    return x

def sin(x,m='d'):
    if m == 'd': x = d2r(x)
    elif m == 'm': x = m2r(x)
    x = np.sin(x)
    z = np.isclose(x,0.)
    if isinstance(z,np.ndarray):
        if np.any(z):
            for i,zi in enumerate(z):
                if zi: x[i] = 0.
    elif z: x = 0.
    return x

def tan(x,m='d'):
    if m == 'd': x = d2r(x)
    elif m == 'm': x = m2r(x)
    x = np.tan(x)
    z = np.isclose(x, 0.)
    if isinstance(z,np.ndarray):
        if np.any(z):
            for i,zi in enumerate(z):
                if zi: x[i] = 0.
    elif z: x = 0.
    return x

def arccos(x,m='d'):
    if m == 'd': x = d2r(x)
    elif m == 'm': x = m2r(x)
    x = np.arccos(x)
    z = np.isclose(x, 0.)
    if isinstance(z,np.ndarray):
        if np.any(z):
            for i,zi in enumerate(z):
                if zi: x[i] = 0.
    elif z: x = 0.
    return x

def arcsin(x,m='d'):
    if m == 'd': x = d2r(x)
    elif m == 'm': x = m2r(x)
    x = np.arcsin(x)
    z = np.isclose(x, 0.)
    if isinstance(z,np.ndarray):
        if np.any(z):
            for i,zi in enumerate(z):
                if zi: x[i] = 0.
    elif z: x = 0.
    return x

def arctan(x,m='d'):
    if m == 'd': x = d2r(x)
    elif m == 'm': x = m2r(x)
    x = np.arctan(x)
    z = np.isclose(x, 0.)
    if isinstance(z,np.ndarray):
        if np.any(z):
            for i,zi in enumerate(z):
                if zi: x[i] = 0.
    elif z: x = 0.
    return x

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

# TODO: the below need work
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
