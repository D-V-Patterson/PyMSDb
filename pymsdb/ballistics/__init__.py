#!/usr/bin/env python
""" ballistics
Copyright (C) 2020 Dale V. Patterson (dale.v.patterson@gmail.com)

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

 TODO: description

Top level directory for ballistics package
Contains:
 TODO

Defines conversion functions/constants
"""

#__name__ = 'ballistics'
__license__ = 'GPLv3'
__version__ = '0.0.1'
__date__ = 'May 2021'
__author__ = 'Dale V. Patterson'
__maintainer__ = 'Dale V. Patterson'
__email__ = 'dale.v.patterson@gmail.com'
__status__ = 'Development'

import numpy as np

# CONSTANTS
# position and velocity are treated as 3d vectors indexed as
IX = 0 # axis index x
IY = 1 # axis index y
IZ = 2 # axis index z

# CONVERSIONS

# degree/radian/minute of angle(moa)
def d2m(d): return d*60.0             # degrees to MOA
def d2r(d): return d*np.pi/180.0      # degrees to radians
def r2d(r): return r*180.0/np.pi      # radians to degrees
def r2m(r): return r*60.0*180.0/np.pi # radians to MOA
def m2d(m): return m/60.0             # MOA to degrees
def m2r(m): return m/60.0*np.pi/180.0 # MOA to radians

# imperial to metric
GR2GM   = 0.0647989  # grains to grams
GR2MG   = 64.7989    # grains to milligrams
GR2KG   = 6.47989e-5 # grains to kilograms
IN2MM   = 25.4       # inches to millimeter
IN2CM   = 2.54       # inches to centimeters

# metric conversions
MM2CM   = 0.1    # mm to cm
MM2M    = 0.001  # mm to m
M2MM    = 1e6    # m to mm

# metric to imperial
msd2isd = 1/703  # sectional density (SD), ballistic coefficent (BC) to imperial
ma2ia   = 1/1550 # cross-sectional area (A) to imperial

# hack trigonometry fcts to return 0. rather than very small number for certain
# degrees i.e. sin of 180 degrees

# all trig fcts expect values in radians but have built-in convert-fromparameter
# where convert-from is one of {'r':radians i.e. do nothing,'d'=degrees.'m'=moa}
# TODO:
#  these are built for singleton cases, in the event of an array being passed
#  as, will have to catch the error and not evaluate
def cos(x,m='r'):
    if m == 'd': x = d2r(x)
    elif m == 'm': x = m2r(x)
    x = np.cos(x)
    if type(x) is not np.ndarray and np.isclose(x,0.): x = 0.
    return x

def sin(x,m='r'):
    if m == 'd': x = d2r(x)
    elif m == 'm': x = m2r(x)
    x = np.sin(x)
    if type(x) is not np.ndarray and np.isclose(x,0.): x = 0.
    return x

def tan(x,m='r'):
    if m == 'd': x = d2r(x)
    elif m == 'm': x = m2r(x)
    x = np.tan(x)
    if type(x) is not np.ndarray and np.isclose(x, 0.): x = 0.
    return x

def arccos(x,m='r'):
    if m == 'd': x = d2r(x)
    elif m == 'm': x = m2r(x)
    x = np.arccos(x)
    if type(x) is not np.ndarray and np.isclose(x,0.): x = 0.
    return x

def arcsin(x,m='r'):
    if m == 'd': x = d2r(x)
    elif m == 'm': x = m2r(x)
    x = np.sin(x)
    if type(x) is not np.ndarray and np.isclose(x,0.): x = 0.
    return x

def arctan(x,m='r'):
    if m == 'd': x = d2r(x)
    elif m == 'm': x = m2r(x)
    x = np.tan(x)
    if type(x) is not np.ndarray and np.isclose(x,0.): x = 0.
    return x
