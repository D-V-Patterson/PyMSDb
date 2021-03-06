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

# CONSTANTS
# position and velocity are treated as 3d vectors indexed as
IX = 0 # axis index x
IY = 1 # axis index y
IZ = 2 # axis index z

# CONVERSIONS

# imperial to metric
GR2GM   = 0.0647989  # grains to grams
GR2MG   = 64.7989    # grains to milligrams
GR2KG   = 6.47989e-5 # grains to kilograms
IN2MM   = 25.4       # inches to millimeter
IN2CM   = 2.54       # inches to centimeters

# caliber to metric
CL2M = 2.54e-4

# metric conversions
MM2M    = 0.001 # mm to m
MS2MMS  = 1e6   # m^2 to mm^2
MMS2MS  = 1e-6  # mm^2 to m^2
MC2MMC  = 1e9   # m^3 to mm^3
MMC2MC  = 1e-9  # mm^3 to m^3
KG2G    = 1000  # kilogram to gram

# metric to imperial
SD2IMP  = 1/703  # sectional density (SD) & ballistic coefficent (BC) to imperial
CSA2IMP = 1/1550 # cross-sectional area (A) to imperial