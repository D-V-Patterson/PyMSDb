#!/usr/bin/env python
""" PyMSDb v 0.1.0: A Python 3.x MSDB tool
Copyright (C) 2020 Dale V. Patterson (dale.v.patterson@gmail.com)

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

Top level directory for PyMDSb program files
 TODO

Requirements
 - Python 3.x https://www.python.org

External Dependencies
 - regex 2.5.93
    https://bitbucket.org/mrabarnett/mrab-regex
    pip3 install regex
 - numpy 1.20.3
    https://numpy.org
    pip3 install numpy
 - scipy 1.6.3
    https://www.scipy.org
    pip3 install scipy
 - matplotlib 3.4.2
    https://matplotlib.org
    pip3 install matplotlib

Usage
 - add a .pth file to /usr/lib/python3/dist-packages with the full path of PyMSDb

DO NOT IMPORT *
"""

#__name__ = 'rkba'
__license__ = 'GPLv3'
__version__ = '0.1.0'
__date__ = 'May 2021'
__author__ = 'Dale V. Patterson'
__maintainer__ = 'Dale V. Patterson'
__email__ = 'dale.v.patterson@gmail.com'
__status__ = 'Development'

class PyMSDbException(Exception):
    def __init__(self,msg): super().__init__(msg)

"""
POINT MASS
from importlib import reload
import numpy as np
import pymsdb.ballistics as bls
import pymsdb.utils as utils
import pymsdb.ballistics.model as model
import pymsdb.ballistics.force as force
import pymsdb.ballistics.atmosphere as atm
import pymsdb.caliber as caliber
import pymsdb.ballistics.bullet as bullet

reload(caliber)
reload(bullet)
cals = caliber.CaliberDS()
cal1=cals.get_caliber('7.62x39')
blt1=cal1.blt_idx(0)
cal2=cals.get_caliber('5.56x45')
blt2=cal2.blt_idx(0)
cal4=cals.get_caliber('9x19')
blt4=cal4.blt_idx(3)

"""
