#!/usr/bin/env python
""" anthropometry
Copyright (C) 2020 Dale V. Patterson (dale.v.patterson@gmail.com)

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

Defines constants for body measurements using below sources
 1. https://multisite.eos.ncsu.edu/www-ergocenter-ncsu-edu/wp-content/uploads/sites/
    18/2016/06/Anthropometric-Detailed-Data-Tables.pdf
 2. https://www.researchgate.net/figure/Chest-Breadth-measurement_fig15_273134617
"""

#__name__ = 'anthropometry'
__license__ = 'GPLv3'
__version__ = '0.0.1'
__date__ = 'May 2021'
__author__ = 'Dale V. Patterson'
__maintainer__ = 'Dale V. Patterson'
__email__ = 'dale.v.patterson@gmail.com'
__status__ = 'Production'

# Anthropometric Data

S_H    = 1.44  # estimated height of shoulder fired weapon (of a male) (m))
W_HM   = 1.13  # estimated height of waist on male (m)
W_HF   = 1.06  # estimated height of waist on female (m)
T_AMH  = 0.384 # target area (shoulder to waist on male) (m)
T_AMB  = 0.290 # target area (breadth on male chest) (m)
T_AFH  = 0.352 # target area (shoulder to waist on female) (m)
T_AFB  = 0.278 # target area (breadth on female chest) (m)