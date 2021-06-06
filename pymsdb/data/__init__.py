#!/usr/bin/env python
""" data
Copyright (C) 2020 Dale V. Patterson (dale.v.patterson@gmail.com)

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

defines paths to data resources IOT find paths regardless of
 download/installation path
"""

#__name__ = 'data'
__license__ = 'GPLv3'
__version__ = '0.0.1'
__date__ = 'May 2021'
__author__ = 'Dale V. Patterson'
__maintainer__ = 'Dale V. Patterson'
__email__ = 'dale.v.patterson@gmail.com'
__status__ = 'Development'

import os

pth_msdb = os.path.join(os.path.dirname(__file__),'msdb.tsv')
pth_firearmspec = os.path.join(os.path.dirname(__file__),'firearm_specs.tsv')
pth_caliber = os.path.join(os.path.dirname(__file__),'caliber.tsv')
pth_states = os.path.join(os.path.dirname(__file__),'states.tsv')
pth_clog = os.path.join(os.path.dirname(__file__),'consistency.tsv')