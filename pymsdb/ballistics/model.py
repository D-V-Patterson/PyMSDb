#!/usr/bin/env python
"""  model.py
Copyright (C) 2021 Dale Patterson (dale.v.patterson@gmail.com)

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

Defines drag models and associated constants/functions
"""

#__name__ = 'model'
__license__ = 'GPLv3'
__version__ = '0.0.1'
__date__ = 'May 2021'
__author__ = 'Dale Patterson'
__maintainer__ = 'Dale Patterson'
__email__ = 'dale.v.patterson@gmail.com'
__status__ = 'Development'

import numpy as np
from scipy.interpolate import interp1d
from pymsdb import PyMSDbException

"""
 TODO:
  1. tangent vs secant ogival volumes/centers of gravity
   a. G1 ogival radius is calculated with tangent ogive formula (are G1 tangent?)
     see Crowell 1996 for conical nose cone definition
"""

# MODELS

"""
list of the drag models (NOTE: only G1 and G7 are currently used)
 Ingalls/G1 flatbase bullet short ogive
 Aberdeen/G2 conical ogive 
 G5 tangent ogive cylinder boattail similar to G7 with short ogival radius
 G6 flat-based secant ogive cylinder
 G7 tangent ogive cylinder with long boattail
 G8 flat based secant ogive cylinder
 GL blunt nose NOT USED 
"""
models = ['Ingalls','G1','G2','G5','G6','G7','G8','GL']
imodels = ['G1','G7'] # implemented models

"""
 defines standard model ratios to define bullet specs
 each model is a dict with keys
  ogv-shp = ogival shape one of {conical,tangent,secant,blunt}
  ogv-len = ratio of ogive-len/bullet-len
  ogv-r = ogival radius
  cyl-len = the ratio of cylinder-len/blt-len
  frs-shp = one of {boattail,none}
  frs-len = the ratio of frustrum-len/blt-len (if present)
  frs-base = diameter of the 'bottom' of the frustum
  measurements are from 
   https://www.alternatewars.com/BBOW/Ballistics/Ext/Drag_Tables.htm 
"""

# measurements in calibers (Cl) 1 Cl = 0.254mm of each standard model bullet
OGV_CON = 0 # conical ogive
OGV_TAN = 1 # tangent ogive
OGV_SEC = 2 # secant ogive
OGV_BLT = 3 # blunt nosed
FRS_NNE = 0 # no frustum
FRS_BTL = 1 # boattail
STD_OAL = 0 # overal length (OABL)
STD_LO  = 1 # length of ogive
STD_LT  = 2 # length of ogive throat (0. if not present) NOTE: incl in ogive len
STD_LC  = 3 # length of cylinder
STD_LF  = 4 # length of frustum (0. if not present)
STD_D   = 5 # diameter
STD_DB  = 6 # diameter of base, at frustum (0. if not present)
STD_MB  = 7 # meplat diameter (0. if not present)
STD_RO  = 8 # radius of ogive
STD_TF  = 9 # taper (degrees) of frustum (0. if not present
G1_SHPE = np.array([3.28,1.32,0.,1.98,0.,1.,1.,0.,2.,0.],np.double)
G2_SHPE = np.array([5.19,2.7,0.41,1.99,0.5,1.,0.8949,0.12,2.34,6.],np.double)
G5_SHPE = np.array([4.29,2.1,0.,1.7,0.49,1.,0.842,0.,6.19,7.5],np.double)
G6_SHPE = np.array([4.81,2.53,0.,2.28,0.,1.,1.,0.,6.99,0.],np.double)
G7_SHPE = np.array([4.23,2.18,0.,1.45,0.6,1.,0.842,0.,10.,7.5],np.double)
G8_SHPE = np.array([3.64,2.18,0.,1.46,0.,1.,1.,0.,10.,0.],np.double)
mdl_specs = {
    'G1':{
        'ogv-shp':OGV_TAN,
        'ogv-len':G1_SHPE[STD_LO]/G1_SHPE[STD_OAL],
        'ogv-r':G1_SHPE[STD_RO]/G1_SHPE[STD_LO],
        'cyl-len':G1_SHPE[STD_LC]/G1_SHPE[STD_OAL],
        'frs-shp':FRS_NNE,
        'frs-len':G1_SHPE[STD_LF]/G1_SHPE[STD_OAL],
        'frs-base':G1_SHPE[STD_DB],
    },
    'G2':{
        'ogv-shp':OGV_CON,
        'ogv-len':G2_SHPE[STD_LO]/G2_SHPE[STD_OAL],
        'ogv-r':G2_SHPE[STD_RO]/G2_SHPE[STD_LO],
        'cyl-len':G2_SHPE[STD_LC]/G2_SHPE[STD_OAL],
        'frs-shp':FRS_BTL,
        'frs-len':G2_SHPE[STD_LF]/G2_SHPE[STD_OAL],
        'frs-base':G2_SHPE[STD_DB],
    },
    'G5':{
        'ogv-shp':OGV_TAN,
        'ogv-len':G5_SHPE[STD_LO]/G5_SHPE[STD_OAL],
        'ogv-r': G5_SHPE[STD_RO] / G5_SHPE[STD_LO],
        'cyl-len':G5_SHPE[STD_LC]/G5_SHPE[STD_OAL],
        'frs-shp':FRS_BTL,
        'frs-len':G5_SHPE[STD_LF]/G5_SHPE[STD_OAL],
        'frs-base':G5_SHPE[STD_DB],
    },
    'G6':{
        'ogv-shp':OGV_SEC,
        'ogv-len':G6_SHPE[STD_LO]/G6_SHPE[STD_OAL],
        'ogv-r':G6_SHPE[STD_RO]/G6_SHPE[STD_LO],
        'cyl-len':G6_SHPE[STD_LC]/G6_SHPE[STD_OAL],
        'frs-shp':FRS_NNE,
        'frs-len':G6_SHPE[STD_LF]/G6_SHPE[STD_OAL],
        'frs-base':G6_SHPE[STD_DB],
    },
    'G7':{
        'ogv-shp':OGV_TAN,
        'ogv-len':G7_SHPE[STD_LO]/G7_SHPE[STD_OAL],
        'ogv-r':G7_SHPE[STD_RO]/G7_SHPE[STD_LO],
        'cyl-len':G7_SHPE[STD_LC]/G7_SHPE[STD_OAL],
        'frs-shp':FRS_BTL,
        'frs-len':G7_SHPE[STD_LF]/G7_SHPE[STD_OAL],
        'frs-base':G7_SHPE[STD_DB],
    },
    'G8':{
        'ogv-shp':OGV_SEC,
        'ogv-len':G8_SHPE[STD_LO]/G8_SHPE[STD_OAL],
        'ogv-r':G8_SHPE[STD_RO]/G8_SHPE[STD_LO],
        'cyl-len':G8_SHPE[STD_LC]/G8_SHPE[STD_OAL],
        'frs-shp':FRS_NNE,
        'frs-len':G8_SHPE[STD_LF]/G8_SHPE[STD_OAL],
        'frs-base':G8_SHPE[STD_DB],
    },
}

#### STD DRAG MODELS COEFFICIENT OF DRAG
# each model Cd is define as two arrays Velocity (mach number) and Cd which
# are interpolated to define the standard model's Cd at given velocity
#  data is from JBM ballistics

def std_cd(blt):
    """
	 determines the coefficient of drag for the 'standard' model bullet
	:param blt: Bullet object
	:return: standard model Cd at given velocity
	"""
    if blt.mdl == 'G1': return g1_std_cd(blt.mach())
    elif blt.mdl == 'G2': return g2_std_cd(blt.mach())
    elif blt.mdl == 'G5': return g5_std_cd(blt.mach())
    elif blt.mdl == 'G6': return g6_std_cd(blt.mach())
    elif blt.mdl == 'G7': return g7_std_cd(blt.mach())
    elif blt.mdl == 'G8': return g8_std_cd(blt.mach())
    else: raise PyMSDbException("Invalid 'drag model' ({})".format(blt.mdl))

#  G1 https://www.jbmballistics.com/ballistics/downloads/text/mcg1.txt
g1_std_cd_v = np.array( # mach number
	[0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.7,
	 0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 1.0,
	 1.025, 1.05, 1.075, 1.1, 1.125, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5,
	 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2.0, 2.05, 2.1, 2.15, 2.2,
	 2.25, 2.3, 2.35, 2.4, 2.45, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4,
	 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0],np.double
)
g1_std_cd_cd = np.array( # Cd
	[0.2629, 0.2558, 0.2487, 0.2413, 0.2344, 0.2278, 0.2214, 0.2155, 0.2104, 0.2061,
	 0.2032, 0.202, 0.2034, 0.2165, 0.223, 0.2313, 0.2417, 0.2546, 0.2706, 0.2901,
	 0.3136, 0.3415, 0.3734, 0.4084, 0.4448, 0.4805, 0.5136, 0.5427, 0.5677, 0.5883,
	 0.6053, 0.6191, 0.6393, 0.6518, 0.6589, 0.6621, 0.6625, 0.6607, 0.6573, 0.6528,
	 0.6474, 0.6413, 0.6347, 0.628, 0.621, 0.6141, 0.6072, 0.6003, 0.5934, 0.5867,
	 0.5804, 0.5743, 0.5685, 0.563, 0.5577, 0.5527, 0.5481, 0.5438, 0.5397, 0.5325,
	 0.5264, 0.5211, 0.5168, 0.5133, 0.5105, 0.5084, 0.5067, 0.5054, 0.504, 0.503,
	 0.5022, 0.5016, 0.501, 0.5006, 0.4998, 0.4995, 0.4992, 0.499, 0.4988],np.double
)
g1_std_cd = interp1d(g1_std_cd_v,g1_std_cd_cd,kind='linear')

# https://www.jbmballistics.com/ballistics/downloads/text/mcg2.txt
g2_std_cd_v = np.array( # mach number
    [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65,
     0.7, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 1.0,
     1.025, 1.05, 1.075, 1.1, 1.125, 1.15, 1.175, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45,
     1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2.0, 2.05, 2.1, 2.15,
     2.2, 2.25, 2.3, 2.35, 2.4, 2.45, 2.5, 2.55, 2.6, 2.65, 2.7, 2.75, 2.8, 2.85,
     2.9, 2.95, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.2, 4.4,
     4.6, 4.8, 5.0],np.double
)
g2_std_cd_cd = np.array( #Cd
    [0.2303, 0.2298, 0.2287, 0.2271, 0.2251, 0.2227, 0.2196, 0.2156, 0.2107, 0.2048,
     0.198, 0.1905, 0.1828, 0.1758, 0.1702, 0.1669, 0.1664, 0.1667, 0.1682, 0.1711,
     0.1761, 0.1831, 0.2004, 0.2589, 0.3492, 0.3983, 0.4075, 0.4103, 0.4114, 0.4106,
     0.4089, 0.4068, 0.4046, 0.4021, 0.3966, 0.3904, 0.3835, 0.3759, 0.3678, 0.3594,
     0.3512, 0.3432, 0.3356, 0.3282, 0.3213, 0.3149, 0.3089, 0.3033, 0.2982, 0.2933,
     0.2889, 0.2846, 0.2806, 0.2768, 0.2731, 0.2696, 0.2663, 0.2632, 0.2602, 0.2572,
     0.2543, 0.2515, 0.2487, 0.246, 0.2433, 0.2408, 0.2382, 0.2357, 0.2333, 0.2309,
     0.2262, 0.2217, 0.2173, 0.2132, 0.2091, 0.2052, 0.2014, 0.1978, 0.1944, 0.1912,
     0.1851, 0.1794, 0.1741, 0.1693, 0.1648],np.double
)
g2_std_cd = interp1d(g2_std_cd_v,g2_std_cd_cd,kind='linear')

# https://www.jbmballistics.com/ballistics/downloads/text/mcg5.txt
g5_std_cd_v = np.array( # mach number
    [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65,
     0.7, 0.75, 0.8, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 1.0, 1.025, 1.05, 1.075,
     1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75,
     1.8, 1.85, 1.9, 1.95, 2.0, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4, 2.45,
     2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9,
     4.0, 4.2, 4.4, 4.6, 4.8, 5.0],np.double
)
g5_std_cd_cd = np.array( # CD
    [0.171, 0.1719, 0.1727, 0.1732, 0.1734, 0.173, 0.1718, 0.1696, 0.1668, 0.1637,
     0.1603, 0.1566, 0.1529, 0.1497, 0.1473, 0.1463, 0.1489, 0.1583, 0.1672, 0.1815,
     0.2051, 0.2413, 0.2884, 0.3379, 0.3785, 0.4032, 0.4147, 0.4201, 0.4278, 0.4338,
     0.4373, 0.4392, 0.4403, 0.4406, 0.4401, 0.4386, 0.4362, 0.4328, 0.4286, 0.4237,
     0.4182, 0.4121, 0.4057, 0.3991, 0.3926, 0.3861, 0.38, 0.3741, 0.3684, 0.363,
     0.3578, 0.3529, 0.3481, 0.3435, 0.3391, 0.3349, 0.3269, 0.3194, 0.3125, 0.306,
     0.2999, 0.2942, 0.2889, 0.2838, 0.279, 0.2745, 0.2703, 0.2662, 0.2624, 0.2588,
     0.2553, 0.2488, 0.2429, 0.2376, 0.2326, 0.228],np.double
)
g5_std_cd = interp1d(g5_std_cd_v,g5_std_cd_cd,kind='linear')

# https://www.jbmballistics.com/ballistics/downloads/text/mcg6.txt
g6_std_cd_v = np.array( # Mach
    [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65,
     0.7, 0.75, 0.8, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 1.0, 1.025, 1.05, 1.075,
     1.1, 1.125, 1.15, 1.175, 1.2, 1.225, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55,
     1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2.0, 2.05, 2.1, 2.15, 2.2, 2.25,
     2.3, 2.35, 2.4, 2.45, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5,
     3.6, 3.7, 3.8, 3.9, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0],np.double
)
g6_std_cd_cd = np.array( # CD
    [0.2617, 0.2553, 0.2491, 0.2432, 0.2376, 0.2324, 0.2278, 0.2238, 0.2205, 0.2177,
     0.2155, 0.2138, 0.2126, 0.2121, 0.2122, 0.2132, 0.2154, 0.2194, 0.2229, 0.2297,
     0.2449, 0.2732, 0.3141, 0.3597, 0.3994, 0.4261, 0.4402, 0.4465, 0.449, 0.4497,
     0.4494, 0.4482, 0.4464, 0.4441, 0.439, 0.4336, 0.4279, 0.4221, 0.4162, 0.4102,
     0.4042, 0.3981, 0.3919, 0.3855, 0.3788, 0.3721, 0.3652, 0.3583, 0.3515, 0.3447,
     0.3381, 0.3314, 0.3249, 0.3185, 0.3122, 0.306, 0.3, 0.2941, 0.2883, 0.2772,
     0.2668, 0.2574, 0.2487, 0.2407, 0.2333, 0.2265, 0.2202, 0.2144, 0.2089, 0.2039,
     0.1991, 0.1947, 0.1905, 0.1866, 0.1794, 0.173, 0.1673, 0.1621, 0.1574],np.double
)
g6_std_cd = interp1d(g5_std_cd_v,g5_std_cd_cd,kind='linear')

# https://www.jbmballistics.com/ballistics/downloads/text/mcg7.txt
g7_std_cd_v = np.array( # mach number
	[0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65,
	 0.7, 0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975,
	 1.0, 1.025, 1.05, 1.075, 1.1, 1.125, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.5,
	 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2.0, 2.05, 2.1, 2.15, 2.2,
	 2.25, 2.3, 2.35, 2.4, 2.45, 2.5, 2.55, 2.6, 2.65, 2.7, 2.75, 2.8, 2.85, 2.9,
	 2.95, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.2, 4.4, 4.6,
	 4.8, 5.0],np.double
)
g7_std_cd_cd = np.array( # Cd
	[0.1198, 0.1197, 0.1196, 0.1194, 0.1193, 0.1194, 0.1194, 0.1194, 0.1193,
	 0.1193, 0.1194, 0.1193, 0.1194, 0.1197, 0.1202, 0.1207, 0.1215, 0.1226,
	 0.1242, 0.1266, 0.1306, 0.1368, 0.1464, 0.166, 0.2054, 0.2993, 0.3803, 0.4015,
	 0.4043, 0.4034, 0.4014, 0.3987, 0.3955, 0.3884, 0.381, 0.3732, 0.3657, 0.358,
	 0.344, 0.3376, 0.3315, 0.326, 0.3209, 0.316, 0.3117, 0.3078, 0.3042, 0.301,
	 0.298, 0.2951, 0.2922, 0.2892, 0.2864, 0.2835, 0.2807, 0.2779, 0.2752, 0.2725,
	 0.2697, 0.267, 0.2643, 0.2615, 0.2588, 0.2561, 0.2533, 0.2506, 0.2479, 0.2451,
	 0.2424, 0.2368, 0.2313, 0.2258, 0.2205, 0.2154, 0.2106, 0.206, 0.2017, 0.1975,
	 0.1935, 0.1861, 0.1793, 0.173, 0.1672, 0.1618],np.double
)
g7_std_cd = interp1d(g7_std_cd_v,g7_std_cd_cd,kind='linear')

# https://www.jbmballistics.com/ballistics/downloads/text/mcg8.txt
g8_std_cd_v = np.array( # Mach
    [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65,
     0.7, 0.75, 0.8, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 1.0, 1.025, 1.05,
     1.075, 1.1, 1.125, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65,
     1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2.0, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35,
     2.4, 2.45, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7,
     3.8, 3.9, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0],np.double
)
g8_std_cd_cd = np.array(
    [0.2105, 0.2105, 0.2104, 0.2104, 0.2103, 0.2103, 0.2103, 0.2103, 0.2103,
     0.2102, 0.2102, 0.2102, 0.2102, 0.2102, 0.2103, 0.2103, 0.2104, 0.2104,
     0.2105, 0.2106, 0.2109, 0.2183, 0.2571, 0.3358, 0.4068, 0.4378, 0.4476,
     0.4493, 0.4477, 0.445, 0.4419, 0.4353, 0.4283, 0.4208, 0.4133, 0.4059,
     0.3986, 0.3915, 0.3845, 0.3777, 0.371, 0.3645, 0.3581, 0.3519, 0.3458,
     0.34, 0.3343, 0.3288, 0.3234, 0.3182, 0.3131, 0.3081, 0.3032, 0.2983,
     0.2937, 0.2891, 0.2845, 0.2802, 0.272, 0.2642, 0.2569, 0.2499, 0.2432,
     0.2368, 0.2308, 0.2251, 0.2197, 0.2147, 0.2101, 0.2058, 0.2019, 0.1983,
     0.195, 0.189, 0.1837, 0.1791, 0.175, 0.1713],np.double
)
g8_std_cd = interp1d(g8_std_cd_v,g8_std_cd_cd,kind='linear')