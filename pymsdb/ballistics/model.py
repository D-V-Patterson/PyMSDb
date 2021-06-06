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
import pymsdb.ballistics.atmosphere as atm
from pymsdb import PyMSDbException

"""
 TODO:
  1. tangent vs secant ogival volumes/centers of gravity
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

def geometry(blt):
    """
     estimates bullet measurements based on bullet's length & std model ratios
    :param blt: Bullet object
    :return: tuple t = (
        ttl-volume (mm^3),
        base-diameter (mm),
        ogive-len (mm),
        cylinder-len (mm),
        frustum-len (mm)
        )
    NOTE;
     bullet must have diameter, bullet length and bullet model configured
    """
   # from the bullet length & diameter, estimate other measurements
    bl = blt.blt_len
    d = blt.d
    mspec = mdl_specs[blt.mdl]

    # estimate lengths of ogive, cylinder, and frustum (if present)
    lo = np.round(mspec['ogv-len'] * bl,3)
    lc = np.round(mspec['cyl-len'] * bl,3)
    lf = np.round(mspec['frs-len'] * bl,3) if mspec['frs-len'] else 0.
    db = np.round(mspec['base-dia'] * d,3)

    # we need the radius of the cylinder and frustum
    rc = d / 2  # cylinder radius
    rf = db / 2 # frustum radius

    # estimate the volume of the bullet by summing the volume of its parts
    # for cylinder and frustum (google equations)
    #  Vc = pi*r^2*h where r=cylinder radius and h = cylinder length
    #  Vf = 1/3 * pi * (h*r1^2 + (r1*r2) + r2^2)
    #   where
    #    h = frustum length
    #    r1 = radius of frustum base
    #    r2 = radius of frustum top
    vc = np.round(np.pi*np.power(rc,2)*lc,3)
    vf = (np.pi*lf)/3 * (np.power(rc,2) + (rc*rf) + np.power(rf,2))

    # the ogival volume is more complicated
    # see http://www.if.sc.usp.br/~projetosulfos/artigos/NoseCone_EQN2.PDF
    # (use wayback machine 11-Apr-2001 capture to download it) pg 12
    # calculate p and p^2
    #  p = (R^2 + L^2)/2R
    #   where
    #    R = radius of the base of the ogive (use rc)
    #    L = length of the ogive (lo)
    #
    # and then for the volume
    #  Vo = pi*(L*p^2) - L^3/3 - (p-R)p^2sin^-1(L/p))
    #  where
    #   p is calculated above,
    #   L = ogive length (lo)
    #   R = ogive base radius (rc)
    # TODO: is this different for secant vice tanget ogive
    p = (np.power(rc,2)+np.power(lo,2)) / (2*rc)
    p2 = np.power(p,2) # calculate once
    vo = np.pi*(lo*p2 - np.power(lo,3)/3 - (p-rc)*p2*np.arcsin(lo/p))

    # add the constituent volumes
    vol = np.round(vo + vc + vf,3)

    # estimate the center of gravity
    # Requires the cog for the cylinder, the frustum and the ogive as well as
    # the weights of each of these
    # CGcylinder = l/2
    #  where
    #   l = length of cylinder
    CGc = lc / 2

    # CGfrustum = l/4 * (rf^2 + 2*(rf*rc) + 3*rc^2)/(rf^2+rf*rc + rc^2))
    #  where
    #   l = length of frustum (or height),
    #   rf = radius of frustum (top of frustum), and
    #   rc = radius of cylinder (bottom of frustum)
    rf2 = np.power(rf,2)
    rc2 = np.power(rc,2)
    CGf = (lf/4) * (rf2 + 2*rf*rc + 3*rc2) / (rf2 + rf*rc + rc2)

    # TODO:
    #  calculating ogival cog requires integration
    #  y(x) = sqrt(r^2 - x^2) + yi 0 <= x <= lo
    #  r = diameter/4 + length^2/diameter
    #  yi = diameter/2 - r
    #  then if the ogive has constant density,
    #   y(x)(x-xcg)dx = 0 and solve for xcg
    # We will estimate using this equation from a forum
    # (http://www.craftsclassic.com/7_928f8ada4e2c3123_1.htm)
    CGo = (0.3479 + 0.0122*np.power(d/lo,2)) * lo

    # estimate the weights (we are using weight and mass interchangeably) -
    # requires the mean density of the bullet: density = mass / volume
    # w(cyl) = pi * rc^2* lc * density
    # w(frs) = 1/3 * pi * lf * (rb^2 + rb*r * r^2) * density
    # w(ogv) = w(bullet) - w(cyl) - w(frs)
    #  where
    #   rc = the radius of the bullet,
    #   lc = length of the cylinder,
    #   lf = length of the frustum, and
    #   rb = radius of the base (top of frustum)
    p = blt.mass / vol
    wc = (np.pi * rc2 * lc) * p
    wf = (1/3 * np.pi * lf * (rf2 + rf * rc + rc2)) * p
    wo = blt.mass - wc - wf

    # using https://www.grc.nasa.gov/www/k-12/rocket/rktcg.html
    cog = np.round((wf*CGf + wc*(CGc + lf) + wo*(CGo + lf + lc))/blt.mass,3)

    return vol,db,cog,lo,lc,lf

"""
 defines standard model ratios to define bullet specs
 each model is a dict with keys
  ogv-len = ratio of ogive-len/bullet-len
  cyl-len = the ratio of cylinder-len/blt-len
  frs-len = the ratio of frustrum-len/blt-len (if present)
  rrs-dia = diameter of the 'top' of the frustum
  measurements are from 
   https://www.alternatewars.com/BBOW/Ballistics/Ext/Drag_Tables.htm
  TODO: use the above to get G2 and other reference numbers 
"""

# measurements in calibers (Cl) 1Cl = 0.254mm of each standard model bullet
STD_OAL = 0 # overal length (OAL)
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
# TODO: add ogive and frustum degrees
mdl_specs = {
    'G1':{
        'ogv-len':np.round(G1_SHPE[STD_LO]/G1_SHPE[STD_OAL],3),
        'cyl-len':np.round(G1_SHPE[STD_LC]/G1_SHPE[STD_OAL],3),
        'frs-len':np.round(G1_SHPE[STD_LF]/G1_SHPE[STD_OAL],3),
        'base-dia':G1_SHPE[STD_DB],
    },
    'G2':{
        'ogv-len':np.round(G2_SHPE[STD_LO]/G2_SHPE[STD_OAL],3),
        'cyl-len':np.round(G2_SHPE[STD_LC]/G2_SHPE[STD_OAL],3),
        'frs-len':np.round(G2_SHPE[STD_LF]/G2_SHPE[STD_OAL],3),
        'base-dia':G2_SHPE[STD_DB],
    },
    'G5':{
        'ogv-len':np.round(G5_SHPE[STD_LO]/G5_SHPE[STD_OAL],3),
        'cyl-len':np.round(G5_SHPE[STD_LC]/G5_SHPE[STD_OAL],3),
        'frs-len':np.round(G5_SHPE[STD_LF]/G5_SHPE[STD_OAL],3),
        'base-dia':G5_SHPE[STD_DB],
    },
    'G6':{
        'ogv-len':np.round(G6_SHPE[STD_LO]/G6_SHPE[STD_OAL],3),
        'cyl-len':np.round(G6_SHPE[STD_LC]/G6_SHPE[STD_OAL],3),
        'frs-len':np.round(G6_SHPE[STD_LF]/G6_SHPE[STD_OAL],3),
        'base-dia':G6_SHPE[STD_DB],
    },
    'G7':{
        'ogv-len':np.round(G7_SHPE[STD_LO]/G7_SHPE[STD_OAL],3),
        'cyl-len':np.round(G7_SHPE[STD_LC]/G7_SHPE[STD_OAL],3),
        'frs-len':np.round(G7_SHPE[STD_LF]/G7_SHPE[STD_OAL],3),
        'base-dia':G7_SHPE[STD_DB],
    },
    'G8':{
        'ogv-len':np.round(G8_SHPE[STD_LO]/G8_SHPE[STD_OAL],3),
        'cyl-len':np.round(G8_SHPE[STD_LC]/G8_SHPE[STD_OAL],3),
        'frs-len':np.round(G8_SHPE[STD_LF]/G8_SHPE[STD_OAL],3),
        'base-dia':G8_SHPE[STD_DB],
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
    if blt.mdl == 'G1': return g1_std_cd(atm.mach(blt.v_i))
    elif blt.mdl == 'G7': return g7_std_cd(atm.mach(blt.v_i))
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

