#!/usr/bin/env python
"""  force.py
Copyright (C) 2021 Dale Patterson (dale.v.patterson@gmail.com)

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

Defines functions related to forces 'opposing' the bullet's trajectory
"""

#__name__ = 'drag'
__license__ = 'GPLv3'
__version__ = '0.0.7'
__date__ = 'May 2021'
__author__ = 'Dale Patterson'
__maintainer__ = 'Dale Patterson'
__email__ = 'dale.v.patterson@gmail.com'
__status__ = 'Development'

import numpy as np
import pymsdb.ballistics.atmosphere as atm

"""
 TODO:
  1. in Fd and mFd add option to compute rho via atmospheric functions (in atm)
   based off of blt's current elevation although I don't think this will play
   much roll
    i. this implies the need to know the shooter's starting elevation as well 
"""

#### GRAVITY

def gravity():
	""" returns a xyz vector of gravity's affect given angle """
	return -np.array([0.,atm.G,0.],np.double)

#### DRAG RELATED

def Fd(blt,rho=atm.RHO):
	"""
	 calculates the force of drag (Fd) given projectiles drag coefficient,
	 cross-sectional area, velocity and air density
	:param blt: Bullet object
	:param rho: air density (kg/m3)
	:return: numpy array of the 3d components of the Fd = Force of Drag
	uses Fd = (-0.5 * Cd * V_x * vv * rho * A) / mass
	 McCoy 1999, Eq 5.3, 5.7
	 where
	  Cd = coefficient of drag,
	  V_x = speed at time i of x-component of velocity
	  vv = xyz components of velocity at time i
	  rho = air density
	  A = cross-sectional area
	  Fd in (N or kg⋅m/s^2)
	"""
	return mFd(blt,rho)*blt.velocity*blt.v_i
	#return -(0.5*blt.Cd()*blt.velocity*blt.v_i*rho*blt.A)/blt.mass

def mFd(blt,rho=atm.RHO):
	"""
	 calculates the force of drag (Fd) magnitude given projectiles drag coefficient,
	 cross-sectional area, velocity and air density
	:param blt: the Bullet object
	:param rho: air density (kg/m3)
	:return: magnitude of the force of drag
	uses Fd = (-0.5 * Cd * rho * A) / mass
	 McCoy 1999, Eq 5.7
	 where
	  Cd = coefficient of drag,
	  rho = air density
	  A = cross-sectional area
	  Fd in (m I think)
	"""
	return -(0.5*blt.Cd()*rho*blt.A)/blt.mass


