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
import pymsdb.utils as utils
import pymsdb.ballistics as bls
import pymsdb.ballistics.atmosphere as atm

"""
 TODO:
  1. in Fd add option to compute rho via atmospheric functions (in atm)
   based off of blt's current elevation although I don't think this will play
   much roll
    i. this implies the need to know the shooter's starting elevation as well 
  2. Forces/Moments:
   b. overturning moment (requires pitching/overturning moment coefficient)
  3. See McCoy eq 9.18, do we need this?
"""

"""
 Forces have the form F<force> where <force> is lowercase letter(s)
 Moments have the form M>moment> where <moment> is lowercase letter(s)
"""

#### GRAVITY

def Fg():
	""" returns a xyz vector of gravity's affect given angle """
	return -np.array([0.,atm.G,0.],np.double)

#### DRAG

def Fd(blt,rho=atm.RHO):
	"""
	 calculates the force of drag (Fd) magnitude given projectiles drag coefficient,
	 cross-sectional area, velocity and air density
	:param blt: the Bullet object
	:param rho: air density (kg/m3)
	:return: magnitude of the force of drag
	uses McCoy 1999, Eq 9.10 - 9.16 (note add - to force of drag)
	 Fd = rho * v * S * Cd / 2 * m
	 where
	   rho = air density (kg/m^3),
	   v = velocity magnitude (m/s),
	   S = cross-sectional area (m^2),
	   Cd = coefficient of drag,
	   m = mass (kg)
	"""
	return -(rho*blt.v_t*blt.A*blt.Cd()) / (2*blt.mass)

#### LIFT

def Fl(blt,rho=atm.RHO):
	"""
	 calculates the force of lift (Fl) magnitude given the air density and the
	 bullet's lift coefficient, cross-sectional area, velocity, mass
	:param blt: Bullet object
	:param rho: air density (kg/m^3)
	:return: force of lift
	Uses McCoy 1999, eq 9.10 - 9.16 pg 190
	 Fl = rho * v * S * Cl / 2 * m
	  where
	   rho = air density (kg/m^3),
	   v = velocity magnitude (m/s),
	   S = cross-sectional area (m^2),
	   Cl = coefficient of lift,
	   m = mass (kg)
	"""
	return (rho*blt.v_t*blt.A*blt.Cl()) / (2*blt.mass)

#### MAGNUS MOMENT

def Mm(blt,rho=atm.RHO):
	"""
	 calculates the Magnus moment magnitude given the air density and the bullet's
	 magnus moment coefficient, velocity, cross-sectional area, diameter and
	 transverse moment of inertia
	:param blt: Bullet object
	:param rho: air density
	:return: the magnus force
	Uses McCoy 1999, eq 9.10 - 9.16 pg 190
	 Mm = rho * S * d^2 * Cmm * p  / 2 * Iy
	  where
	   rho = air density (kg/m^3),
	   S = cross-sectional area (m^2),
	   d = diameter (m),
	   Cmm = magnus moment coefficient,
	   p = spin rate (radians/s), and
	   Iy = transverse moment of inertia
	"""
	_1 = rho*blt.A*np.power(blt.d*bls.MM2M,2)*blt.Cmm()*blt.p_i
	return _1 / (2*blt.Iy())

#### PITCH DAMPING MOMENT

def Mpd(blt,rho=atm.RHO):
	"""
	 calculates the pitch damping moment given the air density and the bullet's
	 pitch damping coefficients, velocity, cross-sectional area, diameter and
	 mass
	:param blt: Bullet object
	:param rho: air density
	:return: pitch damping force
	 Uses McCoy 1999, eq 9.10 - 9.16 pg 190
	  Mpd = rho * v * S * d^2 * (Cmq + Cma) / 2 * Iy
	  where
	   rho = air density (kg/m^3),
	   v = velocity (m/s),
	   S = cross-sectional area (m^2),
	   d = diameter (m),
	   Cmq, Cma = pitch damping coefficients, and
	   Iy = transverse moment of inertia (radians/s)
	"""
	_1 = rho*blt.v_t*blt.A*np.power(blt.d*bls.MM2M,2)*sum(blt.Cmp())
	return _1 / (2*blt.Iy())

def Mo(blt,rho=atm.RHO):
	"""
	 calculates the overturning (aka pitching) moment given the air density and
	 the bullet's pitching moment coefficient, velocity, cross-sectional area,
	 diameter and transverse moment of inertia
	:param blt: Bullet object
	:param rho: air density
	:return: pitch damping force
	 Uses McCoy 1999, eq 9.10 - 9.16 pg 190
	  Mpd = rho * v * S * d * Cma / 2 * Iy
	  where
	   rho = air density (kg/m^3),
	   v = velocity (m/s),
	   S = cross-sectional area (m^2),
	   d = diameter (m),
	   Cma = pitching moment coefficient, and
	   Iy = transverse moment of inertia (radians/s)
	"""
	return (rho*blt.v_t*blt.A*(blt.d*bls.MM2M)*blt.Cmp()[0]) / (2*blt.Iy())

#### CORIOLIS

def Ec(blt,lat=32.,az=90):
	"""
	calculates the coriolis effect on bullet
	:param blt: Bullet object
	:param lat: latitude of firer (degrees)
	:param az: azimuth of fire (degrees) due North = 0, east = 90
	:return: the coriolis acceleration (array)
	from McCoy 1999, eq 8.27 pg 178
	 Ec_x = 2*Omega * (-Vy*cos(L)*sin(AZ) - Vz*sin(L)
	 Ec_y = 2*Omega * (Vx*cos(L)*sin(Az) + Vz*cos(L)*cos(AZ)
	 Ec_z = 2*Omega * (Vx*sin(L) - Vy*cos(L)*cos(AZ))
	 where
	  Vx, Vy, Vz = the bullet's velocity components,
	  L = latitude of firer (degrees), and
	  AZ = azimuth of fire (degrees) due North = 0, east = 90
	"""
	vx,vy,vz = blt.vV_t
	return 2*atm.OMEGA * np.array(
		[
			-vy*utils.sin(lat)*utils.sin(az) - vz*utils.sin(lat),
			vx*utils.cos(lat)*utils.sin(az) + vz*utils.cos(lat)*utils.cos(az),
			vx*utils.sin(lat) - vy*utils.cos(lat)*utils.cos(az)
		],np.double
	)
