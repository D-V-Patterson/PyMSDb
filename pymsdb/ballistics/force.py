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
__version__ = '0.0.8'
__date__ = 'June 2021'
__author__ = 'Dale Patterson'
__maintainer__ = 'Dale Patterson'
__email__ = 'dale.v.patterson@gmail.com'
__status__ = 'Development'

import numpy as np
import pymsdb.utils as utils
import pymsdb.ballistics.atmosphere as atm

"""
 TODO:
  1. in Fd add option to compute rho via atmospheric functions (in atm)
   based off of blt's current elevation although I don't think this will play
   much roll
    i. this implies the need to know the shooter's starting elevation as well 
  2. Forces/Moments:
   b. overturning moment (requires pitching/overturning moment coefficient)
"""

"""
 Forces have the form F<force> where <force> is lowercase letter(s)
 Moments have the form M>moment> where <moment> is lowercase letter(s)
"""

#### GRAVITY

def gravity():
	""" returns a xyz vector of gravity's affect given angle """
	return -np.array([0.,atm.G,0.],np.double)

#### DRAG

def Fd(blt,rho=atm.RHO):
	"""
	 calculates the force of drag (Fd) magnitude given the air density and the
	 bullet's drag coefficient, cross-sectional area, velocity, mass
	:param blt: the Bullet object
	:param rho: air density (kg/m^3)
	:return: magnitude of the force of drag
	uses Fd = (-0.5 * Cd * rho * S) / m
	 McCoy 1999, eq 9.10 - 9.16 pg 190
	 where
	  Cd = coefficient of drag,
	  rho = air density,
	  S = cross-sectional area, and
	  m = mass
	"""
	return (0.5*rho*blt.v_i*blt.A*blt.Cd())/blt.mass

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
	return (0.5*rho*blt.v_i*blt.A*blt.Cl()) / blt.mass

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
	return (0.5*rho*blt.A*np.power(blt.d,2)*blt.Cmm()*blt.p_i) / blt.Iy()

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
	return (0.5*rho*blt.v_i*blt.A*np.power(blt.d,2)*sum(blt.Cmp())) / blt.Iy()

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
	return (0.5*rho*blt.v_i*blt.A*blt.d*blt.Cmp()[0]) / blt.Iy()

#### CORIOLIS

def Ec(blt,lat=32.,az=90):
	"""
	calculates acceleration due to the coriolis effect
	:param blt: Bullet object
	:param lat: latitude of firer (degrees)
	:param az: azimuth of fire (degrees) due North = 0, east = 90
	:return: the coriolis acceleration (array)
	from McCoy 1999, eq 8.27 pg 178 and eq. 9.20 pg 191
	 Ec_x = 2*Omega * (-Vy*cos(L)*sin(AZ) - Vz*sin(L)
	 Ec_y = 2*Omega * (Vx*cos(L)*sin(Az) + Vz*cos(L)*cos(AZ)
	 Ec_z = 2*Omega * (Vx*sin(L) - Vy*cos(L)*cos(AZ))
	 where
	  Vx, Vy, Vz = the bullet's velocity components,
	  L = latitude of firer (degrees), and
	  AZ = azimuth of fire (degrees) due North = 0, east = 90
	"""
	x,y,z = blt.velocity # get axial components
	return 2*atm.OMEGA * np.array(
		[
			-y*utils.cos(lat)*utils.sin(az) - z*utils.sin(lat),
			x*utils.cos(lat)*utils.sin(az) + z*utils.cos(lat)*utils.cos(az),
			x*utils.sin(lat) - y*utils.cos(lat)*utils.cos(az)
		],np.double
	)
