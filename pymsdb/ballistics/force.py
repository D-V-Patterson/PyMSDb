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
  1. in Fd and mFd add option to compute rho via atmospheric functions (in atm)
   based off of blt's current elevation although I don't think this will play
   much roll
    i. this implies the need to know the shooter's starting elevation as well 
  2. Forces:
   a. Magnus effect - 
      i. need p = (Ix/Iy) * (h dot x) don't have h or x but we do have the spin
       rate as calculated from Lahti. Will this suffice
      ii. there is a magnus force coefficient and a magnus momment coefficient,
      McCoy 1999, pg 189 See McCoy 1964, pg 20 a Magnus moment will be produced 
      if the line of action of the magnus force does not pass through the center 
      of mass of projectile
   b. As for Magnus there is a pitch damping force and pitch damping moment where
    Cnq + Cna = pitch damping force coefficient and
    Cmq + Cma = pitch damping moment coefficient 
   c. spin damping force requires spin damping coefficient
   d. rolling moment force requires rolling moment coefficient
"""

#### GRAVITY

def gravity():
	""" returns a xyz vector of gravity's affect given angle """
	return -np.array([0.,atm.G,0.],np.double)

#### DRAG

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

def mFd(blt,rho=atm.RHO):
	"""
	 calculates the force of drag (Fd) magnitude given the air density and the
	 bullet's drag coefficient, cross-sectional area, velocity, mass
	:param blt: the Bullet object
	:param rho: air density (kg/m^3)
	:return: magnitude of the force of drag
	uses Fd = (-0.5 * Cd * rho * S) / m
	 McCoy 1999, Eq 5.7 (see also pg 190)
	 where
	  Cd = coefficient of drag,
	  rho = air density,
	  S = cross-sectional area, and
	  m = mass
	"""
	# TODO:  see mFl and others, conform this (include velocity) or should we
	#  conform below to this and multiply velocity in bullet?
	return (0.5*blt.Cd()*rho*blt.A)/blt.mass

#### LIFT

def mFl(blt,rho=atm.RHO):
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

#### MAGNUS FORCE

def mFm(blt,rho=atm.RHO):
	"""
	 calculates the Magnus Force magnitude given the air density and the bullet's
	 magnus force coefficient, velocity, cross-sectional area, diameter and axial
	 spin rate
	:param blt: Bullet object
	:param rho: air density
	:return: the magnus force
	Uses McCoy 1999, eq 9.10 - 9.16 pg 190
	 Fm = rho * S * d * Cm * p / 2 * m
	  where
	   rho = air density (kg/m^3),
	   S = cross-sectional area (m^2),
	   d = diameter (m),
	   Cn = magnus force coefficient,
	   p = axial spin rate (radians/s), and
	   m = mass
	"""
	# TODO: using the spin rate from Lahti vice McCoy (Ix/Iy) * (h dot x)
	return (0.5*rho*blt.A*blt.d*blt.Cmm()*blt.p_i) / blt.mass

# TODO: magnus moment, will need magnus moment coefficient

#### PITCH DAMPING

# TODO: see Lahit, McCoy the below Cmp() returns the pitch damping moment
#  coefficients and not the pitch damping force coefficient
#def mFpd(blt,rho=atm.RHO):
#	"""
#	 calculates the pitch damping force given the air density and the bullet's
#	 pitch damping coefficients, velocity, cross-sectional area, diameter and
#	 mass
#	:param blt: Bullet object
#	:param rho: air density
#	:return: pitch damping force
#	 Uses McCoy 1999, eq 9.10 - 9.16 pg 190
#	  Fpd = rho * v * S * d * (Cnq + Cna) / 2 * ma
#	  where
#	   rho = air density (kg/m^3),
#	   v = velocity (m/s),
#	   S = cross-sectional area (m^2),
#	   d = diameter (m),
#	   Cnq, Cna = pitch damping coefficients, and
#	   m = mass
#	"""
#	return (0.5*rho*blt.v_i*blt.A*blt.d*sum(blt.Cmp())) / blt.mass


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
