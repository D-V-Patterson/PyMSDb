#!/usr/bin/env python
""" bullet.py
Copyright (C) 2020 Dale V. Patterson (dale.v.patterson@gmail.com)

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

Defines caliber dataset reading/writing and the Caliber class
"""

#__name__ = 'bullet'
__license__ = 'GPLv3'
__version__ = '0.1.7.1'
__date__ = 'June 2021'
__author__ = 'Dale Patterson'
__maintainer__ = 'Dale Patterson'
__email__ = 'dale.v.patterson@gmail.com'
__status__ = 'Development'

import os
import numpy as np
from scipy.interpolate import interp1d
import pymsdb.ballistics as bls
import pymsdb.ballistics.force as force
import pymsdb.ballistics.model as model
import pymsdb.ballistics.atmosphere as atm
import pymsdb.utils as utils
from pymsdb import PyMSDbException
from pymsdb.ballistics.anthropometry import S_H

class BulletException(PyMSDbException):
    def __init__(self,msg): super().__init__(msg)

"""
TODO: 
  3. take a look at Boatright
  5. in addition to tof, elevation, add velocity at range x
   velocity at range x (using k3) is off by 100 m/s (low) for 7.62 but other 
   estimates tof, elevation, zero angle and fall angle compare favorably to 
   trajectory results
  6. See above, also add range when given angle of fire
  8. Cd is currently using method that does not specify ballistic models
  9. Want range based functions vice time (requires integration?) 
 11. continue allowing return of nan in Lahti defined coefficients? or return None
  or throw error?
 18. For ke, momentum do we want to round or not?
 19. round off errors
 20. Find equation for barrel length effect on initial velocity, perhaps Le Duc
 21. Have seen sectional density (SD) defined as the (1) ratio of mass to cross-
  sectional area and (2) ratio of mass to the square of the diameter. Currently
  using 2, should this be kept as is or changed. 
"""

class Bullet(object):
    """
     - Bullet represents a single 'record' (i.e. weight, velocity) of the
     corresponding Caliber Object
     - Bullet geometric properties are in metric units
      d = diameter at widest portion of body (mm)
      db = base diameter/frustum diameter (mm). (If the bullet has no frustum
       this will = d).
      oal = (Overall length), the length of the round from the bullet tip to
       the end of the cartridge (mm)
      oabl = (Overall bullet length), the length of the bullet from the tip to
       the base (mm)
      mass = (using weight/mass interchangeably) (kg)
      volume = (m^3)
      vp = 3d position vector the position w.r.t to the shooter (see below) (m)
      xcg = axial center of gravity, the distance from the base of the bullet
       to the center of gravity on the x-axis (centerline) (mm)
      SD = sectional density (see https://www.chuckhawks.com/sd.htm) (kg/m^2)
       to convert to imperial use balllistics.msd2isd
      A = cross-sectional area ((https://www.chuckhawks.com/frontal_area.htm)
       (m^2) to convert to imperial use ballistics.ma2ia
     - Bullet velocity dependent properties
      BC = ballistic coefficient (kg/m^2) is calculate dynamically using the
       current velocity to convert to imperial use ballistics.msd2isd
      i = form factor (McCoy 1991, eq 6.44 pg 101 (unitless)
      ke = kinetic energy dynamic using current velocity (J, kg*m^2/s^2)
      momentum = (kg*m/s)
     Bullets additionally have a velocity vector which consists of the axial
     components (x,y,z) of velocity and position = (x,y,z)
     Properties will be returned in the standard units of measurements kilogram,
     meter by default

     Adopting the coordinate system of McCoy 1999, Figure 5.1

         y=elevation
         |
         |
         |
         |
         |
         |
        S|_____________x=range  T
        /(0,0,0)
       /
      /
     /
    /
    z=windage
    where S is the shooter, T is the target

    and viewing the shooter top-down as a clock. The z-axis runs from 9 o'clock to
    3 o'clock with the x-axis running 6 o'clock to 12 o'clock

          9                          N = 0
          |                          |
    6-----S-----12 T    W = 270------*------E = 90
          |                          |
          3                          S = 180

    where S = shooter, T = target. We will assume that due North = 0 degrees is
    at 9 o'clock, and due West = 270 degrees is 6 o'clock, etc
    """

    def __init__(self,name,cd,sz,d,clen,blen,w,v,cg,mdl):
        """
         initialize bullet
        :param name: string name of caliber
        :param cd: one of {H,R,S} firearm used to fire bullet
        :param sz: string representing size (i.e. 7.62x39) (not the same as dxclen)
        :param d: diameter in mm
        :param clen: case length mm
        :param blen: bullet length mm
        :param w: weight (gr)
        :param v: muzzle velocity (m/s)
        :param cg: charge (gr)
        :param mdl: drag model to use (currently only G1 or G7) if None set
          based on type of bullet (G1 for handgun vs G7 for long gun)
        :return: a Bullet object

        In addition to initialized parameters, also creates internal variables
         _A = <double> cross-sectional area (m^2)
         _SD = <double> sectional density (kg/m^2)
         _m = <double> mass (kg) by converting w in gr to kg
         _t = <double> 'current' time (s)
         _t = <double> 'current' time (s)
         _vi = <double> 'current' velocity at time _t (m/s)
         _vv = <3D vector> 'current' velocity components (m/s) on the x,y,z axis
         _vp = <3D vector> 'current' position (m) on the x,y,z axis
         _p0 = <double> initial spin rate (radians/s) defined in Lahti 2019 eq 1
         _pi = <double> current spin rate (radians/s) defined in Lahti 2019 eq 2
         _vol = <double> estimated volume of the bullet (mm^3)
         _db = <double> base diameter (for bullets with frustums else db = d)
         _xcg = <double> bullet's x-axis center of gravity (dist. from base) (mm)
         _Ix = <double> bullet's axial moment of inertia (kg-m^2)
         _Iy = <double> bullet's transverse moment of inertia (kg-m^2)
         _ts = <list of tuples> steps for iteration where each
           tuple t = (Time,current velocity,position vector,velocity vector)
        """
        # get passed in variables (assume correct)
        self._name = "{} ({}gr)".format(name,w)
        self._cd = cd
        self._sz = sz
        self._d = np.double(d)
        self._oal = np.double(clen)
        self._oabl = np.double(blen) if blen else None
        self._v0 = np.double(v)
        self._w = w
        self._cg = cg

        # create mass properties that can derived at initiation, that are model
        # dependent and variables for trajectory calculations
        self._m = np.double(self._w*bls.GR2KG)
        self._A = np.pi*np.power(self._d*bls.MM2M/2,2)
        self._SD = self._m / np.power(self._d*bls.MM2M,2)
        self._db = None
        self._vol = None
        self._xcg = None
        self._Ix = None
        self._Iy = None
        self._t = 0.
        self._vi = None # at time t = i
        self._vp = None # at time t = i
        self._vv = None # at time t = i
        self._vp0 = None # at time t = 0
        self._vv0 = None # at time t = 0
        self._p0 = None # spin rate at time t = 0
        self._pi = None # spin rate at time t = i
        self._ts = []

        # reset the bullet and set the model
        self.reset()
        self.mdl = mdl

    def reset(self):
        """ resets bullet to initial state (immediately prior to fire) """
        # set elpased time back and steps to empty
        self._t = np.double(0.)
        self._ts = []

        # position and velocity vectors are set to prior to fire, speeds are
        # reset & spin rate is not a number
        self._vp0 = np.array([0.,0.,0.],np.double)
        self._vv0 = np.array([0.,0.,0.],np.double)
        self._vp = np.array([0.,0.,0.],np.double)
        self._vv = np.array([0.,0.,0.],np.double)
        self._vi = self._v0
        self._pi = self._p0 = np.nan

    @property # initial/muzzle velocity (m/s)
    def name(self): return self._name

    @property # current velocity (m/s)
    def v_0(self): return self._v0

    @property # current velocity (m/s)
    def v_i(self): return self._vi

    @property # initial spin rate (radians/s)
    def p_0(self): return self._p0

    @property # current spin rate (radians/s)
    def p_i(self): return self._pi

    @property
    def mdl(self): return self._mdl

    @mdl.setter
    def mdl(self,mdl):
        """
        set the bullet's drag model to mdl, update the volume and base diameter
        of the bullet
        :param mdl:
        :return:
        """
        if mdl not in model.models:
            raise BulletException(
                "{}: invalid 'drag model' ({})".format(self._name,mdl)
            )
        self._mdl = mdl
        self._mass_props_()

    @property
    def pos(self): return self._vp

    @property
    def velocity(self): return self._vv

    @property
    def x(self): return self._vp[bls.IX] # 'current' horizontal position

    @property
    def x_0(self): return self._vp0[bls.IX] # initial horizontal position

    @property
    def y(self): return self._vp[bls.IY] # 'current' vertical position

    @property
    def y_0(self): return self._vp0[bls.IY]  # initial vertical position

    @property
    def z(self): return self._vp[bls.IZ] # 'current' windage distance

    @property
    def z_0(self): return self._vp[bls.IZ]  # initial windage position

    @property
    def v_x(self): return self._vv[bls.IX] # 'current' x-component of velocity

    @property
    def v_x0(self): return self._vv0[bls.IX]  # initial x-component of velocity

    @property
    def v_y(self): return self._vv[bls.IY] # 'current' y-component of velocity

    @property
    def v_y0(self): return self._vv0[bls.IY]  # initial y-component of velocity

    @property
    def v_z(self): return self._vv[bls.IZ] # 'current' x-component of velocity

    @property
    def v_z0(self): return self._vv0[bls.IZ]  # initial x-component of velocity

    @property
    def grains(self): return self._w

    @property # (kg)
    def mass(self): return self._m

    @property # sectional density (kg/m^2)
    def SD(self): return self._SD

    @property # cross sectional area (m^2)
    def A(self): return self._A

    @property # diameter (m)
    def d(self): return self._d*bls.MM2M

    @property # diameter (mm)
    def d_mm(self): return self._d

    @property # radius (m)
    def r(self): return self.d/2

    @property # diameter of base (diameter of frustum) (m)
    def db(self): return self._db*bls.MM2M

    @property # diameter of base (diameter of frustum) (mm)
    def db_mm(self): return self._db

    @property # radius of base (radius of frustum) (mm)
    def rb(self): return self.db/2

    @property # overal length (including cartridge) (m)
    def oal(self): return self._oal*bls.MM2M

    @property # overal length (includeing cartridge) (mm)
    def oal_mm(self): return self._oal

    @property # overal length of bullet (m)
    def oabl(self): return self._oabl*bls.MM2M

    @property # overal length of bullet (mm)
    def oabl_mm(self): return self._oabl

    @property
    def charge(self): return self._cg

    @property # x-dimension of center of gravity (m)
    def xcg(self): return self._xcg*bls.MM2M

    @property # x-dimension of center of gravity (mm)
    def xcg_mm(self): return self._xcg

    @property # axial moment of inertia
    def Ix(self): return self._Ix

    @property # transvers moment of inertia
    def Iy(self): return self._Iy

    @property # estimated volume (m^3)
    def volume(self): return self._vol

#### OUTPUT

    def print_steps(self,i=0,j=None):
        """
         print the steps defined in the range i to j
        :param i: first step
        :param j: last step use -1 for last element and None to print only at i
        """
        idx = 0 # step counter
        strf = "({}): t={:.2f}. Vel.= {:.1f}. Pos: (x={:.2f},y={:.2f},z={:.2f})" \
               " and Vel: (x={:.3f},y={:.3f},z={:.3f})."
        if j is None:
            print(
                strf.format(
                    idx,
                    self._ts[i][0],
                    self._ts[i][1],
                    self._ts[i][2][bls.IX],
                    self._ts[i][2][bls.IY],
                    self._ts[i][2][bls.IZ],
                    self._ts[i][3][bls.IX],
                    self._ts[i][3][bls.IY],
                    self._ts[i][3][bls.IZ]
                )
            )
            return
        if j == -1: j = len(self._ts)
        for i in range(i,j):
            print(
                strf.format(
                    idx,
                    self._ts[i][0],
                    self._ts[i][1],
                    self._ts[i][2][bls.IX],
                    self._ts[i][2][bls.IY],
                    self._ts[i][2][bls.IZ],
                    self._ts[i][3][bls.IX],
                    self._ts[i][3][bls.IY],
                    self._ts[i][3][bls.IZ],
                )
            )
            idx += 1

#### CALCULABLE ATTRIBUTES
# change w.r.t velocity etc

    def ke(self): return np.round(self._m*np.power(self._vi,2)*0.5,3)

    def momentum(self): return np.round(self._m*self._vi,3)

    def mach(self): return atm.mach(self._vi)

    def BC(self):
        """
         Ballistic Coefficient (BC) can be found in ammunition manufacturer
         labels or via equation. As manufacturer can be missleading and due
         to BC changing w.r.t velocity, will implement as an equation of velocity
         :return: ballistic coefficient in kg/m^2
         see http://www.dexadine.com/qbcunits.html
         BC = m / (i(v,md)*d^2)
          where
           m = mass,
           i(v,md) = Form Factor at velocity v for model md and
           d = diameter

         i = Cd(v) / Cd(m,v)
          where
           Cd(v) = bullet's coefficient of drag at velocity v and
           Cd(m,v) = coefficient of drag for standard model bullet at velocity v

         Note
          See SD, the BC is in kg/m^2 whereas BCs are published by
           ft/in^2 - to convert this BC to imperial, divide by 703
        """
        return self._m/(self.i()*np.power(self.d,2))

    def i(self):
        """
        calculates bullet's form factor
        :return: the form factor
         McCoy 1999 pg 101
         i = Cd(v)/Cd(v,m)
          where
           Cd is the coefficient of drag at velocity v and
           Cd(v,m) is the coefficient of drag for the standard model m at velocity
            v
        """
        try:
            return self.Cd()/model.std_cd(self)
        except TypeError:
            raise BulletException("{} Cd is undefined".format(self._name))

    def Cd(self):
        """
        calculates the drag coefficient (uses hooke's method)
        :return: drag coefficient
        """
        try:
            return self._hooke_cd_()
        except TypeError: # shouldn't happen
            raise BulletException(
                "{} velocity is undefined/invalid".format(self._name)
            )

    def Cl(self):
        """
        calculates the lift coefficient
        :return: lift coefficient
        uses Lahti 2019 eq 8, pg 41
        Cl = Cn - Cd0
         where
          Cn = normal force coefficient slope, and
          Cd0 = zero-yaw drag coefficient
        """
        return self.Cn() - self.Cd0()

    def Cmp(self):
        """
        calculates the pitch damping moment coefficients
        :return: tuple of the pitch damping moment coefficients
          Cma = coefficient due to rate of change of angle of attack alpha, and
          Cmq = coefficient due to transverse angular velocity
        uses Lahti 2019 eq 9, 10 pg 41
         eq 9:  Cma = 2 * (Vo - Sb * (l-xcg) / (sqrt(M)*S*d) (M >= 1 supersonic)
         eq 10: Cma = 2*sqrt(M) * (Vo - Sb * (l-xcg)) / (S*d) (M < 1 subsonic)
          where
           Vo = Volume (m^3),
           S = cross sectional area (A()) (m^2),
           Sb = base area (m^2),
           xcg = location of center of gravity on spin axis (m),
           l = bullet length (m), and
           d = bullet diameter (m)

        This gives us Cma in McCoy eq 10.11, pg 223 to find Cmq we'll use
        Lahti eq 12, 13 pg 41
         Cmq + Cma = -2*Cn * ((l-xcg)/d)^2 (M >= 1 supersonic)
         Cmq + Cma = -Cn * ((l-xcg)/d)^2
          where
           Cn = normal force coefficient,
           l = bullet length (m),
           xcg = center of gravity (m), and
           d = bullet diameter (m)
        TODO: Lahti measures CG from the nose of the bullet, we measure
         from the base so we are not subtracting the xcg from the length
        """
        # calculate frustum base area
        Sb = np.pi*np.power(self.rb,2) # m^2

        # get the mach number (set to 0.5 if less than 0.5)
        M = self.mach()
        if M < 0.5: M = 0.5 # use 0.5 for any mach < 0.5

        # pre-calculated commonalitys of eq 9 & 10 and eq 12 & 13
        _1 = (self._vol - Sb*self.xcg) / (self.A*self.d)
        _2 = -self.Cn()*np.power(self.xcg/self.d,2)

        # factors of eq 9 & 10, and of eq 12 & 13 mach dependent
        if M > 1.:
            fa = 2/np.sqrt(M)
            fq = 2
        else:
            fa = 2*np.sqrt(M)
            fq = 1

        # calculate Cma
        Cma = fa*_1
        Cmq = (fq*_2) - Cma
        return Cma,Cmq

    def Cmm(self):
        """
         calculates the Magnus moment coefficient
        :return: mangus momenet coefficient
         uses Lahti 2019 eq 11 pg 41
          (M/2 - 1) * 1/4*db * (l-xcg)/2
          where
           M = mach number,
           l = bullet length (mm),
           db = diameter of base (mm),
           xcg = center of gravity (mm), and
           d = bullet diameter (mm)
         NOTE: conversions to m are not required as they cancel each other out
         TODO: see Cmp, not subtracting CG from bullet length as we measure
         CG from the base vice nose as Lahti does
        """
        M = self.mach()
        if M > 2.5: M = 2.5
        return (M/2 - 1) * (self._oabl/(4*self._db)) * (self._xcg/self._d)

    def Cn(self):
        """
        calculates the normal force coefficient
        :return: normal force coefficient
        uses Lahti 2019 eqs 3 and 4
         Cn = sqrt(M)*sqrt(l/d)*db/d (when M >= 1, supersonic)
         Cn = sqrt(l/d)*db/d (when M < 1, subsonic)
          where
           M = mach number,
           l = length (mm),
           d = diamter (mm), and
           db = base diameter (mm)
        """
        # get mach and calculate factor (sqrt(M) if M > 1, 1 otherwise)
        M = self.mach()
        f = 1 if M < 1 else np.sqrt(M)
        return f * np.sqrt(self._oabl/self._d) * (self._db/self._d)

    def Cd0(self):
        """
        return the zero yaw coefficient of drag (Cd) based on Lahti 2019 eq 5,6
        :return: Cd
        uses Lahti 2019 eq5, 6m 7 pg 40
         Cd = max-Cd / sqrt(M) for M >= 1 (supersonice)
         Cd = max-Cd / 3 for M < 1 (subsonice
         where
          max-Cd = sqrt(d/l) d=diamaeter, l=length, and
          M = mach number
        """
        M = self.mach()
        _1 = np.sqrt(M) if M >= 1. else 3
        return np.sqrt(self._d/self._oabl) / _1

    def Sg(self,rho=atm.RHO):
        """
        calculates the gyroscopic stability parameter
        :return: the gyroscopic stability parameter
        uses Lahti 2019, eq 14
        Sg = (Ix^2 * p^2) / (2*rho*Iy*S*d*Cma)
        where
         Ix = axial moment of inertia (kg*m^2),
         p = spin rate (radians/s),
         rho = air density (kg/m^3),
         Iy = transverse moment of inertia (kg*m^2),
         S = body reference area (m^2),
         d = body diameter (m),
         v = velocity (m/s^2), and
         Cma = pitching moment coefficient slope
         TODO: see Lahti eq 15 "The classical gyroscopic stability critierion for
          spin-stablized bullets is Sg > 1
          For 7.62x39 we get 7.66e-08
        """
        Cma = self.Cmp()[0]
        return (np.power(self._Ix,2)*np.power(self._pi,2)) / (
                2*rho*self._Iy*self._A*self.d*np.power(self._vi,2)*Cma
        )

    def Sd(self):
        """
        calculates the dynamic stability parameter NOT TO BE CONFUSED WITH
        SD = SECTIONAL DENSITY
        :return: the dynamic stability parameter
        use Lahti 2019, eq 16
         Sd = 2*(Cla + (m*d^2)/(2*Ix)*Cmm) /
              Cla - Cd0 - (m*d^2)/2*Iy * (Cmq - Cma)
        where
         Cla = lift force coefficient slope,
         m = total mass (kg),
         Cmm = magnus moment coefficient slope
         Cd0 = zero-yaw drag coefficient, and
         Cmq, Cma are the pitch damping coefficients
        """
        d2 = np.power(self.d,2) # used multiple times
        _1 = 2 * (self.Cl() + ((self._m*d2)/(2*self._Ix)) * self.Cmm())
        _2 = self.Cl() - self.Cd0() - ((self._m*d2)/(2*self._Iy))*sum(self.Cmp())
        return _1 / _2

    def dynamic_stability(self,rho=atm.RHO):
        """
        determines whether the bullet is dynamically stable Lahit Eqs 17, 18
        :return: True if the bullet is stable, False otherwise
        """
        Sg = self.Sg(rho)
        Sd = self.Sd()
        return Sg >= 1/(Sd*(2-Sd))

#### FIRING METHODS

    def chamber(self,h=0.,theta=0.,phi=0.,tr=254,bl=None):
        """
        sets position, velocity vectors and time to instaneous time of fire
        immediately preceding the bullet leaving the muzzle
        :param h: the initial elevation of the firearm (m)
        :param theta: elevation (y-axis) angle of fire (degrees)
        :param phi: windage (z-axis) angle of fire (degrees)
        :param tr: twist rate (mm) i.e. 1:8 = 8*25.4 default is 1:10" or 1:254mm
        :param bl: barrel length (mm) not implemented
        :return:
         #TODO: add atmopsheric conditions
        """
        # 1) reset the bullet to 'empty/starting' state.
        # 2) Set the initial/current position vectors to (0,h,0).
        # 3) set initial and current velocity vectors as well as initial speed
        # 4) Set initial/current spin rate
        self.reset()
        self._vp0[bls.IY] = h
        self._vp = self._vp0.copy()
        self._set_initial_(theta,phi,tr)

    def zero_angle(self,fs,x):
        """
        estimates zero angle at range x given front sight post height (based
        of the k3 constant McCoy 1999, Sec. 5.7, pg 94.
        :param fs: height of front sight post (mm)
        :param x: range to zero at (m)
        :return: zero angle (degrees)
        NOTE: bullet should be reset/initial state
        """
        # 1) use flat-fire trajectory with constant k to estimate rough za as
        #    start point (McCoy 1999, Sec 5.5 pg 92), creating array of angles
        #    from 0 to estimated za
        # 2) interpolate the results of elevation using generated angles
        # 3) use this to find the angle the puts us at the elevation of the
        #    front sight post
        # 4) in the event that the resulting angle does not give of us a
        #    final destination close to 0, recalculate
        # TODO: can't use fill_value='extrapolate', look at extrap1d
        try:
            phis = np.linspace(0,self._za_k_(fs,x)[2],num=1000)
            zs = interp1d(self.elevation(x,fs*bls.MM2M,phis),phis)
            za = float(zs(fs*bls.MM2M))
        except ValueError: # assuming caused by za = float(zs(fs*bls.MM2M))
            raise BulletException(
                "{}: Interpolation error. Try resetting bullet".format(self._name)
            )
        return za

    def tof(self,x,theta):
        """
         estimates time of flight to range x when fired at angle theta
        :param x: range to estimate tof to (m)
        :param theta: angle of gun (degrees)
        :return: estimated flight time (s)
         uses
          t = x/v_x0 * sqrt(v_x0/v_0) McCoy 1999, eq 5.68 pg. 94
           where
           x = range to estimate (m),
           v_x0 = calculated above (m/s), and
           v_0 = muzzle velocity (m/s)
        """
        vx0 = self._vx_i_(theta)
        return (x/vx0) * np.sqrt(vx0/self._vx_k3_(x,theta))

    def elevation(self,x,h,theta):
        """
        estimates the elevation at range x when fired at angle theta
        :param x: range to estimate y at (m)
        :param h: initial elevation (m)
        :param theta: angle of gun (degrees)
        :return: estimated elevation at range x given theta
        uses
         vx as defined in tof(),
         tof(), and McCoy 1999, eq. 5.70 pg 94
         y = h + x*tan(theta) - 0.5*g*t^2 * (1/3 * (1+2*sqrt(vx/vx0)))
         where
          h = initial elevation (m),
          x = range to estimate at (m),
          theta = angle of fire (degrees),
          g = gravity (m/s),
          t = time of flight (s),
          vx = velocity at range x (m/s) as defined in tof, and
          vx0 = initial x-component of velocity at angle theta (m/s)
        """
        vx = self._vx_k3_(x,theta)
        vx0 = self._vx_i_(theta)
        t = self.tof(x,theta)
        return h + (x*utils.tan(theta)) - (0.5*atm.G*np.power(t,2)) * \
               ((1 + (2*np.sqrt(vx/vx0)))/3)

    def fall_angle(self,x,theta):
        """
         calculates the falling angle of the bullet at distance when fired at
         angle theta
        :param x: range (m)
        :param theta: angle of fire (degrees)
        :return: the angle of fall (degrees)
         uses McCoy 1999, eq. 5.69, pg 94
         phi = arctan(tan(theta) - (G*t)/vx0 * (1/3 *(1 + sqrt(vx0/vx) + vx0/vx))))
         where
          theta = firing angle (degrees),
          G = gravity (m/s),
          t = time of flight at distance x (seconds),
          vx0 = x-component of velocity at time of fire (m/s), and
          vx = x-component of velocity at distance x (m/s)
        """
        vx = self._vx_k3_(x,theta)
        vx0 = self._vx_i_(theta)
        t = self.tof(x,theta)
        tanp = utils.tan(theta) - ((atm.G*t)/vx0) * \
               ((1/3) * (1 + np.sqrt(vx0/vx) + vx0/vx))
        return utils.r2d(utils.arctan(tanp,'r'))

    def trajectory(self,dt=0.01,**kwargs):
        """
         calculates trajectory of bullet at intervals of dt seconds when fired
         from h and stopping when bullet elevelation is below th
        :param dt: time interval for calculations (s)
        :param kwargs: keyword argments
         height: height of firearm at time of fire (m) default is average
          height of a shoulder fired firearm
         elev_angle: elevation angle y-axis (degrees) default is perpendicular
          to tangential earth = 0
         windage_angle: windage angle (degrees) default is 0 ('left' < 0)
         brl_twist: barrel twist rate of firearm (mm) (default is 1:254 = 1:10")
         brl_len: barrel length (NOT IMPLEMENTED)
         poi: point of impact, a 3d position vector [x,y,z] of desired poi (m)
          default is [inf,-inf,0] target height (m) default is 0 such that
          calculation stops if bullet's x-position > poi.x or bullet's y-position
          < poi.y
         wind: wind vector (NOTE IMPLEMENTED)
         maxd: stop at pos > maximum distance (m) default is infinity
        """
        # setup default parameters if not specified
        # parameters associated with firearm at time of fire
        h = kwargs['height'] if 'height' in kwargs else S_H
        theta = kwargs['elev_angle'] if 'elev_angle' in kwargs else 0.
        phi = kwargs['windage_angle'] if 'windage_angle' in kwargs else 0.
        bt = kwargs['brl_twist'] if 'brl_twist' in kwargs else 254
        bl = kwargs['brl_len'] if 'brl_len' in kwargs else None
        # parameters associated with poi
        poi = kwargs['poi'] if 'poi' in kwargs else np.array([np.inf,-np.inf,0.])
        # parameters asscoiated with atmosheric conditions
        vw = kwargs['wind'] if 'wind' in kwargs else np.zeros(3,np.double)

        # don't allow firing 'behind' the firer to the left, right, top, bottom
        if 90 < theta < 270:
            raise BulletException(
                "{}: elevation angle is non-sensical".format(self._name,theta)
            )
        if 90 < phi < 270:
            raise BulletException(
                "{}: windage angle is non-sensical".format(self._name,phi)
            )

        # initialize the bullet, then loop until it meets a specified stop cond.
        # and/or IAW McCoy while the remaining velocity is at least 1/2 of the
        # muzzle velocity i.e v_x >= v_x0/2
        # TODO: have to figure out stopping at poi y-position in cases where
        #  the shooter is firing below the target
        self.chamber(h,theta,phi,bt,bl)
        while self._vp[bls.IX] <= poi[bls.IX] and self._vi >= self._v0/2:
            self._ts.append(
                (self._t,np.linalg.norm(self._vv),tuple(self._vp),tuple(self._vv))
            )
            self._trajectory_(dt,vw)

#### REPORTING

    def range_tbl(self,dpath,inc,h,fa,vw=None):
        """
        constructs a (.tsv) range card
        :param dpath: directory path to write
        :param inc: increment of table rows (m)
        :param h: height of fire
        :param fa: angle of fire
        :param vw: wind vector (not used now)
        """
        # TODO: I don't like this function here in the Bullet class
        # check for existence of directory. If it exists, append the filename
        if not os.path.isdir(dpath):
            raise BulletException(
                "{}: directory {} does not exist".format(self._name,dpath)
            )

        fpath = os.path.join(dpath,"{}m range tabel {}.tsv".format(inc,self._name))

        # run trajectory and get interpolations
        self.trajectory(h=h,theta=fa,vw=vw)
        if not self._ts:
            print("{}: No trajectory for height={}, theta-{}".format(
                self._name,h,fa)
            )
            return
        r2t,r2v,r2pos,_ = self.interpolate()

        # set up max distance
        maxd = int(self._ts[-1][2][0] // inc * inc)

        cnt = 0
        fout = None
        try:
            # open file and write details
            fout = open(fpath,'w')
            fout.write(
                "Bullet: {}\nHeight: {:.2f}\nAngle {:.3f}\n".format(self._name,h,fa)
            )

            # write table header
            hdr = "Range (m)\tTime (s)\tVelocity (m/s)\tHeight (m)\tDrop (cm)\t" \
                  "Windage (m)\tDrift (cm)\n"
            fout.write(hdr)

            # write rows starting at 0 increasing by inc
            for i in range(0,maxd+inc,inc):
                it = r2t(i)
                iv = r2v(i)
                ih,iw = r2pos(i)

                fout.write(
                    "{}\t{:.2f}\t{:.3f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\n".format(
                        i,it,iv,ih,(ih-h)*100,iw,(iw-0)*100
                    )
                )
                cnt += 1

        except ValueError as e:
            # assume due to interpolations
            print(
                "{}: Error executing range card {}\nWrote {} rows".format(
                    self._name,e,cnt
                )
            )
        except OSError as e:
            print("{}: {}. Quitting".format(self._name,e))
        finally:
            if fout: fout.close()

    def interpolate(self):
        """
        for now interpolates data accumulated after firing a trajectory based
        on range (horizontal distance)
        TODO: has to be a better way of doing this
        :return: 4 interpolations range to time, velocity, pos. vector & vel. vector
        """
        try:
            xs,ys,zs = zip(
                *[(x,y,z) for x,y,z in [pos for _,_,pos,_ in self._ts]]
            )
            vxs,vys,vzs = zip(
                *[(x,y,z) for x,y,z in [vel for _,_,_,vel in self._ts]]
            )
            ts, vs = zip(*[(x, y) for x, y, _, _ in self._ts])
            r2t = interp1d(xs, ts, kind='linear')
            r2v = interp1d(xs, vs, kind='linear')
            r2pos = interp1d(xs, (ys, zs), kind='linear')
            r2vel = interp1d(xs, (vxs, vys, vzs), kind='linear')
            return r2t, r2v, r2pos, r2vel
        except ValueError:
            raise BulletException("Run trajectory first")

#### PRIVATE HELPERS

    def _trajectory_(self,t,vw):
        """
         calculates trajectory for time interval t
        :param t: the time step
        :param vw: wind vector
        :return: the new vector of velocity
        """
        # TODO: can remove the return and accept of vv
        # 1. update directional components of velocity,
        # 2. calculate acceleration components,
        # 3. calculate then new directional components of velocity,
        # 4. update position vector (using average velocity), and
        # 5. update elapsed time
        # 6. update directional components of velocity
        # 7. calculate new velocity
        # 8. calculate new spin rate
        va = self._va_(vw)
        vv_i = self._vv + va*t
        self._vp += ((self._vv+vv_i)/2)*t + 0.5*va*np.power(t,2)
        self._t += t
        self._update_vi_(np.linalg.norm(vv_i))
        self._vv = vv_i

    def _set_initial_(self,theta,phi,tr):
        """
        initializes internal velocity and spin rate variables to time t=0
        :param theta: elevation angle (degrees)
        :param phi: windage angel (degrees)
        :param tr: firearm twist rate (mm)
        NOTE: initialization of v_i is reduandant as the speed of the bullet
         at time 0 will always = the muzzle velocity
         For initial spin rate, uses Lahti 2019, eq 1
          p0 = 2*pi*v0/R
           where
            v0 = initial velocity (m/s), and
            R = twist rate
        """
        # TODO: should we convert the twist rate to meters? or leave as mm?
        self._vv0 = self._vcomp_i_(theta,phi)       # init velocity vector
        self._vv = self._vv0.copy()                 # copy to current
        self._p0 = (2*np.pi*self._v0)/(tr*bls.MM2M) # init. spin rate
        self._update_vi_(np.linalg.norm(self._vv))  # init velocity magnitude

    def _update_vi_(self,v):
        """
        sets current velocity (magnitude) to v, updates _pi the spin
        :param v: velocity/speed to set _vi to
        Uses Lahti 2019, eq 2, pg 39
         p = p0 * cbrt(v/v0)
          where
           p0 = initial spin rate (radians/s),
           v = current velocity (m/s), and
           v0 = initial velocity (m/s)
        """
        self._vi = v
        self._pi = self._p0 * np.cbrt(self._vi/self._v0)

    def _vcomp_i_(self,theta=0.,phi=0.):
        """
        calculates the 3d components of current velocity w/ angle of fire theta
        and windage angle phi
        :param theta: angle of fire (degrees)
        :param phi: angle of fire along the z-axis (windage)
        :return: the components of velocity as a array [vx,vy,vz]
        """
        return np.array(
            [self._vx_i_(theta,phi),self._vy_i_(theta,phi),self._vz_i_(phi)],
            np.double
        )

    def _vx_i_(self,theta,phi=0.):
        """
        calculates the x component of current velocity at angle of fire theta
        and windage angle phi
        :param theta: angle of fire (degrees)
        :param phi: angle of fire along the z-axis (windage)
        :return: the x component of current velocity
        """
        return self._vi*utils.cos(theta)*utils.cos(phi)

    def _vy_i_(self,theta,phi=0.):
        """
        calculates the y component of current velocity at angle of fire theta
        and windage angle phi
        :param theta: angle of fire (degrees)
        :param phi: angle of fire along the z-axis (windage)
        :return: the y component of current velocity
        """
        return self._vi*utils.sin(theta)*utils.cos(phi)

    def _vz_i_(self,phi=0):
        """
        calculates the z component of current velocity at windage angle phi
        :param phi: angle of fire along the z-axis (windage)
        :return: the z component of current velocity
        """
        return self._vi*utils.sin(phi)

    def _za_k_(self,fs,x):
        """
        estimates the firing angle required to to zero at x using a constant
        drag coefficient k
        :param fs: front sight height over center of bore (mm)
        :param x: desired zero distance (m)
        :return: vx,t,za
          vx = calculate x-component of velocity at distance x
          t = time of flight to x
          za = angle of gun over tangential horizon to zero at x
         uses McCoy 1999, eq. 5.45, pg 92
         phi = arctan(-fs - 0.5*(g * t*2 * (0.5 + s1 - s2) )) / x)
          where
           fs = front sight height (mm),
           g = gravity (m/s),
           t = time of flight to x (s) see _flat_fire_(),
           x = zero range (m),
           and s1, s2 are defined
            s1 = (v_0/v_x - 1)^-1
            s2 = (v_0/v_x - 1)^-2
            where
             v_0 is muzzle velocity (m/s), and
             v_x is velocity at range x (see _flat_fire_())

         reqquires calculation of k, v_x and t
          v_x = v_0/e^(-k*x) (McCoy 1999 eq 5.46, pg 92)
          where
           v_0 = muzzle velocity (m/s),
           k = contanst drag (_k_), and
           x = desired range (m)

         t = x/v_0 * (v_0/v_x - 1) * ln(v_0/v_X) (McCoy 1999, eq 5.47 pg. 92)
          where x = desired range (m),
          v_0  = muzzle velocity (m/s), and
          v_x = velocity at range x (m/s)

         k = p*A/2*m * Cd
          where
           p = atmosphere pressure (kg/m^3),
           A = cross-sectional density (m^2)
           m = mass (kg)
        """
        k = (atm.RHO * self.A) / (2 * self._m) * self.Cd()
        v_x = self._v0 * np.exp(-k*x)
        t = ((x / self._v0) * (self._v0 / v_x - 1)) / np.log(self._v0/v_x)
        dv = self._v0/v_x
        _1 = np.power(dv-1,-1)
        _2 = np.power(dv-1,-2) * np.log(dv)
        za = utils.arctan(
            (fs*bls.MM2M+(0.5*atm.G*np.power(t,2)) * (0.5+_1-_2))/x,'r'
        )
        return v_x,t,utils.r2d(za)

    def _k3_(self):
        """
        estimates the constant k3 (drag coefficient inversely proportional to
        the square root of the mach #) which approximates hooke's Cd equation
        :return: k3
         uses:
          solves for K3 in Cd = K3 / sqrt(M(v)) McCoy 1999 eq. 5.60, pg 94
          where
           Cd = coefficient of drag at muzzle velocity, and
           M(v) = mach number at muzzle velocity
         and
          k3 = p*A/2*m * K3*sqrt(SOUND) McCoy 1999 eq 5.62, pg 94
          where
           p = air density (kg/m^3),
           A = cross sectional density (m^2),
           m = mass of bullet (kg), and
           K3 = calculated above
         NOTE: not using a fixed a Cd which can lead to some problems
        """
        K3 = self.Cd() * np.sqrt(self.mach())
        return ((atm.RHO*self.A) / (2*self._m)) * K3*np.sqrt(atm.SOUND)

    def _vx_k3_(self,x,theta):
        """
        calculates x-component of velocity at range x fired from theta when
        using k3 drag coefficient
        :param x: range to estimate (m)
        :param theta: angel of fire
        :return:
        vx = (sqrt(v_x0) - 0.5*k3*x)^2 McCoy 1999 eq 5.67 pg. 94
          where
           where v_x0 is the x component of velocity at angle theta (m/s),
           k3 = drag coefficient (see _k3_), and
           x = range to estimate (m)
        """
        return np.power(np.sqrt(self._vx_i_(theta)) - (0.5*self._k3_()*x),2)

    def _va_(self,vw):
        """
         calculates the acceleration vector with wind
        :param vw: wind vector
        :return: the acceleration vector (m/s)
        uses McCoy 1999, eqs 8.3, 8.4 & 8.5
         va_x = mFd * V~ * (v_x - w_x)
         va_y = mFd * V~ * (v_y - w_y) - g
         va_z = mFd * V~ * (v_z - w_z)
         where
          mFd = is magnitude of force,
          v~ is the scalar magnitude of the vector velocity (i.e. vi) see eq 8.6
          v_x, v_y, v_z = the x,y,z -components of velocity at current time,
          w_x, w_y, w_z = the x,y,z -components of the wind velocity, and
          g = gravity

        NOTE: _vv is a velocity vector v(x,y,z) which along with the wind
         velocity vector w(x,y,z) gives us v-w = (v_x-w_x, v_y-w_y, v_z-w_z)
         multipled with the scalars mFd and (vi) and added to the
         gravity vector g(0,-9,8,0)

        TODO:
         1 for this, not sure wind vector should be subtracted as McCoy does
         since the axial components of wind are negative - when they 'oppose'
          the bullet i.e. headwinds and downdrafts. However given a extremely
          strong headwind -200 m/s and adding, the bullet's speed and x-component
          of velocity increased when adding and decreased when subtracting -
          not intuitive to me
          Example for 7.62x39 fired at angle of 0.5 with dt = 0.1
           At t = 0 vv = (730.273,6.373,0) & vi = 730.301
           Let vw = (0,0,0) at t = 0.1, vv =(665.300,4.825,0) & vi = 665.318
           Let vw = (-200,0,0) meaning impossible winds
            i. subtract vv & vw: vv = (647.506,4.825,0) & vi = 647.524
            ii add vv & vw: vv = (683.095,4.825,0) & vi = 683.112
          Meaning the vecotrs must be subtracted
          Because mFd in force.py is negative
        """
        return -force.mFd(self) * self._vi*(self._vv-vw) + force.gravity()

    def _hooke_cd_(self):
        """
        returns the coefficient of drag (Cd) based on Hooke 2015 methodology (Table 1,
        pg. 6 using the mach number (no use/definition of model)
        :return: Cd
        NOTE:
         1. entirely velocity dependent (no accounting for model)
         2. Hooke defines the last equation only up to 2.996
         3. Hooke does not specify if this for a particular model or not and could
            only be for 7.62x51 NATO round
        """
        m = self.mach()
        if 0.0 <= m < 0.681: return 0.0017*np.power(m,2) - 0.0009*m + 0.1391
        elif 0.681 <= m < 0.867: return 0.5443*np.power(m,2) - 0.7482*m + 0.397
        elif 0.867 <= m < 0.934: return 0.0072*np.exp(3.5404*m)
        elif 0.934 <= m < 1.000: return 4e-6*np.exp(11.564*m)
        elif 1.000 <= m < 1.075: return -3.5585*np.power(m,2) + 7.6189*m - 3.6256
        elif 1.075 <= m: return 0.634*np.exp(-0.31*m)
        else:
            raise BulletException(
                "{}: Invalid 'mach number' ({})".format(self.name, m)
            )

    def _mass_props_(self):
        """
         estimates mass properties of the bullet based on bullet's length &
         std model ratios. Estimates:
          1. the diameter of the frustum if present,
          2. the volume of the bullet,
          3. the Center of Gravity (CG) of the bullet, and
          4. the moments of inertia (axial & transverse)
         NOTE: Requires that bullet has a known diameter, length and drag model
        """
        # Given the OABL and drag model, estimate ogival, body & frustum lengths,
        # the frustum base diameter and the ogival radius - throw exception of
        # any model that is not a tanget ogive as these are not calculate yet
        mspec = model.mdl_specs[self._mdl]
        if mspec['ogv-shp'] != model.OGV_TAN:
            raise BulletException(
                "{}: Ogive Volume for ({}) not implemented".format(
                    self.name,mspec['ogv-shp']
                )
            )
        lo = mspec['ogv-len'] * self._oabl
        lc = mspec['cyl-len'] * self._oabl
        lf = mspec['frs-len'] * self._oabl
        self._db = mspec['frs-base'] * self._d
        ro = mspec['ogv-r'] * lo

        # get radii of components (Note for do not confuse the ogival radius ro
        # (above) with the radius of the nose's base which = the body radius
        # convert to meters
        rc = self._d / 2  # cylinder/ogive base radius
        rf = self._db / 2 # frustum base radius

        # estimate the volume, center of gravity and moments of inertia
        # 1. VOLUME estimate bullet volume by summing constituent volumes
        # equations for cylinder & frustum volumes can be googled
        #  Vc = pi*r^2*h where r=cylinder radius and h = cylinder length
        #  Vf = 1/3 * pi * (h*r1^2 + (r1*r2) + r2^2)
        #   where
        #    h = frustum length
        #    r1 = radius of frustum base
        #    r2 = radius of frustum top
        # The ogival volume is more complicated, requiring additional parameters
        #  f, the ogival head radius of curvature and nbl the non-dimensional
        #  ballistic length as defined in Seglates 2019 Fig. 2 and Eq. 1
        #  f = Rogive/R
        #  nbl = sqrt(2*f-1)
        #  where
        #   Rogive = ogival radius
        #   R = radius of  base of ogive
        # From these, the ogive volume in Seglates 2019, Eq. 2 (tangent ogives)
        #  Vogive/R^3 = pi * (f^2 - 1/3nbl^2 - f^2*(f-1)*sin^-1(nbl/f))
        #   where f and nbl are defined above
        vc = np.pi * np.power(rc,2) * lc
        vf = (np.pi*lf)/3 * (np.power(rc,2) + (rc*rf) + np.power(rf,2))
        f = ro / rc
        nbl = np.sqrt(2*f-1)
        vo_r3 = np.pi*( # used later to calculate ogive's center of gravity
                (np.power(f,2) - np.power(nbl,2)/3)*nbl -
                (np.power(f,2)*(f-1)*utils.arcsin(nbl/f,'r'))
        )
        vo = vo_r3 * np.power(rc,3)             # complete eq.of ogive volume
        self._vol = (vo + vc + vf) * bls.MMC2MC # store as cubic m

        # 2. Center of Gravity - center of gravity (CG) of the bullet is measured
        # from the base of the bullet on the centerline such that in the 3d space
        # of the bullet (x,y,z), the cg is (Xcg,0,0)

        # As in the case of volume, CG is calcuated by summing the CG for the
        # consituent compoents and is 'simple' for the body and frustum but more
        # complex for the ogive. However, see
        #  https://www.grc.nasa.gov/www/k-12/rocket/rktcg.html
        # we cannot simply add the consituent CGs but must consider the the weights
        # of each component as well such that the CG of the bullet is defined
        #  CGb*Wb = sum(Di*Wi)
        #  where
        #   CGb = center of gravity of bullet,
        #   Wb = weight of bullet,
        #   Di = distance of component i from a reference line
        #   Wi = weight of component i
        # equation for cylinder at
        #  https://engineersfield.com/centre-gravity-centroid-areas-equilibrium/
        # equation for frustum at
        #  https://archive.lib.msu.edu/crcmath/math/math/c/c581.htm
        # CGcylinder = l / 2
        # CGfrustum = l/4 * (rf^2 + 2*(rf*rc) + 3*rc^2)/(rf^2+rf*rc + rc^2))
        #  where
        #   l = length of cylinder,
        #   rf = radius of frustum base (top), and
        #   rc = radius of cylinder,
        # and the equation for an ogive is found in Seglates 2019 Eq. 3
        # CGogive = pi/Vo/R3 * (
        #     -2/3*(f-1)*(f^3-(f-1)^3) + (1/2*(f^2+(f-1)^2)*nbl^2) - 1/4nbl^4
        #     )
        # where
        #  Vo/R3, f, nbl are defined above in the volume of the ogive
        CGc = lc / 2
        CGf = (lf/4) * (
                (np.power(rf,2) + 2*rf*rc + 3*np.power(rc,2)) /
                (np.power(rf,2) + rf*rc + np.power(rc, 2))
        )
        CGo = (np.pi/vo_r3) * (
                -(2/3 * (f-1) * (np.power(f,3) - np.power(f-1,3))) +
                (0.5 * (np.power(f,2) + np.power(f-1,2)) * np.power(nbl,2)) -
                (np.power(nbl,4)/4)
        )

        # and to find the CG of the bullet,
        #  https://www.grc.nasa.gov/www/k-12/rocket/rktcg.html
        # we need to estimate the weights (using weight & mass interchangeably)
        # requires the mean density of the bullet: density = mass / volume
        # w(cyl) = pi * rc^2* lc * density
        # w(frs) = (1/3 * pi * lf * (rb^2 + rb*r * r^2)) * density
        # w(ogv) = w(bullet) - w(cyl) - w(frs)
        #  where
        #   rc = the radius of the bullet,
        #   lc = length of the cylinder,
        #   lf = length of the frustum, and
        #   rb = radius of the base (top of frustum)
        # NOTE: we estimate density as kg/mm^3 to maintain consistency with
        #  radii and lengths which are in millimeter. As such, need to convert
        #  volume to cubic mm
        p = self._m / (self._vol*bls.MC2MMC)
        mc = (np.pi * np.power(rc,2) * lc) * p
        mf = ((np.pi * lf*(np.power(rf,2) + (rf*rc) + np.power(rf,2)))/3) * p
        mo = self._m - mc - mf
        self._xcg = (mf*CGf + mc*(CGc + lf) + mo*(CGo + lf + lc))/self._m

        # 3. Moments of Inertia use the sum of constituent components Boynton,
        # Wiener 2001. We will use kg-m^2 and need to convert lengths, radii
        rf *= bls.MM2M
        lf *= bls.MM2M
        rc *= bls.MM2M
        lc *= bls.MM2M

        # for the frustum use the equation of moments of inertia of a cone
        # https://en.wikipedia.org/wiki/List_of_moments_of_inertia
        # subtracting Iorg - Itop and Iyorg-Iytop where org = 'original cone'
        # and top = the portion of the 'original cone' cut to make the frustum
        # the equations for the moments of a cone are
        # Ix = 3/10 * m * r^2
        # Iy = m*(3/20*r^2 + 3/80*h^2)
        #  where
        #   m = mass of the cone,
        #   r = radius of the base of the cone, and
        #   h = height of the cone
        # Since we know the radius of the orginal cone = rc, the radius of the
        # top cone = rf and the height of the frustum = lf we will use similar
        # triangles to rewrite the heights of the cones in terms of the frustum's
        # height
        #  rc/horg = rf/(horg-lf)
        #  horg = lf*rc/(rc-rf)
        #   and
        #  rc/(htop+lf) = rf/htop
        #  htop = lf*rf/(rc-rf)
        # Likewise, we will first, rewrite the mass of the top in terms of the
        # mass of the cone then rewrite the mass of the cone in terms of the
        # mass of the frustum see
        #  https://www.quora.com/How-can-I-show-that-the-moment-of-inertia-of-
        #  a-truncated-cone-about-its-axis-with-the-radius-of-its-ends-being-a-
        #  and-b-is-frac-3m-10-frac-a-5-b-5-a-3-b-3
        # for details
        # morg = mf/(1-(rf^3/rc^3))
        #  and
        # mtop = (rf^3/rc^3)*morg
        # plugging the values into the moments equations simplifications can be
        # made such that
        # Ixf = 3/10 * mf/(rc*3-rf^3) * (rc^5 - rf^5)
        #  and
        # Iyf = (mf*rc^3)/(rc^3-rf^3) * (3/20*rc^2 + 3/80*horg^2) -  (Original)
        #       (mf*rf^3)/(rc^3-rf^3) * (3/20*rf^2 + 3/80*htop^2)    (Top)
        Ixf = Iyf = 0.
        if mspec['frs-shp'] != model.FRS_NNE: # no frustum, ignore
            lorg = (lf*rc) / (rc-rf)
            ltop = (lf*rf) / (rc-rf)
            _d = np.power(rc,3) - np.power(rf,3)
            Ixf = .3 * (np.power(mf,2)/(np.power(rc,3) - np.power(rf,3))) * \
                       (np.power(rc,5) - np.power(rf,5))
            Iyorg = ((mf*np.power(rc,3)) / _d) * (
                    (3/20) * np.power(rc,2) + (3/80) * np.power(lorg,2)
            )
            Iytop = ((mf*np.power(rf,3)) / _d) * (
                    (3/20) * np.power(rf,2) + (3/80) * np.power(ltop,2)
            )
            Iyf = Iyorg - Iytop

        # the cylinder Ix and Iy can be found in
        # https://en.wikipedia.org/wiki/List_of_moments_of_inertia
        #  Ix = 1/2 * m * r^2
        #  Iy = 1/12 * m * (3*r^2 + h^2)
        #   where
        #    m = mass of cylinder (kg),
        #    r = radius of cylinder (m), and
        #    h = height of cylinder (m)
        Ixc = 0.5 * mc * np.power(rc,2)
        Iyc = mc/12 * (3*np.power(rc,2) + np.power(lc,2))

        # finally for the ogive
        # Using Seglates 2019 Eq. 5
        # Ixogive/p*R^5 = np.pi*(b0*(f^2*arcsin(nbl/f) - (f-1)*nbl) -
        #                  2*nbl^3*(b1/3 + b2/4 + b3/5))
        # where
        #  f, the ogival head radius of curvature and nbl the non-dimensional
        #  ballistic length are defined above and integrations constants b0 -
        #  b3 are
        # b3 = -1,
        # b2 = 3 + 7/5*b3*f,
        # b1 = -3 + 5/4*b2*f,
        # b0 = 1 + 3/5*b1*f,
        # p = mean density of the bullet, and
        # R = radius of the ogive
        b3 = -1
        b2 = 3 + 7/5*b3*f
        b1 = -3 + 5/4*b2*f
        b0 = 1 + 3/5*b1*f
        Ixo = np.pi * (b0*(np.power(f,2)*utils.arcsin(nbl/f) - (f-1)*nbl) -
                       2*np.power(nbl,3) * (b1/3 + b2/4 + b3/5))
        Ixo *= p*np.power(rc,5)

        # and Seglates 2019 Eq. 6
        # Iyogive/P*R^5 = pi/4 * ((f^2*nbl*(f^2+7/2*(f-1)^2)) +
        #                        1/15*nbl^5 -
        #                        f^2*(f-1)*(5/2*f^2 + 2*(f-1)^2)*arcsin(nbl/f)
        _1 = np.power(f,2) * nbl * (np.power(f,2) + (7/2)*np.power(f-1,2))
        _2 = np.power(nbl,5) / 15
        _3 = np.power(f,2) * (f-1) * (
                ((5/2) * np.power(f,2) + 2*np.power(f-1,2)) * utils.arcsin(nbl/f)
        )
        Iyo = (np.pi/4) * (_1 + _2 - _3)
        Iyo *= p*np.power(rc,5)

        # summing the moments
        self._Ix = Ixf + Ixc + Ixo
        self._Iy = Iyf + Iyc + Iyo