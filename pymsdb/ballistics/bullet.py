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
__version__ = '0.1.7.2'
__date__ = 'July 2021'
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
from pymsdb import PyMSDbException
from pymsdb.ballistics.anthropometry import S_H
from pymsdb.utils import r2d,sin,cos,tan,arctan,arcsin

class BulletException(PyMSDbException):
    def __init__(self,msg): super().__init__(msg)

"""
TODO: 
  8. Cd is currently using method that does not specify ballistic models
  9. Want range based functions vice time (requires integration?) 
 15. add option to set/configure windage angle
 18. For ke, momentum do we want to round or not?
 19. round off errors
 20. Find equation for barrel length effect on initial velocity, perhaps Le Duc
 21. Have seen sectional density (SD) defined as the (1) ratio of mass to cross-
  sectional area and (2) ratio of mass to the square of the diameter. Currently
  using 2, should this be kept as is or changed. 
 22. Figure out what to do with wind vector (see McCoy 1999, eq 9.22 pg 191)
 23. After adding windage angle, should we go back to estimating functions like
  tof and add it or keep as 2d estimating functions
 24. Look at IJCSM1405Gritzapisetal1..pd.pdf sec 9 for aerodynamic jump deflection
  equation
 25. When not specified by user, set barrel twist rate based on handgun or long gun
  type of ammo
 26. No longer using the ch 5 k constant drag coefficients 
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
       to convert to imperial use ballistics.msd2isd
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
         _vol = <double> estimated volume of the bullet (mm^3)
         _db = <double> base diameter (for bullets with frustums else db = d)
         _Ix = <double> bullet's axial moment of inertia (kg-m^2)
         _Iy = <double> bullet's transverse moment of inertia (kg-m^2)
         _t = <double> 'current' time (s)
         _v_t = <double> 'current' velocity magnitude at time _t (m/s)
         _vV_0, _vV_t = <3D vector> velocity components (m/s) on the x,y,z axis
           at initial and time t = self._t
         _vP_0,_vP_t = <3D vector> 'current' position (m) on the x,y,z axis at
           initial and time t = self._t
         _vI, _vJ, _vK = (3D vector>s unit vectors along the bullets axes
          which are used to calculate the angular momentum and orginate at
          (xcg,0,0) (see McCoy 1999, Fig 9.1, sections 9.2 and 9.3)
         _p_0, _p_t = <double> initial and current spin rate (radians/s)
         _xcg = <double> bullet's x-axis center of gravity (dist. from base) (mm)
         _ts = <list of tuples> steps for iteration where each
           tuple t = (Time,current velocity,position vector,velocity vector)

         NOTE: naming conventions
          - time dependent variables are defined with suffix _<time> where time
            is 0=initial or t=current i.e. v_t = current velocity at time
            t = self._t
          - vectors are defined as v<VAR> where the previx 'v' denotes a vector
            and <VAR> is  one or more characters and the first character is
            capitalized i.e. vP_0 is the initial position vector
        """
        # get passed in variables (assume correct) and set variables that can 
        # be defined at initialiation
        self._name = "{} ({}gr)".format(name,w)
        self._cd = cd
        self._sz = sz
        self._d = d
        self._oal = clen
        self._oabl = blen if blen else None
        self._v_0 = v
        self._w = w
        self._cg = cg
        self._m = np.round(self._w*bls.GR2KG,6)
        self._A = np.round(np.pi*np.power(self.d*bls.MM2M,2)/4,6)
        self._SD = np.round(self._m / np.power(self.d*bls.MM2M,2),6)

        # create mass properties that can derived at initiation, properties that
        # are model dependent and variables for trajectory calculations
        self._db = None   # base diameter (frustum base)
        self._vol = None  # volume
        self._xcg = None  # x-axis center of gravity
        self._Ix = None   # axial momement of inertia
        self._Iy = None   # transverse moment of inertia
        self._t = 0.      # current time 
        self._p_0 = None  # initial spin rate
        self._p_t = None  # spin rate at time t = self._t
        self._vP_0 = None # intial position vector
        self._vP_t = None # position vector at time t = self._t
        self._vV_0 = None # initial velocity vector
        self._vV_t = None # velocity vector at time t = self._t
        self._v_t = None  # velocity magnitude at time t = self._t
        self._vI = None   # unit vector along bullet's x-axis
        self._vJ = None   # unit vector along bullet's y-axis
        self._vK = None   # unit vector along bullet's z-axis
        self._ss = []     # steps
        
        # reset the bullet and set the model
        self.reset()
        self.mdl = mdl

    def reset(self):
        """ resets bullet to initial state (immediately prior to fire)  """
        # set elpased time back and steps to empty
        self._t = np.double(0.)
        self._ss = []

        # position and velocity vectors are set to prior to fire, speeds are
        # reset. Bullet axes and  spin rate are undefined
        self._vP_0 = np.array([0.,0.,0.],np.double)
        self._vV_0 = np.array([0.,0.,0.],np.double)
        self._vP_t = np.array([0.,0.,0.],np.double)
        self._vV_t = np.array([0.,0.,0.],np.double)
        self._vI = self._vJ = self._vK = np.nan
        self._v_t = self._v_0
        self._p_t = self._p_0 = np.nan

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
    def name(self): return self._name

    @property
    def charge(self): return self._cg

    @property
    def v_0(self): return self._v_0 # initial/muzzle velocity

    @property
    def v_t(self): return self._v_t # current velocity

    @property
    def p_0(self): return self._p_0 # initial spin rate

    @property
    def p_t(self): return self._p_t # current spin rate

    @property
    def vP_t(self): return self._vP_t

    @property
    def vV_t(self): return self._vV_t

    @property
    def x_0(self): return self._vP_0[bls.IX] # initial horizontal position

    @property
    def x(self): return self._vP_t[bls.IX] # 'current' horizontal position

    @property
    def y_0(self): return self._vP_0[bls.IY]  # initial vertical position

    @property
    def y(self): return self._vP_t[bls.IY] # 'current' vertical position

    @property
    def z_0(self): return self._vP_0[bls.IZ]  # initial windage position

    @property
    def z(self): return self._vP_t[bls.IZ] # 'current' windage distance

    @property
    def weight(self): return self._w

    @property
    def mass(self): return self._m

    @property
    def volume(self): return self._vol

    @property
    def SD(self): return self._SD

    @property
    def A(self): return self._A

    @property
    def d(self): return self._d

    @property
    def db(self): return self._db

    @property
    def oal(self): return self._oal

    @property
    def oabl(self): return self._oabl

    @property
    def xcg(self): return self._xcg # x-dimension of center of gravity

#### OUTPUT

    def print_steps(self,i=0,j=None):
        """
         print the steps defined in the range i to j
        :param i: first step
        :param j: last step use -1 for last element and None to print only at i
        """
        strf = "At time t={:.2f}. Vel.= {:.3f}. Pos: (x={:.3f},y={:.3f},z={:.3f})" \
               " and Vel: (x={:.3f},y={:.3f},z={:.3f})."
        if j is None:
            print(
                strf.format(
                    self._ss[i][0],
                    self._ss[i][1],
                    self._ss[i][2][bls.IX],
                    self._ss[i][2][bls.IY],
                    self._ss[i][2][bls.IZ],
                    self._ss[i][3][bls.IX],
                    self._ss[i][3][bls.IY],
                    self._ss[i][3][bls.IZ]
                )
            )
            return
        if j == -1: j = len(self._ss)
        for i in range(i,j):
            print(
                strf.format(
                    self._ss[i][0],
                    self._ss[i][1],
                    self._ss[i][2][bls.IX],
                    self._ss[i][2][bls.IY],
                    self._ss[i][2][bls.IZ],
                    self._ss[i][3][bls.IX],
                    self._ss[i][3][bls.IY],
                    self._ss[i][3][bls.IZ],
                )
            )

#### CALCULABLE ATTRIBUTES
# change w.r.t velocity etc

    def ke(self): return np.round(self._m*np.power(self._v_t,2)*0.5,6)

    def momentum(self): return np.round(self._m*self._v_t,6)

    def mach(self): return np.round(atm.mach(self._v_t),6)

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
        return np.round(self._m/(self.i()*np.power(self._d*bls.MM2M,2)),6)

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
            return np.round(self.Cd()/model.std_cd(self),6)
        except TypeError:
            raise BulletException("{} Cd is undefined".format(self._name))

    def Cd(self):
        """
        calculates the drag coefficient (uses hooke's method)
        :return: drag coefficient
        """
        try:
            return np.round(self._hooke_cd_(),6)
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
        return np.round(self.Cn() - self.Cd0(),6)

    def Cmp(self):
        """
        calculates the pitch damping moment coefficient
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
        # since A is in m^2 and volume is in m^3 convert, xcg, d to meters
        xcg = self._xcg*bls.MM2M                     # m
        d = self._d*bls.MM2M                         # m
        Sb = np.pi*np.power(self._db*bls.MM2M,2) / 4 # m^2

        # get the mach number (set to 0.5 if less than 0.5)
        M = self.mach()
        if M < 0.5: M = 0.5 # use 0.5 for any mach < 0.5

        # pre-calculated commonalitys of eq 9 & 10 and eq 12 & 13
        _1 = (self._vol - Sb*xcg) / (self.A*d)
        _2 = -self.Cn()*np.power(xcg/d,2)

        # factors of eq 9 & 10, and of eq 12 & 13 mach dependent
        if M > 1.:
            fa = 2/np.sqrt(M)
            fq = 2
        else:
            fa = 2*np.sqrt(M)
            fq = 1

        # calculate Cma
        Cma = np.round(fa*_1,6)
        Cmq = np.round(fq*_2 - Cma,6)
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
        return np.round(
            (M/2 - 1) * (self._oabl/(4*self._db)) * (self._xcg/self._d),6
        )

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
        return np.round(f*np.sqrt(self._oabl/self._d)*(self._db/self._d),6)

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
        return np.round(np.sqrt(self._d/self._oabl)/_1,6)

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
        _1 = np.power(self._Ix,2)*np.power(self._p_t,2)
        _2 = 2*rho*self._Iy*self._A*(self._d*bls.MM2M)*np.power(self._v_t,2)*Cma
        return np.round(_1/_2,6)

    def Sd(self):
        """
        calculates the dynamic stability parameter
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
        d2 = np.power(self._d*bls.MM2M,2) # used multiple times
        _1 = 2 * (self.Cl() + ((self._m*d2)/(2*self._Ix)) * self.Cmm())
        _2 = self.Cl() - self.Cd0() - ((self._m*d2)/(2*self._Iy))*sum(self.Cmp())
        return np.round(_1/_2,6)

    def dynamic_stability(self,rho=atm.RHO):
        """
        determines whether the bullet is dynamically stable Lahit Eqs 17, 18
        :return: True if the bullet is stable, False otherwise
        """
        Sg = self.Sg(rho)
        Sd = self.Sd()
        return Sg >= 1/(Sd*(2-Sd))

#### FIRING METHODS

    """
     IAW McCoy the current firing methods below should only be used out to the
     range at which the remaining velocity is at least 1/2 of the muzzle 
     velocity i.e v_x >= v_x0/2
    """

    def chamber(self,**kwargs):
        """
        based on McCoy 1999 section 9.3, sets initial conditions for trajectory
        calculations:
         initial velocity vector
         initial spin rate
         bullet's axial position vector
         bullet's world position vector
        for time = 0 immediately at the moment of fire
        :param kwargs: keywords arguments as follows
         height: height of firearm at time of fire (m) default is average height
          of a shoulder fired firearm
         elev_angle elevation angle i.e. y-axis (degrees) default is tangetial
         to earth = 0
         windage_angle: windage angle i.e. z-axis (degrees) default is 0
         pitch_angle: the pitch angle at gun muzzle
         yaw_angle: the yaw angle at gun muzzle
         brl_length: (NOT USED YET) the length of the barrel (mm)
         brl_twist: the twist rate (mm) of the barrel default is dependent on
          the type of bullet handgun vs long gun
        NOTE:
         For spin rate, uses McCoy 1999, eq 9.31
          p0 = 2*pi*v0/(n*d)
           where
            v0 = initial velocity (m/s),
            n = twist rate (mm), and
            d = diameter (mm)
        """
        # set up variables to default unless specified in kwargs
        h = kwargs['height'] if 'height' in kwargs else S_H
        phi = kwargs['elev_angle'] if 'elev_angle' in kwargs else 0.
        theta = kwargs['windage_angle'] if 'windage_angle' in kwargs else 0.
        alpha = kwargs['pitch_angle'] if 'pitch_angle' in kwargs else 0.
        beta = kwargs['yaw_angle'] if 'yaw_angle' in kwargs else 0.
        bt = kwargs['brl_twist'] if 'brl_twist' in kwargs else 254
        #bl = kwargs['brl_len'] if 'brl_len' in kwargs else None

        # 1) reset the bullet to 'empty/starting' state.
        # 2) Set the initial and current position vectors to (0,h,0).
        # 3) set initial and current velocity vectors and get current velocity
        #  magnitude
        # 4) Set initial and current spin rate
        # 5) set x,y,z unit vectors need pitch angle and yaw angle
        self.reset()
        self._vP_0[bls.IY] = h
        self._vP_t = self._vP_0.copy()
        self._vV_0 = self._vv_3d_(phi,theta)
        self._vV_t = self._vV_0.copy()
        self._v_t = np.linalg.norm(self._vV_t)
        self._p_0 = (2*np.pi*self._v_0)/(self._d*(bt*bls.MM2M))
        self._p_t = self._p_0
        self._vI,self._vJ,self._vK = self._vb_3d_(phi,theta,alpha,beta)

    def zero_angle(self,fs,x):
        """
        estimates zero angle at range x given front sight post height (based
        of the k3 constant McCoy 1999, Sec. 5.7, pg 94.
        :param fs: height of front sight post (mm)
        :param x: range to zero at (m)
        :return: zero angle (degrees)
        NOTE: bullet should be reset/initial state
        """
        # create array of angles from 0 to 45 and interpolate the results of
        # elevation using generated angles. use the interpolation to find the
        # angle the puts us at the elevation of the front sight post
        # TODO:
        #  - can't use fill_value='extrapolate', look at extrap1d
        #  - minimize the number of angles, thetas
        try:
            fs *= bls.MM2M
            #thetas = np.linspace(0,45.,num=1000)
            thetas = np.linspace(0,1.)
            zs = interp1d(self.elevation(x,fs,thetas),thetas)
            return np.round(zs(fs),6)
        except ValueError: # assuming caused by za = float(zs(fs*bls.MM2M))
            raise BulletException(
                "{}: Interpolation error. Try resetting bullet".format(self._name)
            )

    def tof(self,x,phi):
        """
         estimates time of flight to range x when fired at angle phi
        :param x: range to estimate tof to (m)
        :param phi: angle of gun (degrees)
        :return: estimated flight time (s)
         uses
          t = x/v_x0 * sqrt(v_x0/v_0) McCoy 1999, eq 5.68 pg. 94
           where
           x = range to estimate (m),
           v_x0 = calculated above (m/s), and
           v_0 = muzzle velocity (m/s)
        """
        vx0 = self._vx_t_(phi)
        return np.round((x/vx0) * np.sqrt(vx0/self._vx_k3_(x,phi)),6)

    def elevation(self,x,h,phi):
        """
        estimates the elevation at range x when fired at angle phi
        :param x: range to estimate y at (m)
        :param h: initial elevation (m)
        :param phi: angle of gun (degrees)
        :return: estimated elevation at range x given phi
        uses
         vx as defined in tof(),
         tof(), and McCoy 1999, eq. 5.70 pg 94
         y = h + x*tan(phi) - 0.5*g*t^2 * (1/3 * (1+2*sqrt(vx/vx0)))
         where
          h = initial elevation (m),
          x = range to estimate at (m),
          phi = angle of fire (degrees),
          g = gravity (m/s),
          t = time of flight (s),
          vx = velocity at range x (m/s) as defined in tof, and
          vx0 = initial x-component of velocity at angle phi (m/s)
        """
        vx = self._vx_k3_(x,phi)
        vx0 = self._vx_t_(phi)
        t = self.tof(x,phi)
        return np.round(
            h + (x*tan(phi)) -
            (0.5*atm.G*np.power(t,2)) * ((1 + (2*np.sqrt(vx/vx0)))/3),6
        )

    def fall_angle(self,x,phi):
        """
         calculates the falling angle of the bullet at distance when fired at
         angle phi
        :param x: range (m)
        :param phi: angle of fire (degrees)
        :return: the angle of fall (degrees)
         uses McCoy 1999, eq. 5.69, pg 94
         theta = arctan(tan(phi) - (G*t)/vx0 * (1/3 *(1 + sqrt(vx0/vx) + vx0/vx))))
         where
          phi = firing angle (degrees),
          G = gravity (m/s),
          t = time of flight at distance x (seconds),
          vx0 = x-component of velocity at time of fire (m/s), and
          vx = x-component of velocity at distance x (m/s)
        """
        vx = self._vx_k3_(x,phi)
        vx0 = self._vx_t_(phi)
        t = self.tof(x,phi)
        tanp = tan(phi) - ((atm.G*t)/vx0) * \
               ((1/3) * (1 + np.sqrt(vx0/vx) + vx0/vx))
        return np.round(r2d(arctan(tanp,'r')),6)

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
         wind: wind vector (NOT IMPLEMENTED)
         maxd: stop at pos > maximum distance (m) default is infinity
        """
        # setup default parameters if not specified
        # parameters associated with firearm at time of fire
        h = kwargs['height'] if 'height' in kwargs else S_H
        phi = kwargs['elev_angle'] if 'elev_angle' in kwargs else 0.
        theta = kwargs['windage_angle'] if 'windage_angle' in kwargs else 0.
        bt = kwargs['brl_twist'] if 'brl_twist' in kwargs else 254
        bl = kwargs['brl_len'] if 'brl_len' in kwargs else None
        # parameters associated with poi
        poi = np.array([np.inf,-np.inf,0.])
        if 'poi' in kwargs and kwargs['poi'] is not None: poi = kwargs['poi']
        # parameters asscoiated with atmosheric conditions
        vw = np.zeros(3,np.double)
        if 'wind' in kwargs and kwargs['wind'] is not None: vw = kwargs['wind']

        # don't allow firing 'behind' the firer to the left, right, top, bottom
        if 90 < phi < 270:
            raise BulletException(
                "{}: elevation angle is non-sensical".format(self._name,phi)
            )
        if 90 < theta < 270:
            raise BulletException(
                "{}: windage angle is non-sensical".format(self._name,theta)
            )

        # initialize the bullet, then loop until it meets a specified stop cond.
        # and/or IAW McCoy while the remaining velocity is at least 1/2 of the
        # muzzle velocity i.e v_x >= v_x0/2
        # TODO: have to figure out stopping at poi y-position in cases where
        #  the shooter is firing below the target
        self.chamber(
            height=h,elev_angle=phi,windage_angle=theta,brl_twist=bt,brl_len=bl
        )
        while self._vP_t[bls.IX] <= poi[bls.IX] and self._v_t >= self._v_0/2:
            self._ss.append(
                (
                    self._t,np.linalg.norm(self._vV_t),
                    tuple(self._vP_t),
                    tuple(self._vV_t)
                )
            )
            self._trajectory_(dt,vw)

#### REPORTING

    def range_tbl(self,dpath,inc,h,fa,vw=None):
        """
        constructs a (.tsv) range card
        :param dpath: directory path to write
        :param inc: increment of table rows (m)
        :param h: height of fire (m)
        :param fa: angle of fire (degrees)
        :param vw: wind vector (not used now)
        """
        # TODO: I don't like this function here in the Bullet class
        # check for existence of directory. If it exists, append the filename
        if not os.path.isdir(dpath):
            raise BulletException(
                "{}: directory {} does not exist".format(self._name,dpath)
            )

        fpath = os.path.join(dpath,"{}m range-tbl {}.tsv".format(inc,self._name))

        # run trajectory and get interpolations
        self.trajectory(height=h,elev_angle=fa,wind=vw)
        if not self._ss:
            print("{}: No trajectory for height={}, phi-{}".format(
                self._name,h,fa)
            )
            return
        r2t,r2v,r2pos,_ = self.interpolate()

        # set up max distance
        maxd = int(self._ss[-1][2][0] // inc * inc)

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
                *[(x,y,z) for x,y,z in [pos for _,_,pos,_ in self._ss]]
            )
            vxs,vys,vzs = zip(
                *[(x,y,z) for x,y,z in [vel for _,_,_,vel in self._ss]]
            )
            ts, vs = zip(*[(x, y) for x, y, _, _ in self._ss])
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
        Uses
         Lahti 2019, eq 2, pg 39
         p = p0 * cbrt(v/v0)
          where
           p0 = initial spin rate (radians/s),
           v = current velocity (m/s), and
           v0 = initial velocity (m/s)
        TODO:
         update & annotate calculation of vpi
        """
        # 1. update directional components of velocity,
        # 2. calculate acceleration components,
        # 3. calculate then new directional components of velocity,
        # 4. update position vector (using average velocity), and
        # 5. update elapsed time
        # 6. update directional components of velocity
        # 7. calculate new velocity
        # 8. calculate new spin rate
        vA = self._va_(vw)
        vV_t = self._vV_t + vA * t
        self._vP_t += ((self._vV_t + vV_t) / 2) * t + 0.5 * vA * np.power(t,2)
        self._t += t
        self._v_t = np.linalg.norm(vV_t)
        self._p_t = self._p_0 * np.cbrt(self._v_t / self._v_0)
        self._vV_t = vV_t

    def _vv_3d_(self,phi=0.,theta=0.):
        """
        calculates the 3d components of current velocity w/ angle of fire phi
        and windage angle theta
        :param phi: angle of fire (degrees)
        :param theta: angle of fire along the z-axis (windage)
        :return: the components of velocity as a array [vx,vy,vz]
        """
        return np.array(
            [self._vx_t_(phi,theta),self._vy_t_(phi,theta),self._vz_t_(theta)],
            np.double
        )

    def _vb_3d_(self,phi=0,theta=0.,alpha=0.,beta=0.):
        """
         calculates the three unit vectors x, y and z of the bullet's axes
        :param phi: elevation angle (degrees)
        :param theta: windage angle (degrees)
        :param alpha: pitch angle
        :param beta: yaw angle
        :return: x, y and z
        uses McCoy 1999, eq 9.23 (x), 9.24 (y), 9.25 (z) and 9.26 (Q)
         The x,y,z unit vectors form a triad with z perpendicular to x and lying
         in the horizontal plane. y is thend defined as the cross product
          y = z X x
        """
        # TODO: for z unit vector eq 9.25, is the 2nd component 0? looks like
        #  an 'O'
        # get the reciprical of the square root of Q
        Q = 1/np.sqrt(self._Q_(phi,theta,alpha,beta))

        # calculate bullet's x unit vector
        bx = np.array(
            [
                cos(phi+alpha)*cos(theta+beta),
                sin(phi+alpha)*cos(theta+beta),
                sin(theta+beta)
            ],np.double
        )

        # and y unit vector
        by = Q * np.array(
            [
                -np.power(cos(theta),2)*sin(phi,alpha)*cos(phi+alpha),
                np.power(cos(theta+beta),2)*np.power(cos(phi+alpha),2) +
                    np.power(sin(theta+beta),2),
                -sin(theta+beta)*cos(theta+beta)*sin(phi+alpha)
            ],np.double
        )

        # finally the z unit vector
        bz = Q * np.array(
            [-sin(theta+beta),0,cos(theta+beta)*cos(phi+alpha)],np.double
        )

        # return
        return bx,by,bz

    def _vx_t_(self,phi,theta=0.):
        """
        calculates the x component of current velocity at angle of fire phi
        and windage angle theta
        :param phi: angle of fire (degrees)
        :param theta: angle of fire along the z-axis (windage)
        :return: the x component of current velocity
        """
        return self._v_t*cos(phi)*cos(theta)

    def _vy_t_(self,phi,theta=0.):
        """
        calculates the y component of current velocity at angle of fire phi
        and windage angle theta
        :param phi: angle of fire (degrees)
        :param theta: angle of fire along the z-axis (windage)
        :return: the y component of current velocity
        """
        return self._v_t*sin(phi)*cos(theta)

    def _vz_t_(self,theta=0):
        """
        calculates the z component of current velocity at windage angle theta
        :param theta: angle of fire along the z-axis (windage)
        :return: the z component of current velocity
        """
        return self._v_t*sin(theta)

    def _za_k_(self,fs,x):
        """
        estimates the firing angle required to to zero at x using a constant
        drag coefficient k
        :param fs: front sight height over center of bore (m)
        :param x: desired zero distance (m)
        :return: vx,t,za
          vx = calculate x-component of velocity at distance x
          t = time of flight to x
          za = angle of gun over tangential horizon to zero at x
         uses McCoy 1999, eq. 5.45, pg 92
         theta = arctan(-fs - 0.5*(g * t*2 * (0.5 + s1 - s2) )) / x)
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
        v_x = self._v_0 * np.exp(-k * x)
        t = ((x / self._v_0) * (self._v_0 / v_x - 1)) / np.log(self._v_0 / v_x)
        dv = self._v_0 / v_x
        _1 = np.power(dv - 1, -1)
        _2 = np.power(dv - 1, -2) * np.log(dv)
        za = arctan((fs * bls.MM2M + (0.5 * atm.G * np.power(t, 2)) * (0.5 + _1 - _2)) / x, 'r')
        return v_x, t, r2d(za)

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

    def _vx_k3_(self,x,phi):
        """
        calculates x-component of velocity at range x fired from phi when
        using k3 drag coefficient
        :param x: range to estimate (m)
        :param phi: angel of fire
        :return:
        vx = (sqrt(v_x0) - 0.5*k3*x)^2 McCoy 1999 eq 5.67 pg. 94
          where
           where v_x0 is the x component of velocity at angle phi (m/s),
           k3 = drag coefficient (see _k3_), and
           x = range to estimate (m)
        """
        return np.power(np.sqrt(self._vx_t_(phi)) - (0.5*self._k3_()*x),2)

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
          Meaning the vectors must be subtracted
          Because mFd in force.py is negative
        """
        return force.Fd(self) * (self._vV_t-vw) + force.Fg()

    def _Q_(self,phi,theta,alpha,beta):
        """
         implements McCoy 1999, eq 9.26 pg 192 for Q
        :param phi: elevation angle
        :param theta: windage angle
        :param alpha: pitch angle
        :param beta: yaw angle
        :return: Q
         Q = sin^2(phi + beta) + cos^2(phi+beta)*cos^2(theta,alpha)
         where
          phi = elevation angle (degrees),
          theta = windage angle (degrees),
          alpha = pitch angle (degrees), and
          beta = yaw angle (degrees)
        """
        return  np.power(sin(theta+beta),2) + (
                np.power(cos(theta+beta),2) *
                np.power(cos(phi+alpha),2)
        )

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
         Given the OABL & drag model, estimates ogival, body & frustum lengths,
         the frustum base diameter and the ogival radius
         std model ratios. Estimates:
          1. the diameter of the frustum if present,
          2. the volume of the bullet,
          3. the Center of Gravity (CG) of the bullet, and
          4. the moments of inertia (axial & transverse)
         NOTE: Requires that bullet has a known diameter, length and drag model
        """
        # throw exception for non-tangent ogive models until calculated
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
                (np.power(f,2)*(f-1)*arcsin(nbl/f,'r'))
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
        # rc *= bls.MM2M
        # lc *= bls.MM2M
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
        Ixo = np.pi * (b0*(np.power(f,2)*arcsin(nbl/f) - (f-1)*nbl) -
                       2*np.power(nbl,3) * (b1/3 + b2/4 + b3/5))
        Ixo *= p*np.power(rc,5)

        # and Seglates 2019 Eq. 6
        # Iyogive/P*R^5 = pi/4 * ((f^2*nbl*(f^2+7/2*(f-1)^2)) +
        #                        1/15*nbl^5 -
        #                        f^2*(f-1)*(5/2*f^2 + 2*(f-1)^2)*arcsin(nbl/f)
        _1 = np.power(f,2) * nbl * (np.power(f,2) + (7/2)*np.power(f-1,2))
        _2 = np.power(nbl,5) / 15
        _3 = np.power(f,2) * (f-1) * (
                ((5/2) * np.power(f,2) + 2*np.power(f-1,2)) * arcsin(nbl/f)
        )
        Iyo = (np.pi/4) * (_1 + _2 - _3)
        Iyo *= p*np.power(rc,5)

        # sum the moments
        self._Ix = Ixf + Ixc + Ixo
        self._Iy = Iyf + Iyc + Iyo
