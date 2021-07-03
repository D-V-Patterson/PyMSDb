#!/usr/bin/env python
""" atmosphere.py
Copyright (C) 2021 Dale Patterson (dale.v.patterson@@mail.com)

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

Defines atmospheric factors
"""

#__name__ = 'atmosphere'
__license__ = 'GPLv3'
__version__ = '0.0.3'
__date__ = 'May 2021'
__author__ = 'Dale Patterson'
__maintainer__ = 'Dale Patterson'
__email__ = 'dale.v.patterson@@mail.com'
__status__ = 'Development'

import numpy as np
import random
from scipy.interpolate import interp1d
import pymsdb.ballistics as bls
from pymsdb import PyMSDbException

"""
 TODO:
  2. calculate density altitude
  3. calculate mach given atmospheric parameters as well given the interpolations
   a2c, t2c, p2c defined below
"""

# CONSTANTS

# Air density from : see http://www.wind101.net/air-density/air-density-calculator.htm
# atmospheric found in ICAO pages
# Coriolis effect is apx 2*pi / 24*3600 (radians/sec)
G      = 9.80665   # acceleration of gravity (m/s)
BARP   = 101.325   # barametric pressure (kPa)
ALT    = 0.0       # sea level (m)
TEMP   = 15        # temperature (C)
HUM    = 50        # relative humidity (%)
RHO    = 1.221     # air density at sea level and 20 celsius (kg/m^3)
ABSZ   = 273.15    # absolute zero = -ABSZ (K)
SOUND  = 340.292   # speed of sound (m/s) at above parames (alt.,temp.,pressure)
OMEGA  = 7.722e-05 # angular velocity of earth for calculating coriolis (radians)

#### SPEED OF SOUND

def mach(vs): # TODO: add atmospheric conditions as an option (or at least y_i)
    """ returns the mach number of velocity(s) vs vs/Speed of Sound """
    return vs/SOUND

def atmospherics(y):
    """
    returns nominal atmopsheric values at elevation y
    :param y: elevation (m) above/below sea level = 0
    :return: a tuple t = (T,P,D,C) where
     T = tempature (C)
     P = barametric pressure (kPa)
     D = air density (kg/m^3)
     C = Speed of Sound (m/s)
    """
    return float(_isa_stdt_(y)),\
           float(_isa_stdp_(y)),\
           float(_isa_stdd_(y)),\
           float(_isa_stdc_(y))

#### WIND RELATED
"""
     Adopting the coordinate system of McCoy 1999, Figure 5.1, 7.1 for wind

        y=elevation
        |
        |
        |
        |
        |
        |
       S|_____________x=range  T (12 o'clock)
       /(0,0,0)
      /
     /
    /
   /
   z=windage (3 o'clock)
  
  and viewing the shooter top-down as a clock
  
        9
        |
  6-----S-----12 T
        |
        3 
        
 where S is the shooter, T is the target  
   
   For wind vectors W = (v_x,v_y,v_z) let S = (s_x,s_y,s_z) = (0,0,0) be the 
   shooter's position and T = (t_x,t_y,t_z) be the target's position then
    v_x (the x-component of the wind velocity) is positive for tailwinds blowing 
     from Sx to Tx (or 6 o'clock to 12 o'clock) and negative for headwinds 
     blowing from Sx to Tx (or 12 o'clock to 6 o'clock)
    v_y (the y-component of velocity) is positive when blowing from Sy to y > Sy
     and negative when blowing from y > Sy to Sy
    v_z (z-component of the wind velocity) is positive when blowing 9 o'clock
     to 3 o'clock and negative otherwise 
 
 From http://gyre.umeoce.maine.edu/data/gomoos/buoy/php/variable_description.php?
 variable=wind_2_speed
     
 TODO:
  wind should also have a position vector
"""

# From http://gyre.umeoce.maine.edu/data/gomoos/buoy/php/variable_description.php?
#  variable=wind_2_speed
# NOTE: their ranges do meet at defined borders i.e. light is 0.5 to 1.5 and
# the next, gentle is 2 to 3. We use the bottom end to define the categories
# so that for example, a fresh gale is (WFGLE = 9) >= 20.5 and < 24.0, a strong
# gale
# all winds speeds are in m/s
def wspd(l,h):
    """
    generate a random wind speed between l and h
    :param l: lower defined wind speed code below
    :param h: upper defined wind speed code below
    :return: a random wind speed (m/s) rounded to 3 decimals
    """
    if l < WNONE or l > WHRCN:
        raise PyMSDbException("Invalid 'lower category".format(l))
    elif h < WNONE or h > WHRCN:
        raise PyMSDbException("Invalid 'higher category".format(l))
    return np.round(random.uniform(wind_speeds[l],wind_speeds[h]),3)

WNONE =  0 # no wind
WCALM =  1 # calm
WLAIR =  2 # light air
WLBRZ =  3 # light breeze (this is roughly 10mph winds)
WGBRZ =  4 # gentle breeze
WMBRZ =  5 # moderate breeze
WFBRZ =  6 # fresh breeze
WSBRZ =  7 # strong breeze
WMGLE =  8 # moderate gale
WFGLE =  9 # fresh gale
WSGLE = 10 # strong gale
WWGLE = 11 # whole gale
WSTRM = 12 # storm
WHRCN = 13 # hurricane
wind_speeds = np.array(
    [0.0,0.5,2.0,3.5,5.5,8.5,11.0,14.0,17.0,20.5,24.0,28.0,32.0,],np.double
)

def crosswind(wz):
    """
     generates a constant crosswind w_x = w_y = 0
    :param wz: wind speed of z-component
    :return: returns a crosswind wind vector at velocity wz
    """
    return np.array([0,0,wz],np.double)

def rangewind(wx):
    """
     generates a constant rangewind w_y = w_z = 0
    :param wx: wind speed of z-component
    :return: returns a rangewind wind vector at velocity wx
    """
    return np.array([wx,0,0],np.double)

def deflection(vw,x,t,vv):
    """
     calculates crosswind deflection (crosswind acting constantly everywhere
     on the pm-trajectory)
    :param vw: wind vector such that x = y = 0 and z != 0 (m/s)
    :param x: range to calculate deflection at (m)
    :param t: time of flight (s)
    :param vv: the initial velocity vector (m/s)
    :return: deflection caused by crosswind
    uses McCoy 7.27
     Z = w_z*(t - x/v_xo)
     where
     w_z = z-component of wind velocity (m/s),
     t = time of flight to x (s),
     x = range (m), and
     v_x0 = initial x-component of bullet's velocity
    """
    return vw[bls.IZ] * (t - x/vv[bls.IX])

def deflection_r():
    """
    deflection at range r
    :return:
    uses McCoy 1999, eq 7.28 pg 160
    z(r) = w_zi * [t(r) - t(x_i) - (r-x_i)/v_xi]
    where
     w_zi = z-component of wind velocity at time i (m/s),
     t(r) = time for bullet to reach distance r (seconds),
     t(x_i) = time for bullet to reach distance x_i (seconds),
     r = final distance (m),
     x_i = distance (from shooter) at time i (seconds), and
     v_xi = x-component of bullet velocity at time i
    NOTE:
     McCoy uses a table
     see also eq 7.30 which calculates th crosswind weighting factor f_w_zi
    """
    # pass

### ATMOSPHERIC ENVIRONMENTAL FACTORS

# TODO: may not add anything over using interpolated values from ISA below
def temperature(y,t0):
    """
    estimates the temperature at elevation y given the temperature at the
    firing elevation
    :param y: elevation (m)
    :param t0: temperature at firing location y = 0
    :return: estimated temparture at elevation y (Celsius)
    uses McCoy 1999, eq 8.10 pg 166
     t_y = (t_0 + 273.15)e^(K*y) - 273.15
     where
      t_0 = temperature at firing elevation (C),
      273.15 = Kelving to Celsius constant,
      K = temperature-altitude decay factor, and
      y = elevation to estimate (m)
      since the constant K is expressed in feet we'll use a simpler formula
       t_o - 0.0065*y
        where y = (y - y_f)
    """
    # K = 6.858e-6 + 2.776e-11
    #return (t0 + ABSZ)*np.exp(K*y) - ABSZ
    # NOTE: that the below eq. given altitude of 1000m is approximate to
    # _isa_stdt_(1000)
    return t0 - 0.0065*y # (-0.0065 is the temp. lapse rate in meters)

### MISCELLANEOUS

def vpw(t):
    """
    calculates the vapor pressure of water at temperature t
    :param t: tempature (C)
    :return: Vpw (kPa)
     Buck 1996 (Buck Research Panel)
     Vpw = 6.1121e^(18.678-t/234.5)*t/(257.14+t))
    """
    return 0.61121*np.exp((18.678-(t/234.5))*(t/(257.14-t)))

def stationp(p,a):
    """
    calculates the station barometric pressure
    :param p: barometric pressure (kPa)
    :param a: altitude above sea level (m)
    :return: station barometric pressure (kPa)
    https://panthaera.com/
     the-effects-of-altitude-and-barometric-pressure-on-a-bullet-during-flight/
     SP = BP - (A/9)
      where
       BP = barometric pressure (hPa),
       a = altitude (m)
       NOTE: we convert the given pressure to hPa (x 10) then convert the
        resulting pressure back to kPa (/ 0.1)
    """
    return ((p*10) - (a/9))*0.1

#### THE BELOW COME FROM gnu-ballistics

# TODO:not really necessary for current application
def BCa(BC,a,p,t):
    """
    adjusted BC based on atmoshperic conditions
    :param BC: bullet's ballistic coefficient
    :param a: altitude (m) above sea level
    :param p: barometric pressure (kPa)
    :param t: temperature (C)
    :return: correctd ballistic coefficient (unitless)
     NOTE tried BCa = BC*(Af*(1.0+Tf-Pf)*Hf) from
      https://www.sierrabullets.com/exterior-ballistics/
        5-1-effects-of-altitude-and-atmospheric-conditions/
      but results were way off (possibly error on my part)
      currently using https://riflebarrels.com/air-density-and-bullet-performance/
       BCa = BC * Tf * Pf * Af
      where
       BC = ballistic coefficient,
       Tf = the temperature adjustment factor,
       Pf = the bar. pressure adjustment factor, and
       Af = the altitude adjustment factor
      this equation does not consider humidity
    """
    return BC*_Tf_(t,a)*_Pf_(p,a)*_Af_(a)

#### PRIVATE HELPERS

def _Tf_(t,a):
    """
     temperature adjustment factor
    :param t: temperature (C)
    :param a: altitude (m)
    :return: temp. adjustment factor
     Tf = 1 + (T-Tstd)/(Tstd+FP)
      where
       T is the tempature,
       Tstd is the standard tempature at altitude,
       FP is the freezing point of water (K)
     https://riflebarrels.com/air-density-and-bullet-performance/
    """
    return 1+(t-_isa_stdt_(a))/(_isa_stdt_(a)+ABSZ)

def _Pf_(p,a):
    """
     barometric pressure adjustment factor
    :param p: barametric pressure (kPa)
    :param a: altitude (m)
    :return: pressure adjusstment factor
     Pf = 1 - (P-Pstd)/Pstd
      where
       P is the barometric pressure and
       Pstd is the baometric pressure at altitude
       https://riflebarrels.com/air-density-and-bullet-performance/
    """
    return 1-(p-_isa_stdp_(a))/_isa_stdp_(a)

# CURRENTLY NOT USED
#def _Hf_(h,p,t):
#    """
#     relative humidity adjustment factor
#    :param h: relative humidity (0.0 < h <= 1.0)
#    :param p: barametric pressure (kPa)
#    :param t: temperature (C)
#    :return: rel. humidity adjustment factor
#     Hf = 0.9950*(P/(P*0.3783*H*Vpw))
#      where
#       P is the barometric pressure,
#       H is the relative humidity
#       Vpw is the vapor pressure of water at temperature
#       https://www.sierrabullets.com/exterior-ballistics/
#        5-1-effects-of-altitude-and-atmospheric-conditions/
#    """
#    return 0.9950*(p/(p*0.3783*h*vpw(t)))

def _Af_(a):
    """
    altitude adjustment factor
    :param a: altitude (m) above/below sea level = 0
    :return: alt. adjustment factor
    Af = 1/e^(h*a)
     where h = 0.00003158 and
     a = altitude
     https://riflebarrels.com/air-density-and-bullet-performance/
    """
    # TODO is this for imperial units? converting the constant h to meter
    #  the results are close
    return 1/np.exp(0.00003158*a)

# International Standard Atmosphere
# (https://www.engineeringtoolbox.com/international-standard-atmosphere-d_985.html)
# columns are Altitude (m), Temperature (C), Pressure (kPa), Density (kg/m3),
#  Kin. Viscosity (m2/s), Thermal Conductivity (W/m K) & Speed of Sound (m/s)
# Each column is a seperate numpy array. Somewhat Overkill as we could stop at
# 3000 - 5000 m (Everest is 8.8 km) and denali is highest pt in US at 6.2K
# TODO: should we try to pickle these and open as used?
_as_ = np.array( # altitude
    [
        -2000.0, -1500.0, -1000.0, -500.0, 0.0, 500.0, 1000.0, 1500.0, 2000.0,
        2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0, 5500.0, 6000.0, 6500.0,
        7000.0, 7500.0, 8000.0, 8500.0, 9000.0, 9500.0, 10000.0, 10500.0, 11000.0,
        11500.0, 12000.0, 12500.0, 13000.0, 13500.0, 14000.0, 14500.0, 15000.0,
        15500.0, 16000.0, 16500.0, 17000.0, 17500.0, 18000.0, 18500.0, 19000.0,
        19500.0, 20000.0, 22000.0, 24000.0, 26000.0, 28000.0, 30000.0
    ],np.double
)
_ts_ = np.array( # temperature
    [
        28.05, 24.75, 21.55, 18.25, 15.0, 11.75, 8.55, 5.25, 2.05, -1.25, -4.45,
        -7.75, -10.95, -14.25, -17.45, -20.75, -23.95, -27.25, -30.45, -33.65,
        -36.95, -40.15, -43.45, -46.65, -49.85, -53.15, -56.35, -56.45, -56.45,
        -56.45, -56.45, -56.45, -56.45, -56.45, -56.45, -56.45, -56.45, -56.45,
        -56.45, -56.45, -56.45, -56.45, -56.45, -56.45, -56.45, -54.55, -52.55,
        -50.65, -48.65, -46.65
    ],np.double
)
_ps_ = np.array( # pressure
    [
        127.78, 120.7, 113.93, 107.48, 101.325, 95.46, 89.88, 84.56, 79.5, 74.69,
        70.12, 65.78, 61.66, 57.75, 54.05, 50.54, 47.22, 44.08, 41.11, 38.3, 35.65,
        33.15, 30.8, 28.58, 26.5, 24.54, 22.7, 20.98, 19.4, 17.93, 16.58, 15.33,
        14.17, 13.1, 12.11, 11.2, 10.35, 9.572, 8.85, 8.182, 7.565, 6.995, 6.467,
        5.98, 5.529, 4.047, 2.972, 2.188, 1.616, 1.197
    ],np.double
)
_ds_ = np.array( # density
    [
        1.4782075, 1.411445, 1.34701, 1.2849025, 1.225, 1.1673025, 1.1116875,
        1.058155, 1.0065825, 0.95697, 0.9093175, 0.86338, 0.8194025, 0.7770175,
        0.73647, 0.697515, 0.6601525, 0.62426, 0.5900825, 0.5572525, 0.52577,
        0.4957575, 0.4670925, 0.4396525, 0.41356, 0.38857, 0.364805, 0.3374875,
        0.311885, 0.288365, 0.26656, 0.24647, 0.22785, 0.2107, 0.194775, 0.180075,
        0.1664775, 0.15386, 0.142345, 0.131565, 0.1216425, 0.1124795, 0.10399025,
        0.0961625, 0.0889105, 0.0645085, 0.046942, 0.03426325, 0.02507575,
        0.01841175
    ],np.double
)
_kvs_ = np.array( # kinematic viscosity # TODO: don't think this is necessary
    [
        1.253, 1.301, 1.352, 1.405, 1.461, 1.52, 1.581, 1.646, 1.715, 1.787,
        1.863, 1.943, 2.028, 2.117, 2.211, 2.311, 2.416, 2.528, 2.646, 2.771,
        2.904, 3.046, 3.196, 3.355, 3.525, 3.706, 3.899, 4.213, 4.557, 4.93,
        5.333, 5.768, 6.239, 6.749, 7.3, 7.895, 8.54, 9.237, 9.99, 10.805, 11.686,
        12.639, 13.67, 14.784, 15.989, 22.201, 30.743, 42.439, 58.405, 80.134
    ],np.double
)
_tcs_ = np.array( # thermal conductivity # TODO: don't think this is necessary either
    [
        2.636, 2.611, 2.585, 2.56, 2.534, 2.509, 2.483, 2.457, 2.431, 2.405, 2.379,
        2.353, 2.327, 2.301, 2.275, 2.248, 2.222, 2.195, 2.169, 2.142, 2.115,
        2.088, 2.061, 2.034, 2.007, 1.98, 1.953, 1.952, 1.952, 1.952, 1.952, 1.952,
        1.952, 1.952, 1.952, 1.952, 1.952, 1.952, 1.952, 1.952, 1.952, 1.952, 1.952,
        1.952, 1.952, 1.968, 1.985, 2.001, 2.018, 2.034
    ],np.double
)
_cs_ = np.array( # Speed of Sound
    [
        347.9, 346.0, 344.1, 342.2, 340.3, 338.4, 336.4, 334.5, 332.5, 330.6,
        328.6, 326.6, 324.6, 322.6, 320.5, 318.5, 316.5, 314.4, 312.3, 310.2,
        308.1, 306.0, 303.8, 301.7, 299.8, 297.4, 295.2, 295.1, 295.1, 295.1,
        295.1, 295.1, 295.1, 295.1, 295.1, 295.1, 295.1, 295.1, 295.1, 295.1,
        295.1, 295.1, 295.1, 295.1, 295.1, 296.4, 297.7, 299.1, 300.4, 301.7
    ],np.double
)
_isa_stdt_ = interp1d(_as_,_ts_,kind='linear',fill_value='extrapolate')
_isa_stdp_ = interp1d(_as_,_ps_,kind='linear',fill_value='extrapolate')
_isa_stdd_ = interp1d(_as_,_ds_,kind='linear',fill_value='extrapolate')
#_isa_stdkv_ = interp1d(_as_,_kvs_,kind='linear',fill_value='extrapolate')
#_isa_stdtc_ = interp1d(_as_,_tcs_,kind='linear',fill_value='extrapolate')
_isa_stdc_ = interp1d(_as_,_cs_,kind='linear',fill_value='extrapolate')

