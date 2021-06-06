#!/usr/bin/env python
""" taxonomy.py
Copyright (C) 2020 Dale V. Patterson (dale.v.patterson@gmail.com)

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

Defines classification of firearms
"""

from pymsdb import PyMSDbException

#__name__ = 'taxonomy'
__license__ = 'GPLv3'
__version__ = '0.0.1'
__date__ = 'April 2021'
__author__ = 'Dale Patterson'
__maintainer__ = 'Dale Patterson'
__email__ = 'dale.v.patterson@gmail.com'
__status__ = 'Development'

"""
For purpose of rifle sub-classification
 Barrel Length
  Pistol < 10
  Carbine 16 - 20
  Rifle >= 20
  Could add mid (14-20)?
  
  How to do a tactical rifle that is a carbine i.e. Ruger AR-556 (16.1" barrel)?
"""
# firearm type levels
#  0 firearm
#    1 long gun
#      2 shotgun
#       3 SBS
#      2 rifle
#       3 carbine
#       3 SBR
#        4 ar-15 rifle (5.56)
#        4 ar-10 rifle (7.62)
#        4 AKM rifle
#   1 handgun
#      2 revolver
#      2 derringer
#      2 pistol
#       3 Pocket
#       3 Sub Compact
#       3 Compact
#       3 Full
#       3 LBP # TODO: need a better name
#        4 SMG
#        4 ar-15 pistol

FM_LVL0    =  ['Firearm']
FM_LVL1    = ['Handgun','Long gun']
FM_LVL2_HG = ['Derringer','Revolver','Pistol']
FM_LVL2_LG = ['Shotgun','Rifle']
FM_LVL2    = FM_LVL2_HG + FM_LVL2_LG
FM_LVL3_PS = ['Pocket','Sub-compact','Compact','Full','Large']
FM_LVL3_RF = ['Carbine','SBR']
FM_LVL3    = FM_LVL3_PS + FM_LVL3_RF

def parent_of(f):
    # return the parent of f
    if f == 'Firearm': return None
    if f in FM_LVL1: return 'Firearm'
    elif f in FM_LVL2:
        if f in FM_LVL2_HG: return 'Handgun'
        elif f in FM_LVL2_LG: return 'Long gun'
    elif f in FM_LVL3:
        if f in FM_LVL3_PS: return 'Pistol'
        elif f in FM_LVL3_RF: return 'Rifle'
    else:
        raise PyMSDbException("Invalid Arg 'firearm type' ({})".format(f))

"""
Barrel lengths from current set of known firearms in msdb
 pistol (with SMG)    min = 2.75, max = 10.2
 pistol (exclude SMG) min = 2.75, max = 5.75)
   SMGS
   Tec-9 Mini 2.99	9.49	9.49	2.71  (both Intratec & Navegar)
   M-11       5.08	9.76	20.9	3.5
   Uzi        10.2	17	    25      7.7
   M-10       5.75	10.7	21.57   6.26
   Uzi CB     16.14	26	    26	    7.7
 rifle min = 16.0 max = 24 
 shotgun min = 14 max = 28
 
OAL from current set of known firearms in msdb
 pistol (with SMG) min = 5.16, max = 25.0
 pistol (w/o SMG) min = 5.16, max = 9
 rifle min = 26, max = 43.63
 shotgun min = 26.3, max = 48
"""

# LENGTH CONSTANTS
PST_PKT_BRL = 3
PST_SUB_BRL = 3.5
PST_CMP_BRL = 4.5
PST_FLL_BRL = 6.2  # full size pistol, anything over for now is LBP
HG_BRL      = 10.2 # split point between handgun and long gun
SBR_BRL     = 16   # short barrel rifle barrel length
SBR_OAL     = 26   # short barrel rifle overall length
RFL_CRB_BRL = 20   # carbine barrel length
SBS_BRL     = 18   # short barrel shotgun barrel length
SBS_OAL     = 26   # short barrel shotgun overall length

# NFA Fire Actions
NFA_FA = ['FA','SA-M']

# below are based on the type
firearm_root = 'Firearm'
firearm_lvl_2 = ['Long gun','Handgun']
longgun_lvl_3 = ['Shotgun','Rifle']
handgun_lvl_3 = ['Revolver','Pistol','Derringer']
rifle_lvl_4 = ['Tactical Rifle','Tactical Carbine','Carbine']
tac_rifle_lvl_5 = ['AKM Rifle','AR-15 Rifle','AR-10 Rifle']
pistol_lvl_4 = ['Tactical Pistol']
tac_pistol_lvl_5 = ['SMG','AR-15 Pistol']
firearm_types = [firearm_root] + \
            firearm_lvl_2 + \
            longgun_lvl_3 + \
            rifle_lvl_4 + \
            tac_rifle_lvl_5 + \
            handgun_lvl_3 + \
            pistol_lvl_4 + \
            tac_pistol_lvl_5

# tactical firearm categories and caliber
tac_firearm = ['AR-15 Pistol','SMG','AKM Rifle','AR-15 Rifle','AR-10 Rifle']
tac_firearm_cal = ['.223','.308','5.56x45','7.62x39','5.45x39','6.8mm SPC']
tac_firearm_fa = ['SA','SA-M','SF','FA']

# level 2 categories

def firearm_lvl_2(wtype):
    if is_handgun(wtype): return 'Handgun'
    elif is_longgun(wtype): return 'Long gun'
    else: return None

def is_handgun(wtype):
    return wtype in tac_pistol_lvl_5 + pistol_lvl_4 + handgun_lvl_3 + ['Handgun']

def is_longgun(wtype):
    return wtype in tac_rifle_lvl_5 + rifle_lvl_4 + longgun_lvl_3 + ['Long gun']

# level 3 categories

def firearm_lvl_3(wtype):
    if is_shotgun(wtype): return 'Shotgun'
    elif is_rifle(wtype): return 'Rifle'
    elif is_pistol(wtype): return 'Pistol'
    elif is_revolver(wtype): return 'Revolver'
    elif is_derringer(wtype): return 'Derringer'
    else: return None

def is_shotgun(wtype): return wtype == 'Shotgun'

def is_rifle(wtype): return wtype in tac_rifle_lvl_5 + rifle_lvl_4 + ['Rifle']

def is_revolver(wtype): return wtype == 'Revolver'

def is_derringer(wtype): return wtype == 'Derringer'

def is_pistol(wtype): return wtype in pistol_lvl_4 + tac_pistol_lvl_5 + ['Pistol']

# level 4 categories

def firearm_lvl_4(wtype):
    if is_carbine(wtype): return 'Carbine'
    elif is_tac_rifle(wtype): return 'Tactical Rifle'
    elif is_tac_pistol(wtype): return 'Tactical Pistol'
    else: return None

def is_carbine(wtype): return wtype == 'Carbine'

def is_tac_rifle(wtype): return wtype in tac_rifle_lvl_5 + ['Tactical Rifle']

def is_tac_pistol(wtype): return wtype in tac_pistol_lvl_5 + ['Tactical Pistol']

# level 5 categories

def firearm_lvl_5(wtype):
    if wtype in tac_pistol_lvl_5: return wtype
    elif wtype in tac_rifle_lvl_5: return wtype
    else: return None

def is_tac_firearm(wtype): return is_tac_rifle(wtype) or is_tac_pistol(wtype)

def lcd(wtypes):
    # remove duplicates, then check if all types are the same
    ws = list(set(wtypes))
    if len(ws) == 1: return ws[0]

    # check for 5th level commonness
    ws5 = [firearm_lvl_5(w) for w in ws]
    if len(list(set(ws5))) == 1 and ws5[0] is not None: return ws5[0]

    # check for 4th level commonness
    ws4 = [firearm_lvl_4(w) for w in ws]
    if len(list(set(ws4))) == 1 and ws4[0] is not None: return ws4[0]

    # check for 3rd level commonness
    ws3 = [firearm_lvl_3(w) for w in ws]
    if len(list(set(ws3))) == 1 and ws3[0] is not None: return ws3[0]

    # check for 2nd level commonness
    ws2 = [firearm_lvl_2(w) for w in ws]
    if len(list(set(ws2))) == 1 and ws2[0] is not None: return ws2[0]

    return 'Firearm'

def firearm_cat_lvl(wtype):
    if wtype == 'Firearm': return 1
    elif wtype == 'Long gun' or wtype == 'Handgun': return 2
    elif wtype in ['Shotgun','Rifle','Revolver','Derringer','Pistol']: return 3
    elif wtype in ['Carbine','Tactical Rifle','Tactical Pistol']: return 4
    elif wtype in ['SMG','AR-15 Rifle','AR-10 Rifle','AR-15 Pistol']: return 5
    return -1
