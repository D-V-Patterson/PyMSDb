#!/usr/bin/env python
""" firearmspec.py
Copyright (C) 2020 Dale V. Patterson (dale.v.patterson@gmail.com)

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

Defines the FirearmSpec class - a wrapper around a firearm specification dict

a FirearmSpec is defined as:
 {
    MANUF
	FULL MODEL
	MODEL
	VAR/CONF/GEN/SER
	CATEGORY
	TYPE
	SUB-TYPE
	FIRE_MECH
	CALIBER
	BRL-LEN
	MIN-LEN
	MAX-LEN
	WEIGHT (lbs)
	TGR-PULL
	STOCK
	GRIP
	ACTION
	TWIST
	FSH
	RSH
	SR
	CYCLIC ROF
	EFF ROF
	EFF-RNG-AREA (m)
	EFF-RNG-PT (m)
	MOA
	FEED-SYS
	STD-CAP
	NOTES
}
"""

#__name__ = 'firearmsepc'
__license__ = 'GPLv3'
__version__ = '0.0.2'
__date__ = 'May 2021'
__author__ = 'Dale Patterson'
__maintainer__ = 'Dale Patterson'
__email__ = 'dale.v.patterson@gmail.com'
__status__ = 'Development'

import csv
from pymsdb.data import pth_firearmspec
from pymsdb import PyMSDbException

# constants for file headers
FS_MANUF        =  0
FS_FULL_MDL     =  1
FS_MODEL        =  2
FS_SERIES       =  3
FS_CATEGORY     =  4
FS_TYPE         =  5
FS_SUB_TYPE     =  6
FS_FIRE_MECH    =  7
FS_CALIBER      =  8
FS_BRL_LEN      =  9
FS_LLEN         = 10
FS_HLEN         = 11
FS_WEIGHT       = 12
FS_TPULL        = 13
FS_STOCK        = 14
FS_GRIP         = 15
FS_ACTION       = 16
FS_TWIST        = 17
FS_FSH          = 18
FS_RSH          = 19
FS_SR           = 20
FS_CYCLIC_ROF   = 21
FS_EFF_ROF      = 22
FS_EFF_RNG_AREA = 23
FS_EFF_RNG_PT   = 24
FS_MOA          = 25
FS_FEED_SYS     = 26
FS_FEED_LOC     = 27
FS_STD_CAP      = 28
FS_CLONE        = 29
FS_NOTES        = 30

# fire action types
fire_acts = [
    'BA',   # bolt-action
    'LA',   # lever action, applies only to riles i.e. repeater
    'SA',   # semi-automatic (applies to rifle, pistol, shotgun as does below)
    'SA-M', # semi-automatic, modified in some way to exceed SA & achieve FA or near FA
    'FA',   # full-automatic
    'SAO',  # single action only applies to some revolvers, pistols
    'PA',   # applies only to shotguns
    'SS',   # single shot a derringer or in some cases a shotgun
    'ML',   # muzzle loader
    'SF',   # select fire
]

# allowable stock
stocks = [
    'N/A',         # i.e. handguns
    'None',        # stock does not come standard with firearm (different from N/A)
    'Brace',       # ar pistols
    'Folding',     # stock folds onto side of firearm
    'Fixed',       # stock is part of the 'body' of the firearm i.e. hunting rifles
    'Collapsible', # stock can be adjusted i.e. AR-15s to varying lengths
]

# allowable grip
grips = [
    'HG', # catchall for a handgun (revolver, pistol) grip
    'TH', # thumbhole, fixed stock connected to grip with hole/section for the thumb
    'A2', # catchall for pistol grip on a rifle
    'MC', # montecarlo, the 'heel' is low while the 'comb' is high
    'BH', # birdshead grip
]

# Feed-System is a colon separated field of the form:
#  feed-location : feed unload process : feed load process
#   where feed unload process and feed reload process can be a comma seperated list
# no spaces allowed in fields
# allowable feed locations
feed_sys_loc = [
    'mag-well',       # empty magazine well where a detachable magazine can be inserted
    'mag-well-rot',   # empty magazine well using a detachable rotary mag
    'mag-well-floor', # box mag, fits inside the firearm, covered by a hinged floor
    'cylinder',       # revolvers have a cylinder with x # of chambers for ammo
    'int-box',        # box magazine internal to the firearm, part of it protruding
    'int-blind',      # internal magazine with a hinged floor
    'int-tube',       # internal tubular magazine
    'barrel',         # ammo is inserted directly into the barrel
]

# allowable unloads NOTE: most firearms can be 'unloaded' via cycling until all rounds are
#  are ejected - this is not specified/required to be listed except when necessary
feed_sys_unload = [
    'btn-release',     # push button release drops detachable magazine
    'lvr-release',     # pull lever release drops detachable magazine
    'drop-floor',      # internal box magazines that can be 'dropped', releasing all ammo
    'hinge-floor',     # same as above but blind mag, hinged floor is dropped, releasing all ammo
    'det-hinge-floor', # detachable mag, floor is dropped, mag is removed
    'rel-sng-btm',     # no floor plate, button/lever releases catch on ammo one at time, firearm is upside down
    'break-open',      # the firearm is 'break-opened', operator manually inserts ammo
    'eject-rod-sng',   # SAO revolvers (eject ammo one by one for each chamber)
    'eject-rod-all',   # DS revolvers (latch swings open cylinder, eject-rod drops all ammo
    'rm-tube-plngr',   # remove tube plunger, tilt weapon downward to drop ammo (SA w/ charging handle)
    'cycle-action',    # repeated cycle of action ejects all rounds (lever action)
]

# allowable_load
feed_sys_load = [
    'insert-lock',     # insert a detachable magazine with self-locking mechanism on firearm
    'insert-floor',    # detachable mag is attached to floor and floor is closed
    'manual-top',      # manually insert ammo by ammo into internal mag (box or blind) from top of firearm
    'manual-btm',      # manually insert ammo by ammo into internal mag/tube, firearm is upside down
    'manual-chmbr',    # manually insert ammo by ammo into individual chambers
    'manual-brl',      # manually insert ammo by ammo into barrel
    'manual-feed-rmp', # manually insert ammo by ammo into internal tube generally thru side of weapon
    'spd-loader',      # revolvers, use a speed loader to load all bullets
    'stripper-clip',   # a stripper clip inserted into mag from top of firearm
]

# feed locations (w.r.t to shooting hand grip)
feed_loc = [
    'forward',  # forward of shooting hand grip
    'inside',   # i.e. pistols inside the grip
    'behind',   # behind the shooting hand grip
    'internal', # tube, blind mags inside weapon, generally under barrel
]

fpath = 'data/firearm_specs.tsv'

"""
 # TODO:
   1. define values for stock, grip, action, feed system
   2. check each for nfa/fawb compliance (for example the COLT SP1 used during FAWB was legal
   3. verify for AR-10 and AR-15 if they can fire both calibers
   5. add oal as max of min and max lengths
   6. once done debugging, catch errors in read, print message and continue
   7. add handle, grip as list
   8. figure out terminology for grips, stocks
  11. need to break up where necessary by caliber, will have to insitute a key based 
   on model and caliber as well
  12. add front sight heigts, see if we can find data on rear sight height and
   sight radius
"""

# TODO: currently Colt Competition and beyond is undefined
#  No longer in msdb: 'Israeli Military Industries Uzi', 'Glock 19 Gen4', 'Glock 20 SF',
#  'Savage Arms Mark II', 'Browning BPS', 'DPMS Panther Classic 16', 'Glock 21 SF']

class FirearmSpecException(PyMSDbException):
    def __init__(self,msg): super().__init__(msg)

class FirearmSpec(object):
    def __init__(self):
        self._ks    = []  # ordered list of keys in firearm spec dict (_fs)
        self._fs    = {}  # firearm spec data
        self.read()

    def read(self):
        # reset data
        self._ks = []
        self._fs = {}

        # read file
        fin = None
        try:
            # open file, pass to csv reader and discard headers
            fin = open(pth_firearmspec)
            csvr = csv.reader(fin,delimiter='\t')
            _ = csvr.__next__()

            # iterate the data and save to dict, keys
            for row in csvr:
                # construct key should be no duplicates but make sure
                key = row[FS_MANUF] + ' ' + row[FS_FULL_MDL]
                if key in self._ks:
                    print("Found duplicate key {}, skipping...".format(key))
                    continue

                # read in row and set key lists if read
                if self._process_(key,row):
                    try:
                        self._ks.append(key)
                    except ValueError as e:
                        print("{}, skipping".format(e))
                else:
                    print("Error reading {}, skipping...".format(key))
        except IOError as e:
            print("Error reading {}: {}".format(fpath,e))
        finally:
            if fin: fin.close()

    def keys(self): return self._ks

    def get_key(self,manuf,model):
        # first try keys
        try:
            return self._ks[self._ks.index(manuf + ' ' + model)]
        except (TypeError,ValueError):
            pass
        return None

    def get_rec(self,key): return self._fs[key].copy()

    def get_spec(self,key,spec): return self._fs[key][spec]

#### PRIVATE HELPERS

    def _process_(self,key,row):
        # process row (for now, just reads in data
        self._fs[key] = {
            'manuf':row[FS_MANUF],
            'full-mdl':row[FS_FULL_MDL],
            'model':row[FS_MODEL],
            'series':row[FS_SERIES],
            'category':row[FS_CATEGORY],
            'type':row[FS_TYPE],
            'sub-type':row[FS_SUB_TYPE],
            'fire-mech':row[FS_FIRE_MECH],         # TODO: check for validity
            'caliber':row[FS_CALIBER].split(', '), # TODO check for validty
            'barrel-len':row[FS_BRL_LEN],
            'firearm-min-len':row[FS_LLEN],
            'firearm-max-len':row[FS_HLEN],
            'weight':row[FS_WEIGHT],
            'tgr-pull':row[FS_TPULL],
            'stock':row[FS_STOCK].split(', '),     # TODO: once defined, check for validity
            'grip':row[FS_GRIP].split(', '),       # TODO: once defined, check for validity
            'action':row[FS_ACTION].split(', '),   # TODO: once defined, check for validity
            'twist':row[FS_TWIST],                 # TODO: check for correct pattern
            'fsh':row[FS_FSH],                     # TODO: validate
            'rsh':row[FS_RSH],                     # TODO: validate
            'sr':row[FS_SR],                       # TODO: validate
            'cyclic-rof':row[FS_CYCLIC_ROF],
            'effective-rof':row[FS_EFF_ROF],
            'effective-rng-area':row[FS_EFF_RNG_AREA],
            'effective-rng-pt':row[FS_EFF_RNG_PT],
            'moa':row[FS_MOA],
            'feed-system':row[FS_FEED_SYS],
            'feed-loc':row[FS_FEED_LOC],
            'std-capacity':row[FS_STD_CAP],
            'clone':row[FS_CLONE],
            'notes':row[FS_NOTES],
        }

        # type cast as necessary
        try:
            if self._fs[key]['barrel-len']: self._fs[key]['barrel-len'] = float(self._fs[key]['barrel-len'])
            if self._fs[key]['firearm-min-len']: self._fs[key]['firearm-min-len'] = float(self._fs[key]['firearm-min-len'])
            if self._fs[key]['firearm-max-len']: self._fs[key]['firearm-max-len'] = float(self._fs[key]['firearm-max-len'])
            if self._fs[key]['weight']: self._fs[key]['weight'] = float(self._fs[key]['weight'])
            if self._fs[key]['fsh']: self._fs[key]['fsh'] = float(self._fs[key]['fsh'])
            if self._fs[key]['rsh']: self._fs[key]['rsh'] = float(self._fs[key]['rsh'])
            if self._fs[key]['sr']: self._fs[key]['sr'] = float(self._fs[key]['sr'])
            if self._fs[key]['cyclic-rof']: self._fs[key]['cyclic-rof'] = float(self._fs[key]['cyclic-rof'])
            if self._fs[key]['effective-rof']: self._fs[key]['effective-rof'] = float(self._fs[key]['effective-rof'])
            if self._fs[key]['effective-rng-area']: self._fs[key]['effective-rng-area'] = float(self._fs[key]['effective-rng-area'])
            if self._fs[key]['effective-rng-pt']: self._fs[key]['effective-rng-pt'] = float(self._fs[key]['effective-rng-pt'])
            if self._fs[key]['moa']: self._fs[key]['moa'] = float(self._fs[key]['moa'])
        except ValueError as e:
            raise FirearmSpecException("({}) invalid field {}".format(
                key,str(e).split(': ')[1])
            )

        # check stock for allowable values, set to empty list if nothing there
        if self._fs[key]['stock'] == ['']: self._fs[key]['stock'] = []
        else:
            for s in self._fs[key]['stock']:
                if s not in stocks:
                    raise FirearmSpecException(
                        "Invalid stock code ({}) for {}".format(s,key)
                    )

        # check grip for allowable values, set to empty list if nothing there
        if self._fs[key]['grip'] == ['']: self._fs[key]['grip'] = []
        else:
            for g in self._fs[key]['grip']:
                if g not in grips:
                    raise FirearmSpecException(
                        "Invalid grip code ({}) for {}".format(g,key)
                    )

        # check feed system for validy should be 3 colon seperated codes
        if self._fs[key]['feed-system']:
            try:
                # first split the feed system by colon then split each term by comma
                loc,u,l = self._fs[key]['feed-system'].split(':')
                locs = loc.split(',')
                us = u.split(',')
                ls = l.split(',')
                for loc in locs:
                    if loc not in feed_sys_loc:
                        raise FirearmSpecException(
                            "Invalid code ({}) for feed location for {}".format(loc,key)
                        )
                for u in us:
                    if u not in feed_sys_unload:
                        raise FirearmSpecException(
                            "Invalid code ({}) for feed unload for {}".format(u,key)
                        )
                for l in ls:
                    if l not in feed_sys_load:
                        raise FirearmSpecException(
                            "Invalid code ({}) for feed load for {}".format(l,key)
                        )
                self._fs[key]['feed-system'] = {
                    'location':locs,
                    'unload':us,
                    'load':ls,
                }
            except ValueError:
                raise FirearmSpecException(
                    "Incorrectly formated feed-system ({}) for {}".format(
                        self._fs[key]['feed-system'],key
                    )
                )
            except TypeError as e:
                raise FirearmSpecException(e)

        # check feed location for validity
        if self._fs[key]['feed-loc'] and not self._fs[key]['feed-loc'] in feed_loc:
            raise FirearmSpecException("Invalid feed location ({}) for {}".format(
                key,self._fs[key]['feed-loc']
                )
            )

        # convert and check std-cap for validity
        if self._fs[key]['std-capacity']:
            try:
                self._fs[key]['std-capacity'] = [int(s) for s in self._fs[key]['std-capacity'].split(', ')]
            except (ValueError,TypeError) as e:
                raise FirearmSpecException(
                    "Failed to process {} data for 'std-capacity' {}".format(key,e)
                )

        # convert and check trigger-pull for validity
        if self._fs[key]['tgr-pull']:
            try:
                self._fs[key]['tgr-pull'] = [float(s) for s in self._fs[key]['tgr-pull'].split('-')]
            except (ValueError,TypeError) as e:
                raise FirearmSpecException(
                    "Failed to process {} data for 'trigger-pull' {}".format(key,e)
                )

        return True