#!/usr/bin/env python
""" caliber.py
Copyright (C) 2020 Dale V. Patterson (dale.v.patterson@gmail.com)

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

Defines caliber dataset reading/writing and the Caliber class
"""

#__name__ = 'caliber'
__license__ = 'GPLv3'
__version__ = '0.1.5'
__date__ = 'May 2021'
__author__ = 'Dale Patterson'
__maintainer__ = 'Dale Patterson'
__email__ = 'dale.v.patterson@gmail.com'
__status__ = 'Development'

import csv
import regex as re
import numpy as np
import pymsdb.ballistics as bls
import pymsdb.ballistics.bullet as bullet
from pymsdb import PyMSDbException
from pymsdb.data import pth_caliber

# FILE COLUMN CONSTANTS
CAL_NAME   =  0
CAL_FA     =  1
CAL_SZ     =  2
CAL_DIA    =  3
CAL_CLEN   =  4
CAL_BLEN   =  5
CAL_WGT    =  6
CAL_MV     =  7
CAL_ME     =  8
CAL_MMTN   =  9
CAL_CHRG   = 10

# regex
re_blt_sz = re.compile(r"(\d+(?:\.\d+)*)[×x](\d+(?:\.\d+)*)")
re_cal_shotgun = re.compile(r"^([\d\.])+-(gauge|bore)$")
re_cal_diabylen = re.compile(r"^(\d\d*(?:\.\d+)?[x×]\d+)mmR?(?: ([\w–\/ ]+))?$")
re_cal_dia = re.compile(r"^([\d\.]+)(?:mmR?)? ([\w&\- ]+)$")
re_cal_diabyweight = re.compile(r"^([\d\.]+(?:mm)?-\d+)(?: ([\w&\- ]+))?$")
re_cal_misc = re.compile(r"^([\d\.]+\/\d+)(?: ([\w&\-– ]+))?$") # .577/450

"""
 TODO:
  1. re-add momentum can be calculated from mass*velocity
  3. shotgun mv, me is all estimated/notional
  6. get bullet length (needed for Lahti Cd calculations as well as others we 
   might need BUT should this be in the caliber spreadsheet or elsewhere to 
   minimize clutter
"""

class CaliberException(PyMSDbException):
    def __init__(self,msg): super().__init__(msg)

class CaliberDS(object):
    """ The Caliber Dataset """
    def __init__(self,fpath=pth_caliber,update=False):
        """
        initialize CaliberDS
        :param fpath: path of the dataset to read, if None, creates an empty db
        """
        self._ks = []
        self._cs = {}
        self._short_names = []
        self._family = []
        if fpath: self.read(fpath,update)

    def read(self,fpath=pth_caliber,update=False):
        """
        reads the caliber dataset at path
        :param fpath: path of caliber dataset
        :param update: if True calculate momentum from mass/velocity
        """
        self._ks = []
        self._cs = {}
        self._short_names = []
        self._family = []
        fin = None
        try:
            fin = open(fpath)
            csvr = csv.reader(fin,delimiter='\t')
            _ = csvr.__next__() # ignore headers

            # read contents
            for row in csvr:
                key = None
                try:
                    key = row[CAL_NAME]
                    self._cs[key] = Caliber(key,row,'G1',update)
                    self._ks.append(key)
                    self._short_names.append(self._cs[key].short_name)
                    self._family.append(self._cs[key].family)
                except CaliberException as e:
                    raise CaliberException("({}) {}".format(key,e))
        except IOError as e:
            raise CaliberException("Error reading caliber: ds {}".format(e))
        finally:
            if fin: fin.close()

    def write(self,fpath=pth_caliber):
        """ write caliber dataset to path """
        fout = None
        try:
            fout = open(fpath,'w')
            hdr = "NAME\tTYPE\tSIZE (mm)\tBASE DIAMETER (mm)\tCASE LEN (mm)\t" \
                  "BLT LEN (mm)\tWEIGHT (gr)\tMV (m/s)\tME (N⋅m)\t" \
                  "Momentum (kg⋅m/s)\tCHARGE (gr)\n"
            fout.write(hdr)
            for key in self._ks: fout.write(self._cs[key].__str__())
        except IOError as e:
            print('Error writing file',e)
        finally:
            if fout: fout.close()

    def __getitem__(self,param): return self._cs[self.find_name(param)]
    def __setitem__(self,key,value):
        raise CaliberException("Item assignment not supported")
    def __delitem__(self,key):
        raise CaliberException("Item deletion not supported")

    def keys(self): return self._ks

    def find_name(self,cal):
        # check keys first, then short name, then family
        if cal in self._ks: return self._ks[self._ks.index(cal)]
        elif cal in self._short_names: return self._ks[self._short_names.index(cal)]
        elif cal in self._family: return self._ks[self._family.index(cal)]
        return None

    def get_caliber(self,cal): return self._cs[self.find_name(cal)]

class Caliber(object):
    def __init__(self,key,row,mdl='G1',update=False):
        """
        Create a Caliber object with name key, and data in row
        :param key: name of caliber
        :param row: data
        :param update: if true will update calculable data i.e. muzzle energy,
         momentum
        """
        # initialize record (dict of values read in from roww)
        self._rs = {}
        self._name = key
        self._mdl = mdl

        # the below are used for calculating trajectory
        self._process_row_(row,update)

    def __str__(self):
        """
        :return: string row of internal data
        """
        return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
            self._name,
            self._rs['firearm-cd'],
            self._rs['size'],
            self._rs['d'] if self._rs['d'] else '',
            self._rs['case-len'] if self._rs['case-len'] else '',
            self._rs['blt-len'] if self._rs['blt-len'] else '',
            ", ".join([str(m) for m in self._rs['weight']]),
            ", ".join([str(mv) for mv in self._rs['mv']]),
            ", ".join([str(me) for me in self._rs['me']]),
            ", ".join([str(m) for m in self._rs['momentum']]),
            self._rs['charge'],
        )

    @property
    def v_0(self): return self._rs['mv'] # initial/muzzle velocity

    def blt_idx(self,idx):
        """
         returns a Bullet object representing the Caliber specifics at idx
        :param idx: index to get
        :return: Bullet object
        """
        if idx < 0 or idx > len(self._rs['weight']):
            raise CaliberException("{}: {} is out of range".format(self._name,idx))
        return bullet.Bullet(
            self._name,
            self._rs['firearm-cd'],
            self._rs['size'],
            self._rs['d'],
            self._rs['case-len'],
            self._rs['blt-len'],
            self._rs['weight'][idx],
            self._rs['mv'][idx],
            self._rs['charge'],
            self._mdl
        )

    def blt_gr(self,gr):
        """
         returns a Bullet object representing the Caliber specifics by
         grain
        :param gr: grain of bullet to get
        :return: Bullet object
        """
        for i,g in enumerate(self._rs['weight']):
            if g == gr:
                return bullet.Bullet(
                    self._name,
                    self._rs['firearm-cd'],
                    self._rs['size'],
                    self._rs['d'],
                    self._rs['case-len'],
                    self._rs['blt-len'],
                    self._rs['weight'][i],
                    self._rs['mv'][i],
                    self._rs['charge'],
                    self._mdl
                )
        return None


#### PROPERTIES

    @property
    def name(self): return self._name

    @property
    def short_name(self): return self._short_name_()[0]

    @property
    def family(self): return self._short_name_()[1]

    @property
    def mdl(self): return self._mdl

    @property
    def firearm_code(self): return self._rs['firearm-cd']

    @property # TODO: rename this possibly caliber?
    def size(self): return self._rs['size']

    @property
    def diameter(self): return self._rs['d']

    @property
    def case_length(self): return self._rs['case-len']

    @property
    def blt_length(self): return self._rs['blt-len']

    @property
    def cnt(self): return len(self._rs['weight'])

    @property
    def weight(self): return self._rs['weight']

    @property
    def muzzle_velocity(self): return self._rs['mv']

    @property
    def muzzle_energy(self): return self._rs['me']

    @property
    def momentum(self): return self._rs['momentum']

    @property
    def charge(self): return self._rs['charge']

#### PRIVATE HELPERS

    def _process_row_(self,row,update=False):
        """
         process row of data updating if specified and set model based parameters
        :param row: list of texts
        :param update: updates momentum, energy from mass and velocity
        :return:
        """
        #TODO: add update functionality
        # read in values from row
        self._rs = {
            'firearm-cd':row[CAL_FA],
            'size':row[CAL_SZ],
            'd':row[CAL_DIA],
            'case-len':row[CAL_CLEN],
            'blt-len':row[CAL_BLEN],
            'weight':row[CAL_WGT].split(', '),
            'mv':row[CAL_MV].split(', '),
            'me':row[CAL_ME].split(', '),
            'momentum':row[CAL_MMTN].split(', '),
            'charge':row[CAL_CHRG],
        }

        # convert diameter
        if self._rs['d']:
            try:
                self._rs['d'] = float(self._rs['d'])
            except ValueError:
                raise CaliberException(
                    "Invalid 'diameter' ({})".format(self._rs['d'])
                )
        else: self._rs['d'] = None

        # convert case length
        if self._rs['case-len']:
            try:
                self._rs['case-len'] = float(self._rs['case-len'])
            except ValueError:
                raise CaliberException(
                    "Invalid 'case length' ({})".format(self._rs['case-len'])
                )
        else: # use notional for shotgun
            if self._rs['firearm-cd'] == 'S': self._rs['case-len'] = 0.5
            else: self._rs['case-len'] = None

        # convert bullet length
        if self._rs['blt-len']:
            try:
                self._rs['blt-len'] = float(self._rs['blt-len'])
            except ValueError:
                raise CaliberException(
                    "Invalid 'bullet length' ({})".format(self._rs['bullet-len'])
                )
        else: # use notional for shotgun
            if self._rs['firearm-cd'] == 'S': self._rs['blt-len'] = 0.1
            else: self._rs['blt-len'] = None

        # weight is a list of integers
        try:
            # create numpy array out of mass(s) in grains, then convert to kgs
            self._rs['weight'] = np.array(self._rs['weight'],np.int)
        except ValueError:
            if self._rs['weight'] == ['']: self._rs['weight'] = np.array([],np.int)
            else:
                raise CaliberException(
                    "Invalid 'weight(s)' ({})".format(self._rs['weight'])
                )

        # muzzle velocity is a list of floats
        try:
            self._rs['mv'] = np.array(self._rs['mv'],np.double)
        except ValueError:
            if self._rs['mv'] == ['']: self._rs['mv'] = np.array([],np.double)
            else:
                raise CaliberException(
                    "Invalid 'muzzle velocity' ({})".format(self._rs['mv'])
                )

        # muzzle energy is a list of floats
        try:
            self._rs['me'] = np.array(self._rs['me'],np.double)
        except ValueError:
            if self._rs['me'] == ['']: self._rs['me'] = np.array([],np.int)
            else:
                raise CaliberException(
                    "Invalid 'muzzle energy' ({})".format(self._rs['me'])
                )

        # momentum is list of floats
        if update:
            # new calibers added or old ones updated, overwrite any current
            # momentum with calculated momentums
            self._calc_momentum_()
        else:
            try:
                self._rs['momentum'] = np.array(self._rs['momentum'],np.double)
            except ValueError:
                if self._rs['momentum'] == ['']:
                    self._rs['momentum'] = np.array([],np.int)
                else:
                    raise CaliberException(
                        "Invalid 'momentum' ({})".format(self._rs['momentum'])
                    )

        # & charge is a simple float
        try:
            if self._rs['charge']: self._rs['charge'] = float(self._rs['charge'])
        except ValueError:
            raise CaliberException(
                "Invalid 'charge' ({})".format(self._rs['charge'])
            )

    def _short_name_(self):
        """
        :return: tuple t = (short name,family) where
          short name removes mm(R) and family does not include mm or text name
         i.e. given 5.56x45mm NATO
          short name = 5.56x45 NATO
          family = 5.56x45
        """
        # TODO: is this still necessary
        # try diameterxlength
        try:
            sz,txt = re_cal_diabylen.search(self._name).groups()
            return "{}{}".format(sz, ' ' + txt if txt else ''),sz
        except AttributeError:
            pass

        # then diameter-weight
        try:
            sz,_ = re_cal_diabyweight.search(self._name).groups()
            return self._name,sz
        except AttributeError:
            pass

        # then diameter
        try:
            sz,txt = re_cal_dia.search(self._name).groups()
            sn = sz + ' {}'.format(txt) if txt else ''
            fam = sz if sz.startswith('.') else sz + 'mm'
            return sn,fam
        except AttributeError:
            pass

        # finally misc
        try:
            sz,_ = re_cal_misc.search(self._name).groups()
            return self._name,sz
        except AttributeError:
            pass

        # nothing else to try
        return None,None

#### UPDATE CALCULABLE DATA

    # TODO: the below have to be redone after recent changes
    def _calc_me_(self):
        """
        helper function to calculate muzzle energy when caliber is added/updated
        muzzle energy(KE) = 1/2*m*v^2
        where
         m = mass and
         v = velocity
        result is N⋅m
        """
        try:
            self._rs['me'] = np.round(
                self._rs['weight']*bls.GR2KG*np.power(self._rs['mv'],2)*0.5,3
            )
        except IndexError:
            pass

    def _calc_momentum_(self):
        """
         helper function to calculate momentum when caliber is added/updated
         momentum = m*v
         where
          m = mass and
          v = velocity
        result is in kg⋅m/s
        """
        try:
            self._rs['momentum'] = np.round(
                self._rs['weight']*bls.GR2KG*self._rs['mv'],3
            )
        except IndexError:
            pass