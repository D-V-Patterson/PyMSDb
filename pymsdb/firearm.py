#!/usr/bin/env python
""" firearm.py
Copyright (C) 2020 Dale V. Patterson (dale.v.patterson@gmail.com)

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

Defines the Firearm class - a wrapper around a firearm dict for msdb incidents

a Firearm is defined as:
 {
  model: (string) model i.e. SKS
  type: (string) as defined in firearm_type
  manuf: (string) maker of firearm
  fire-act: (string) as defined in firearmspec.fire_acts
  caliber: (string) imperial or metric caliber of bullet
  cap: (integer) number of rounds in initial magazine
  ttl-mags: (list) of dicts of the form {op,num,cap}
  attachments: (list) of strings
  rnds-fired: (dict) of the form {op,num}
  notes: (string)
  std-mag: (list) of numbers
  mags-changed: (dict) of the form {op,num}
  legality: (string) '/' separated codes as defined in firearm_legal_c2l
  acquired; (string) when the firearms was bought/stole/borrowed
}
 with defined_cat as derived from db and taxonomy
 with desc as derived from firearm dict
 with capacities as derived from firearm dict
 with killed, wounded, malfunction, primary found (if present) in notes field
 with ammo dict as found in caliber db
 with specs dict as found in firearmspec db
 with category, type, sub_type as defined in taxonomy for lvl 1, lvl 2 & lvl 3 respectively
"""

#__name__ = 'firearm'
__license__ = 'GPLv3'
__version__ = '0.0.7'
__date__ = 'May 2021'
__author__ = 'Dale Patterson'
__maintainer__ = 'Dale Patterson'
__email__ = 'dale.v.patterson@gmail.com'
__status__ = 'Development'

import datetime
import regex as re
import pymsdb.taxonomy as tax
import pymsdb.firearmspec as firearmspec
from pymsdb import PyMSDbException
from pymsdb.utils import is_int

#### CONSTANTS
# TODO: at a minimum, need:
#   FS (front sight height),
#   RS (rear sight height), and
#   SR (sight radius)
#  specs are not avaiable for each firearm so estimate it by classification into
#   Pistols
#   Revolvers
#   Shotguns
#   Rifles
#   Carbines
#   and possibly tactical (which generally have higher sights
#  although preferred would be to have the three measurements for every firearm
#  NOTE: this is off a little
#
#           /b|
#          /  | FS
#   = SR  /   |
#        /    |
#       /     |
#      /a_____|
#    RS   = RS-FS
#
# Given SR (the distance between FS and RS) is the hyponteneuse and since RS is
# higher than the FS, the opposite side (of b) = RS-FS we can use the formula
#  b = arcsin((RS-FS)/SR) (in degrees)
# to find the angle b that the barrel is 'tilted' IOT to have horizontal LoS
#  from RS to FS to target
BETA_HG = 0.871
BETA_LG = 0.427
# BETA_HG and BETA_LG are the estimated angles of the bore (for handguns & long
# guns respectively) if both the front sight and rear sight are aligned
# horizontally with the target (i.e. LoS)
# TODO: the below are off from bullet.zero_angle
# Glock 19 GEN 4
#  SR = 152mm, FS = 4.191mm, RS = 6.502 b = 0.871
# M4 Carbine (NOTE: M4 has high sights
#  SR = 381mm, FS = 63.5mm, RS = 66.04mm, b = 0.427
#  NOTE: AK has 378mm SR
# M16
#  SR = 500mm

# Weapon Legality
#Legal/Purchase From/Purchase By/Form 4473/NCIC/NFA/State Permit/GOV Failure
# Legal - was firearm legal in at time/location of incident
# NCIC in 1998, Form 4473 in 1968, NFA in 1934
# U/U/U/U/U/N/N/U

# categories of firearm legality
firearm_lgl_cat = [
    'legal','prch-from','prch-by','form-4473','ncic','nfa','state','gov-err'
]

# firearm legaility code to labels
firearm_lgl_c2l = {
    'legal':{'L':'legal','S':'Stolen','B':'Borrowed','I':'Illegal possession',
             'P':'Prohibited conviction','A':'Adjudicated mental','G':'Ghost',
             'T':'Illegal Transaction','C':'Post Convicted','D':'Post Diagnosis',
             'R':'Unregulated Sale','U':'Unknown'},
    'prch-from':{'F':'FFL','T':'Online transfer to FFL','O':'Online FFL Unknown',
                 'L':'Unlicensed','P':'Private','V':'Victim','N':'FFL not required',
                 'U':'Unknown'},
    'prch-by':{'O':'Owner','F':'Friend','P':'Partner','R':'Relative',
               'S':'Straw-buyer','E':'Employer','U':'Unknown'},
    'form-4473':{'A':'Accurate','L':'Lied','N':'Not required','B':'Bypassed',
               'U':'Unknown'},
    'ncic':{'P':'Passed','E':'Time-limit Expired','N':'Not required',
            'G':'Govt failure','B':'Bypassed','S':'Private seller conducted',
            'U':'Unknown'},
    'nfa':{'T':'Tax Stamp','N':'Not required','I':'Illegal not registered',
           'M':'Illegal Modification','B':'Bypassed','U':'Unknown'},
    'state':{'P':'Permitted','B':'Banned in SOU','N':'Not required',
              'R':'Revoked','E':'Expired','X':'Unregistered','D':'Prohibited',
              'C':'Banned component','O':'Other','U':'Unknown'},
    'gov-err':{'Y':'Other Government failure','N':'No','U':'Unknown'}
}

# allowed values for each firearm legality category
firearm_lgl_legal = [k for k in firearm_lgl_c2l['legal']]
firearm_lgl_prch_from = [k for k in firearm_lgl_c2l['prch-from']]
firearm_lgl_prch_by = [k for k in firearm_lgl_c2l['prch-by']]
firearm_lgl_fm_4473 = [k for k in firearm_lgl_c2l['form-4473']]
firearm_lgl_ncic = [k for k in firearm_lgl_c2l['ncic']]
firearm_lgl_nfa = [k for k in firearm_lgl_c2l['nfa']]
firearm_lgl_st_permit = [k for k in firearm_lgl_c2l['state']]
firearm_lgl_gov = [k for k in firearm_lgl_c2l['gov-err']]

# firearm legality regex
re_firearm_lgl = re.compile(
    r"^([LSBIPAGTCDRU])\/" # legal
    r"([FTOLPVNU])\/"      # purchase from
    r"([OFPRSEU])\/"       # purchase by
    r"([ALNBU])\/"         # Form 4473
    r"([PEBGBUNS])\/"      # NCIC
    r"([TNIMBU])\/"        # NFA
    r"([PBNREXDCOU])\/"    # State/Other Permit
    r"([YNU])$"            # gov
)

#### REGEX patterns for validation of record
# caliber regex
re_shotgun_cal = re.compile(r"^\d+-gauge$")
re_euro_cal = re.compile(r"^\d\d*(?:\.\d+)?x\d+$")
re_us_cal = re.compile(r"^\.\d\d\d*(?:-\d\d\d?)?(?: [\w& ]+)?$")
re_dia_manuf = re.compile(r"^[\d\.]+mm[Rr]? \w+$")

# magazine related
re_mag_cap = re.compile(r"^(\d+)\-rnd$")
re_ttl_mags = re.compile(r"^([<>\?])?(\d+)x (?:(\d+|\?))-rnd$")
re_est_num = re.compile(r"([<>\?])?(\d+)")

#### REGEX patterns for notes section
re_pri_firearm = re.compile(r"Primary(?:-(\w+))?\.")
re_firearm_vic = re.compile(r"Killed (\d+)(?:, Wounded (\d+))?\.")
re_firearm_malfunction = re.compile(r"[Ww]eapon malfunctioned(?: (.+?))?\.")
re_firearm_serno = re.compile(r"SN ([\w\-]+)\.")
re_awb_compliant = re.compile(
    r"(State|Federal) AWB compliant (grandfathered|feature [^\.]+)\."
)

"""
# TODO:
 1. how to translate Kum & GO SKS to get specs and any other
 2. save 'name/key' of caliber and firearmspec
 3. add functionality to read firearm spec db w/out going
  thru msdb
"""

# simple firearm exception
class FirearmException(PyMSDbException):
    def __init__(self,msg): super().__init__(msg)

class Firearm(object):
    """ a wrapper around a firearm dict """
    def __init__(
            self,cs,fs,model,cat,manuf,fa,cal,cap,ttlmags,attachs,fired,notes,
            std,changed,legal,acq,acqst
    ):
        """
        initializes firearm with given paremters where a firearm is defined as
        :param cs: CaliberDS Object
        :param fs: Firearm Object
        :param model: (string) model i.e. SKS
        :param cat: (string) as defined in taxonomy
        :param manuf: (string) maker of firearm
        :param fa: (string) fire action as defined in taxonomy
        :param cal: (string) imperial/metric caliber of bullets
        :param cap: (integer) number of rounds in initial magazine
        :param ttlmags: (string) comma separated magazine descriptions
        :param attachs: (string) comma separated list of attachments
        :param fired: (string) rounds fired
        :param notes: (string) notes on specific firearm
        :param std:  (string) standard magazine capacity (may be comma seperated list)
        :param changed: (string) magazines changed
        :param legal: (string) '/' seperated legality codes as defined above
        :param acq: (string) date firearm acquired in MM/DD/YYYY format
        :param acqst: (string) 2-char state abbreviation
        """
        # initialize internal data to default firearms
        self._r = None           # internal dict as read from firearm row fields
        self._serno = None       # as found in notes field
        self._killed = None      # as found in notes field
        self._wounded = None     # as found in notes field
        self._malfunction = None # as found in notes field
        self._primary = None     # as found in notes field
        self._defined_cat = ''   # as derived from db an taxonomy
        self._capacities = []    # as derived from cap and mags
        self._desc = ''          # as derived from manuf and model
        self._blt = None         # as found in caliber db
        self._specs = {}         # as found in firearmspec db
        self._category = None    # as defined through firearmspecs and taxonomy LVL 1
        self._family = None      # as defined through firearmspecs and taxonomy LVL 2
        self._group = None       # as defined through firearmspecs and taxonomy LVL 3
        self.reset()

        # validate the parameters for internal dict, process the internal dict, using
        # external db to derive additional data and classify the firearm
        self._validate_(
            cs,fs,model,cat,manuf,fa,cal,cap,ttlmags,attachs,fired,notes,
            std,changed,legal,acq,acqst
        )
        self._process_()
        self._classify_()

    def reset(self):
        # resets all internal data
        # empty firearm dict
        self._r = {
            'model':'',
            'type':'Firearm',
            'manuf':'',
            'fire-act':'',
            'caliber':'',
            'capacity':'',
            'ttl-mags':[],
            'attachments':[],
            'rnds-fired':{},
            'notes':'',
            'std-mag':[],
            'mags-changed':None,
            'legality':{
                'legal':'U',
                'prch-from':'U',
                'prch-by':'U',
                'form-4473':'U',
                'ncic':'U',
                'nfa':'N',
                'state':'N',
                'gov-err':'U',
            },
            'acquired':None,
            'acq-state':'XX',
        }

        # reset internal variables
        self._category = self._family = self._group = None
        self._serno = None
        self._killed = None
        self._wounded = None
        self._malfunction = {}
        self._defined_cat = ''
        self._primary = None
        self._desc = ''
        self._capacities = []
        self._blt = None
        self._specs = {}

    def __getitem__(self,param): return self._r[param]
    def __setitem__(self,key,value):
        raise FirearmException("Item assignment not supported")
    def __delitem__(self,key):
        raise FirearmException("Item deletion not supported")

    def write(self,fout):
        # writes self to file object fout
        fout.write(
            "\t\t\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                self._r['model'],
                self._r['type'],
                self._r['manuf'],
                self._r['fire-act'],
                self._r['caliber'],
                "{}-rnd".format(self._r['capacity']) if self._r['capacity'] else '',
                self._format_ttl_mags_(),
                ", ".join(self._r['attachments']),
                self._format_rnds_fired_(),
                self._r['notes'],
                self._format_std_mag_(),
                self._format_mag_change_(),
                "{}/{}/{}/{}/{}/{}/{}/{}".format(
                    self._r['legality']['legal'],
                    self._r['legality']['prch-from'],
                    self._r['legality']['prch-by'],
                    self._r['legality']['form-4473'],
                    self._r['legality']['ncic'],
                    self._r['legality']['nfa'],
                    self._r['legality']['state'],
                    self._r['legality']['gov-err'],
                ),
                self._r['acquired'].strftime('%m/%d/%Y') if self._r['acquired'] else '',
                self._r['acq-state'] if self._r['acq-state'] != 'XX' else '',
            )
        )

#### properties

    @property
    def model(self): return self._r['model']

    @property
    def type(self): return self._r['type']

    @property
    def manufacturer(self): return self._r['manuf']

    @property
    def caliber(self): return self._r['caliber']

    @property
    def diameter(self): return self._r['ammo']['diameter']

    @property
    def case_length(self): return self._r['ammo']['case-len']

    @property
    def muzzle_velocity(self):
        try:
            return self._r['ammo']['muzzle-velocity']
        except KeyError:
            return None

    @property
    def muzzle_energy(self):
        try:
            return self._r['ammo']['muzzle-energy']
        except KeyError:
            return None

    @property
    def capacity(self): return self._r['capacity']

    @property  # note this only a total count, for
    def ttl_mags_cnt(self):
        tm = [mag['num'] for mag in self._r['ttl-mags']]
        if tm: return sum(tm)
        if self._r['capacity']: return 1
        return None

    @property
    def capacities(self): return self._capacities

    @property
    def attachments(self): return self._r['attachments']

    @property
    def rnds_fired(self):
        return self._r['rnds-fired']['num'] if self._r['rnds-fired'] else None

    @property
    def notes(self): return self._r['notes']

    @property
    def legality(self): return firearm_lgl_c2l[self._r['legality']['legal']]

    @property
    def purchased_from(self): return firearm_lgl_c2l[self._r['legality']['prch-from']]

    @property
    def purchased_by(self): return firearm_lgl_c2l[self._r['legality']['prch-by']]

    @property
    def form_4473(self): return firearm_lgl_c2l[self._r['legality']['form-4773']]

    @property
    def ncic(self): return firearm_lgl_c2l[self._r['legality']['ncic']]

    @property
    def nfa(self): return firearm_lgl_c2l[self._r['legality']['nfa']]

    @property
    def state(self): return firearm_lgl_c2l[self._r['legality']['state']]

    @property
    def govt_error(self): return firearm_lgl_c2l[self._r['legality']['gov-err']]

    @property
    def date_acquired(self): return self._r['acquired']

    @property
    def state_acquired(self): return self._r['acq-state']

    @property
    def std_mag_cap(self): return self._r['std-mag']

    @property
    def mags_changed(self):
        return self._r['mags-changed']['num'] if self._r['mags-changed'] else None

    @property
    def malfunction(self):
        try:
            return self._malfunction['type']
        except (KeyError,TypeError):
            return None

    @property
    def killed(self): return self._killed

    @property
    def wounded(self): return self._wounded

    @property
    def fatalities(self): return self._killed,self._wounded

    @property
    def description(self): return self._desc

    @property # this is category as defined by law or manufacturer
    def defined_category(self): return self._defined_cat

    @property
    def serno(self): return self._serno

    # category, type and sub-type are defined by specifications

    @property
    def category(self): return self._category

    @property
    def family(self): return self._family

    @property
    def group(self): return self._group

    # TODO: add properties for firearm specs

    @property
    def fsh(self): return self._specs['fsh'] if self._specs['fsh'] else None

    @property
    def rsh(self): return self._specs['rsh'] if self._specs['rsh'] else None

    @property
    def sr(self): return self._specs['sr'] if self._specs['sr'] else None

#### ATTRIBUTES (DERIVED PROPERTIES)

    def is_primary(self): return self._primary

    def full_name(self):
        # returns the full name = manuf + model if both exist, none otherwise
        if self._r['manuf'] and self._r['model']:
            return "{} {}".format(self._r['manuf'],self._r['model'])
        else: return None

    def max_capacity(self):
        try:
            return max(self._capacities)
        except ValueError: # empty list returned
            return None

    def magazines(self):
        # a list of tuples t = (num,cap) (disregards any operator)
        if self._r['capacity'] and not self._r['ttl-mags']:
            return [(1,self._r['capacity'])]
        else: return [(m['num'],m['cap']) for m in self._r['ttl-mags']]

    def is_tac_firearm(self):
        """ is defined an 'assualt weapon' """
        if tax.is_tac_firearm(self._r['type']): return 'True'
        elif self._r['type'] in ['Rifle','Pistol'] and self._r['caliber'] in tax.tac_firearm_cal:
            if self._r['fire-act'] in tax.tac_firearm_fa: return 'Unknown'
        return 'False'

    def used_hcm(self):
        """ return 'True', 'False', 'Unknown' regard a hcm used in firearm """
        try:
            return 'True' if self.max_capacity() > 10 else 'False'
        except (TypeError,ValueError):
            if self._r['type'] in ['Shotgun','Revolver','Derringer']: return 'False'
        return 'Unknown'

    def used_extended(self):
        """
         return 'Not-listed if std is unknown, 'unknown' if mags not known,
         false or true otherwise
        """
        try:
            std = max(self._r['std'])
            mcap = self.max_capacity()
            if mcap is None: return 'Unknown'
            elif mcap > std: return 'True'
            else: return 'False'
        except ValueError:
            return 'Not-listed'

#### private helper functions

    def _validate_(self,co,fo,model,cat,manuf,fa,cal,cap,ttlmags,attachs,fired,notes,std,changed,legal,acq,acqst):
        """ validates parameters (expects internal data to be set) """
        # accept all params that don't require validation
        self._r['notes'] = notes

        # model has no checks #TODO: how to break out model and gen etc
        self._r['model'] = model

        # firearm type must be a defined/allowed type
        if cat in tax.firearm_types: self._r['type'] = cat
        else: raise FirearmException('Weapon Type ({}) not recognized'.format(cat))

        # make sure manufacturer is not also listed in model
        if manuf:
            if model and manuf in model:
                raise FirearmException('Remove manufacturer ({}) from model'.format(manuf))
            else: self._r['manuf'] = manuf

        # fire-action has defined allowed value
        if fa and not fa in firearmspec.fire_acts:
            raise FirearmException("Invalid 'fire action' ({})".format(fa))
        else: self._r['fire-act'] = fa

        # caliber is checked against predefined expressions
        if cal: # TODO: drop the below and use from caliber.py
            match = False
            if re_shotgun_cal.match(cal): match = True
            elif re_euro_cal.match(cal): match = True
            elif re_us_cal.match(cal): match = True
            elif re_dia_manuf.match(cal): match = True
            if not match: raise FirearmException("Invalid 'caliber' ({})".format(cal))
            self._r['caliber'] = cal

        # capacity will be of the form <NUM>-rnd (even for shotguns)
        if cap:
            try:
                self._r['capacity'] = int(re_mag_cap.search(cap).groups()[0])
            except (AttributeError,ValueError,IndexError):
                raise FirearmException("Invalid 'capaciity' ({})".format(cap))

        # std-mag may be list or single term (terms will have the form <NUM>-rnd)
        if std:
            for mag in std.split(', '):
                try:
                    self._r['std-mag'].append(
                        int(re_mag_cap.search(mag).groups()[0])
                    )
                except (AttributeError,ValueError,IndexError):
                    raise FirearmException("Invalid 'std-mag' ({})".format(cap))

        # magazines is difficult due to sometimes unspecific reporting - will be
        # a list or single term of the form defined above for std-mag
        if ttlmags:
            for mag in ttlmags.split(', '):
                try:
                    m = re_ttl_mags.search(mag).groups()
                    self._r['ttl-mags'].append(
                        {
                            'op':m[0] if m[0] else '=',
                            'num':int(m[1]) if is_int(m[1]) else '?',
                            'cap':int(m[2]) if is_int(m[2]) else '?',
                        }
                    )
                except AttributeError:
                    raise FirearmException(
                        "Invalid 'magazine(s)'' ({})".format(ttlmags)
                    )

        # no check is performed on attachments, just split out individual terms
        self._r['attachments'] = [at for at in attachs.split(', ') if at]

        # rounds fired fits defined criteria
        if fired:
            try:
                m = re_est_num.search(fired).groups()
                self._r['rnds-fired'] = {'op':m[0] if m[0] else '=','num':int(m[1])}
            except (AttributeError,ValueError):
                raise FirearmException("Invalid 'rnds-fired' ({})".format(fired))

        # mags changed
        if changed:
            try:
                op,num = re_est_num.search(changed).groups()
                self._r['mags-changed'] = {'op':op if op else '=','num':int(num)}
            except AttributeError:
                raise FirearmException("Invalid 'mc' ({})".format(changed))

        # legality must be of proper form
        try:
            l,pf,pb,f4,ncic,nfa,st,gov = re_firearm_lgl.match(legal).groups()
            self._r['legality'] = {
                'legal':l,
                'prch-from':pf,
                'prch-by':pb,
                'form-4473':f4,
                'ncic':ncic,
                'nfa':nfa,
                'state':st,
                'gov-err':gov,
            }
        except AttributeError:
            raise FirearmException("Invalid 'weapon legality' ((})".format(legal))

        # acquired date must be mm/dd/yyyy
        # TODO: must be greater than incident date
        if acq:
            try:
                m,d,y = acq.split('/')
                self._r['acquired'] = datetime.date(int(y),int(m),int(d))
            except (ValueError,SyntaxError):
                raise FirearmException(
                    "Invalid 'acquired-date' ({})".format(self._r['acquired'])
                )

        # state must be a valid two-char abbreviation
        # TODO: how to confirm state abbr exists
        if acqst and len(acqst) == 2: self._r['acq-state'] = acqst

        # link to caliber and firearmspec dbs

        # caliber
        try:
            self._blt = co.get_caliber(self._r['caliber'])
        except KeyError:
            if not self._r['caliber']: pass
            else:
                raise FirearmException("'Caliber' ({}) not found in DB".format(
                    self._r['caliber'])
                )

        # firearmspec
        try:
            self._specs = fo.get_rec(fo.get_key(self._r['manuf'],self._r['model']))
        except KeyError:
            # if there is a manuf and model raise an error
            if self._r['manuf'] and self._r['model']:
                raise FirearmException("'Firearm ({} {}) not found in DB".format(
                    self._r['manuf'],self._r['model'])
                )
            self._specs = None

    def _process_(self):
        # parse notes sections for addtional firearm details
        try:
            mf = re_firearm_malfunction.search(self._r['notes']).groups()[0]
            self._malfunction = {'type':mf if mf else 'unk'}
        except AttributeError:
            self._malfunction = {}

        # determine firearm desc
        if self._r['manuf'] and self._r['model']:
            self._desc = "{} {}".format(self._r['manuf'], self._r['model'])
        elif self._r['type'] and self._r['model'] and not self._r['manuf']:
            self._desc = "{} {}".format(self._r['model'],self._r['type'])
        elif not self._r['manuf'] and not self._r['model'] and self._r['type']:
            self._desc = "Unknown {}".format(self._r['type'])
        elif self._r['manuf'] and not self._r['model']:
            self._desc = "Unknown {} {}".format(self._r['manuf'],self._r['type'])

        # determine capacities
        if self._r['capacity'] and not self._r['ttl-mags']:
            self._capacities = [self._r['capacity']]
        else:
            self._capacities = [
                m['cap'] for m in self._r['ttl-mags'] if is_int(m['cap'])
            ]

        # determine firearm category
        # TODO: have to do something with this
        wt = self._r['type']
        if wt == 'Firearm': self._defined_cat = 'Unknown'
        elif tax.firearm_lvl_4(wt) or tax.firearm_lvl_5(wt):
            if tax.is_tac_rifle(wt): self._defined_cat = 'Assault Weapon'
            elif tax.is_tac_pistol(wt): self._defined_cat = 'Assault Weapon'
            elif wt == 'Carbine': self._defined_cat = 'Carbine'
        elif tax.firearm_lvl_3(wt): self._defined_cat = wt
        elif tax.firearm_lvl_2(wt): self._defined_cat = wt

        # determine capacities

        # check notes for killed, wounded
        try:
            k,w = re_firearm_vic.search(self._r['notes']).groups()
            self._killed = int(k)
            self._wounded = int(w) if w else 0
        except AttributeError:
            pass

        # check notes for serial number
        try:
            self._serno = re_firearm_serno.search(self._r['notes']).groups()[0]
        except AttributeError:
            pass

        # check notes for primary
        try:
            if re_pri_firearm.search(self._r['notes']).groups()[0]: self._primary = True
        except (AttributeError,IndexError):
            self._primary = False

    def _classify_(self):
        # check first if we retreived firearm specs for this gun
        if self._specs is None:
            # TODO: should we fall back on type as defined in internal dict
            self._category = self._family = self._group = None
            return

        # grab specs we need
        try:
            brl = self._specs['barrel-len']
            feed = self._specs['feed-system']['location']
            cat = self._specs['category']
            oal = self._specs['firearm-max-len']
        except KeyError as e:
            print(type(self._specs))
            raise FirearmException("'Firearm ({} {}) missing '{}' in DB".format(
                self._r['manuf'], self._r['model'],e.args[0])
            )

        # start with level 1 handgun/Long gun
        if brl < tax.HG_BRL: self._category = 'Handgun'
        else: self._category = 'Long gun'

        # break out into level two
        if self._category == 'Handgun':
            # revolvers have cylinder feed, derringers have barrel, all others are pistols
            # only pistols are broke down into groups
            if 'cyclinder' in feed: self._family = 'Revolver'
            elif 'barrel' in feed: self._family = 'Derringer'
            else:
                self._family = 'Pistol'
                if brl < tax.PST_PKT_BRL: self._group = 'Pocket'
                elif brl < tax.PST_SUB_BRL: self._group = 'Subcompact'
                elif brl < tax.PST_CMP_BRL: self._group = 'Compact'
                elif brl < tax.PST_FLL_BRL: self._group = 'Full'
                else: self._group = "LBP" # TODO: find a better name for this
        else:
            # for now, have to rely on 'typed' category to identify shotguns
            if cat == 'Shotgun':
                self._family = 'Shotgun'
                try:
                    if brl < tax.SBS_BRL or oal < tax.SBS_OAL: self._group = 'SBS'
                except TypeError: # have an undefined length
                    self._group = 'Undefined'
            else:
                self._family = 'Rifle'
                try:
                    if brl < tax.SBR_BRL or oal < tax.SBS_OAL: self._group = 'SBR'
                    elif brl < tax.RFL_CRB_BRL: self._group = 'Carbine'
                except TypeError: # have an undefined length
                    self._group = 'Undefined'

#### writing helpers

    def _format_ttl_mags_(self):
        if not self._r['ttl-mags']: return ''
        tms = []
        for m in self._r['ttl-mags']:
            tms.append(
                "{}{}x {}-rnd".format(m['op'] if m['op'] != '=' else '',m['num'],m['cap'])
            )
        return ', '.join(tms)

    def _format_rnds_fired_(self):
        if not self._r['rnds-fired']: return ''
        op = '' if self._r['rnds-fired']['op'] == '=' else self._r['rnds-fired']['op']
        num = self._r['rnds-fired']['num']
        return "{}{}".format(op,num)

    def _format_std_mag_(self):
        if not self._r['std-mag']: return ''
        return ', '.join(['{}-rnd'.format(m) for m in self._r['std-mag']])

    def _format_mag_change_(self):
        if not self._r['mags-changed']: return ''
        op = '' if self._r['mags-changed']['op'] == '=' else self._r['mags-changed']['op']
        num = self._r['mags-changed']['num']
        return "{}{}".format(op,num)

