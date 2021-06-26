#!/usr/bin/env python
""" msdb.py
Copyright (C) 2020 Dale V. Patterson (dale.v.patterson@gmail.com)

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

Defines the MSDB class - a wrapper around a the msdb.csv contents
"""

#__name__ = 'msdb'
__license__ = 'GPLv3'
__version__ = '0.2.3'
__date__ = 'April 2021'
__author__ = 'Dale Patterson'
__maintainer__ = 'Dale Patterson'
__email__ = 'dale.v.patterson@gmail.com'
__status__ = 'Development'

import csv
import datetime
import regex as re
import numpy as np
from calendar import month_name as i2m
from operator import itemgetter
import pymsdb.states as states
import pymsdb.taxonomy as taxonomy
import pymsdb.firearm as firearm
import pymsdb.caliber as caliber
import pymsdb.firearmspec as firearmspec
from pymsdb import PyMSDbException
from pymsdb.utils import is_int
from pymsdb.data import pth_msdb, pth_clog

"""
MSDB is a dict of dicts 
 {case_id,incident}
where
 case_id is a unique label of the incident,
 incident is a dict of the form:
 {
  shooter: (list of dicts) {last,first,middle}
  shooter-data = (list of dicts) {age,gender,race,religion,education,military}
   where:
    age = age in years
    gender = one of {M=male, F=female}
    race = one of {W=white, B=black, H=hispanic, A=asian,
                   E=Middle Eastern, N=Native, O=other,U=Uunknown}
    religion = {
      C=christion,M=muslim,J=jewish,A=atheist,B=buddhist,N=none,O=other,U=Unknown
    }
    ed-lvl = {
      0=Unknown,1=less than highschool,2=high school,3=some college/CC/Cert,
      4=BS/BA,5=Higher
    }
    military = {
      A=Army, N=Navy, F=Air Force, M=Marines, C=Coast Guard, G=National Guard, 
      X=None,U=Unknown
    }
    leo-knowledge = {U=Unknown, N=No, Y=Yes}
    mental-health = {
      U=Unknown,N=None,D=Diagnosed,S=Signs/Symptoms,
      P=Post-mortem/Post Incident  Diagnosis
    }
    substance-use = { (history
     U=Unknown,N=None,A=Alcohol,D=Drugs,P=Psychiatric Meds,
     C=Combination of drugs/alcohol/prescribed,
    }
  location: (dict) {city,state}
  date: (datetime) the date of the incident
  victim: (string) of the form K/W/T where K = number killed, W = number wounded,
   T = total victims
  leo: (string) same as above
  notes: (string) freeform notes of incident
  report: (integer) one of {
      0:no report,
      1:partial report, (limited)
      2;full report}  (includes firearm, caliber, magazine # & capacity, rounds fired
      from police or govt
  toa: (float) time in minutes of incident upto point where shooter was occupied
   by police, dead or left the scene
  cke: (list of dicts) {cke-code,cke-party} for each shooter
   where cke_code is oneof 
   {
     C:capture,K:kill,S:surrender,E:escape (escapes will have the form E-Code
   }
   and cke_party is one of None,Victim,Police,Bystander,Self,Guard
  scene: (string) comma seperated list description of scene, each item starts with
   inside or outside
  venue: (string) basic venue locaiton i.e. school, factory etc
  type: (character) as defined in shooting_types
  trigger: (character) as defined in tgr_codes
  motive: a list of motive codes as in defined in motive codes
  target: (string) two character specifying target people and target location
  ttl-fa: (integer) count of all firearms on primary scene/used regardless if
   used or not. For example Albert Wong used a shotgun to kill himself but did
   not use it to kill victims 
  firearms: a list of Firearms
 }
"""

# constants
# for incident
INC_NUM =  0
INC_LBL =  1
INC_SHT =  2
INC_SDT =  3
INC_LOC =  4
INC_DTG =  5
INC_VIC =  6
INC_LEO =  7
INC_NTE =  8
INC_RPT =  9
INC_TOA = 10
INC_CKE = 11
INC_SCN = 12
INC_VEN = 13
INC_TYP = 14
INC_TGR = 15
INC_MTV = 16
INC_TGT = 17
INC_ALL = 18

# for firearms
FA_KEY =  2
FA_MDL =  3
FA_TYP =  4
FA_MAN =  5
FA_FAC =  6
FA_CAL =  7
FA_CAP =  8
FA_MAG =  9
FA_ATH = 10
FA_RND = 11
FA_NTE = 12
FA_STD = 13
FA_MCH = 14
FA_LEG = 15
FA_ACQ = 16
FA_AST = 17

# FAWB START/STOP
FAWB_START = datetime.date(1994,9,13)
FAWB_STOP  = datetime.date(2004,9,13)

# DATA START/STOP
DATA_START = datetime.date(1982,1,1)
DATA_STOP = datetime.date.today()

# DEFAULT HIST START/STOP
DEF_START = datetime.date(1900,1,1)
DEF_STOP  = datetime.date.today()+ datetime.timedelta(days=1)

# CRITERIA 'FLAGS'
# strict scene = >= 4 deaths, with one scene >= 4 deaths
# loose  scene = >= 4 deaths, or 3 deaths and >= 5 total victims in at
#  least one scene
# strict incident >= 4 deaths, total victim count is used regardless of
#  scene counts
# loose incident >= 4 deaths, or 3 deaths and >= 5 total victims across
#  the incident
CRIT_NONE         = 0
CRIT_STRICT_SCENE = 1
CRIT_LOOSE_SCENE  = 2
CRIT_STRICT_INC   = 3
CRIT_LOOSE_INC    = 4
CRIT_KLAREVES     = 5

# TODO: order these in order of appearance in the db file
# shooting type is two character set of codes, the first defines the positon relative
# to the kill zone, the second defines the action of the shooter
shooting_type_pos = ['P','I','A']
shooting_type_pos_c2lshooting_types_c2l = {
    'P':'Fixed',  # a fixed, secured, (generally elevated) position overlooking the kill zone
    'I':'Inside', # shooter enters the kill zone
    'A':'Ambush', # shooter sets up hasty position outside the kill zone
}
shooting_type_act = ['R','S']
shooting_type_act_c2l = {
    'R':'Roaming',    # shooter is roaming
    'S':'Stationary', # shooter stays in one spot
}

# allowed initiated codes
# we have come up with the following codes.
# R: (Reactionary) Shooter responded to a trigger within a short time. May have
#  had firearm on person or in vehicle etc.
# P: (Planned) Shooter had a plan, site, targets, time etc
# I: (Instigated) Shooter had a basic intent/idea and was waiting for a trigger
#  to use as a catalyst
# O: (Opportunity) no plan or intention but an opportunity presented itself
# U: (Unknown) cannot determine
tgr_codes = ['R','P','I','O','U']
tgr_c2l = {
    'R':'Reactionary','P':'Planned','I':'Instigated','O':'Opportunity','U':'Unknown'
}

# allowed motive codes
# A = Authority i.e Law Enforcment, Officials
# C = Cultural
# D = Domestic
# E = Education
# F = Finacial
# G = Gender
# I = Terrorism/Extermism Ideology
# L = Legal
# M = Mentally Impaired, Drugs
# N = Notability/Fame
# O = Orientation
# P = Personal
# R = Racial
# S = Religious
# T = Terrorism/Extremism Affiliation
# U = Unknown
# V = Violent Ideology/Tendency
# W = Employment/Work
# Y = Mentally Impaired, Psychosis
mtv_codes = [
    'A','C','D','E','F','G','I','M','N','O','P','R','S','T','L','U','V','W','Y'
]
mtv_c2l = {
    'A':'Authority','C':'Cultural','D':'Domestic','E':'Education','F':'Finacial',
    'G':'Gender','I':'Extemist Ideology','L':'Legal','M':'Drugs','N':'Notability',
    'O':'Orientation','P':'Personal','R':'Race','S':'Religion',
    'T':'Extremist Affiliation','U':'Unknown','V':'Violence','W':'Work','Y':'Psychosis',
}

# allowed target codes (two-character code <people><location>)
# People:
#  G = Group: a group for example coworkers, classmates etc associated with
#    shooter or a demographic associated with shooter's ideology such as race
#  S = Specific: specific targets (by-name list)
#  R = Random: generally in relation to violent ideology, shooters intent is to
#   kill without regard to who to dies
#  U = Unknown: cannot determine
# Location:
#  C = Chosen: location has no association/ties w/ shooter, selected for
#   potential. Generally in relation to violent ideology
#  S = Specific: location is associated with shooter such as a school or bus.
#   or a demographic i.e. Mosque
#  R = Random: generally in relation to an opportunistic or reactionary event
#  U = Unknown: cannot determine
tgt_c2l_tgt = {'G':'Group','S':'Specific','R':'Random','U':'Unknown'}
tgt_c2l_loc = {'C':'Chosen','S':'Specific','R':'Random','U':'Unknown'}

# shooter demographics is a string with a default/empty form of
# age/gender/race/religion/education/military/criminal rec/mental health/subst. use
# 0/M/U/U/0/U/U/U/U

# categories of demographics
sht_demo_cat = [
    'age','gender','race','religion','ed-lvl','military',
    'leo-knowledge','mental-health','substance-use'
]

# demographic code to labels
sht_demo_c2l = {
    'gender':{'M':'Male','F':'Female'},
    'race':{'W':'White','B':'Black','H':'Hispanic','A':'Asian',
            'E':'Middle Eastern','N':'Native','O':'Other','U':'Unknown'},
    'religion':{'C':'Christian','M':'Muslim','J':'Jewish','B':'Buddhist',
                'A':'Atheist','N':'None','O':'Other','U':'Unknown'},
    'ed-lvl':{'0':'Unknown','1':'Some HS','2':'HS/GED',
              '3':'Some College','4':'Bachelors','5':'Higher'},
    'military':{'A':'Army','N':'Navy','F':'Air Force','M':'Marines',
                'C':'Coast Guard','G':'National Guard','X':'None','U':'Unknown'},
    'leo-knowledge':{'U':'Unknown','N':'None','Y':'Police Knowledge'},
    'mental-health':{'N':'None','U':'Unknown','D':'Pre-Diagnosed','S':'Signs/Symptons',
                     'P':'Post-Diagnosed'},
    'substance-use':{'U':'Unknown','N':'None','A':'Alcohol','D':'Drugs',
                     'P':"Pyschiatric",'C':'Combination'}
}

# allowed values for each category
# sht_demo_age = 0 to 100
sht_demo_gender = ['M','F']
sht_demo_race = [k for k in sht_demo_c2l['race']]
sht_demo_rel = [k for k in sht_demo_c2l['religion']]
sht_demo_ed_lvl = [k for k in sht_demo_c2l['ed-lvl']]
sht_demo_mil = [k for k in sht_demo_c2l['military']]
sht_demo_cri_rec = [k for k in sht_demo_c2l['leo-knowledge']]
sht_demo_men_hlt = [k for k in sht_demo_c2l['mental-health']]
sht_demo_subs_use = [k for k in sht_demo_c2l['substance-use']]

# shooter data regex =
re_shooter_data = re.compile(
    r"^(\d+)\/"       # age
    r"([MF])\/"       # gender
    r"([WBHAENOU])\/" # race
    r"([CMJBANOU])\/" # religion
    r"([0-5])\/"      # education level
    r"([ANFMCGXU])\/" # military
    r"([UNY])\/"      # criminal record
    r"([UNDSP])\/"    # mental health
    r"([UNADPC])$"    # substance use
)

# regex for notes
re_shooter_suicide = re.compile(r"[Ss]hooter killed self prior to police arrival\.")
re_subdued_attempt = re.compile(
    r"(?:(\d+)x )?([Vv]ictim|[Bb]ystander) intervention (attempted|succeeded) "
    r"(during|after) incident\."
)
re_ttl_rnds_fired = re.compile(r"[Ee]xpended (\d+)x rounds(?: of (.+?))?\.")
re_ttl_mags_known = re.compile(r"(\d+)x (?:(\d+|\?))-rnd(?: (.+?))? magazines\.")
re_mag_change = re.compile(r"(\d+)x magazine changes\.")
re_ppl_present = re.compile(r"(?:[Ee]stimated )?(\d+) people on site\.")
re_mult_scenes = re.compile(r"Multiple scenes \d+\/\d+(?:, \d+\/\d+)*\.")
re_mult_scenes_tally = re.compile(r"(\d+)\/(\d+)")
re_acq_state = re.compile(r"(ACQ-State [^\.]+\.)")
re_acq_state_term = re.compile(r"([A-Z][A-Z]:(?:\d+|\?))")
re_acq_date = re.compile(r"ACQ-Date ([^\.]+)\.")
re_acq_date_term = re.compile(r"(\d\d\/\d\d\/\d\d\d\d)")

# fct filter
# three types of terms 1) dates, 2) numeric, 3) strings
filter_date = ['date']
filter_numeric = [
    'year','month','dow','victims','killed','wounded','toa','ttl_firearm','ppl_present',
    'min_capacity','max_capacity','capacity','ttl_mags','rnds_fired','mags_changed',
    'min_acq_delta','max_acq_delta',
]
filter_string = [
    'state','venue','trigger','motive','target','fa_type','fa_manuf','fa_fireact',
    'caliber','attachment','legal','acq_state'
]
# regular expressions
re_filter_term = re.compile(re.compile(r"([^|]+)(?:\s|\s)?")) # find '|' separated terms
# (1) date
re_filter_date_range = re.compile(
    r"(!)?([\[\(])(\d\d)\/(\d\d)\/(\d\d\d\d)-(\d\d)\/(\d\d)\/(\d\d\d\d)([\]\)])"
)
re_filter_date_term = re.compile(r"^(>|<|>=|<=|!=)?(\d\d)\/(\d\d)\/(\d\d\d\d)$")
# (2) numeric can be a range or single term (preceded by comparison operators)
re_filter_num_range = re.compile(r"(!)?([\[\(])(\d+)-(\d+)([\]\)])")
re_filter_num_term = re.compile(r"^(>|<|>=|<=|!=)?(\d+)$")
# (3) strings
re_filter_str_range = re.compile(r"in \[([^\)]+)\]^$")
re_filter_str_term = re.compile(r"^(!=)?(.+)$")

# date conversions
i2dow =  [None,'Monday','Tuesday','Wednesday','Thursday','Friday','Saturday','Sunday']

"""
 # TODO: 
  7. Add something about 'HCM' based on a by type/subtype criteria
  9. How to look at the deaths low to high associated with AT firearm vs non-AT 
   firearm i.e. Waffle House Shooting had an AR-15 but only 6 total victims (4 
   dead, 2 wounded) is on the low-end vs Virginia Tech shooting which used 
   pistols and was on the high end
  19. Figure something out for primary notation for two-shooter incidents
  20. Add support to process acquistion dates in incident notes. (see Bruce Pardo)
  25. Add support in MAGAZINE(S) column to specify unknown number of mags (see Camden)
   i.e. ?x 8-rnd
  26. need to look at rounds fired, primary_firearm for validity 
  27. For primary firearm
    'USPS Shooting', both firerarms were 1911A
  28. Finish caliber reading/processing
  29. Per capita firearms in the us outnumber people
  30. fatality rate = killed / total shot
  32. need regex for acquisition dates that cannot be tied to specific weapons
   - regex done but need functionality to capture it. to use regex do
    ad = re_acq_date.search(test).groups()[0]
    dates = re_acq_date_term.findall(ad)
  33. add a 'flatten' functionality that will write the msdb to tsv/csv without
   subrows for weapons
"""

# data manip functions for hists/pivot

def key_sort(s): return float('inf') if type(s) == type('') else s

def hist_cnt(hs,rnull=False,sort=False,sortby=1):
    """
    returns the hist hs as a dict of key->cnt or list of tuples (key,cnt) i.e
    reduces the hist's list of incidents to a count
    if rnull is True, values of 0 are removed
    if sort is True, the hist is returned as a list sorted by the count
    if sortby == 0; will sort by the key otherwise if 1 will sort by the cnt
    """
    h = {k:len(hs[k]) for k in hs if hs[k]} if rnull else {k:len(hs[k]) for k in hs}
    if sort: return sorted([(k,h[k]) for k in h],key=lambda x: key_sort(x[sortby]))
    else: return h

def hist_sort(hs):
    """ sorts the hist in ascending order based on the number of incidents """
    return sorted([(k,hs[k]) for k in hs],key=lambda x: len(x[1]))

def hist_rank(hs):
    """ sorts the hist in ascending order based on the key """
    return sorted([(k,hs[k]) for k in hs],key=lambda x: key_sort(x[0]))

def pvt_cnt(pvt,sort=False):
    """
    returns the pivot as a dict of dicts of the form {fct1:{fct2:cnt}} reducing
    the orginal pivot's list of incidents to a count
    :param pvt: the orginal pivot
    :param sort: if True, return pivot as a sorted list of tuples
    :return: the reduced pivot
    """
    cpvt = {}
    if sort:
        for p in pvt:
            cpvt[p] = sorted([(s,len(pvt[p][s])) for s in pvt[p]],key=itemgetter(1))
    else:
        for p in pvt: cpvt[p] = {s:len(pvt[p][s]) for s in pvt[p]}
    return cpvt

def pvt_print(pvt):
    """ pretty prints a pivot """
    for p in pvt:
        print(p)
        for s in pvt[p]:
            print("\t({}) {}: {}".format(len(pvt[p][s]),s,pvt[p][s]))

#

#### HELPER CLASSES TODO: better in utils? same as fcts above?

class BucketException(PyMSDbException):
    def __init__(self,msg): super().__init__(msg)
class EpochException(PyMSDbException):
    def __init__(self,msg): super().__init__(msg)
class MSDBException(PyMSDbException):
    def __init__(self,msg): super().__init__(msg)

class Bucket(object):
    """ A class for the msdb fct buckets() """
    def __init__(self,left,right,data=None):
        self._data = data if data else []
        try:
            _ = left > 0 and right > 0
        except TypeError:
            raise BucketException(
                "Invalid Arg 'left/right' ({}/{}".format(left,right)
            )
        self._left = left
        self._right = right

    @property
    def label(self): return "{} to {}".format(self._left,self._right)

    @property
    def count(self): return len(self._data)

    @property
    def left(self): return self._left

    @property
    def right(self): return self._right

    def get_data(self): return self._data
    def set_data(self,d):
        if not type(d) == type([]):
            raise BucketException("Invalid Arg 'data' ({})".format(d))
        self._data = d
    def del_data(self): self._data = []
    data = property(get_data,set_data,del_data)

class Epoch(object):
    """ A class for the msdb fct epochs (similar to above) """
    def __init__(self,start,stop,data=None):
        self._data = data
        if not type(start) == datetime.date:
            raise EpochException("Invalid Arg 'start'".format(start))
        if not type(stop) == datetime.date:
            raise EpochException("Invalid Arg 'stop'".format(stop))
        self._start = start
        self._stop = stop

    @property
    def label(self):
        return "{} to {}".format(
            self.start.strftime("%m/%d/%Y"),self.stop.strftime("%m/%d/%Y")
        )

    @property
    def start(self): return self._start

    @property
    def stop(self): return self._stop

    @property
    def delta(self): return self.stop - self.start

    def get_data(self): return self._data
    def set_data(self,d):
        if not type(d) == type({}):
            raise EpochException("Invalid Arg 'data' ({})".format(d))
        self._data = d
    def del_data(self): self._data = {}
    data = property(get_data,set_data,del_data)

class MSDB(object):
    """ a more manageable wrapper around the msdb dictionary """
    def __init__(self,init=True):
        # initiate to none and read in states file
        self._rs = self._ks = self._ks = None
        try:
            self._ss = states.read()
        except IOError as e:
            raise MSDBException(
                'Unable to read state population file: {}'.format(e)
            )

        if init: self.read()
        else: self.reset()

    def reset(self):
        self._rs = {} # incident records dict
        self._ks = [] # key list

    def read(self,fpath=pth_msdb):
        # reset variables
        fin = None
        self.reset()
        try:
            # open a file object for the csv reader & read in header row
            fin = open(fpath)
            csvr = csv.reader(fin,delimiter='\t')
            _ = csvr.__next__()

            # create caliber and firearmspec objects
            try:
                cs = caliber.CaliberDS()
                fs = firearmspec.FirearmSpec()
            except caliber.CaliberException as e:
                print("Error reading caliber db: ({}). quitting...".format(e))
                return
            except firearmspec.FirearmSpecException as e:
                print("Error reading firearmspec db: ({}). quitting...".format(e))
                return

            # read contents and instantiate records db
            for row in csvr:
                try:
                    if 'MODEL' in row[FA_MDL]: continue
                    elif is_int(row[INC_NUM]): self._validate_inc_rec_(row)
                    else: self._validate_firearm_rec_(cs,fs,row)
                except firearm.FirearmException as e:
                    print("({}) invalid firearm: {}".format(self._ks[-1],e))
                    break

        except IOError as e:
            print("Error reading {}: {}".format(fpath,e))
        finally:
            if fin: fin.close()

    def write(self,fpath=pth_msdb):
        # set up headers
        hdr = "NUM\tCASE\tSHOOTER\tSHT DATA\tLOCATION\tDTG\tVICTIM\tLEO\tNOTES\t"
        hdr += "RPT\tTOA\tCKE\tSCENE\tVENUE\tTYPE\tTGR\tMOTIVE\tTGT\tTTL FA\n"
        fa_hdr = "\t\tFIREARMS\tMODEL\tTYPE\tMANUF\tFIRE ACT\tCALIBER\tCAPACITY\t"
        fa_hdr += "MAGAZINE(S)\tATTACHMENTS\tRNDS\tNOTES\tSTD\tMC\tLEG\tACQ DT\tACQ ST\n"

        # open file for writing
        fout = None
        try:
            # open the file for writing
            fout = open(fpath,'w')

            # write the main header
            fout.write(hdr)

            # write the records
            for i,key in enumerate(self._ks):
                # incident rec
                rec = self._rs[key]
                fout.write(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                        i+1,
                        key,
                        self._write_shooter_(rec['shooter']),
                        self._write_shooter_data_(rec['shooter-data']),
                        "{}, {}".format(rec['location']['city'],rec['location']['state']),
                        self._write_dtg_(rec['dtg']),
                        "{}/{}/{}".format(rec['killed'],rec['wounded'],rec['total']),
                        "{}/{}/{}".format(rec['leo-killed'],rec['leo-wounded'],rec['leo-total']),
                        rec['notes'],
                        rec['report'],
                        "{:.2f}".format(rec['toa']) if rec['toa'] else '',
                        self._write_cke_cause_(rec['cke']),
                        ", ".join(rec['scene']),
                        rec['venue'],
                        '{}{}'.format(rec['type']['pos'],rec['type']['act']) if rec['type'] else '',
                        rec['trigger'],
                        ", ".join(rec['motive']),
                        rec['tgt'],
                        rec['ttl-fa'],
                    )
                )

                # firearms records
                fout.write(fa_hdr)
                for f in rec['firearms']: f.write(fout)
        except IOError as e:
            print(e)
        finally:
            if fout: fout.close()

    def __getitem__(self,key): return self._rs[key]
    def __setitem__(self,key,value):
        raise MSDBException("Item assignment not supported")
    def __delitem__(self,key):
        raise MSDBException("Item deletion not supported")

    def keys(self,exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        return [k for k in self._ks if self._keep_(k,exclude,start,stop)]

#### INDIVIDUAL RECORDS

    def meets_criteria(self,key,criteria=CRIT_STRICT_SCENE):
        """
         Determines whether the incident at key meets 'Mass Shooting' criteria.
         There are multiple definitions. The Mother Jones DB used the criteria of
         4 or more killed upto 2013. Afterwhich IAW the Investigative Assistance
         for Violent Crimes Act of 2012 (Obama) definition of "three or more".
         This presents problems as the comparison of counts across the 2012/2013
         will be inaccurately reflected. Another issues with Mother Jones is they
         include the shooter in the victim count. (We have removed 8 cases where
         the total killed count (exclusive the shooter) is less than four and the
         total victim count is <= 4.) We want to use a compromise between USA
         Today (4 or more killed) and Mother Jones (3 or more victims) and use
         the following definition:
              4 or more killed excluding the shooter
           or 5 or more victims with at least 3 dead.
         But to get fair/accurate comparison, we need to find all cases prior to
         2013 that meet this criteria and were not entered in the MJ DB. Until
         then, this function will determine if the given case meets the orginal
         standard of 4 or more non-shooter victims.
        :param key: the incident to check
        :param criteria: as defined in CRIT_*
        :return: True or false
        """
        tk = self._rs[key]['killed']
        tw = self._rs[key]['wounded']

        if criteria == CRIT_NONE: return True
        elif criteria == CRIT_KLAREVES:
            if tk >= 6: return True
            else: return False
        elif criteria == CRIT_STRICT_SCENE:
            # check in notes for multiple scenes and require at least one to have 4 deaths
            if re_mult_scenes.search(self._rs[key]['notes']):
                for k,_ in re_mult_scenes_tally.findall(self._rs[key]['notes']):
                    if int(k) >= 4: return True
                return False
            else: return tk >= 4
        elif criteria == CRIT_LOOSE_SCENE:
            # check in notes for multiple scenes and require at least one to have 4 deaths
            if re_mult_scenes.search(self._rs[key]['notes']):
                for k,w in re_mult_scenes_tally.findall(self._rs[key]['notes']):
                    if int(k) >= 4 or (int(k) >= 3 and int(w) >= 2): return True
                return False
            else: return tk >= 4 or (tk >= 3 and tw >= 2)
        elif criteria == CRIT_STRICT_INC:
            return tk >= 4
        elif criteria == CRIT_LOOSE_INC:
            return tk >= 4 or (tk >= 3 and tw >= 2)
        else:
            raise MSDBException("Invalid Arg 'criteria' ({})".format(criteria))

    def idx2case(self,i,adjust=True):
        return self._ks[i-1] if adjust else self._ks[i]

    def case2idx(self,c,adjust=True):
        return self._ks.index(c)+1 if adjust else self._ks.index(c)

#### CHECKS/CONSISTENCY

    def consistency_check(self,fpath=pth_clog):
        """
        :param fpath' path of log file
        analyzes firearms by comparing with caliber and firearmspec for consistency
        and validity. writes a log file of results
        """
        # create caliber and firearmspec objects
        oC = caliber.CaliberDS()
        oF = firearmspec.FirearmSpec()

        inc = [] # incomplete firearm, cannot match specs (missing model and/or manuf
        nf = []  # firearms not found in firearmspec
        nc = []

        for key in self.keys():
            for i,f in enumerate(self._rs[key]['firearms']):
                # get manuf and model
                manuf,model = f.manufacturer,f.model

                # check for incomplete manuf/model
                if not manuf:
                    inc.append((key,i,"missing manufacturer"))
                    continue
                elif not model:
                    inc.append((key,i,"missing model"))
                    continue

                # get the firearm rec for the weapon
                fkey = oF.get_key(manuf,model)
                if not fkey:
                    nf.append(
                        (key,i,"no matching firearm specification record for {} {}".format(manuf,model))
                    )
                    continue

                # check caliber (but continue processing if not found)
                ckey = oC.find_name(f.caliber)
                if not ckey:
                    nc.append(
                        (key,i,"no matching caliber record for ({})".format(f.caliber))
                    )

        fout = None
        try:
            fout = open(fpath,'w')

            # start with incomplete records
            fout.write("## INCOMPLETE FIREARMS RECORDS ##\n")
            for i,(key,widx,msg) in enumerate(inc):
                fout.write("{}. {}: firearm #{} {}\n".format(i,key,widx+1,msg))
            if not inc: fout.write("No incomplete firearm records\n")
            fout.write("\n\n")

            # do firearm spec not found
            fout.write("## UNDEFINED FIREARM SPECIFICATION ##\n")
            for i,(key,widx,msg) in enumerate(nf):
                fout.write("{}. {}: firearm #{} {}\n".format(i,key,widx+1,msg))
            if not nf: fout.write("No undefined firearm specifications\n")
            fout.write("\n\n")

            # and caliber not found
            fout.write("## UNDEFINED CALIBER ##\n")
            for i,(key,widx,msg) in enumerate(nc):
                fout.write("{}. {}: firearm #{} {}\n".format(i,key,widx+1,msg))
            if not nc: fout.write("No undefined caliber\n")
            fout.write("\n\n")

        except IOError:
            print('Failed to write consistency.log')
        finally:
            if fout: fout.close()

#### BASIC COUNTS

    def num_incidents(self,exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        return len([k for _,k in self._enumerate_(exclude,start,stop)])

    def num_shooters(self,exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        return sum([len(self._rs[k]['shooter']) for _,k in self._enumerate_(exclude,start,stop)])

    def num_firearms(self,exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        return sum([len(self._rs[k]['firearms']) for _,k in self._enumerate_(exclude,start,stop)])

    def num_victims(self,column='killed',exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        return sum(self._rs[k][column] for _,k in self._enumerate_(exclude,start,stop))

#### NOTES SEARCH

    def shooter_ended_incident(self,key):
        """ shooter killed self, stopped or left the scene """
        # Note for cause, we assume if there is more than one shooter, both shooters
        # actions were the same
        if re_shooter_suicide.match(self._rs[key]['notes']): return True
        if self._rs[key]['cke'][0]['code'].endswith('S'): return True
        if self._rs[key]['cke'][0]['code'].startswith('E'):
            # have to make sure victim or bystander did not force shooter toleave
            for party in ['Victim','Bystander']:
                ss,_ = self.subdue_attempt(key,party)
                for s in ss:
                    if s[1] == 'during': return False
            return True
        return False

    def civilian_ended_incident(self,key,party='Victim'):
        """ victime subdued/killed or forced shooter to leave """
        if not party in ['Bystander','Victim']:
            raise MSDBException("Invalid Arg 'party ({})".format(party))
        for a in re_subdued_attempt.findall(self._rs[key]['notes']):
            if a[1] == party and a[2] == 'succeeded' and a[3] == 'during':
                return True
        c,p = self._rs[key]['cke'][0]['code'],self._rs[key]['cke'][0]['party']
        if not c.startswith('E') and p == party: return True
        return False

    def subdue_attempt(self,key,party):
        # there should only be one match in each notes based on collected incidents
        # but just in case
        if not party in ['Bystander','Victim']:
            raise MSDBException("Invalid Arg 'party' ({})".format(party))
        ss = [] # successes
        fs = [] # failures
        for a in re_subdued_attempt.findall(self._rs[key]['notes']):
            if a[1] != party: continue
            if a[2] == 'succeeded': ss.append((int(a[0]) if a[0] else 1,a[3]))
            elif a[2] == 'attempted': fs.append((int(a[0]) if a[0] else 1,a[3]))
        return ss,fs

    def ttl_rnds_fired(self,key):
        # check the firearms data first, then the notes
        try:
            return sum([f.rnds_fired for f in self._rs[key]['firearms']])
        except TypeError:
            rf = sum([int(m[0]) for m in re_ttl_rnds_fired.findall(self._rs[key]['notes'])])
            if rf: return rf
        return None

    def firearm_malfunction(self,key):
        # returns count of malfunctioning firearms in case
        return len([f for f in self._rs[key]['firearms'] if f.malfunction])

    def ttl_mags(self,key):
        # check event notes first
        tm = sum([int(m[0]) for m in re_ttl_mags_known.findall(self._rs[key]['notes'])])
        if tm: return tm
        tm = sum([f.ttl_mags_cnt for f in self._rs[key]['firearms'] if f.ttl_mags_cnt])
        if tm: return tm
        return None

    def mags_changed(self,key):
        # same as rnds fired, check firearms first, then the notes
        try:
            return sum([f.mags_changed for f in self._rs[key]['firearms']])
        except TypeError:
            rf = sum([int(m) for m in re_mag_change.findall(self._rs[key]['notes'])])
            if rf: return rf
        return None

    def people_present(self,key):
        m = re_ppl_present.search(self._rs[key]['notes'])
        if m: return int(m.groups()[0])
        else: return None

    def acq_state(self,key):
        # only looks in note section
        try:
            ast = re_acq_state.search(self._rs[key]['notes']).group()
            ast = [x.split(':') for x in re_acq_state_term.findall(ast)]
            return [(x[0],int(x[1]) if is_int(x[1]) else 1) for x in ast]
        except AttributeError:
            return None

#### WEAPON RELATED

    def firearms(self,key): return [(i,f) for i,f in enumerate(self._rs[key]['firearms'])]

    def get_firearm(self,key,widx): return self._rs[key]['firearms'][widx]

    def same_category(self,key):
        """
        Returns the firearm category if all firearms in the incident were of the same
        category otherwise returns None
        """
        wc = list(set(f.category for f in self._rs[key]['firearms']))
        if len(wc) > 1: return None
        else: return wc[0]

    def same_caliber(self,key):
        wc = [f.caliber for f in self._rs[key]['firearms'] if f.caliber]
        for i,c in enumerate(wc):
            if c == '.233' or c == '5.56x45': wc[i] = '.233/5.45x45'
            elif c == '.308' or c == '7.62x39': wc[i] = '.308/7.62x39'
        wc = list(set(wc))
        if not wc or len(wc) > 1: return None
        else: return wc[0]

    def majority_rnds_fired(self,key,thresh=0.85):
        """
         returns the firearm firing the majority of rounds
         NOTE: only looks at incidents where the number rounds fired is known for
          each firearm
        """
        try:
            rfs = [f.rnds_fired for f in self._rs[key]['firearms']]
            ttl = sum(rfs)
            for i,rf in enumerate(rfs):
                if rf/ttl > thresh: return i
        except TypeError:
            return None
        return None

    def primary_firearm(self,key,thresh=0.85):
        """
        returns index of the primary firearm (if exists) for incident defined as:
         a) the only firearm used
         b) firearm identified as primary
         c) the firearm that fired a majority of rounds (as defined in thresh)
         d) the first firearm if all firearms are the same (manuf, model)
         e) the firearm category if all firearms are of the same category
          i the firearm that fired the most rounds if category is 4 and 5
          ii the firearm with highest/known mag cap if category level 3 and same caliber

        """
        # a) the only firearm used
        if len(self._rs[key]['firearms']) == 1: return 0

        # local variable setup
        # check notes for a total rounds fired
        trf = self.ttl_rnds_fired(key)
        if not trf: trf = None
        rfs = []
        caps = []
        same_cat = self.same_category(key)
        same_cal = self.same_caliber(key)
        same_firearm = True
        first_firearm = {
            # TODO: have to handle case like SKS, Glock that may have a generation
            'manuf':self._rs[key]['firearms'][0].manufacturer,
            'model':self._rs[key]['firearms'][0].model,
        }

        # iterate each firearm
        for i,f in enumerate(self._rs[key]['firearms']):
            # TODO: have to handle multiple shooters with identified primary firearms
            if f.is_primary(): return i
            rfs.append(f.rnds_fired)
            caps.append(f.capacity if f.capacity else float('-inf'))
            manuf,model = f.manufacturer,f.model
            if manuf != first_firearm['manuf'] or model != first_firearm['model']:
                same_firearm = False

        # c) the firearm that fired a majority of rounds (as defined in thresh)
        if trf and None not in rfs:
            for i,rf in enumerate(rfs):
                if rf/trf > thresh: return i

        # d) the first firearm if all firearms are the same
        # TODO: return the firearm with the highest mag cap
        if same_firearm: return 0

        # e) firearms all share the same category
        if same_cat:
            # i return the firearm that fired the most rounds if known
            if taxonomy.firearm_cat_lvl(same_cat) >= 4:
                if not None in rfs: return rfs.index(max(rfs))
            elif taxonomy.firearm_cat_lvl(same_cat) == 3 and same_cal:
                return caps.index(max(caps))

        return None

    def acquired_dates(self,key):
        """ returns a list of acquistion dates if known """
        return [f.date_acquired for f in self._rs[key]['firearms']]

    def firearm_victims(self,key):
        """
         returns list of tuples t = (widx,killed,wounded) (if unknown, will be
         set to None
        """
        # easy if only one firearm
        if len(self._rs[key]['firearms']) == 1:
            return [(0,self._rs[key]['killed'],self._rs[key]['wounded'])]

        # if more than one firearm, have to check for annotations in firearm notes
        vs = []
        for i,f in enumerate(self._rs[key]['firearms']):
            k,w = f.fatalities
            if k is not None: vs.append((i,k,w))
            else: vs.append((i,None,None))
        return vs

    def ttl_firearms(self,key): return self._rs[key]['ttl-fa']

    def firearms_used(self,key): return len(self._rs[key]['firearms'])

#### COUNTS BY INCIDENT

    def firearms_per_incident(self,exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        """
        returns a dict of the form case # ->(used,total)
        :param exclude; if True will not tally records that do not meet criteria
        :param start: the first date to tally
        :param stop: the last date to tally
        :return: firearm per incident count
        """
        return {
            i+1:(len(self._rs[k]['firearms']),self._rs[k]['ttl-fa']) for i,k in self._enumerate_(exclude,start,stop)
        }

    def fawb_per_incident(self,exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        """
        returns a dict of the form case # -> total fawb violations
        :param exclude; if True will not tally records that do not meet criteria
        :param start: the first date to tally
        :param stop: the last date to tally
        :return: assault type firearm per incident count
        """
        at = self.tac_per_incident(exclude,start,stop)
        hcm = self.hcm_per_incident(exclude,start,stop)
        return {i:at[i]+hcm[i] for i in at}

    def tac_per_incident(self,exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        """
        returns a dict of the form case # -> total # AT Weapons
        :param exclude; if True will not tally records that do not meet criteria
        :param start: the first date to tally
        :param stop: the last date to tally
        :return: assault type firearm per incident count
        """
        tac = {}
        for i,key in self._enumerate_(exclude,start,stop):
            tn = 0
            for f in self._rs[key]['firearms']:
                if taxonomy.is_tac_firearm(f.type): tn += 1
            tac[i+1] = tn
        return tac

    def hcm_per_incident(self,exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        """
        returns a dict of the form case # -> total # hcm
        :param exclude; if True will not tally records that do not meet criteria
        :param start: the first date to tally
        :param stop: the last date to tally
        :return: hcm incident count
        """
        hcm = {}
        for i,key in self._enumerate_(exclude,start,stop):
            cnt = 0
            for f in self._rs[key]['firearms']:
                for n,c in f.magazines():
                    if is_int(c):
                        if not is_int(n): n = 1 # set unknowns to 1
                        if c > 10: cnt += n
            if not hcm: # check notes if nothing was found in firearms
                for num,cap,_ in re_ttl_mags_known.findall(self._rs[key]['notes']):
                    if is_int(cap) and int(cap) > 10: cnt += int(num)
            hcm[i+1] = cnt
        return hcm

    def extended_per_incident(self,exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        """
        returns a dict of the form case # -> total # exctended mags
        :param exclude; if True will not tally records that do not meet criteria
        :param start: the first date to tally
        :param stop: the last date to tally
        :return: over standard mag during incident count
        """
        # NOTE: cannot use notes section here as over-std applies to each firearm
        ost = {}
        for i,key in self._enumerate_(exclude,start,stop):
            cnt = 0
            for f in self._rs[key]['firearms']:
                std = f.std_mag_cap
                if not std: continue
                else: std = max(std)
                for n,c in f.magazines():
                    if is_int(c):
                        if not is_int(n): n = 1 # set unknowns to 1
                        if c > std: cnt += n
            ost[i+1] = cnt
        return ost

    # TODO: make this a by key function??
    def fatality_rate(self,exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        """
        returns a dict of the form case # -> fatality rate where fatality rate =
         # killed / # shot
        :param exclude; if True will not tally records that do not meet criteria
        :param start: the first date to tally
        :param stop: the last date to tally
        :return: fatality rate incident
        """
        return {
            i+1:self._rs[k]['killed']/self._rs[k]['total'] for i,k in self._enumerate_(exclude,start,stop)
        }

#### WEAPONS BY INCIDENT

    def capacities(self,key):
        """
        returns a list of the unique magazine capacities used in the incident
        NOTE: this does not consider mag capacity on a firearm by firearm basis
        NOTE: a list of returned capacities does not imply all firearm capacities
         are known
        """
        # check notes first and return if found
        caps = []
        for _,cap,_ in re_ttl_mags_known.findall(self._rs[key]['notes']):
            try:
                caps.append(int(cap))
            except ValueError:
                pass
        if caps: return list(set(caps))

        # then ttl mags, firearm capacity
        for f in self._rs[key]['firearms']: caps.extend(f.capacities())

        return list(set(caps))

    def tac_firearm_used(self,key):
        """ returns 'True', 'False', 'Unknown' for at firearm used in incident """
        ret = 'False'
        for f in self._rs[key]['firearms']:
            tw = f.is_tac_firearm()
            if tw == 'True': return 'True'
            elif tw == 'Unknown': ret = 'Unknown'
        return ret

    def hcm_used(self,key):
        """ returns 'True', 'False', 'Unknown' for hcm used in incident """
        # Three posibilities:
        #  a) one or more firearms has unknown capacity magazine
        #  b) one or more firearms has capacity of 10 or more
        #  c) all firearms have capacity of 10 or less

        # TODO: have to ensure notes covers all firearms
        # checks notes ASSUMING all firearms are covered
        caps = []
        unk = False
        for _,cap,_ in re_ttl_mags_known.findall(self._rs[key]['notes']):
            if is_int(cap): caps.append(int(cap))
            elif cap == '?': unk = True
        if caps:
            if unk: # return True if max is over 10, otherwise unknown
                if max(caps) > 10: return 'True'
                else: return 'Unknown'
            else: # return true if max is over 10, otherwise false
                if max(caps) > 10: return 'True'
                else: return 'False'

        # check each firearm
        mags = [f.used_hcm() for f in self._rs[key]['firearms']]

        # check for an unknown
        if 'True' in mags: return 'True'
        elif 'Unknown' in mags: return 'Unknown'
        else: return 'False'

    def over_std_used(self,key):
        """
         returns 'True', 'False', 'Unknown' or 'Not-listed' for at firearm used
         in incident
        """
        overs = [f.used_over_std() for f in self._rs[key['firearms']]]
        if 'True' in overs: return 'True'
        elif 'Not-listed' in overs: return 'Not-listed'
        elif 'Unknown' in overs: return 'Unknown'
        else: return 'False'

    def acq_states(self,key):
        """ returns a list of tuples t = (State,Count) """
        ss = self.acq_state(key)
        if ss:
            # set ast and append unknowns
            ast = ss
            ast.append(('XX',len(self._rs[key]['firearms'])-sum([w[1] for w in ast])))
        else:
            ss = {}
            for st in [f.state_acquired for f in self._rs[key]['firearms']]:
                if st in ss: ss[st] += 1
                else: ss[st] = 1
            ast = [(st,ss[st]) for st in ss]
        return ast

#### HISTOGRAMS

    def hist_state(self,exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        """
        returns hist of ms by state through period if defined in start, stop
        :param exclude; if True will not tally records that do not meet criteria
        :param start: the first date to tally
        :param stop: the last date to tally
        :return: hist of states during start-stop
        """
        hs = {}
        for i,key in self._enumerate_(exclude,start,stop):
            v = self._rs[key]['location']['state']
            if v in hs: hs[v].append(i+1)
            else: hs[v] = [i+1]
        return hs

    def hist_msr(self,by='2019',exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        """
        returns hist mass shooting rate (msr) for each state (found in msdb)
        :param by: what population to use
        :param exclude; if True will not tally records that do not meet criteria
        :param start: the first date to tally
        :param stop: the last date to tally
        :return: hist of states during start-stop
        see https://oag.ca.gov/sites/all/files/agweb/pdfs/cjsc/stats/computational_formulas.pdf
        """
        # validate by
        col = 'pop-'
        if by in ['1990','2000','2010','2019']: col += by[2:]
        else: raise MSDBException("Invalid Arg 'by' ({})".format(by))

        # get the state hist & return the adjusted count
        h = self.hist_state(exclude,start,stop)
        return {k:round(len(h[k])/self._ss[k][col]*100000,4) for k in h}

    def hist_msr2(self,exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        """
        returns hist mass shooting rate (msr) for each state (found in msdb)
        :param exclude; if True will not tally records that do not meet criteria
        :param start: the first date to tally
        :param stop: the last date to tally
        :return: hist of states during start-stop
        see https://oag.ca.gov/sites/all/files/agweb/pdfs/cjsc/stats/computational_formulas.pdf
        """
        d1990 = datetime.date(1990,12,31)
        d2000 = datetime.date(2000,12,31)
        d2010 = datetime.date(2010,12,31)
        d2019 = datetime.date(2019,12,31)
        ps = [d1990,d2000,d2010,d2019]
        hs = {d1990:{},d2000:{},d2010:{},d2019:{}}
        for i,key in self._enumerate_(exclude,start,stop):
            # get incident date and determine where to place the record
            doi = self._rs[key]['dtg']['date']
            d = d2019
            if doi <= d1990: d = d1990
            elif doi <= d2000: d = d2000
            elif doi <= d2010: d = d2010
            elif doi <= d2019: d = d2019

            # get the state and add to the right subdict
            v = self._rs[key]['location']['state']
            if v in hs[d]: hs[d][v].append(i+1)
            else: hs[d][v] = [i+1]

        # determine the msr for each date period
        msr = {}
        for p in ps:
            for s in hs[p]:
                pi = len(hs[p][s]) / self._ss[s]['pop-' + str(p.year)[2:]] * 100000
                if s in msr: msr[s] += pi
                else: msr[s] = pi
        return {s:round(msr[s],4) for s in msr}

    def hist_date(self,by='year',exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        """
        returns a dict of dates according to by
        :param by: oneof {'year','month','dow'=(day of week)}
        :param exclude; if True will not tally records that do not meet criteria
        :param start: the first date to tally
        :param stop: the last date to tally
        :return: dict of date based count
        """
        if not by in ['year','month','dow']:
            raise MSDBException("Invalid Arg 'by' ({})".format(by))
        hd = {}
        for i,key in self._enumerate_(exclude,start,stop):
            dt = self._rs[key]['dtg']['date']
            if by == 'year': dt = dt.year
            elif by == 'month': dt = i2m[dt.month]
            else: dt = i2dow[dt.isoweekday()]
            if dt in hd: hd[dt].append(i+1)
            else: hd[dt] = [i+1]
        return hd

    def hist_elapsed(self,exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        """
        returns hist of ms by elapsed time between shootings in days
        :param exclude; if True will not tally records that do not meet criteria
        :param start: the first date to tally
        :param stop: the last date to tally
        :return: elapsed time hist
        """
        he = {}
        for i,e in self.chrono_elapsed(exclude,start,stop):
            if e in he: he[e].append(i)
            else: he[e] = [i]
        return he

    def chrono_elapsed(self,exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        """
        returns list of elapsed times between shootings
        :param exclude; if True will not tally records that do not meet criteria
        :param start: the first date to tally
        :param stop: the last date to tally
        :return: chronological list of tuples t = (index,elapsed time from previous index)
        NOTE:
         1. the elapsed time for index = i is the difference between i and i-1
         2. the elapsed time for the first record is 0
        """
        es = []
        last = self.num_incidents() - 1
        for i,key in self._enumerate_(exclude,start,stop):
            if i == last: es.append((i+1,0)) # last record has no preceding
            else:
                es.append(
                    (i+1,
                     abs((self._rs[key]['dtg']['date'] - self._rs[self._ks[i+1]]['dtg']['date']).days)
                     )
                )
        return es[::-1] # return from earliest to latest event

    def hist_demographics(self,by='age',exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        """
        returns hist of ms by shooter demographics
        :param by: one of demographic category found in sht_demo_cat
        :param exclude; if True will not tally records that do not meet criteria
        :param start: the first date to tally
        :param stop: the last date to tally
        :return: dict of shooter demographics
        """
        if not by in sht_demo_cat:
            raise MSDBException("Invalid Arg 'by' ({})".format(by))
        hd = {}
        for i,key in self._enumerate_(exclude,start,stop):
            for sd in self._rs[key]['shooter-data']:
                d = sd[by]
                try:
                    d = sht_demo_c2l[by][d]
                except KeyError:
                    pass
                if d in hd: hd[d].append(i+1)
                else: hd[d] = [i+1]
        return hd

    def hist_victims(self,column='total',vtype='all',exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        """
        returns a dict of victims under column
        :param column: = oneof {'killed','wounded','total'}
        :param vtype: = oneof {'all','civilian','leo']
        :param exclude; if True will not tally records that do not meet criteria
        :param start: the first date to tally
        :param stop: the last date to tally
        :return: dict of victims->count
        """
        if not column in ['killed','wounded','total']:
            raise MSDBException("Invalid Arg 'column' ({})".format(column))
        if not vtype in ['all','civilian','leo']:
            raise MSDBException("Invalid Arg 'vtype' ({})".format(vtype))

        hv = {}
        for i,key in self._enumerate_(exclude,start,stop):
            # find the specified value
            if vtype == 'all': v = self._rs[key][column]
            elif vtype == 'leo': v = self._rs[key]['leo-'+column]
            else: v = self._rs[key][column] - self._rs[key]['leo-'+column]
            # add to histogram
            if v in hv: hv[v].append(i+1)
            else: hv[v] = [i+1]
        return hv

    def hist_report(self,exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        """
        returns hist of ms by report type
        :param exclude; if True will not tally records that do not meet criteria
        :param start: the first date to tally
        :param stop: the last date to tally
        :return: dict of official reports
        """
        hr = {0:[],1:[],2:[]}
        for i,key in self._enumerate_(exclude,start,stop):
            hr[self._rs[key]['report']].append(i+1)
        return {'full':hr[2],'partial':hr[1],'none':hr[0]}

    def hist_toa(self,exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        """
        returns hist of cke-time by interval
        :param exclude; if True will not tally records that do not meet criteria
        :param start: the first date to tally
        :param stop: the last date to tally
        :return: the cke time histograms
        """
        ht = {}
        for i,key in self._enumerate_(exclude,start,stop):
            ct = self._rs[key]['toa']
            if ct is None: ct = 'Unknown'
            if ct in ht: ht[ct].append(i+1)
            else: ht[ct] = [i+1]
        return ht

    def hist_cke(self,pri=True,exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        """
        returns hist of cke causes as a dict
        :param pri: if True will differentiate E-K from K otherwise will only
        look at the 'K'
        :param exclude; if True will not tally records that do not meet criteria
        :param start: the first date to tally
        :param stop: the last date to tally
        :return; the cke cause histogram
        """
        hc = {}
        for i,key in self._enumerate_(exclude,start,stop):
            # for now all multi-shooter incidents have the same cause
            cause = self._rs[key]['cke'][0]
            c,p = cause['code'],cause['party']
            if not pri and '-' in c: c = c.split('-')[1]
            cause = ':'.join((c,p))
            if cause in hc: hc[cause].append(i+1)
            else: hc[cause] = [i+1]
        return hc

    def hist_stopped_by(self,exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        """
        returns a hist of the party stopping incident, one of
        {
         'shooter' = the shooter ended it by either leaving, suicide (prior to leo
          intervention) or just stopped
         'leo' = law enforcement officer
         'victim' = victim(s) subdued shooter(s) until leo arrival or actions led
           to victim's death/capture
         'bystander' = bystander(s) subdued shooter(s0 until leo arrival or actions
           lead to shooter's death/capture
        }
        :param exclude: if True will not tally records that do not meet criteria
        :param start: the first date to tally
        :param stop: the last date to tally
        :return: hist of responsible parties
        NOTE: If for victime or bystander, the intervention occurred after the
         shooter left, the shooter will be tallied as the 'ending' party
        """
        hs = {'shooter':[],'leo':[],'victim':[],'bystander':[]}
        for i,key in self._enumerate_(exclude,start,stop):
            if self.shooter_ended_incident(key): hs['shooter'].append(i+1)
            elif self.civilian_ended_incident(key,'Victim'): hs['victim'].append(i+1)
            elif self.civilian_ended_incident(key,'Bystander'): hs['bystander'].append(i+1)
            else: hs['leo'].append(i+1)
        return hs

    def hist_intervention_attempt(self,exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        """
        returns matrix hist of intervention attempts
        :param exclude; if True will not tally records that do not meet criteria
        :param start: the first date to tally
        :param stop: the last date to tally
        :return; the intervention histogram
        """
        ia = {
            'Victim':{
                'Success':{'During':[],'After':[]},
                'Failure':{'During':[],'After':[]},
            },
            'Bystander':{
                'Success':{'During':[],'After':[]},
                'Failure':{'During':[],'After':[]},
            },
        }
        for i,key in self._enumerate_(exclude,start,stop):
            for party in ['Victim','Bystander']:
                ss,fs = self.subdue_attempt(key,party)
                for s in ss: ia[party]['Success'][s[1].capitalize()].append(i+1)
                for f in fs: ia[party]['Failure'][f[1].capitalize()].append(i+1)
        return ia

    def hist_scene(self,exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        """
        returns hist of scenes (only inside, outside) as a dict
        :param exclude; if True will not tally records that do not meet criteria
        :param start: the first date to tally
        :param stop: the last date to tally
        :return; the scene histogram
        """
        hs = {'inside':[],'outside':[]}
        for i,key in self._enumerate_(exclude,start,stop):
            for scene in self._rs[key]['scene']: hs[scene.split(' ')[0]].append(i+1)
        return hs

    def hist_venue(self,exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        """
        returns hist of venue
        :param exclude; if True will not tally records that do not meet criteria
        :param start: the first date to tally
        :param stop: the last date to tally
        :return; the venue histogram
        """
        hv = {}
        for i,key in self._enumerate_(exclude,start,stop):
            v = self._rs[key]['venue']
            if v in hv: hv[v].append(i+1)
            else: hv[v] = [i+1]
        return hv

    def hist_trigger(self,exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        """
        returns hist of trigger codes
        :param exclude; if True will not tally records that do not meet criteria
        :param start: the first date to tally
        :param stop: the last date to tally
        :return; the initiation code histogram
        """
        hs = {tc:[] for tc in tgr_codes}
        for i,key in self._enumerate_(exclude,start,stop):
            hs[self._rs[key]['trigger']].append(i+1)
        return {tgr_c2l[tc]:hs[tc] for tc in hs if hs[tc]}

    def hist_motives(self,exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        """
        returns hist of initiation codes
        :param exclude; if True will not tally records that do not meet criteria
        :param start: the first date to tally
        :param stop: the last date to tally
        :return; the initiation code histogram
        """
        hm = {mv:[] for mv in mtv_codes}
        for i,key in self._enumerate_(exclude,start,stop):
            for mv in self._rs[key]['motive']: hm[mv].append(i+1)
        return {mtv_c2l[mv]:hm[mv] for mv in hm if hm[mv]}

    def hist_target(self,exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        """
        returns hist of target tuples
        :param exclude; if True will not tally records that do not meet criteria
        :param start: the first date to tally
        :param stop: the last date to tally
        :return; the target tuple histogram
        """
        ht = {}
        for i,key in self._enumerate_(exclude,start,stop):
            tt = self._rs[key]['tgt']
            tt = "{}:{}".format(tgt_c2l_tgt[tt[0]],tgt_c2l_loc[tt[1]])
            if tt in ht: ht[tt].append(i+1)
            else: ht[tt] = [i+1]
        return ht

    def hist_used_firearms(self,exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        """
        returns hist of # of firearms of used
        :param exclude; if True will not tally records that do not meet criteria
        :param start: the first date to tally
        :param stop: the last date to tally
        :return; the # firearms used histogram
        """
        hw = {}
        for i,key in self._enumerate_(exclude,start,stop):
            l = len(self._rs[key]['firearms'])
            if l in hw: hw[l].append(i+1)
            else: hw[l] = [i+1]
        return hw

    def hist_all_firearms(self,exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        """
        returns hist of # of total firearms
        :param exclude; if True will not tally records that do not meet criteria
        :param start: the first date to tally
        :param stop: the last date to tally
        :return; the all firearms histogram
        """
        ha = {}
        for i,key in self._enumerate_(exclude,start,stop):
            a = self._rs[key]['ttl-fa']
            if a in ha: ha[a].append(i+1)
            else: ha[a] = [i+1]
        return ha

    def hist_people_present(self,exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        """
        returns hist of # of total firearms
        :param exclude: if True will not tally records that do not meet criteria
        :param start: the first date to tally
        :param stop: the last date to tally
        :return; the all firearms histogram
        """
        hp = {}
        for i,key in self._enumerate_(exclude,start,stop):
            pp = self.people_present(key)
            if not pp: pp = 'Unknown'
            if pp in hp: hp[pp].append(i+1)
            else: hp[pp] = [i+1]
        return hp

    def hist_tac_firearm(self,exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        """
        returns a dict of tactical firearms used in incident
        :param exclude: if True will not tally records that do not meet criteria
        :param start: the first date to tally
        :param stop: the last date to tally
        :return; the tactical firearm used histogram
        """
        ht = {'True':[],'False':[],'Unknown':[]}
        for i,key in self._enumerate_(exclude,start,stop):
            ht[self.tac_firearm_used(key)].append(i+1)
        return ht

    def hist_hcm(self,exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        """
        returns a dict of hcm used in incident
        :param exclude; if True will not tally records that do not meet criteria
        :param start: the first date to tally
        :param stop: the last date to tally
        :return; the hcm used histogram
        """
        ht = {'True':[],'False':[],'Unknown':[]}
        for i,key in self._enumerate_(exclude,start,stop):
            ht[self.hcm_used(key)].append(i+1)
        return ht

    def hist_over_std(self,exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        """
        returns a dict of over standard magazines used in incident
        :param exclude; if True will not tally records that do not meet criteria
        :param start: the first date to tally
        :param stop: the last date to tally
        :return; the over std magazine used histogram
        """
        ht = {'True':[],'False':[],'Unknown':[],'Not-listed':[]}
        for i,key in self._enumerate_(exclude,start,stop):
            ht[self.over_std_used(key)].append(i+1)
        return ht

    def hist_fawb(self,exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        """
        returns a dict of a fawb ban firearm/magazine used in incident
        :param exclude; if True will not tally records that do not meet criteria
        :param start: the first date to tally
        :param stop: the last date to tally
        :return; the fawb ban firearm/magazine used histogram
        """
        ht = {'True':[],'False':[],'Unknown':[]}
        for i,key in self._enumerate_(exclude,start,stop):
            t = self.tac_firearm_used(key)
            h = self.hcm_used(key)
            if t == 'True' or h == 'True': ht['True'].append(i+1)
            elif t == 'False' and h == 'False': ht['False'].append(i+1)
            else: ht['Unknown'].append(i+1)
        return ht

    def hist_firearm(self,per_case=True,exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        """
        returns a dict of firearms (manuf model)
        :param per_case: if True will only consider a firearm once per case
        :param exclude; if True will not tally records that do not meet criteria
        :param start: the first date to tally
        :param stop: the last date to tally
        :return; the firearm histogram
        NOTE:
          if per_case is True, for each firearm, will return the incidents the
            firearm was used
          if per_case is False, for each firearm, will return the total number of
            times the firearm was used acroos all incidents
          For example in the Planned Parenthood shooting, the shooter used four
          SKS's if per_case is True, SKS will only be counted once, otherwise it
          will be counted four times
        """
        hw = {}
        for i,key in self._enumerate_(exclude,start,stop):
            recorded = []
            for f in self._rs[key]['firearms']:
                wname = f.name
                if per_case:
                    if wname in recorded: continue
                    recorded.append(wname)
                if wname in hw: hw[wname].append(i+1)
                else: hw[wname] = [i+1]
        return hw

    def hist_firearm_type(self,per_case=True,exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        """
        returns a hist of firearm types as a dict
        :param per_case: if True will only consider a firearm once per case
        :param exclude; if True will not tally records that do not meet criteria
        :param start: the first date to tally
        :param stop: the last date to tally
        :return; the firearm type histogram
        NOTE:
          if per_case is True, for each firearm type, will return the incidents the
            firearm type was used
          if per_case is False, for each firearm type, will return the total number
            of times the firearm type was used acroos all incidents
        """
        hw = {wt:[] for wt in taxonomy.firearm_types}
        for i,key in self._enumerate_(exclude,start,stop):
            recorded = []
            for f in self._rs[key]['firearms']:
                wt = f.type
                if per_case:
                    if wt in recorded: continue
                    recorded.append(wt)
                hw[wt].append(i+1)
        return hw

    def hist_firearm_cat(self,per_case=True,exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        """
        returns a hist of firearm categories as a dict
        :param per_case: if True will only consider a firearm once per case
        :param exclude; if True will not tally records that do not meet criteria
        :param start: the first date to tally
        :param stop: the last date to tally
        :return; the firearm category histogram
        NOTE:
          if per_case is True, for each firearm category, will return the incidents the
            firearm type was used
          if per_case is False, for each firearm type, will return the total number
            of times the firearm category was used acroos all incidents
        """
        # TODO: allow for category level to be passed
        hc = {}
        for i,key in self._enumerate_(exclude,start,stop):
            recorded = []
            for f in self._rs[key]['firearms']:
                wc = f.category
                if per_case:
                    if wc in recorded: continue
                    recorded.append(wc)
                if wc in hc: hc[wc].append(i+1)
                else: hc[wc] = [i+1]
        return hc

    def hist_fire_action(self,per_case=True,exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        """
        returns a hist of the fire action as a dict
        :param per_case: if True will only consider a firearm once per case
        :param exclude; if True will not tally records that do not meet criteria
        :param start: the first date to tally
        :param stop: the last date to tally
        :return; the fire action histogram
        NOTE:
          if per_case is True, for fire action type, will return the incidents the
            fire action type was used
          if per_case is False, for each fire action type, will return the total
            number of times the fire action type was used acroos all incidents
        """
        hf = {f:[] for f in firearmspec.fire_acts+['Unknown']}
        for i, key in self._enumerate_(exclude, start, stop):
            recorded = []
            for f in self._rs[key]['firearms']:
                fa = f.fire_act if f.fire_act else 'Unknown'
                if per_case:
                    if fa in recorded: continue
                    recorded.append(fa)
                hf[fa].append(i+1)
        return hf

    def hist_caliber(self,per_case=True,exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        """
        returns a hist of caliber and counts
        :param per_case: if True will only consider a firearm once per case
        :param exclude; if True will not tally records that do not meet criteria
        :param start: the first date to tally
        :param stop: the last date to tally
        :return; the caliber histogram
        NOTE:
          if per_case is True, for each caliber, will return the incidents the
            caliber was used
          if per_case is False, for each caliber, will return the total number of
            times the caliber was used acroos all incidents
        """
        hc = {}
        for i,key in self._enumerate_(exclude,start,stop):
            recorded = []
            for f in self._rs[key]['firearms']:
                cal =  f.caliber if f.caliber else 'Unknown'
                if per_case:
                    if cal in recorded: continue
                    recorded.append(cal)
                if cal in hc: hc[cal].append(i+1)
                else: hc[cal] = [i+1]
        return hc

    def hist_capacity(self,per_case=True,exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        """
        returns a hist of capacity
        :param per_case: if True will only consider a firearm once per case
        :param exclude; if True will not tally records that do not meet criteria
        :param start: the first date to tally
        :param stop: the last date to tally
        :return; the initial mag capacity histogram
        NOTE:
          if per_case is True, for each capacity, will return the incidents the
            capacity was used
          if per_case is False, for each caliber, will return the total number of
            times the caliber was used acroos all incidents
        NOTE:
          this is not a 1 to 1 histogram see below
        """
        hc = {}
        for i,key in self._enumerate_(exclude,start,stop):
            recorded = []
            for f in self._rs[key]['firearms']:
                cap = f.capacity if f.capacity else 'Unknown'
                if per_case:
                    if cap in recorded: continue
                    recorded.append(cap)
                if cap in hc: hc[cap].append(i+1)
                else: hc[cap] = [i+1]
        return hc

    def hist_highest_capacity(self,exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        """
        returns a hist of the highest capacity used in each incident
        :param exclude; if True will not tally records that do not meet criteria
        :param start: the first date to tally
        :param stop: the last date to tally
        :return; the highest capacity histogram
        """
        hc = {}
        for i,key in self._enumerate_(exclude,start,stop):
            try:
                cap = max(self.capacities(key))
            except ValueError:
                cap = 'Unknown'
            if cap in hc: hc[cap].append(i+1)
            else: hc[cap] = [i+1]
        return hc

    def hist_ttl_mags(self,exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        """
        returns a hist of total mags used per incident
        :param exclude; if True will not tally records that do not meet criteria
        :param start: the first date to tally
        :param stop: the last date to tally
        :return; the total mags histogram
        """
        # NOTE: for this we ignore the operator, treating each one as '='
        ht = {}
        for i,key in self._enumerate_(exclude,start,stop):
            tm = self.ttl_mags(key)
            if tm is None: tm = 'Unknown'
            if tm in ht: ht[tm].append(i+1)
            else: ht[tm] = [i+1]
        return ht

    def hist_attachments(self,per_case=True,exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        """
        returns a hist of attachments
        :param per_case: if True will only consider a firearm once per case
        :param exclude; if True will not tally records that do not meet criteria
        :param start: the first date to tally
        :param stop: the last date to tally
        :return; the attachments histogram
        """
        ha = {}
        for i,key in self._enumerate_(exclude,start,stop):
            recorded = []
            for f in self._rs[key]['firearms']:
                if not f.attachments:
                    if per_case:
                        if 'None' in recorded: continue
                        recorded.append('None')
                    if 'None' in ha: ha['None'].append(i+1)
                    else: ha['None'] = [i+1]
                for atc in f.attachments:
                    if per_case:
                        if atc in recorded: continue
                        recorded.append(atc)
                    if atc in ha: ha[atc].append(i+1)
                    else: ha[atc] = [i+1]
        return ha

    def hist_rnds_fired(self,exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        """
        returns a hist of rounds fired per incident
        :param exclude; if True will not tally records that do not meet criteria
        :param start: the first date to tally
        :param stop: the last date to tally
        :return; the rnds fired per incident histogram
        """
        # NOTE: not considering the operator here, only the number
        hr = {}
        for i,key in self._enumerate_(exclude,start,stop):
            rf = self.ttl_rnds_fired(key)
            if rf is None: rf = 'Unknown'
            if rf in hr: hr[rf].append(i+1)
            else: hr[rf] = [i+1]
        return hr

    def hist_mags_changed(self,exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        """
        returns a hist of mags changed
        :param exclude; if True will not tally records that do not meet criteria
        :param start: the first date to tally
        :param stop: the last date to tally
        :return; the mags changed histogram
        """
        # NOTE: not considering the operator here, only the number
        hm = {}
        for i,key in self._enumerate_(exclude,start,stop):
            mc = self.mags_changed(key)
            if mc is None: mc = 'Unknown'
            if mc in hm: hm[mc].append(i+1)
            else: hm[mc] = [i+1]
        return hm

    def hist_legality(self,by='legal',exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        """
        returns a hist of firearms in ms by legality category
        :param by; one of hist code categories found in firearm
        :param exclude; if True will not tally records that do not meet criteria
        :param start: the first date to tally
        :param stop: the last date to tally
        :return; the firearm legality histogram
        """
        if not by in firearm.firearm_lgl_cat:
            raise MSDBException("Invalid Arg 'by' ({})".format(by))
        hl = {}
        for i,key in self._enumerate_(exclude,start,stop):
            for f in self._rs[key]['firearms']:
                cd = f['legality'][by]
                lbl = firearm.firearm_lgl_c2l[by][cd]
                if lbl in hl: hl[lbl].append(i+1)
                else: hl[lbl] = [i+1]
        return hl

    def hist_acq_delta(self,exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        """
        returns a hist of days between firearm acquisition and incident
        :param exclude; if True will not tally records that do not meet criteria
        :param start: the first date to tally
        :param stop: the last date to tally
        :return; the acquisition delta histogram
        """
        hd = {'Unknown':[]}
        for i,key in self._enumerate_(exclude,start,stop):
            d = self._rs[key]['dtg']['date']
            for f in self._rs[key]['firearms']:
                try:
                    delta = (d - f.date_acquired).days
                    if delta in hd: hd[delta].append(i+1)
                    else: hd[delta] = [i+1]
                except TypeError:
                    hd['Unknown'].append(i+1)
        return hd

    def hist_acq_state(self,exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        """
        returns a hist of states each firearm was acquired in
        :param exclude; if True will not tally records that do not meet criteria
        :param start: the first date to tally
        :param stop: the last date to tally
        :return; the state of acquisition histogram
        """
        ha = {}
        for i,key in self._enumerate_(exclude,start,stop):
            for s,n in self.acq_states(key):
                if s in ha: ha[s].extend([i+1]*n)
                else: ha[s] = [i+1]*n
        return ha

    def hist_sou2soa(self,exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        """
        returns a hist of state of use to state state of acquistion for each firearm
        :param exclude; if True will not tally records that do not meet criteria
        :param start: the first date to tally
        :param stop: the last date to tally
        :return; the state of use to state of acquisition histogram
        """
        hs = {'Unknown':[],'In-state':[],'Out-of-state':[]}
        for i,key in self._enumerate_(exclude,start,stop):
            sou = self._rs[key]['location']['state']
            for s,n in self.acq_states(key):
                if s == 'XX': s = 'Unknown'
                elif s != sou: s = 'Out-of-state'
                else: s = 'In-state'
                hs[s].extend([i+1]*n)
        return hs

    def hist_outofstate(self,exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        """
        returns a hist of out of state firearm used i.e the incident involved at
        least one out of state gun
        :param exclude; if True will not tally records that do not meet criteria
        :param start: the first date to tally
        :param stop: the last date to tally
        :return; the out of state firearm use histogram
        """
        hs = {'Unknown':[],'Yes':[],'No':[]}
        for i,key in self._enumerate_(exclude,start,stop):
            soa = self.acq_states(key)
            if len(soa) == 1 and soa[0][0] == 'XX': hs['Unknown'].append(i+1)
            else:
                sou = self._rs[key]['location']['state']
                found = False
                for s,_ in soa:
                    if s != 'XX' and s != sou:
                        hs['Yes'].append(i+1)
                        found = True
                        break
                if not found: hs['No'].append(i+1)
        return hs

    def hist_exported(self,exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        """
        returns a hist of states that exported firearms in ms
        at least one out of state gun
        :param exclude; if True will not tally records that do not meet criteria
        :param start: the first date to tally
        :param stop: the last date to tally
        :return; the exported state histogram
        """
        he = {}
        for i,key in self._enumerate_(exclude,start,stop):
            sou = self._rs[key]['location']['state']
            for s,n in self.acq_states(key):
                if s != 'XX' and s != sou:
                    if s in he: he[s].extend([i+1]*n)
                    else: he[s] = [i+1]*n
        return he

#### ADVANCED FUNCTIONS

    def pivot(self,fct1,fct2,**kwargs):
        """
        returns a pivot from fct1 to fct 2
        :param fct1: (string) the function to pivot on
        :param fct2: (string) the function to pivot to
        :param kwargs: list of args to pass to the functions
                       args to pass to fct1 will be preceded by 'f1_'
                       args to pass to fct2 will be preceded by 'f2_'
                       args to pass to both functions will be unprefixed
        :return: the pivot
        Example: to pivot from race to motives use
        pivot('hist_demographics',
              'hist_motives',
              f1_by='race',
              exclude=msdb.CRIT_STRICT_SCENE,
              start=msdb.FAWB_START,
              stop=msdb.FAWB_STOP)
        """
        # datetime object is janky passed through eval, do work around here
        start = DEF_START
        stop = DEF_STOP

        _=self # get rid of pycharm warning
        # evaluate the functions
        try:
            fct1 = eval('self.'+fct1)
            fct2 = eval('self.'+fct2)
        except AttributeError:
            raise

        # compile the arguments and execute the functions
        args1 = {}
        args2 = {}
        for k,v in kwargs.items():
            if k.startswith('f1_'): args1[k.replace('f1_','')] = v
            elif k.startswith('f2_'): args2[k.replace('f2_','')] = v
            else:
                if k == 'start': start = v
                elif k == 'stop': stop = v
                else:
                    args1[k] = v
                    args2[k] = v
        h1 = fct1(start=start,stop=stop,**args1)
        h2 = fct2(start=start,stop=stop,**args2)

        # create the pivot from fct1 to fct2
        pvt = {i:{} for i in h1}
        for i in h1:
            for idx in h1[i]:
                for j in h2:
                    if idx in h2[j]:
                        if j in pvt[i]: pvt[i][j].append(idx)
                        else: pvt[i][j] = [idx]
        return pvt

    def mpivot(self,*args,**kwargs):
        """
        returns a pivot from from fct1 to fctn as defined in kwargs
        :param args: list of hist functions names to pivot on
        :param kwargs: list of args to pass to the functions
                       args to pass to fct1 will be preceded by 'f1_'
                       args to pass to fct2 will be preceded by 'f2_'
                       args to pass to fctn will be preceded by 'fn_'
                       args to pass to both functions will be unprefixed
        :return: ???
        """
        _ = self  # get rid of pycharm warning
        re_fct_pre = re.compile(r"^f(\d+)_(.+)$")
        # evaluate the functions and execute the functions
        fcts = []
        for fct in args: fcts.append(eval('self.'+fct))

        # compile the arguments
        args = [{} for _ in fcts]  # individual fct args
        for k,v in kwargs.items():
            try:
                f,a = re_fct_pre.match(k).groups()
                args[int(f) - 1][a] = v
            except AttributeError:
                for arg in args: arg[k] = v
        hs = [f(**args[i]) for i,f in enumerate(fcts)]

        # create the pvt from the first function to the others
        h1,hs = hs[0],hs[1:]
        pvt = {i:{} for i in h1}
        for i in h1:
            for _ in h1[i]: pass
            # TODO: not sure how useful this will be,
        return pvt

    def buckets(self,first,last,sz,fct,**kwargs):
        """
        creates buckets from the histogram .
        :param first: (integer) the first bin
        :param last: (integer) the last bin
        :param sz: (integer) the bin size
        :param fct: the histogram function
        :param kwargs: the arguments to pass to fct
        :return: a tuple buckets,extra where
         buckets is list of Buckets
          and
         extra is a dict of the form key->list where key is a non-numeric value
        Note: only works for histogram functions where the key is numeric i.e.
         hist_rnds_fired
        """
        _ = self  # get rid of pycharm warning
        # check function for existence
        try:
            fct = eval('self.' + fct)
        except AttributeError:
            raise

        # check if the bucket size is a float or integer & create bins accordingly
        if type(sz) == type(0.1): bins = [(i,i+sz) for i in np.arange(first,last,sz)]
        else: bins = [(i,i+sz) for i in range(first,last,sz)]
        bins = [(float('-inf'),first)]+bins+[(last,float('inf'))]

        # initialize the buckets
        bkts = [Bucket(b[0],b[1]) for b in bins]
        abkts = {} # non-numeric buckets i.e. 'Unknown'

        # execute and assign to buckets
        hist = fct(**kwargs)
        for k in hist:
            for bkt in bkts:
                try:
                    if bkt.left <= k < bkt.right:
                        bkt.data.extend(hist[k])
                        break
                except TypeError:
                    if k in abkts: abkts[k].extend(hist[k])
                    else: abkts[k] = hist[k]
                    break

        # return the buckets and left over buckets
        return bkts,abkts

    def epochs(self,start,stop,delta,fct,**kwargs):
        """
        creates buckets of histogram function fct based on time periods, executing
        given function over those periods with the given keyword arguments
        :param start: first bucket start
        :param stop: last buck stop
        :param delta: the time period
         must be a dict of the form
           {'time-period':tp,
            'amount':n}
             where tp is one of {year,month} and amount is a whole number
            TODO: may add days but for now it is too 'specific'
        :param fct: (string) the function to execute for each bucket
        :param kwargs: keyword arguments to pass to fct
        :return: a list of dicts of the form {period-start,period-stop,data}
        NOTE: fct must have arguments start and stop
        """
        _ = self  # get rid of pycharm warning
        # check function for existence & evaluate it
        try:
            fct = eval('self.'+fct)
        except AttributeError:
            raise

        # check delta for validity
        try:
            tp,n = delta['time-period'],delta['amount']
            if tp not in ['year', 'month']:
                raise MSDBException(
                    "Invalid Arg 'delta.time-period' ({})".format(tp)
                )
            elif not is_int(n) or n <= 0:
                raise MSDBException("Invalid Arg 'delta.amount' ({})".format(n))
            elif tp == 'month' and n > 11:
                raise MSDBException(
                    "Invalid Arg 'delta.amount' ({}) for tp=month".format(n)
                )
            elif tp == 'day' and n > 30:
                raise MSDBException(
                    "Invalid Arg 'delta.amount' ({}) for tp=day".format(n)
                )
        except KeyError:
            raise MSDBException("Invalid Arg 'delta' ({})".format(delta))

        epochs = []
        bp = start
        while bp < stop:
            # calculate the new end period
            y,m,d = bp.year,bp.month,bp.day
            if tp == 'year': y += n
            elif tp == 'month':
                mod,rem = (m+1) // 12,(m+1) % 12
                if mod == 1 and rem == 1:
                    m = 1
                    y += 1
                else: m += 1

            try:
                ep = datetime.date(y,m,d)
                if ep > stop: ep = stop
            except ValueError:
                d -= 1  # we got an invalid date i.e. february 29, so move back 1 day
                ep = datetime.date(y,m,d)

            # execute the function, them set beginning to end
            epochs.append(Epoch(bp,ep,fct(**kwargs,start=bp,stop=ep)))
            bp = ep

        return epochs

    def ms_rate(self,by='2019',sort=False,exclude=CRIT_LOOSE_INC,start=DEF_START,stop=DEF_STOP):
        """
        returns mass shooting rate for each state
        see https://oag.ca.gov/sites/all/files/agweb/pdfs/cjsc/stats/computational_formulas.pdf
        :param by: what year to use
        :param sort: if True returns a list of sorted tuples by ms rate otherwise
         returns a dict of state->rate
        :param exclude; if True will not tally records that do not meet criteria
        :param start: the first date to tally
        :param stop: the last date to tally
        :return: mass shooting rate per state
        """
        # validate by
        col = {
            '1990':'pop-90','2000':'pop-00','2010':'pop-10','2019':'pop-19'
        }
        try:
            by = col[by]
        except KeyError:
            raise MSDBException("Invalid Arg 'by' ({})".format(by))

        # get the by state hist
        h = self.hist_state(exclude,start,stop)

        # return the adjusted count
        if sort:
            return sorted(
                [(k,round(len(h[k])/self._ss[k][by]*100000,4)) for k in h],
                key=itemgetter(1)
            )
        else:
            return {k:round(len(h[k])/self._ss[k][by]*100000,4) for k in h}

#### PRIVATE FUNCTIONS ####

    def _process_(self): pass

    def _enumerate_(self,e,s,f):
        """
        private helper fct for that returns list of tuples (i=index,k=key)
        :param e: if True does not consider records failing criteria
        :param s: start of period to consider
        :param s: stop of period to consider
        :return: True if record should be considered, False otherwise
        """
        return [(i,k) for i,k in enumerate(self._ks) if self._keep_(k,e,s,f)]

    def _keep_(self,k,e,s,f):
        """
        private helper fct for hist record
        :param k: key of record
        :param e: if True does not consider records failing criteria
        :param s: start of period to consider
        :param f: stop of period to consider
        :return: True if record should be considered, False otherwise
        """
        if s <= self._rs[k]['dtg']['date'] <= f: return self.meets_criteria(k,e)
        return False

    def _validate_inc_rec_(self,r):
        """ formats and validates the incident record """
        # append label to key list and set initial values for record
        try:
            key = r[INC_LBL]
            rec = {
                'shooter':r[INC_SHT],
                'shooter-data':r[INC_SDT],
                'location':r[INC_LOC],
                'dtg':r[INC_DTG],
                'killed':r[INC_VIC],
                'leo-killed':r[INC_LEO],
                'notes':r[INC_NTE],
                'report':r[INC_RPT],
                'toa':r[INC_TOA],
                'cke':r[INC_CKE],
                'scene':r[INC_SCN],
                'venue':r[INC_VEN],
                'type':r[INC_TYP],
                'trigger':r[INC_TGR],
                'motive':r[INC_MTV],
                'tgt':r[INC_TGT],
                'ttl-fa':r[INC_ALL],
                'firearms':[]
            }
        except IndexError: # should not get here unless file is malformed
            raise MSDBException("Missing fields in DB")

        # format as necessary, validating in process
        # parse out shooter(s)
        shooter = rec['shooter'].split(' & ')
        rec['shooter'] = []
        for s in shooter:
            try:
                last,rem = s.split(', ')
            except ValueError:
                raise MSDBException(
                    "Malformed 'name' ({}) in record {}".format(s,key)
                )
            if not rem.endswith('.'):
                first = rem
                m = None
            else:
                first = rem[:-3]
                m = rem[-2]
            rec['shooter'].append({'last':last,'first':first,'middle':m})

        # make a list for shooter data, parse out and validate
        sd = None
        try:
            # split and validate
            sd = rec['shooter-data'].split(' & ')
            if len(sd) != len(rec['shooter']):
                raise MSDBException(
                    "Malformed 'shooter data' ({}) in record {}".format(
                        key,rec['shooter-data']
                    )
                )

            rec['shooter-data'] = []
            for d in sd:
                a,g,r,l,e,m,c,h,s = re_shooter_data.match(d).groups()
                rec['shooter-data'].append(
                    {'age':int(a),
                     'gender':g,
                     'race':r,
                     'religion':l,
                     'ed-lvl':e,
                     'military':m,
                     'leo-knowledge':c,
                     'mental-health':h,
                     'substance-use':s}
                )
        except AttributeError:
            raise MSDBException(
                "Malformed 'shooter data' ({}) in record {}".format(sd,key)
            )

        # location should be city, state
        try:
            city,state = rec['location'].split(', ')
            if not state in self._ss:
                raise MSDBException(
                    "Invalid 'state' ({}) in record {}".format(state,key)
                )
            rec['location'] = {'city':city,'state':state}
        except ValueError:
            raise MSDBException(
                "Invalid 'location' ({}) in record {})".format(rec['location'],key)
            )

        # DTG should be mm/dd/yyyy and optional HHMM preceded by a space
        try:
            dtg = rec['dtg'].split(' ') # split into date and time

            # check date first, is required
            m,d,y = dtg[0].split('/')
            rec['dtg'] = {'date':datetime.date(int(y),int(m),int(d)),'time':None}

            # then time if present
            if len(dtg) == 2:
                hr,mi = dtg[1][0:2],dtg[1][2:]
                # noinspection PyTypeChecker
                rec['dtg']['time'] = datetime.time(int(hr),int(mi))
        except (ValueError,IndexError):
            raise MSDBException(
                "Invalid 'DTG'({}) in record {}".format(rec['dtg'],key)
            )

        # for both victims and leos, split the tallies and set appropriate dict entries
        try:
            k,w,t = rec['killed'].split('/')
            rec['killed'] = int(k)
            rec['wounded'] = int(w)
            rec['total'] = rec['killed'] + rec['wounded']
        except ValueError:
            raise MSDBException(
                "Invalid 'victims' ({}) in record {}".format(rec['killed'],key)
            )

        try:
            k,w,t = rec['leo-killed'].split('/')
            rec['leo-killed'] = int(k)
            rec['leo-wounded'] = int(w)
            rec['leo-total'] = rec['leo-killed'] + rec['leo-wounded']
        except ValueError:
            raise MSDBException(
                "Invalid 'leo-victims' ({}) in record {}".format(
                    rec['leo-killed'],key
                )
            )

        # strip leading/trailing whitespace and remove newlines from notes
        rec['notes'] = rec['notes'].strip()
        rec['notes'] = rec['notes'].replace('\n','')

        # report should be 1, 2 or 3
        try:
            rec['report'] = int(rec['report'])
            if not rec['report'] in [0,1,2]: raise ValueError
        except ValueError:
            raise MSDBException(
                "Invalid 'report' ({}) in record {}".format(rec['report'],key)
            )

        # cke time should be float
        try:
            rec['toa'] = float(rec['toa']) if rec['toa'] else None
        except ValueError:
            raise MSDBException(
                "Invalid 'CKE time' ({}) in record {}".format(rec['toa'],key)
            )

        # for cke cause, we may have two or more depending on # of shooters
        cke = rec['cke'].split(' & ')
        if len(cke) != len(rec['shooter']):
            raise MSDBException(
                "Invalid 'CKE Cause' length ({}) in rcord {}".format(
                    len(cke), key
                )
            )

        rec['cke'] = []
        for c in cke:
            try:
                c,p = c.split(':')
                rec['cke'].append({'code':c,'party':p})
            except ValueError:
                raise MSDBException(
                    "Invalaid 'CKE cause'({}) in record {}".format(cke,key)
                )

        # scene should be a comma separated list with each element beginning with
        # inside or outside
        sc = rec['scene'].lower().split(', ')
        rec['scene'] = []
        for s in sc:
            if s.startswith('inside') or s.startswith('outside'):
                rec['scene'].append(s)
            else:
                raise MSDBException(
                    "Invalid 'scene' ({}) in record {}".format(s,key)
                )

        # TODO: validate venue

        # type must be in defined codes
        if rec['type']:
            if len(rec['type']) != 2:
                raise MSDBException(
                    "Invalid 'shooting type' ({}) in record {}".format(
                        rec['type'],key
                    )
                )
            elif not rec['type'][0] in shooting_type_pos:
                raise MSDBException(
                    "Invalid 'shooting type.position' ({}) in record {}".format(
                        rec['type'],key
                    )
                )
            elif not rec['type'][1] in shooting_type_act:
                raise MSDBException(
                    "Invalid 'shooting type.action' ({}) in record {}".format(
                        rec['type'],key
                    )
                )
            rec['type'] = {'pos':rec['type'][0],'act':rec['type'][1]}
        else: rec['type'] = None

        # trigger must be as defined
        if not rec['trigger'] in tgr_codes:
            raise MSDBException(
                "Invalid 'init-trigger'{}) in record {}".format(rec['trigger'],key)
            )

        # ensure motive codes are valid
        mcs = rec['motive'].split(', ')
        if mcs == ['']: rec['motive'] = []
        else:
            for mc in mcs:
                if mc not in mtv_codes:
                    raise MSDBException(
                        "Invalid 'motive code' ({}) in record {}".format(mc,key)
                    )
            rec['motive'] = mcs

        # ensure target codes are valid
        try:
            p,l = rec['tgt']
            if p not in tgt_c2l_tgt:
                raise MSDBException(
                    "Invalid 'target code' ({}) in record {}. "
                    "'People' ({}) not recognized.".format(rec['tgt'],key,p)
                )
            if l not in tgt_c2l_loc:
                raise MSDBException(
                    "Invalid 'target code' ({}) in record {}. "
                    "'Location' ({}) not recognized.".format(rec['tgt'],key,l)
                )
        except ValueError: # target is not two-character
            raise MSDBException(
                "Invalid 'target code' ({}) in record {}".format(rec['tgt'],key)
            )

        # all/ttl firearms will be an integer (greater than 0)
        try:
            rec['ttl-fa'] = int(rec['ttl-fa'])
            if rec['ttl-fa'] < 1: raise ValueError
        except ValueError:
            raise MSDBException(
                        "Invalid 'TTL FA' ({}) in record {}".format(
                            key,rec['ttl-fa']
                        )
                    )

        # all went well add key to the key list and the inc. record to records
        self._ks.append(key)
        self._rs[key] = rec

    def _validate_firearm_rec_(self,c,f,r):
        """
         creates a Firearm Object from row
         :param c: Caliber Object
         :param f: FirearmSpec Object
         :param r: the row of date to read
        """
        key = self._ks[-1] # get the last appended key in the key list
        try:
            self._rs[key]['firearms'].append(
                firearm.Firearm(
                    c,
                    f,
                    r[FA_MDL],
                    r[FA_TYP],
                    r[FA_MAN],
                    r[FA_FAC],
                    r[FA_CAL],
                    r[FA_CAP],
                    r[FA_MAG],
                    r[FA_ATH].lower(),
                    r[FA_RND],
                    r[FA_NTE],
                    r[FA_STD],
                    r[FA_MCH],
                    r[FA_LEG],
                    r[FA_ACQ],
                    r[FA_AST],
                )
            )
            # TODO: have to check that acq-date is not after inc date
        except IndexError:
            raise MSDBException("Invalid 'firearm data' in record {}".format(key))
        except firearm.FirearmException as e:
            raise MSDBException(
                "Firearm creation failed: {} in reocrd {}".format(e.args[0],key)
            )

    @staticmethod
    def _write_shooter_(ss):
        out = ""
        for s in ss:
            if out: out += " & "
            out += "{}, {}".format(s['last'],s['first'])
            if s['middle']: out += " {}.".format(s['middle'])
        return out

    @staticmethod
    def _write_shooter_data_(sd):
        out = ""
        for s in sd:
            if out: out += " & "
            out += "{}/{}/{}/{}/{}/{}/{}/{}/{}".format(
                s['age'],
                s['gender'],
                s['race'],
                s['religion'],
                s['ed-lvl'],
                s['military'],
                s['leo-knowledge'],
                s['mental-health'],
                s['substance-use']
            )
        return out

    @staticmethod
    def _write_dtg_(dtg):
        out = dtg['date'].strftime('%m/%d/%Y')
        if dtg['time']: out += " " + dtg['time'].strftime('%H%M')
        return out

    @staticmethod
    def _write_cke_cause_(cs):
        return " & ".join(["{}:{}".format(c['code'],c['party']) for c in cs])
