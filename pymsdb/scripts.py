#!/usr/bin/env python
""" scripts.py
Copyright (C) 2020 Dale V. Patterson (dale.v.patterson@gmail.com)

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

Defines funtions outside of msdb scope to work with ms
"""

#__name__ = 'scripts'
__license__ = 'GPLv3'
__version__ = '0.0.2'
__date__ = 'April 2021'
__author__ = 'Dale Patterson'
__maintainer__ = 'Dale Patterson'
__email__ = 'dale.v.patterson@gmail.com'
__status__ = 'Development'

import datetime
import numpy as np
from scipy.stats import pearsonr,spearmanr,kendalltau
import pymsdb.msdb as msdb
import pymsdb.firearm as firearm
import pymsdb.draw as draw
import pymsdb.utils as utils
import pymsdb.taxonomy as tax
from pymsdb import PyMSDbException

# decade data from grant duwe
duwe = [
    (1900,0),(1910,2),(1920,2),(1930,9),(1940,8),(1950,1),(1960,6),(1970,13),
    (1980,32),(1990,42),(2000,28),(2010,14)
]
def plot_duwe():
    # TODO: fix axis especially y-axis
    #draw.plot_ver_bar(duwe)
    draw.plot_time_series(duwe)

# values for weapon specification
wpn_spec = [
    'max-cap',         # maximum magazine capacity
    'muzzle-energy',   # bullet muzzle energy
    'muzzle-velocity', # bullet muzzle velocity
]

# values for victim type
vic_spec = [
    'killed',   # victims killed
    'wounded',  # victims wounded
    'total',    # total victims (killed/wounded) due to shooter gunfire
    'fat-rate', # fatality rate (killed / total shot)
]

####
## weapon to fatalities and correlations
####

def firearm2vic(ms,known=False,exclude=0,start=None,stop=None):
    """
    calculates victims by weapon for each based on:
     1. weapon is only weapon
     2. individual weapons have annotated victims
     3. weapon is identified as primary
    events are not included if one of the three above is not met
    :param ms: the msdb object
    :param known: if set only returns weapons from 1 and 2
    :param exclude: msdb criteria to use
    :param start: the start date
    :param stop: the stop date
    :return: a list of tuples t = (key,weapon id,killed,wounded)
    """
    # set default start/stop dates if not set
    if not start: start = datetime.date(1949,1,1)
    if not stop: stop = datetime.date.today() + datetime.timedelta(1)

    # start first with known weapon to fatalities (either the lone weapon used or
    # individual weapons in an event with known fatalities, If not found, try
    # primary weapon for that event
    fs = []
    for key in ms.keys():
        # 1. only 1 weapon used?
        if ms.firearms_used(key) == 1:
            fs.append((key,0,ms[key]['killed'],ms[key]['wounded']))
            continue

        # 2. ind. weapons have known fatalities
        ind = False
        for fidx,farm in ms.firearms(key):
            k,w = farm.fatalities
            if k is None: break
            else:
                ind = True
                fs.append((key,fidx,k,w))
        if ind: continue

        # 3. finally, check for primary weapon
        if not known:
            try:
                pw = ms.primary_firearm(key)
                if pw is not None:
                    fs.append((key,pw,ms[key]['killed'],ms[key]['wounded']))
            except TypeError:
                print(key)
                raise
    return fs

def wpnspec2vic(ms,spec='max-cap',vtype='killed',known=False,exclude=0,start=None,stop=None):
    """
    finds records where both a weapon can be tied to a specification (i.e. magazine cap
     or caliber) and to fatalities
    :param ms: the msdb object
    :param spec: as defined in wpn_spec
    :param vtype: as defined in vic_spec
    :param known: if set
    :param exclude: msdb criteria to use
    :param start: the start date
    :param stop: the stop date
    :return: a list of tuples = (case id, wpn id, spec, vic value)

    NOTE: this will return different results based on the spcification.
     if max-cap is set and capacities are not known, it won't be included
     if me or mv is set and and caliber is not known, it won't be included
    """
    #TODO: figure out how to add weapon type (need a way to sort weapons)
    #set default start/stop dates if not set
    if not start: start = datetime.date(1949,1,1)
    if not stop: stop = datetime.date.today() + datetime.timedelta(1)

    ts = []
    for key,widx,k,w in firearm2vic(ms,known,exclude,start,stop):
        # get event id, firearm object and set up vtype value
        wpn = ms.get_wpn(key,widx)
        i = ms.case2idx(key)
        val = None
        if vtype == 'killed': val = k
        elif vtype == 'wounded': val = w
        elif vtype == 'total': val = k+w
        elif vtype == 'fat-rate': val = round(k / (k+w),2)
        else: raise PyMSDbException("Invalid Arg 'vtype' ({})".format(vtype))

        if spec == 'max-cap':
            mc = wpn.max_capacity()
            if mc: ts.append((i,widx,mc,val))
        elif spec == 'muzzle-energy' or spec == 'muzzle-velocity':
            # with caliber related, have to make sure there is a caliber and the
            # caliber has a known muzzle-energy or muzzle-velocity
            wcal = wpn.caliber
            if wcal:
                try:
                    sval = None
                    if spec == 'muzzle-velocity': sval = wpn.blt_muzzle_velocity
                    else: sval = wpn.blt_muzzle_energy
                    if sval: ts.append((i,widx,sval,val))
                except KeyError: # cal cannot be found
                    pass
        else: raise PyMSDbException("Invalid Arg 'spec' ({})".format(spec))
    return ts

def wpnspec2vic_correlation(ts,fct='pearson'):
    """
    calculates corrleation of capacity to killed using various functions
    :param ts: data list of tuples t = (case id,wpn id, spec value ,vic value)
    :param fct: oneof
     pearson (both variables x and y should be normally distributed
     spearman
     kendall
    :return: correlation coefficient, p-value
    """
    r = p = None
    xs = [x for _,_,x,_ in ts]
    ys = [y for _,_,_,y in ts]
    if fct == 'pearson': r,p = pearsonr(xs,ys)
    elif fct == 'spearman': r,p = spearmanr(xs,ys)
    elif fct == 'kendall': r,p = kendalltau(xs,ys)
    else: raise PyMSDbException("Invalid Arg 'fct' ({})".format(fct))
    return r,p

def draw_wpnspec2vic(ms,spec='max-cap',vtype='killed',known=False,lobf=False,exclude=0,start=None,stop=None):
    """
    draw a scatter plot of weapon spec to victim
    :param ms: the msdb object
    :param spec: specification to correlate can be one of
      'max-cap': the maximum capacity of the weapon
      'me': muzzle energy
      'mv': muzzle velocity
    :param vtype: oneof 'killed','wounded','total','fr' = fatality rate
    :param known: checks only weapons with known fataliities
    :param lobf: draw a line of best fit
    :param exclude: msdb criteria to use
    :param start: the start date
    :param stop: the stop date
    """
    # get the spec 2 vic and extract vic count
    pts = [(s,v) for _,_,s,v in wpnspec2vic(ms,spec,vtype,known,exclude,start,stop)]
    ss = [s for s,_ in pts] # spec values
    vs = [v for _,v in pts] # vtype values

    # setup axis tick marks and lables
    xts = []
    yts = []
    if spec == 'max-cap':
        xts = [(x,str(x)) for x in range(5,max(ss)+5,5)]
    elif spec == 'muzzle-energy' or spec == 'muzzle-velocity':
        hstep = 200 if spec == 'muzzle-velocity' else 500
        #low = int(min(ss) // hstep) * hstep
        high = (int(round(max(ss),0)) + hstep) // hstep * hstep
        xts = [(x,str(x)) for x in range(0,high+5,hstep)]
    if vtype == 'fat-rate':
        yts = [(y,str(y)) for y in [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]]
    else:
        vstep = 5 if (max(vs) - min(vs)) / len(vs) < 1 else 50
        yts = [(y,str(y)) for y in range(vstep,max(vs)+5,vstep)]

    # setup labels
    ttl = ylbl = xlbl = None

    # xlabel and first part of title is spec
    if spec == 'max-cap':
        xlbl = "Maximum Magazine Capacity"
        ttl = "Magazine Capacity to "
    elif spec == 'muzzle-energy':
        xlbl = 'Muzzle Energy (N-m)'
        ttl = "Muzzle Energy to to "
    elif spec == 'muzzle-velocity':
        xlbl = 'Muzzle Velocity (m/s)'
        ttl = "Muzzle Velocity to "

    # ylabel and last part of title is vtype
    if vtype == 'killed':
        ylbl = "Victims Killed"
        ttl += ylbl
    elif vtype == 'wounded':
        ylbl = "Victims Wounded"
        ttl += ylbl
    elif vtype == 'total':
        ylbl = "Total Victims"
        ttl += ylbl
    elif vtype == 'fat-rate':
        ylbl = 'Fatality Rate'
        ttl += ylbl

    # plot line of best fit
    if lobf: lobf = utils.lobf(pts)
    else: lobf = []

    # TODO: add text box with kendall tau

    # plot the scatterplot
    draw.plot_scatter(pts,ttl,xlbl,ylbl,xts,yts,lobf)

####
## events over time
####

DB_DATE_START = datetime.date(1965,9,14) # skips 1949
DB_DATE_STOP = datetime.date.today()     # stop today

def time_series_incidents(ms,exclude=0,start=None,stop=None,amt=1):
    # skips the early 1949 incident
    if not start: start = datetime.date(1966,1,1)
    if not stop: stop = datetime.date.today()
    es = ms.epochs(
        exclude=exclude,
        start=start,
        stop=stop,
        delta={'time-period':'year','amount':amt},
        fct='incidents'
    )
    return [(e.start.year,e.data if e.data else 0) for e in es]

def time_series_victims(ms,vtype='killed',exclude=0,start=None,stop=None,amt=1):
    # skips the early 1949 incident
    if not start: start = datetime.date(1966,1,1)
    if not stop: stop = datetime.date.today()
    es = ms.epochs(
        exclude=exclude,
        start=start,
        stop=stop,
        delta={'time-period':'year','amount':amt},
        fct='victims',
        column=vtype
    )
    return [(e.start.year,e.data if e.data else 0) for e in es]

def time_series_vicperincident(ms,exclude=0,start=None,stop=None,amt=1):
    if not start: start = datetime.date(1966,1,1)
    if not stop: stop = datetime.date.today()
    esi = time_series_incidents(ms,start,stop,amt,exclude)
    esv = time_series_victims(ms,'killed',start,stop,amt,exclude)
    rs = []
    for i,e in enumerate(esi):
        rs.append((e[0],esv[i][1]/e[1] if e[1] else 0))
    return rs

def draw_time_series(ms,tstype,clbl=None,exclude=0,start=None,stop=None,amt=1,hlfawb=False,trend=0):
    """
    plot a time series from start to none
    :param ms: the mass shooting db
    :param tstype: the type of time series
    :param clbl: overwrite default lables if set should be a list of strings
     [Title,X-Axis,Y-Axis] where any can be None or defined
    :param exclude: msdb criteria to use
    :param start: the start date
    :param stop: the stop date
    :param amt: # of years for each period
    :param hlfawb: highlight the fawb period as gray box
    :param trend: the period for the trend line, none if 0
    """
    if not tstype in ['incidents','killed','wounded','total','kpi']:
        raise PyMSDbException("Invalid Arg 'tstype' ({})".format(tstype))

    # add default start stop dates if not set
    if not start: start = datetime.date(1966,1,1)
    if not stop: stop = datetime.date.today()

    # get data and set labels
    ts = ttl = xlbl = ylbl = hl = es = None
    if clbl:
        if clbl[0]: ttl = clbl[0]
        if clbl[1]: xlbl = clbl[1]
        if clbl[2]: ylbl = clbl[2]

    if not xlbl: xlbl = "Year"
    if tstype == 'incidents':
        ts = time_series_incidents(ms,start,stop,amt,exclude)
        if not ttl: ttl = "Mass Shooting Events"
        if not ylbl: ylbl = "Number of Events"
        yts = [(x,str(x)) for x in range(0,max([y for _,y in ts])+amt)]
    elif tstype == 'kpi':
        ts = time_series_vicperincident(ms,start,stop,amt,exclude)
        if not ttl: ttl = "Mass Shooting Fatalities per Event"
        if not ylbl: ylbl = "Fatalities per Event"
        yts = [(x,str(x)) for x in range(0,int(round(max([r for _,r in ts]),0)),amt)]
    else:
        ts = time_series_victims(ms,tstype,start,stop,amt,exclude)
        if not ttl: ttl = "Mass Shooting Victims"
        if not ylbl:
            if tstype == 'killed': ylbl = "Number of Fatalities"
            elif tstype == 'wounded': ylbl = "Number of Injured"
            else: ylbl = "Total Number of Victims"
        yts = [(x,str(x)) for x in range(0,max([y for _,y in ts])+10,10)]

    # set up axis tick marks and values
    xts = [(x,str(x)) for x,_ in ts]

    # setup highlight fawb box if set - Unlike bar charts, we are using the x-axis
    #  values, namely the year
    if hlfawb: hl = [1994,2004,0,1,'FAWB']

    # get trend line if set
    if trend == 0: trend = []
    else: es = utils.ema([v for _,v in ts],trend)

    draw.plot_time_series(ts,ttl,xlbl,ylbl,xts,yts,hl,es)

def draw_snapshot(ms,tstype,cttl='',exclude=0,start=None,stop=None,amt=1,hlfawb=False):
    # add default start stop dates if not set
    if not start: start = datetime.date(1966, 1, 1)
    if not stop: stop = datetime.date.today()

    ts = ttl = ylbl = hl = None
    xlbl = "Period"
    if tstype == 'incidents':
        ts = time_series_incidents(ms,start,stop,amt,exclude)
        if cttl is not None: ttl = "Mass Shooting Events" if cttl == '' else cttl
        ylbl = "Total Events"
    elif tstype == 'kpi':
        ts = time_series_vicperincident(ms,start,stop,amt,exclude)
        if cttl is not None: ttl = "Mass Shooting Fatalities" if cttl == '' else cttl
        ylbl = "Fatalities per Event"

    # rewrite the ts list as tuples t = ("time range",value) before drawing
    # have to check if the last period was full period
    if amt > 1: ps = [("{}-{}".format(t,t+amt),v) for t,v in ts]
    else: ps = [(str(t),v) for t,v in ts]
    if int(ts[-1][0]) + amt > datetime.date.today().year:
        ps[-1] = ("{} - {}".format(ts[-1][0],datetime.date.today().year),ps[-1][1])

    # find the index of the fawb start and subtract .5, do the same for fawb start
    if hlfawb:
        b = e = 0
        for i,t in enumerate(ts):
            if t[0] == 1994: b = i - 0.5
            elif t[0] == 2004: e = i -0.5
        hl = [b,e,0,1,'FAWB']

    # and draw it
    draw.plot_ver_bar(ps,ttl,xlbl,ylbl,hl)

####
## SCP / event frequency analysis
####
# TODO: standard deviation is too large to be useful
def scp_parameters(ms,l,start=None,stop=None,exclude=0):
    """
    caluclates scp parameters from ms on elapsed time
    :param ms: the mass shooting db
    :param l: lambda
    :param start: the start date
    :param stop: the stop date
    :param exclude: msdb criteria to use
    returns elapsed times,mean,std,var,range
    """
    # get elapsed times for specified time period and criteria, extract times
    if start is None: start = msdb.DEF_START
    if stop is None: stop = msdb.DEF_STOP
    es = ms.chrono_elapsed(exclude,start,stop)
    ts = np.array([t for _,t in es if not np.isnan(t)])

    ewma = utils.ewma(ts,l)
    m = np.average(ewma)
    std = np.std(ewma)
    ucl = m + 3 * std * np.sqrt(l/(2-l))
    lcl = m - 3 * std * np.sqrt(l/(2-l))
    return ewma,m,std,ucl,lcl

    # calculate parameters
    #m = np.average(ts)
    #s = np.std(ts)
    #v = np.var(ts)
    #r = max(ts) - min(ts)

    #return es,m,s,v,r

####
## FREQUENCY OF SHOOTINGS
####

# A CUMULATIVE MEAN
def draw_freq2time_series(ms,clbl=None,exclude=0,start=None,stop=None,amt=1,hlfawb=False):
    """
    plot a time series from start to none
    :param ms: the mass shooting db
    :param clbl: overwrite default lables if set should be a list of strings
     [Title,X-Axis,Y-Axis] where any can be None or defined
    :param exclude: msdb criteria to use
    :param start: the start date
    :param stop: the stop date
    :param amt: # of years for each period
    :param hlfawb: highlight the fawb period as gray box
    """
    # add default start stop dates if not set
    if not start: start = datetime.date(1966,1,1)
    if not stop: stop = datetime.date.today()

    # set labels
    ts = ttl = xlbl = ylbl = hl = es = None
    if clbl:
        if clbl[0]: ttl = clbl[0]
        if clbl[1]: xlbl = clbl[1]
        if clbl[2]: ylbl = clbl[2]

    # get elapsed times (frequency)
    #ts = time_series_incidents(ms,exclude,start,stop,1)
    fs = ms.chrono_elapsed(exclude,start,stop)

    # calculate the cumulative moving average
    cma = [0]*len(fs)
    for i,(_,et) in enumerate(fs):
        if i == 0: cma[i] = 0
        else: cma[i] = (et + i * cma[i-1]) / (i+1)
    es = [(i,x[1]) for i,x in enumerate(fs)]
    cma = [(i,x) for i,x in enumerate(cma)]

    # configure the x-axis ticks to show the first every 5 years
    years = [str(x) for x in range(1970,2025,5)]
    xts = []
    for i,y in [
        (i,str(ms[ms.idx2case(idx)]['date'].year)) for i,(idx,_) in enumerate(fs)
    ]:
        if not years: break
        elif y < years[0]: continue
        elif y == years[0]:
            xts.append((i,y))
            years.pop(0)

    # draw it
    draw.plot_time_series_with_freq(es,cma,xticks=xts)

    return es,cma

#def draw_comparison(ms):
#    #### DEFUNCT ####
#    # get the dates for each db and tally for years
#    tvp,stf,mj,ts = tally.get_dates(ms)
#    tvp_n,stf_n,mj_n,ts_n = tally.tally(tvp,stf,mj,ts)
#
#    hl = ['1994','2004',0,1,'FAWB']
#
#    # set up axis tick marks and values
#    xts = [(x,str(x)) for x,_ in ts_n]
#    yts = [(x,str(x)) for x in range(0,max([y for _,y in stf_n])+1)]
#
#    draw.plot_db_comparison(
#        tvp_n,stf_n,mj_n,ts_n,
#        'MS Databases: Comparison','Year','Events per Year',xts,yts,
#        hl
#    )

# test
# find weapon category associated with killed/injured
# looking at level 3 and level 4 only

TYPE_CAT_2   = 2 # Level 2 only handgun and long gun
TYPE_CAT_3   = 3 # Level 3 only shotgun,rifle,revolver,derringer,pistol
TYPE_CAT_3_4 = 4 # Level 4 where possible (i.e. carbine, tact-wpn), level 3 otherwise
def draw_lethality_by_type(ms,lvl=TYPE_CAT_3,exclude=0,start=None,stop=None):
    # add default start stop dates if not set, set up lvl for wpn_type_lethality
    if not start: start = datetime.date(1949,1,1)
    if not stop: stop = datetime.date.today()

    tlvl = 2 if lvl == TYPE_CAT_2 else 3

    # get lethality by types
    ds,nds = wpn_type_lethality(ms,tlvl,exclude,start,stop)

    # cast to appropriate level
    tds = None
    if lvl == TYPE_CAT_2:
        tds = [(key,idx,tax.firearm_lvl_2(t),k,w) for key,idx,t,k,w in ds]
    elif lvl == TYPE_CAT_3:
        tds = [(key,idx,tax.firearm_lvl_3(t),k,w) for key,idx,t,k,w in ds]
    else:
        tds = [(key,idx,_wpn_lvl_3_4_(t),k,w) for key,idx,t,k,w in ds]

    ts = {}
    for key,_,t,k,w in tds:
        if t in ts: ts[t].append((key,k,w,k+w))
        else: ts[t] = [(key,k,w,k+w)]
    return ts
#for t in ts:
#    ks = [k for k,_,_ in ts[t]]
#    print(t,len(ts[t]),np.average(ks),np.std(ks),np.var(ks),max(ks)-min(ks))

def _wpn_lvl_3_4_(wtype):
    # try level 4 first, otherwise level 3
    wlvl = tax.firearm_lvl_4(wtype)
    try:
        if wlvl.startswith('Tactical'): wlvl = 'Assault Weapon'
    except AttributeError:
        wlvl = tax.firearm_lvl_3(wtype)
    return wlvl

def wpn_type_lethality(ms,lvl=3,exclude=0,start=None,stop=None):
    # add default start stop dates if not set
    if not start: start = datetime.date(1949,1,1)
    if not stop: stop = datetime.date.today()

    ds = []  # determinable (can put weapon type to killed, injured
    nds = [] # non-determinable (cannot put weapon type to killed, injured

    for key in ms.keys():
        if ms.wpns_used(key) == 1:
            # for single weapon incidents, check if the category is known at the
            #  specified category level
            wtype = ms[key]['weapons'][0]['type']
            if tax.firearm_cat_lvl(wtype) >= lvl:
                ds.append((key,0,wtype,ms[key]['killed'],ms[key]['wounded']))
            else: nds.append(key)
        else:
            # for multiple weapon incidents determine if each weapon has k/w listed
            found = False
            for i,wpn in enumerate(ms[key]['weapons']):
                try:
                    k,w = firearm.re_firearm_vic.search(wpn['notes']).groups()
                    wtype = wpn['type']
                    if tax.firearm_cat_lvl(wtype) >= lvl:
                        ds.append((key,i,wtype,int(k),int(w) if w else 0))
                        found = True
                except AttributeError:
                    pass

            # or if the weapons are of the same category at the specified category level
            if not found:
                lcd = tax.lcd([wpn['type'] for wpn in ms[key]['weapons']])
                if tax.firearm_cat_lvl(lcd) >= lvl:
                    ds.append((key,'ALL',lcd,ms[key]['killed'],ms[key]['wounded']))
                    found = True

            if not found: nds.append(key)
    return ds,nds