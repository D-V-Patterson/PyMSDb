#!/usr/bin/env python
""" states.py
Copyright (C) 2020 Dale V. Patterson (dale.v.patterson@gmail.com)

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

Defines funtions to read/write the State population db
"""

#__name__ = 'states'
__license__ = 'GPLv3'
__version__ = '0.0.2'
__date__ = 'March 2021'
__author__ = 'Dale Patterson'
__maintainer__ = 'Dale Patterson'
__email__ = 'dale.v.patterson@gmail.com'
__status__ = 'Development'

import csv
from operator import itemgetter
from pymsdb import PyMSDbException
from pymsdb.data import pth_states

ST_NAME =  0
ST_ABBR =  1
ST_GLS  =  2
ST_P90  =  3
ST_PI90 =  4
ST_P00  =  5
ST_PI00 =  6
ST_P10  =  7
ST_PI10 =  8
ST_P19  =  9
ST_PI19 = 10
ST_DLTA = 11
ST_DENS = 12
ST_PIAV = 13

fpath = 'data/states.tsv'

def read():
    """ reads in state pop db and returns a dcit """
    ss = {}
    fin = None
    try:
        fin = open(pth_states)
        csvr = csv.reader(fin, delimiter='\t')
        _ = csvr.__next__()
        for row in csvr:
            ss[row[ST_ABBR]] = {
                'name':row[ST_NAME],
                'gls':int(row[ST_GLS]),
                'pop-90':int(row[ST_P90]),
                'pi-90':float(row[ST_PI90]),
                'pop-00':int(row[ST_P00]),
                'pi-00':float(row[ST_PI00]),
                'pop-10':int(row[ST_P10]),
                'pi-10':float(row[ST_PI10]),
                'pop-19':int(row[ST_P19]),
                'pi-19':float(row[ST_PI19]),
                'pop-delta':float(row[ST_DLTA]),
                'pop-dens':float(row[ST_DENS]),
                'pi-avg':float(row[ST_PIAV]),
            }
    except IOError:
        raise
    finally:
        if fin: fin.close()
    return ss

def ranked(ss,by='gls'):
    """
    returns a list of keys ranked by specified criteria
    where by is one of {'gls','1990','2000','2010','2019'}
     (gls) = gifford gun law state ranking
    """
    # check by for validity
    if by != 'gls': by = 'pop-' + by
    if not by in ['gls','pop-1990','pop-2000','pop-2010','pop-2019']:
        raise PyMSDbException("Invalid Arg 'by' ({})".format(by))

    # return keys ranked by criteria
    return [x for x,_ in sorted([(k,ss[k][by]) for k in ss],key=itemgetter(1))]

def write(ss):
    """ writes states population to file """
    fout = None
    hdr = "NAME\tABBRV\tPOP 1990\tPOP IDX 1990\tPOP 2000\tPOP IDX 2000\tPOP 2010\t"
    hdr += "POP IDX 2010\tPOP 2019\tPOP IDX 2019\tPOP DELTA\tPOP DENS\tPI AVG\n"
    try:
        fout = open(fpath,'w')
        fout.write(hdr)
        for key in ss:
            fout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                ss[key]['name'],
                key,
                ss[key]['gls'],
                ss[key]['pop-90'],
                ss[key]['pi-90'],
                ss[key]['pop-00'],
                ss[key]['pi-00'],
                ss[key]['pop-10'],
                ss[key]['pi-10'],
                ss[key]['pop-19'],
                ss[key]['pi-19'],
                ss[key]['pop-delta'],
                ss[key]['pop-dens'],
                ss[key]['pi-avg']
            ))
    except IOError:
        pass
    finally:
        if fout: fout.close()

def calc_index(ss):
    """ calculate pop index for each given year """
    t90 = sum([ss[k]['pop-90'] for k in ss])
    t00 = sum([ss[k]['pop-00'] for k in ss])
    t10 = sum([ss[k]['pop-10'] for k in ss])
    t19 = sum([ss[k]['pop-19'] for k in ss])

    for k in ss:
        ss[k]['pi-90'] = round(ss[k]['pop-90']/t90,4)
        ss[k]['pi-00'] = round(ss[k]['pop-90']/t00,4)
        ss[k]['pi-10'] = round(ss[k]['pop-90']/t10,4)
        ss[k]['pi-19'] = round(ss[k]['pop-90']/t19,4)
        ss[k]['pi-avg'] = round(
            sum([ss[k]['pi-90']+ss[k]['pi-00']+ss[k]['pi-10']+ss[k]['pi-19']])/4,4
        )