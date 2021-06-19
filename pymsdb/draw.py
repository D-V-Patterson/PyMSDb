#!/usr/bin/env python
""" draw.py
Copyright (C) 2020 Dale V. Patterson (dale.v.patterson@gmail.com)

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

Defines funtions to plot data
"""

#__name__ = 'draw'
__license__ = 'GPLv3'
__version__ = '0.0.1'
__date__ = 'March 2021'
__author__ = 'Dale Patterson'
__maintainer__ = 'Dale Patterson'
__email__ = 'dale.v.patterson@gmail.com'
__status__ = 'Development'

import numpy as np
import matplotlib.pyplot as plt
import pymsdb.ballistics as bls
#from matplotlib.patches import RegularPolygon
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes
#from mpl_toolkits.axes_grid1 import make_axes_locatable
#from matplotlib.projections import get_projection_class
#from mpl_toolkits.mplot3d import Axes3D
#from mpl_toolkits.mplot3d import proj3d
#import matplotlib.cm as cm
#import matplotlib.colors as colors
#import pymsdb.ballistics.force as force

def plot_hor_bar(ts):
    """ plots ts = [(x,y),(x,y)....] as hor bar chart """
    ys = range(len(ts))
    xs = [x for _,x in ts]
    ls = [l for l,_ in ts]

    fig,ax = plt.subplots()
    plt.tick_params(axis='x',which='both',bottom=False,top=False)
    plt.tick_params(axis='y',which='both',right=False,left=False)
    ax.barh(ys,xs,1.,align='center')
    ax.set_yticks(ys)
    ax.set_yticklabels(ls, verticalalignment='center')
    ax.axis('tight')
    #plt.tight_layout()
    plt.show()

def plot_ver_bar(ts,ttl='',xlbl='',ylbl='',hl=None):
    """
    :param ts: the data ts = [(x,l),(x,l)....] as hor bar chart where x is the data value
     and l is the label
    :param ttl: title if any of the chart
    :param xlbl: label if any of the x-axis
    :param ylbl: label if any of the y-axis
    :param hl: coordiantes if any of a highlight box to draw
    """
    xs = range(len(ts))
    ys = [v for _,v in ts]
    ls = [l for l,_ in ts]
    xangle = 90 if len(xs) > 10 else 0 # set x-axis label orientation

    fig,ax = plt.subplots()
    ax.bar(xs,ys,width=1.,align='center')
    if hl:
        ax.axvspan(
            hl[0],hl[1],hl[2],hl[3],color='gray',alpha=0.5,label=hl[4] if len(hl)==5 else ''
        )
    # TODO: add text label at top of span

    # format graph
    plt.tick_params(axis='x',which='both',bottom=False,top=False)
    plt.tick_params(axis='y',which='both',right=False,left=False)
    if ttl: plt.title(ttl)
    if xlbl: plt.xlabel(xlbl)
    if ylbl: plt.ylabel(ylbl)
    ax.set_yticks(ys)
    ax.set_xticks(xs)
    ax.set_xticklabels(ls,rotation=xangle)
    if hl and len(hl) == 5: ax.legend(loc='upper left')
    ax.axis('tight')
    plt.tight_layout()

    # draw it
    plt.show()

def plot_histogram(xs,n=5):
    """ plots the data in xs as a histogram of n bins """
    plt.hist(xs,n,facecolor='blue')
    plt.show()

def plot_scatter(ts,ttl='',xlbl='',ylbl='',xticks=None,yticks=None,lobf=None):
    """
     plots the data in ts as a scatter plot
    :param ts: data, ts = [(x1,y1),(x2,y2),...,(xn,yn)]
    :param ttl: title of the plot
    :param xlbl: the x label of the plot
    :param ylbl: the y label of the plot
    :param xticks: a list of tuples t = (pt,label) for the x axis
    :param yticks: a list of tuples t = (pt,label) for the y axis
    :param lobf: list of tuples t = (x,y) of line to plot
    """
    # set up data points
    x = [x for x,_ in ts]
    y = [y for _,y in ts]

    # plot it, add line of best fit if set
    fig,ax = plt.subplots()
    plt.scatter(x,y)
    if lobf: plt.plot([x for x,_ in lobf],[y for _,y in lobf],'-k')

    # format plot
    plt.ylim([min(y)-0.1,max(y)+0.1])
    plt.xlim([min(x)-0.1,max(x)+0.1])
    if ttl: plt.title(ttl)
    if xlbl: plt.xlabel(xlbl)
    if ylbl: plt.ylabel(ylbl)
    if xticks:
        ax.set_xticks([xt for xt,_ in xticks])
        ax.set_xticklabels([xl for _,xl in xticks])
    if yticks:
        ax.set_yticks([yt for yt,_ in yticks])
        ax.set_yticklabels([yl for _,yl in yticks])
        ax.yaxis.grid(linestyle='-') # TODO: add this as fct parameter
    plt.tick_params(axis='both',which='both',top=False,right=False)

    # and show it
    plt.show()

def plot_db_comparison(tvp,stf,mj,ts,ttl='',xlbl='',ylbl='',xticks=None,yticks=None,hl=None):
    """
    :param tvp: violence project tally
    :param stf: standford project tally
    :param mj:  mother jones tally
    :param ts: MSDB
    :param ttl: title of the plot
    :param xlbl: the x label of the plot
    :param ylbl: the y label of the plot
    :param xticks: a list of tuples t = (pt,label) for the x axis
    :param yticks: a list of tuples t = (pt,label) for the y axis
    :param hl: coordiantes if any of a highlight box to draw
    """
    # draw plots and highlight box is set
    fig,ax = plt.subplots()
    plt.plot([x for x,_ in tvp],[y for _,y in tvp],label='The Violence Project',color='red')
    plt.plot([x for x,_ in stf],[y for _,y in stf],label='Standford MSA',color='darkblue')
    plt.plot([x for x,_ in mj],[y for _,y in mj],label='Mother Jones',color='green')
    plt.plot([x for x,_ in ts],[y for _,y in ts],label='MSDB',color='black')
    if hl:
        ax.axvspan(
            hl[0],hl[1],hl[2],hl[3],color='gray',alpha=0.5,label=hl[4] if len(hl)==5 else ''
        )

    # format plot
    if ttl: plt.title(ttl)
    if xlbl: plt.xlabel(xlbl)
    if ylbl: plt.ylabel(ylbl)
    if xticks:
        ax.set_xticks([xt for xt, _ in xticks])
        ax.set_xticklabels([xl for _, xl in xticks], rotation=90)
    if yticks:
        ax.set_yticks([yt for yt, _ in yticks])
        ax.set_yticklabels([yl for _, yl in yticks])
        ax.yaxis.grid(linestyle='-')
    plt.tick_params(axis='both', which='both', top=False, bottom=False, left=False, right=False)
    ax.margins(x=0.01, y=0.01)
    plt.legend(loc='upper left')
    plt.show()

def plot_time_series(ts,ttl='',xlbl='',ylbl='',xticks=None,yticks=None,hl=None,trend=None):
    """
     plots the data in ts as a time series
    :param ts: data, ts = [(x1,y1),(x2,y2),...,(xn,yn)] (where x is a time period)
    :param ttl: title of the plot
    :param xlbl: the x label of the plot
    :param ylbl: the y label of the plot
    :param xticks: a list of tuples t = (pt,label) for the x axis
    :param yticks: a list of tuples t = (pt,label) for the y axis
    :param hl: coordiantes if any of a highlight box to draw
    :param trend: trend values if any (array of y-axis values)
    """
    # plot the time-series, and if set the highlight box and trend line
    fig,ax = plt.subplots()
    plt.plot([x for x,_ in ts],[y for _,y in ts],'-o')
    if hl:
        ax.axvspan(
            hl[0],hl[1],hl[2],hl[3],color='gray',alpha=0.5,label=hl[4] if len(hl)==5 else ''
        )
    if trend is not None: plt.plot([x for x,_ in ts],trend,'-k',label='Trend')

    # format plot
    if ttl: plt.title(ttl)
    if xlbl: plt.xlabel(xlbl)
    if ylbl: plt.ylabel(ylbl)
    if xticks:
        ax.set_xticks([xt for xt,_ in xticks])
        ax.set_xticklabels([xl for _,xl in xticks],rotation=90)
    if yticks:
        ax.set_yticks([yt for yt,_ in yticks])
        ax.set_yticklabels([yl for _,yl in yticks])
        ax.yaxis.grid(linestyle='-')
    plt.tick_params(axis='both',which='both',top=False,bottom=False,left=False,right=False)
    ax.margins(x=0.01,y=0.01)
    if hl and len(hl) == 5: plt.legend(loc='upper left')

    plt.show()

def plot_time_series_with_freq(fs,cs,ttl='',xlbl='',ylbl='',xticks=None,yticks=None,hl=None):
    """
     plots the data in ts as a time series
    :param fs: data, fs = [(x1,y1),(x2,y2),...,(xn,yn)] (where x is a time period)
    :param cs: data, Cs = [f1,f2,f3, ... ,fn] fi is the rolling mean of frequencies for xi
    :param ttl: title of the plot
    :param xlbl: the x label of the plot
    :param ylbl: the y label of the plot
    :param xticks: a list of tuples t = (pt,label) for the x axis
    :param yticks: a list of tuples t = (pt,label) for the y axis
    :param hl: coordiantes if any of a highlight box to draw
    """
    fig,ax = plt.subplots()
    plt.plot([x for x,_ in fs],[y for _,y in fs],'-o',label='Elapsed Time')
    plt.plot([x for x,_ in cs],[y for _,y in cs],label='CMA',color='red',linewidth=5)

    # configure plot
    if ttl: plt.title(ttl)
    if xlbl: plt.xlabel(xlbl)
    if ylbl: plt.ylabel(ylbl)
    plt.legend(loc='upper right')
    if xticks:
        ax.set_xticks([xt for xt,_ in xticks])
        ax.set_xticklabels([xl for _,xl in xticks],rotation=90)
    if yticks: pass # TODO:
    if hl: pass # TODO:
    plt.tick_params(axis='both',which='both',top=False,bottom=False,left=False,right=False)
    plt.show()

#### RELATED TO TRAJECTORYS

def plot_trajectory2d(steps):
    """
    plots the positions along the trajectory
    :param steps: a list of bullet steps where x and y are at idx 2,0 and 2,1
     respectively
    """
    # TODO: use zip and tuples to get xs and ys at the same time
    xs = np.round([x for x,_,_ in [pos for _,_,pos,_ in steps]],3)
    ys = np.round([y for _,y,_ in [pos for _,_,pos,_ in steps]],3)

    plt.plot(xs,ys)
    plt.show()

def plot_trajectory3d(steps):
    """
     plots a 3-d trajectory using the z-position as well
    :param steps: list of trajectory steps
    :return:
    """
    # TODO: draw x-y as line plot, z on a polar plot or some other top-down
    # setup up plotting
    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')

    # get the x, y and z positions
    xs = np.round([x for x,_,_ in [pos for _,_,pos,_ in steps]],3)
    ys = np.round([y for _,y,_ in [pos for _,_,pos,_ in steps]],3)
    zs = np.round([z for _,_,z in [pos for _,_,pos,_ in steps]],3)

    ax.plot(xs,ys,zs)
    plt.show()

def plot_coefficient(blt,param='Cd'):
    """
    plot the cd of blt from 0 to initial velocity
    :param blt: Bullet Object
    :param param: coefficient to plot
    :return:
    """
    xs = np.arange(0,blt.v_0,0.1)
    ys = []

    # NOTE: in order to do this, we have to set the 'current' velocity of the
    # bullet via access to private attribute
    v = blt.v_i
    for i in xs:
        blt._vi = i
        if param == 'BC-i':
            ys.append(blt.BC()*bls.msd2isd)
            plt.title("{}: BC imperial".format(blt.name))
            plt.ylabel("BC")
        elif param == 'BC-m':
            ys.append(blt.BC())
            plt.title("{}: BC metric".format(blt.name))
            plt.ylabel("BC")
        elif param == 'Cd':
            ys.append(blt.Cd())
            plt.title("{}: {}".format(blt.name,"$C_{D}$"))
            plt.ylabel("$C_{D}$")
        elif param == "Cd0":
            ys.append(blt.Cd0())
            plt.title("{}: {}".format(blt.name,"$C_{D0}$"))
            plt.ylabel("$C_{D0}$")
        elif param == 'Cl':
            ys.append(blt.Cl())
            plt.title("{}: {}".format(blt.name, "$C_{L}$"))
            plt.ylabel("$C_{L}$")
        elif param == 'Cn':
            ys.append(blt.Cn())
            plt.title("{}: {}".format(blt.name, "$C_{N}$"))
            plt.ylabel("$C_{N}$")
        elif param == 'Cma':
            ys.append(blt.Cmp()[0])
            plt.title("{}: {}".format(blt.name, "$C_{M_{a}}$"))
            plt.ylabel("$C_{M_{a}}$")
        elif param == 'Cmq':
            ys.append(blt.Cmp()[0])
            plt.title("{}: {}".format(blt.name, "$C_{M_{q}}$"))
            plt.ylabel("$C_{M_{q}}$")
        elif param == 'Cmm':
            ys.append(blt.Cn())
            plt.title("{}: {}".format(blt.name, "$C_{M}$"))
            plt.ylabel("$C_{M}$")
    blt._vi = v

    plt.xlabel("Velocity (m/s)")
    plt.plot(xs,ys)
    plt.show()
