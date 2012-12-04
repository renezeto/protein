# -*- coding: utf-8 -*-
import pylab, numpy, sys, math
from pylab import *
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Slider, Button, RadioButtons

t_steps = 18

nATP = ['natp']*(t_steps)
for n in range(0,t_steps):
    if n < 10:
        nATP[n] = 'natp00'+repr(n)+'.dat'
    else:
        nATP[n] = 'natp0'+repr(n)+'.dat'

ne = ['ne']*(t_steps)
for n in range(0,t_steps):
    if n < 9:
        ne[n] = 'ne00'+repr(n)+'.dat'
    else:
        ne[n] = 'ne0'+repr(n)+'.dat'
      
def maxnum(A,ylen):
    Z = [0]*ylen
    for d in range(ylen):
        Z[d] = max(A[d])
    maxval = max(Z)
    return maxval

def minnum(A,ylen):
    Z = [0]*ylen
    for d in range(ylen):
        Z[d] = min(A[d])
    minval = min(Z)
    return minval

def animate_circles(B):
    pylab.ion()
    for k in range(0,len(B)):
        A = numpy.loadtxt(B[k])
        ylen = A.shape[0]
        zlen = A.shape[1]
        ax = subplot(111, autoscale_on=True, aspect='equal')
        ax.set_xbound(0,zlen)
        ax.set_ybound(0,ylen)
        pylab.ylabel('Y axis position')
        pylab.xlabel('Z axis position')
        pylab.title('relative density at time: '+repr(5*(k+1))+'s')
        for i in range(zlen):
            for j in range(ylen):
                if (A[j][i] < minnum(A,ylen)):
                    rad = 0
                else:
                    rad = .5*(A[j][i]/maxnum(A,ylen))
                cir = plt.Circle((i,j), radius= rad,  fc='y')
                ax.add_patch(cir)
        pylab.draw()           
        clf()
    print('graph done')
    close()
    
animate_circles(ne)
animate_circles(nATP)

'''
def interactive_circles(B):
    k=0
    axfreq = axes([0.1, 0.025, 0.1, 0.04], axisbg='r')
    kslide = Slider(axfreq, 'Time', 0, 10, valinit=1) 

    def update(val):
        clf()
        k = round(kslide.val)
        draw()
    kslide.on_changed(update)
    
    A = numpy.loadtxt(B[k])
    ylen = A.shape[0]
    zlen = A.shape[1]
    ax = subplot(111, autoscale_on=True, aspect='equal')
    subplots_adjust(left=0.25, bottom=0.25)
    t = arange(0.0, 1.0, 0.001)
    ax.set_xbound(0,zlen)
    ax.set_ybound(0,ylen)
    pylab.ylabel('Y axis position')
    pylab.xlabel('Z axis position')
    pylab.title('relative density at time: '+repr(5*(k+1))+'s')
    for i in range(zlen):
        for j in range(ylen):
            if (A[j][i] < minnum(A,ylen)):
                rad = 0
            else:
                rad = .5*(A[j][i]/maxnum(A,ylen))
            cir = plt.Circle((i,j), radius= rad,  fc='y')
            ax.add_patch(cir)
    show()
'''
