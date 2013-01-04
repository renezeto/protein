import pylab, numpy, sys, math
from pylab import *
import matplotlib.pyplot as plt
import numpy as np
import time
import sys

t_steps = 18
cellshape = sys.argv[1]
dimA = sys.argv[2]
dimB = sys.argv[3]
dimC =sys.argv[4]
dimD =sys.argv[5]

nATP = ['natp']*(t_steps)
for n in range(0,t_steps): #change '+sys.argv[1]+' to '+sys.argv[1]+' on diff comp
    if n < 10:
        nATP[n] = 'shape-'+sys.argv[1]+'/natp-'+sys.argv[1]+'-'+sys.argv[2]+'-'+sys.argv[3]+'-'+sys.argv[4]+'-'+sys.argv[5]+'-00'+repr(n)+'.dat'
    else:
        nATP[n] = 'shape-'+sys.argv[1]+'/natp-'+sys.argv[1]+'-'+sys.argv[2]+'-'+sys.argv[3]+'-'+sys.argv[4]+'-'+sys.argv[5]+'-0'+repr(n)+'.dat'

ne = ['ne']*(t_steps)
for n in range(0,t_steps):
    if n < 10:
        ne[n] = 'shape-'+sys.argv[1]+'/ne-'+sys.argv[1]+'-'+sys.argv[2]+'-'+sys.argv[3]+'-'+sys.argv[4]+'-'+sys.argv[5]+'-00'+repr(n)+'.dat'
    else:
        ne[n] = 'shape-'+sys.argv[1]+'/ne-'+sys.argv[1]+'-'+sys.argv[2]+'-'+sys.argv[3]+'-'+sys.argv[4]+'-'+sys.argv[5]+'-0'+repr(n)+'.dat'

nADP = ['nADP']*(t_steps)
for n in range(0,t_steps):
  if n < 10:
    nADP[n] = 'shape-'+sys.argv[1]+'/nadp-'+sys.argv[1]+'-'+sys.argv[2]+'-'+sys.argv[3]+'-'+sys.argv[4]+'-'+sys.argv[5]+'-00'+repr(n)+'.dat'
  else:
    nADP[n] = 'shape-'+sys.argv[1]+'/nadp-'+sys.argv[1]+'-'+sys.argv[2]+'-'+sys.argv[3]+'-'+sys.argv[4]+'-'+sys.argv[5]+'-0'+repr(n)+'.dat'

nd = ['nd']*(t_steps)
for n in range(0,t_steps):
  if n < 10:
    nd[n] = 'shape-'+sys.argv[1]+'/nd-'+sys.argv[1]+'-'+sys.argv[2]+'-'+sys.argv[3]+'-'+sys.argv[4]+'-'+sys.argv[5]+'-00'+repr(n)+'.dat'
  else:
    nd[n] = 'shape-'+sys.argv[1]+'/nd-'+sys.argv[1]+'-'+sys.argv[2]+'-'+sys.argv[3]+'-'+sys.argv[4]+'-'+sys.argv[5]+'-0'+repr(n)+'.dat'

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

def timemax(A,ylen):
    Z = [0]*t_steps
    for d in range(t_steps):
        Z[d] = max(maxnum(A,ylen))
    maxval = max(Z)
    return maxval
    

def contourplt(B):
    pylab.ion()
    for k in range(0,len(B)):
        A = numpy.loadtxt(B[k])
        zlen = np.arange(0,A.shape[1],1)
        ylen = np.arange(0,A.shape[0],1)
        Z, Y = np.meshgrid(zlen, ylen)
        pylab.ylabel('Y axis position')
        pylab.xlabel('Z axis position')
        pylab.title('density at time: '+repr(5*(k+1))+'s')
        CS = plt.contourf(Z, Y, A, 100, cmap=plt.cm.jet,origin='lower')
        cbar = plt.colorbar(CS)
        plt.clim(0,400)
        cbar.ax.set_ylabel('density')
        dt = 1
        time.sleep(dt)
        plt.axes().set_aspect('equal')
        pylab.draw()           
        clf()
    print('graph done')
    close()
    
contourplt(nATP)
contourplt(nADP)
contourplt(ne)
contourplt(nd)


