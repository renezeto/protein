import pylab, numpy, sys, math
from pylab import *
import matplotlib.pyplot as plt
import numpy as np
import time
import sys

t_steps = 18
cellshape = "p"#sys.argv[1]
dimA = 4#sys.argv[2]
dimB = 1#sys.argv[3]
dimC =1#sys.argv[4]
dimD =1#sys.argv[5]

nATP = ['natp']*(t_steps)
for n in range(0,t_steps):
    if n < 10:
        nATP[n] = 'natp-'+repr(cellshape)+'-0'+repr(dimA)+'-0'+repr(dimB)+'-0'+repr(dimC)+'-0'+repr(dimD)+'-00'+repr(n)+'.dat'
    else:
        nATP[n] = 'natp-'+repr(cellshape)+'-0'+repr(dimA)+'-0'+repr(dimB)+'-0'+repr(dimC)+'-0'+repr(dimD)+'-0'+repr(n)+'.dat'

ne = ['ne']*(t_steps)
for n in range(0,t_steps):
    if n < 10:
        ne[n] = 'ne-'+repr(cellshape)+'-0'+repr(dimA)+'-0'+repr(dimB)+'-0'+repr(dimC)+'-0'+repr(dimD)+'-00'+repr(n)+'.dat'
    else:
        ne[n] = 'ne-'+repr(cellshape)+'-0'+repr(dimA)+'-0'+repr(dimB)+'-0'+repr(dimC)+'-0'+repr(dimD)+'-0'+repr(n)+'.dat'

nADP = ['nADP']*(t_steps)
for n in range(0,t_steps):
  if n < 10:
    nADP[n] = 'nadp-'+repr(cellshape)+'-0'+repr(dimA)+'-0'+repr(dimB)+'-0'+repr(dimC)+'-0'+repr(dimD)+'-00'+repr(n)+'.dat'
  else:
    nADP[n] = 'nadp-'+repr(cellshape)+'-0'+repr(dimA)+'-0'+repr(dimB)+'-0'+repr(dimC)+'-0'+repr(dimD)+'-0'+repr(n)+'.dat'

nd = ['nd']*(t_steps)
for n in range(0,t_steps):
  if n < 10:
    ne[n] = 'nd-'+repr(cellshape)+'-0'+repr(dimA)+'-0'+repr(dimB)+'-0'+repr(dimC)+'-0'+repr(dimD)+'-00'+repr(n)+'.dat'
  else:
    ne[n] = 'nd-'+repr(cellshape)+'-0'+repr(dimA)+'-0'+repr(dimB)+'-0'+repr(dimC)+'-0'+repr(dimD)+'-0'+repr(n)+'.dat'

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
        CS = plt.contourf(Z, Y, A, 30, cmap=plt.cm.jet,origin='lower')
        cbar = plt.colorbar(CS)
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


'''
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
        if (A[j][i] < 1):#minnum(A,ylen)):
          rad = 0
        else:
          rad = .5*(A[j][i]/20)#maxnum(A,ylen))
        cir = plt.Circle((i,j), radius= rad,  fc='y')
        ax.add_patch(cir)
    pylab.draw()
    clf()
  print('graph done')
  close()

'''
