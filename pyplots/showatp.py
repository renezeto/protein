from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import sys
import time
import file_loader as load
#
f_shape = sys.argv[1]
f_param1 = sys.argv[2]
f_param2 = sys.argv[3]
f_param3 = sys.argv[4]
f_param4 = sys.argv[5]
f_param5 = sys.argv[6]

natp = load.data(protein="natp")
ne = load.data(protein="ne")
nadp = load.data(protein="nadp")
nd = load.data(protein="nd")

def maxnum(page):
    Z = [0. for i in range(page.shape[0])]
    for i in range(page.shape[0]):
        Z[i] = max(page[i])
    maxval = max(Z)
    return maxval

def minnum(page):
    Z = [0. for i in range(page.shape[0])]
    for i in range(page.shape[0]):
        Z[i] = min(page[i])
        if Z[i] == 0:
            Z[i] = 50000000 #change this
    minval = min(Z)
    return minval

def timemax(protein):
    Z = [0. for i in range(protein.tsteps)]
    for i in (range(protein.tsteps-1)):
        Z[i+1] = maxnum(protein.dataset[i+1])
    maxval = max(Z)
    return maxval

def timemin(protein):
    Z = [0. for i in range(protein.tsteps)]
    for i in (range(protein.tsteps-1)):
        Z[i+1] = minnum(protein.dataset[i+1])
    minval = min(Z)
    return minval

def contourplt(protein):
    plt.ion()
    maxval = timemax(protein)
    minval = timemin(protein)
    for k in range(len(protein.dataset)):
        page = protein.dataset[k]
        Z, Y = np.meshgrid(np.arange(0,protein.datashape[1],1), np.arange(0,protein.datashape[0],1))
        plt.axes().set_aspect('equal', 'datalim')
        plt.ylabel('Y axis position')
        plt.xlabel('Z axis position')
        plt.title('density at time: '+repr(5*(k+1))+'s') #not correct
        CS = plt.contourf(Z, Y, page, cmap=plt.cm.jet,origin='lower',levels=np.arange(minval,maxval+10,5))
        cbar = plt.colorbar(CS)
        plt.clim(minval,maxval)
        cbar.ax.set_ylabel('density')
        plt.axes().set_aspect('equal')
        plt.draw()
        plt.clf()
    close()

contourplt(natp)
contourplt(ne)
contourplt(nadp)
contourplt(nd)
