import pylab, numpy, sys, math
from pylab import *
import matplotlib.pyplot as plt
import numpy as np
import time
import sys
import glob

f_shape = sys.argv[1]
f_param1 = sys.argv[2]
f_param2 = sys.argv[3]
f_param3 = sys.argv[4]
f_param4 = sys.argv[5]
#f_density = sys.argv[6]
f_param5 = sys.argv[6]

#for fn in glob.iglob('../data/shape-'+f_shape+'/natp-'+f_shape+'-'+f_param1+'-'+f_param2+'-'+f_param3+'-'+f_param4+'-'+f_param5+'*.dat'):
#        dat_filenames.append(fn)

dat_filenames = []
for fn in glob.iglob('./data/shape-'+f_shape+'/natp-'+f_shape+'-'+f_param1+'-'+f_param2+'-'+f_param3+'-'+f_param4+'-'+f_param5+'*.dat'):
        dat_filenames.append(fn)
dat_filenames = sorted(dat_filenames)
print "Filenames are: "
print f_param5
for i in range(4):
        dat_filenames.pop(0)
t_steps = len(dat_filenames)
data_natp_set = array([np.loadtxt(dat_filenames[i]) for i in range(t_steps)])

dx = .05 #microns
data = data_natp_set[0] #retrieves data format.
data_shape = [data.shape[n] for n in range(len(data.shape))]
data_size = [data.shape[n]*dx for n in range(len(data.shape))]
axis = [arange(0,data_size[n],dx) for n in range(len(data.shape))]

def maxnum(page):
    Z = [0. for i in range(data.shape[0])]
    for i in range(data.shape[0]):
        Z[i] = max(page[i])
    maxval = max(Z)
    return maxval

def minnum(page):
    Z = [0. for i in range(data.shape[0])]
    for i in range(data.shape[0]):
        Z[i] = min(page[i])
        if Z[i] == 0:
            Z[i] = 50000000
    minval = min(Z)
    return minval

def timemax(dataset):
    Z = [0. for i in range(t_steps)]
    for i in (range(t_steps-1)):
        page = dataset[i+1]
        Z[i+1] = maxnum(page)
    maxval = max(Z)
    return maxval

def timemin(dataset):
    Z = [0. for i in range(t_steps)]
    for i in (range(t_steps-1)):
        page = dataset[i+1]
        Z[i+1] = minnum(page)
    minval = min(Z)
    return minval

def contourplt(dataset):
    pylab.ion()
    minval = timemin(dataset)
    maxval = timemax(dataset)
    for k in range(len(dataset)):
        page = dataset[k]
        Z, Y = np.meshgrid(np.arange(0,data.shape[1],1), np.arange(0,data.shape[0],1))
        pylab.ylabel('Y axis position')
        pylab.xlabel('Z axis position')
        pylab.title('density at time: '+repr(5*(k+1))+'s')
        CS = plt.contourf(Z, Y, page, cmap=plt.cm.jet,origin='lower',levels=np.arange(minval,maxval+10,5))
        cbar = plt.colorbar(CS)
        plt.clim(minval,maxval)
        cbar.ax.set_ylabel('density')
        dt = 1
        time.sleep(dt)
        plt.axes().set_aspect('equal')
        pylab.draw()
        clf()
    print('graph done')
    close()

contourplt(data_natp_set)

## old file importing - has been replaced with a better method
# nATP = ['natp']*(t_steps)
# for n in range(0,t_steps): #change '+sys.argv[1]+' to '+sys.argv[1]+' on diff comp
#     if n < 10:
#         nATP[n] = 'shape-'+sys.argv[1]+'/natp-'+sys.argv[1]+'-'+sys.argv[2]+'-'+sys.argv[3]+'-'+sys.argv[4]+'-'+sys.argv[5]+'-'+sys.argv[6]+'-00'+repr(n)+'.dat'
#     else:
#         nATP[n] = 'shape-'+sys.argv[1]+'/natp-'+sys.argv[1]+'-'+sys.argv[2]+'-'+sys.argv[3]+'-'+sys.argv[4]+'-'+sys.argv[5]+'-'+sys.argv[6]+'-0'+repr(n)+'.dat'

# ne = ['ne']*(t_steps)
# for n in range(0,t_steps):
#     if n < 10:
#         ne[n] = 'shape-'+sys.argv[1]+'/ne-'+sys.argv[1]+'-'+sys.argv[2]+'-'+sys.argv[3]+'-'+sys.argv[4]+'-'+sys.argv[5]+'-'+sys.argv[6]+'-00'+repr(n)+'.dat'
#     else:
#         ne[n] = 'shape-'+sys.argv[1]+'/ne-'+sys.argv[1]+'-'+sys.argv[2]+'-'+sys.argv[3]+'-'+sys.argv[4]+'-'+sys.argv[5]+'-'+sys.argv[6]+'-0'+repr(n)+'.dat'

# nADP = ['nADP']*(t_steps)
# for n in range(0,t_steps):
#   if n < 10:
#     nADP[n] = 'shape-'+sys.argv[1]+'/nadp-'+sys.argv[1]+'-'+sys.argv[2]+'-'+sys.argv[3]+'-'+sys.argv[4]+'-'+sys.argv[5]+'-'+sys.argv[6]+'-00'+repr(n)+'.dat'
#   else:
#     nADP[n] = 'shape-'+sys.argv[1]+'/nadp-'+sys.argv[1]+'-'+sys.argv[2]+'-'+sys.argv[3]+'-'+sys.argv[4]+'-'+sys.argv[5]+'-'+sys.argv[6]+'-0'+repr(n)+'.dat'

# nd = ['nd']*(t_steps)
# for n in range(0,t_steps):
#   if n < 10:
#     nd[n] = 'shape-'+sys.argv[1]+'/nd-'+sys.argv[1]+'-'+sys.argv[2]+'-'+sys.argv[3]+'-'+sys.argv[4]+'-'+sys.argv[5]+'-'+sys.argv[6]+'-00'+repr(n)+'.dat'
#   else:
#     nd[n] = 'shape-'+sys.argv[1]+'/nd-'+sys.argv[1]+'-'+sys.argv[2]+'-'+sys.argv[3]+'-'+sys.argv[4]+'-'+sys.argv[5]+'-'+sys.argv[6]+'-0'+repr(n)+'.dat'

# cellshape = sys.argv[1]
# dimA = sys.argv[2]
# dimB = sys.argv[3]
# dimC =sys.argv[4]
# dimD =sys.argv[5]
