from __future__ import division
from numpy import *
from matplotlib import *
from pylab import *
import sys
import glob

f_shape = sys.argv[1]
f_param1 = sys.argv[2]
f_param2 = sys.argv[3]
f_param3 = sys.argv[4]
f_param4 = sys.argv[5]
f_seed = sys.argv[6]

dat_filenames = []
for fn in glob.iglob('shape-'+f_shape+'/natp-'+f_shape+'-'+f_param1+'-'+f_param2+'-'+f_param3+'-'+f_param4+'-'+f_seed+'*.dat'):
        dat_filenames.append(fn)
dat_filenames = sorted(dat_filenames)
dat_filenames.pop(0)
t_steps = len(dat_filenames)
data_natp_set = array([np.loadtxt(dat_filenames[i]) for i in range(t_steps)])

dx = .05 #microns
data = data_natp_set[0] #retrieves data format.
data_shape = [data.shape[n] for n in range(len(data.shape))]
data_size = [data.shape[n]*dx for n in range(len(data.shape))]
axis = [arange(0,data_size[n],dx) for n in range(len(data.shape))]

#between tstep-1 and tstep, track the location of the maximum. 

X, Y = meshgrid(arange(0,10,.1),arange(0,10,.1))
Z = -(X - 5*(X+1)/(X+1))**2 - (Y - 5*(Y+1)/(Y+1))**2 #+ 1*(X+1)/(X+1)

def globalmax(page):
    lmax = []
    for i in range(1,len(page)-1):
        index_pairs = list(enumerate(page[i]))
        for j in range(1,len(page[0])-1):
            if (page[i][j] > page[i-1][j-1]) and (page[i][j] > page[i-1][j]) and (page[i][j] > page[i-1][j+1]) and (page[i][j] > page[i][j-1]) and (page[i][j] > page[i][j+1]) and (page[i][j] > page[i+1][j-1]) and (page[i][j] > page[i+1][j]) and (page[i][j] > page[i+1][j+1]):
                lmax += [(i,index_pairs[j][0],index_pairs[j][1])]
    if lmax == []:
        return "empty"
    gmax = [0.]
    for i in range(len(lmax)):
        if (lmax[i][2] > max(gmax)):
            gmax = [lmax[i]]
        elif (lmax[i][2] == max(gmax)):
            gmax += [lmax[i]]
    if gmax[0] == 0.:
        gmax.pop(0)
    return gmax

def maxtracker(dataset):
    locs = zeros_like(dataset[0])
    for t in range(t_steps):
        if globalmax(dataset[t]) != "empty":
            for i in range(len(globalmax(dataset[t]))):
                x, y, m = globalmax(dataset[t])[i]
                locs[x][y] = m
    return locs

def maxvectorplot(dataset):
    positions = []
    displacements = []
    for t in range(t_steps):
        if globalmax(dataset[t]) != "empty":
            positions += [(globalmax(dataset[t])[0][0],globalmax(dataset[t])[0][1])]
    for i in range(1,len(positions)):
        displacements += [(positions[i-1][0]-positions[i][0],positions[i-1][1]-positions[i][0])]
    return 0

f = maxvectorplot(data_natp_set)
