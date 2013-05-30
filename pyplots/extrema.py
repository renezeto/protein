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
for fn in glob.iglob('../data/shape-'+f_shape+'/natp-'+f_shape+'-'+f_param1+'-'+f_param2+'-'+f_param3+'-'+f_param4+'-'+f_seed+'*.dat'):
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

cell_membrane = np.loadtxt('membrane.dat') #update membrane.dat filename syntax


def unzip_membrane(memdat):
    x = []
    y = []
    for i in range(len(memdat)):
        for j in range(len(memdat[i])):
            if memdat[i][j] != 0:
                x += [i]
                y += [j]
    return x,y

def localmax(page):
    lmax = []
    for i in range(1,len(page)-1):
        index_pairs = list(enumerate(page[i]))
        for j in range(1,len(page[0])-1):
            if (page[i][j] > page[i-1][j-1]) and (page[i][j] > page[i-1][j]) and (page[i][j] > page[i-1][j+1]) and (page[i][j] > page[i][j-1]) and (page[i][j] > page[i][j+1]) and (page[i][j] > page[i+1][j-1]) and (page[i][j] > page[i+1][j]) and (page[i][j] > page[i+1][j+1]):
                lmax += [(i,index_pairs[j][0],index_pairs[j][1])]
    if (lmax == []):
        return "empty"
    return lmax

def globalmax(page):
    gmax = [0.]
    lmax = localmax(page)
    if lmax != "empty":
        for i in range(len(lmax)):
            if (lmax[i][2] > max(gmax)):
                gmax = [lmax[i]]
            elif (lmax[i][2] == max(gmax)):
                gmax += [lmax[i]]
                if gmax[0] == 0.:
                    gmax.pop(0)
    else: 
        return "empty"
    return gmax

def maxtracker(dataset):
    locs = zeros_like(dataset[0])
    for t in range(1,t_steps-1):
        if globalmax(dataset[t]) != "empty":
            for i in range(len(globalmax(dataset[t]))):
                if (globalmax(dataset[t])[i] > globalmax(dataset[t+1])[i]) and (globalmax(dataset[t])[i] > globalmax(dataset[t-1])[i]):
                    x, y, m = globalmax(dataset[t])[i] #[i] accounting for possibility of multiple maxima
                    locs[x][y] = m
    return locs

def maxvectorplot(dataset):
    positions = []
    displacements = []
    for t in range(1,t_steps-1):
        if globalmax(dataset[t]) != "empty":
            for i in range(len(globalmax(dataset[t]))):
                if (globalmax(dataset[t])[i] > globalmax(dataset[t+1])[i]) and (globalmax(dataset[t])[i] > globalmax(dataset[t-1])[i]):
                    positions += [[globalmax(dataset[t])[i][0],globalmax(dataset[t])[i][1]]]
    for i in range(1,len(positions)):
        displacements += [[positions[i-1][0],positions[i-1][1],positions[i][0]-positions[i-1][0],positions[i][1]-positions[i-1][1]]] #[tail_x, tail_y, head_x, head_y]
    X,Y,U,V = zip(*displacements)
    plt.figure()
    ax = plt.gca()
    ax.quiver(X,Y,U,V,angles='xy',scale_units='xy',scale=1,linewidths=.05)
<<<<<<< HEAD:extrema.py
    #scatter(unzip_membrane(cell_membrane)[0],unzip_membrane(cell_membrane)[1], marker='.')
=======
#    scatter(unzip_membrane(cell_membrane)[0],unzip_membrane(cell_membrane)[1], marker='.')
>>>>>>> 541a87a02008f4ac1b4e28dc7bdf163be85b54d5:pyplots/extrema.py
    plt.xlim((0,data_shape[0]))
    plt.ylim((0,data_shape[1]))
    show()
    return 0

maxvectorplot(data_natp_set)
