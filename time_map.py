from __future__ import division
from matplotlib.pyplot import *
from numpy import *
from time import *
from random import *

dx = .05 #microns
data = zeros((40,40,40),dtype='float') 
data_shape = [data.shape[n] for n in range(len(data.shape))]
data_size = [data.shape[n]*dx for n in range(len(data.shape))]
axis = [arange(0,data_size[n],dx) for n in range(len(data.shape))]
tsteps = 20 
data_set = [0 for t in range(tsteps)]

#randomize some data
for t in range(tsteps):
    for n in range(data.shape[0]):
        for m in range(data.shape[1]):
            for o in range(data.shape[2]):
                data[n][m][o] = random() #loadtxt('filename.dat',delimiter=',')
    data_set[t] = data

def average_location(page):
    tsum = zeros_like(data_set[0][page])
    for t in range(tsteps):
        tsum += data_set[t][page]
    return tsum/tsteps

contourf(axis[1],axis[2],average_location(2),500)
show()
