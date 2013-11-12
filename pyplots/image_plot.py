from __future__ import division
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import sys
import time
import file_loader as load
import Image

f_shape = sys.argv[1]
f_param1 = sys.argv[2]
f_param2 = sys.argv[3]
f_param3 = sys.argv[4]
f_param4 = sys.argv[5]
f_param5 = sys.argv[6]
start_time = float(sys.argv[7])
end_time = float(sys.argv[8])


dump_time_step = 2.5

#computes the global maximum over a set of two dimensional arrays (stored as files)
def timemax(protein):
    return max([v.max() for v in protein])

#computes the global minimum over a set of two dimensional arrays (stored as files)
def timemin(protein):
    return min([v.min() for v in protein])

proteinList = ["nADP","nATP","ND","NDE","nE"]
proteinLabels = ["MinD:ADP (ctyo)",
                 "MinD:ATP (cyto)",
                 "MinD:ATP (mem)",
                 "MinE:MinD:ATP",
                 "MinE (cyto)"]

proteins = [0]*len(proteinList)

for i in range(len(proteinList)):
    proteins[i] = load.data(proteinList[i],start_time,end_time)

plt.figure(figsize=(10,4))
numtimes = int(end_time/dump_time_step)- int(start_time/dump_time_step)
numproteins = len(proteins)
dZ = proteins[0].datashape[1]*1.05
dY = proteins[0].datashape[0]*1.1
kvals = range(int(end_time/dump_time_step)-int(start_time/dump_time_step))

#axes = plt.subplot(1, 1, 1, axisbg='black')

for i in range(len(proteins)):
    if (proteinList[i]=="ND" or proteinList[i]=="NDE"):
        maxval = timemax(proteins[i].dataset)/4
    else:
        maxval = timemax(proteins[i].dataset)
    minval = timemin(proteins[i].dataset[int(start_time/dump_time_step):])
    nlevels = 10
    mylevels = np.linspace(minval,(1+1.0/nlevels)*maxval,nlevels)
    Z, Y = np.meshgrid(np.arange(0,proteins[i].datashape[1],1), np.arange(0,proteins[i].datashape[0],1))
#generate a sequence of .png's for each file (printed time step). these will be used to create a gif.
    for k in kvals:
        page = proteins[i].dataset[k]
        page[page>maxval] = maxval
        plt.contourf(Y+k*dY, Z+i*dZ, page, cmap=plt.cm.hot_r, levels=mylevels)
plt.axes().get_yaxis().set_ticks([(i+0.5)*dZ for i in range(len(proteinList))])
plt.axes().get_yaxis().set_ticklabels(proteinLabels)
plt.axes().get_xaxis().set_ticks([(0.5+k)*dY for k in kvals[::4]])
plt.axes().get_xaxis().set_ticklabels([int(k*dump_time_step) for k in kvals[::4]])
plt.axes().set_aspect('equal')
plt.xlabel('time (s)')
plt.tight_layout()
plt.savefig(load.print_string("image-plot",""))
