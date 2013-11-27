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


dump_time_step = 0.5

#computes the global maximum over a set of two dimensional arrays (stored as files)
def timemax(protein):
    return max([v.max() for v in protein])

#computes the global minimum over a set of two dimensional arrays (stored as files)
def timemin(protein):
    return min([v.min() for v in protein])

proteinList = ["nADP","nATP","ND","NDE","nE"]
proteinLabels = ["MinD:ADP (cyto)",
                 "MinD:ATP (cyto)",
                 "MinD:ATP (mem)",
                 "MinE:MinD:ATP (mem)",
                 "MinE (cyto)"]

proteins = [0]*len(proteinList)

for i in range(len(proteinList)):
    proteins[i] = load.data(proteinList[i],start_time,end_time)

plt.figure(figsize=(9,3.5))
numtimes = int(end_time/dump_time_step)- int(start_time/dump_time_step)
numproteins = len(proteins)
skip_times = 2 # only plot every skip_times of the snapshots
dZ = proteins[0].datashape[1]*1.05
dY = proteins[0].datashape[0]*1.1/skip_times
kvals = range(0, int(end_time/dump_time_step)-int(start_time/dump_time_step), skip_times)

#axes = plt.subplot(1, 1, 1, axisbg='black')

for i in range(len(proteins)):
    if (proteinList[i]=="ND" or proteinList[i]=="NDE"):
        maxval = timemax(proteins[i].dataset)/4
    else:
        maxval = timemax(proteins[i].dataset)
    #minval = timemin(proteins[i].dataset[int(start_time/dump_time_step):])
    nlevels = 20
    mylevels = np.linspace(0,(1+1.0/nlevels)*maxval,nlevels)
    Z, Y = np.meshgrid(np.arange(0,proteins[i].dataset[0].shape[1],1),
                       np.arange(0,proteins[i].dataset[0].shape[0],1))
#generate a sequence of .png's for each file (printed time step). these will be used to create a gif.
    for k in kvals:
        page = proteins[i].dataset[k]
        page[page>maxval] = maxval
        plt.contourf(Y+k*dY, proteins[0].datashape[1]-Z+i*dZ, page, cmap=plt.cm.hot_r, levels=mylevels)
plt.axes().get_yaxis().set_ticks([(i+0.5)*dZ for i in range(len(proteinList))])
plt.axes().get_yaxis().set_ticklabels(proteinLabels)
plt.axes().get_xaxis().set_ticks([(0.5+k)*dY for k in kvals[::int(2.5*skip_times)]])
plt.axes().get_xaxis().set_ticklabels([int(k*dump_time_step) for k in kvals[::int(2.5*skip_times)]])
plt.axes().set_aspect('equal')
plt.axes().get_yaxis().set_ticks_position('none')
plt.axes().get_xaxis().set_ticks_position('none')
plt.xlabel('time (s)')

# I got the colorbar bit here from
# http://matplotlib.org/users/tight_layout_guide.html, which explains
# the difficulties with tight_layout() and colorbar positioning.
from mpl_toolkits.axes_grid1 import make_axes_locatable
divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", "3%", pad="1%")
cbar = plt.colorbar(cax=cax,ticks=[0,maxval])
cbar.ax.set_yticklabels(['0', 'max'])

#plt.colorbar(ticks=[])
plt.tight_layout()
plt.savefig(load.print_string("image-plot",""))
