from __future__ import division
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import file_loader as load
import sys

# WIP

for protein in load.proteinList:
    file = np.loadtxt("./data/shape-%s/%s%s%stime-map-%s-%s-%s-%s-%s-%s-%s.dat"%(load.f_shape,load.debug_str,load.hires_str,load.slice_str,protein,load.f_shape,load.f_param1,load.f_param2,load.f_param3,load.f_param4,load.f_param5))
    plt.figure()
    plt.pcolor(file)
    plt.xlim(0,file.shape[1])
    plt.ylim(0,file.shape[0])
    plt.xlabel("Z grid position")
    plt.ylabel("Y grid position")
    plt.title("Time averaged view of %s"%(protein))
    plt.colorbar()
    plt.savefig(load.print_string("time-map",protein))
