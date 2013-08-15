from __future__ import division
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import file_loader_dump as load
import sys

# WIP

for protein in load.proteinList:
    file = np.loadtxt("./data/shape-%s/time-map-%s-%s-%s-%s-%s-%s.dat"%(load.f_shape,protein,load.f_param1,load.f_param2,load.f_param3,load.f_param4,load.f_param5))
    plt.figure()
    plt.pcolor(file)
    plt.xlim(0,file.shape[1])
    plt.ylim(0,file.shape[0])
    plt.xlabel("Z grid position")
    plt.ylabel("Y grid position")
    plt.title("Time averaged view of %s"%(protein))
    plt.colorbar()
    plt.savefig("./data/shape-%s/plots/time-map-%s-%s-%s-%s-%s-%s.pdf"%(load.f_shape,protein,load.f_param1,load.f_param2,load.f_param3,load.f_param4,load.f_param5))
