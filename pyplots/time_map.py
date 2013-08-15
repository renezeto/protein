from __future__ import division
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import file_loader as load
import sys
##


for protein in load.proteinList:
    plt.figure()
    contourf(np.loadtxt("./data/shape-%s/time-map-%s-%s-%s-%s-%s-%s.dat"%(load.f_shape,protein,load.f_param1,load.f_param2,load.f_param3,load.f_param4,load.f_param5)))
    plt.savefig("./data/shape-%s/plots/time-map-%s-%s-%s-%s-%s-%s.pdf"%(load.f_shape,protein,load.f_param1,load.f_param2,load.f_param3,load.f_param4,load.f_param5))
