from __future__ import division
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
import pylab
import file_loader as load

data = np.loadtxt("./data/shape-%s/membrane_files/%s%s%ssections-%s-%s-%s-%s-%s-%s.dat"%(load.f_shape,load.debug_str,load.hires_str,load.slice_str,load.f_shape,load.f_param1,load.f_param2,load.f_param3,load.f_param4,load.f_param5))

plt.figure()
plt.pcolormesh(data)
plt.savefig(load.print_string("sections_plot",""))
