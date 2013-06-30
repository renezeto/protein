from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import file_loader as load
import sys

natp = load.data(protein="natp")
ne = load.data(protein="ne")
nadp = load.data(protein="nadp")
nd = load.data(protein="nd")

f_shape = sys.argv[1]
f_param1 = sys.argv[2]
f_param2 = sys.argv[3]
f_param3 = sys.argv[4]
f_param4 = sys.argv[5]
f_param5 = sys.argv[6]

def average_location(dataset):
    tsum = np.zeros_like(dataset[0])
    for t in range(natp.tsteps):
            tsum += dataset[t]
    return tsum/(natp.tsteps)

for p in [natp, ne, nadp, nd]:
    plt.figure()
    plt.contourf(p.axes[1],p.axes[0],average_location(p.dataset),500) #why are axes backwards?
    plt.colorbar()
    plt.savefig('./data/shape-'+f_shape+'/plots/time_map-'+str(p.protein)+'-'+f_shape+'-'+f_param1+'-'+f_param2+'-'+f_param3+'-'+f_param4+'-'+f_param5+'.pdf')
    print "time_map plot generated: " + str(p.protein)
