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

plt.figure(1)
plt.subplot(221)
plt.contourf(natp.axes[1],natp.axes[0],average_location(natp.dataset),500)

plt.subplot(222)
plt.contourf(ne.axes[1],ne.axes[0],average_location(ne.dataset),500)

plt.subplot(223)
plt.contourf(nadp.axes[1],nadp.axes[0],average_location(nadp.dataset),500)

plt.subplot(224)
plt.contourf(nd.axes[1],nd.axes[0],average_location(nd.dataset),500)
#plt.show(1)
plt.savefig('./data/shape-'+f_shape+'/plots/time_map_compare-'+f_shape+'-'+f_param1+'-'+f_param2+'-'+f_param3+'-'+f_param4+'-'+f_param5+'.pdf')
print "time_map plot generated: subplot comparison"

for p in [natp, ne, nadp, nd]:
    plt.figure()
    plt.contourf(p.axes[1],p.axes[0],average_location(p.dataset),500)
    plt.axes().set_aspect('equal', 'datalim')
    plt.colorbar()
    plt.savefig('./data/shape-'+f_shape+'/plots/time_map-'+str(p.protein)+'-'+f_shape+'-'+f_param1+'-'+f_param2+'-'+f_param3+'-'+f_param4+'-'+f_param5+'.pdf')
    print "time_map plot generated: " + str(p.protein)
