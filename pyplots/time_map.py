from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import file_loader as load
import sys

nflE = load.data(protein="nflE")
nflD = load.data(protein="nflD")
nATP = load.data(protein="nATP")
nE = load.data(protein="nE")
nADP = load.data(protein="nADP")
Nd = load.data(protein="Nd")

f_shape = sys.argv[1]
f_param1 = sys.argv[2]
f_param2 = sys.argv[3]
f_param3 = sys.argv[4]
f_param4 = sys.argv[5]
f_param5 = sys.argv[6]

def average_location(dataset):
    tsum = np.zeros_like(dataset[0])
    for t in range(nATP.tsteps):
            tsum += dataset[t]
    return tsum/(nATP.tsteps)

plt.figure(1)
plt.subplot(221)
plt.contourf(nATP.axes[1],nATP.axes[0],average_location(nATP.dataset),500)

plt.subplot(222)
plt.contourf(nE.axes[1],nE.axes[0],average_location(nE.dataset),500)

plt.subplot(223)
plt.contourf(nADP.axes[1],nADP.axes[0],average_location(nADP.dataset),500)

plt.subplot(224)
plt.contourf(Nd.axes[1],Nd.axes[0],average_location(Nd.dataset),500)
#plt.show(1)
plt.savefig('./data/shape-'+f_shape+'/plots/time_map_compare-'+f_shape+'-'+f_param1+'-'+f_param2+'-'+f_param3+'-'+f_param4+'-'+f_param5+'.pdf')
print "time_map plot generated: subplot comparison"

for p in [nflE, nflD, nATP, nE, nADP, Nd]:
    plt.figure()
    plt.contourf(p.axes[1],p.axes[0],average_location(p.dataset),500)
    plt.axes().set_aspect('equal', 'datalim')
    plt.colorbar()
    plt.savefig('./data/shape-'+f_shape+'/plots/time_map-'+str(p.protein)+'-'+f_shape+'-'+f_param1+'-'+f_param2+'-'+f_param3+'-'+f_param4+'-'+f_param5+'.pdf')
    print "time_map plot generated: " + str(p.protein)
