from __future__ import division
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import file_loader as load
import sys

#create protein data objects (see file_loader.py)
nflE = load.data(protein="nflE")
nflD = load.data(protein="nflD")
nATP = load.data(protein="nATP")
nE = load.data(protein="nE")
nADP = load.data(protein="nADP")
Nd = load.data(protein="Nd")

f_shape = sys.argv[1]
f_param1 = int(10*(float(sys.argv[2])))
f_param2 = int(10*(float(sys.argv[3])))
f_param3 = int(10*(float(sys.argv[4])))
f_param4 = int(10*(float(sys.argv[5])))
f_param5 = int(10*(float(sys.argv[6])))

#compute time average (removing first 10% happens in file loading)
def average_location(dataset):
    tsum = np.zeros_like(dataset[0])
    for t in range(nATP.tsteps):
            tsum += dataset[t]
    return tsum/(nATP.tsteps)

print './data/shape-'+f_shape+'/plots/'+load.hires_str+load.m_str+'time-map-compare-'+f_shape+'-'+str(f_param1)+'-'+str(f_param2)+'-'+str(f_param3)+'-'+str(f_param4)+'-'+str(f_param5)+'.pdf'
#create comparison subplots and print them to one file. 
plt.figure(1)
plt.subplot(221)
plt.contourf(nATP.axes[1],nATP.axes[0],average_location(nATP.dataset),500)
plt.title('nATP', fontsize=20)

plt.subplot(222)
plt.contourf(nE.axes[1],nE.axes[0],average_location(nE.dataset),500)
plt.title('nE', fontsize=20)

plt.subplot(223)
plt.contourf(nADP.axes[1],nADP.axes[0],average_location(nADP.dataset),500)
plt.title('nADP', fontsize=20)

plt.subplot(224)
plt.contourf(Nd.axes[1],Nd.axes[0],average_location(Nd.dataset),500)
plt.title('Nd', fontsize=20) #load.hires_str and load.m_str are empty strings unless -hires and -slice are in the sys.argv list. see file_loader.py
plt.savefig('./data/shape-'+f_shape+'/plots/'+load.hires_str+load.m_str+'time-map-compare-'+f_shape+'-'+str(f_param1)+'-'+str(f_param2)+'-'+str(f_param3)+'-'+str(f_param4)+'-'+str(f_param5)+'.pdf')
print "time_map plot generated: subplot comparison"

#plot each of the data set time maps individually.
for p in [nflE, nflD, nATP, nE, nADP, Nd]:
    plt.figure()
    plt.contourf(p.axes[1],p.axes[0],average_location(p.dataset),500)
    plt.axes().set_aspect('equal', 'datalim')
    plt.colorbar()
    plt.savefig('./data/shape-'+f_shape+'/plots/'+load.hires_str+load.m_str+'time-map-'+str(p.protein)+'-'+f_shape+'-'+str(f_param1)+'-'+str(f_param2)+'-'+str(f_param3)+'-'+str(f_param4)+'-'+str(f_param5)+'.pdf')
    print "time_map plot generated: " + str(p.protein)
