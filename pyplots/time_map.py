from __future__ import division
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import file_loader as load
import sys

#create protein data objects (see file_loader.py)
NflE = load.data(protein="NflE")
NflD = load.data(protein="NflD")
nATP = load.data(protein="nATP")
nE = load.data(protein="nE")
nADP = load.data(protein="nADP")
Nd = load.data(protein="Nd")
Nde = load.data(protein="Nde")

#compute time average (removing first 10% happens in file loading)
def average_location(dataset):
    tsum = np.zeros_like(dataset[0])
    for t in range(nATP.tsteps):
            tsum += dataset[t]
    return tsum/(nATP.tsteps)

#plot each of the data set time maps individually.
for p in [NflE, NflD, nATP, nE, nADP, Nd]:
    plt.figure()
    plt.contourf(p.axes[0],p.axes[1],average_location(p.dataset),500)
    plt.axes().set_aspect('equal', 'datalim')
    plt.colorbar()
    plt.savefig(load.print_string("time-map",p))
