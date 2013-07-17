from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import sys
import file_loader as load

#params integer-ized to deal with latex issues
f_shape = sys.argv[1]
f_param1 = int(10*(float(sys.argv[2])))
f_param2 = int(10*(float(sys.argv[3])))
f_param3 = int(10*(float(sys.argv[4])))
f_param4 = int(10*(float(sys.argv[5])))
f_param5 = int(10*(float(sys.argv[6])))

#create the protein data objects (see file_loader.py)
nflE = load.data(protein="nflE")
nflD = load.data(protein="nflD")
nATP = load.data(protein="nATP")
nE = load.data(protein="nE")
nADP = load.data(protein="nADP")
Nd = load.data(protein="Nd")

#store the membrane info in a way that can be scatter plotted (update later with mem_f print)
def unzip_membrane(memdat):
    x = []
    y = []
    for i in range(len(memdat)):
        for j in range(len(memdat[i])):
            if memdat[i][j] != 0:
                x += [i]
                y += [j]
    return x,y

#function to compute the local maxima of a 2 dimensional array (used for local maxima in space)
def localmax(page):
    lmax = []
    for i in range(1,len(page)-1):
        index_pairs = list(enumerate(page[i]))
        for j in range(1,len(page[0])-1):
            if (page[i][j] > page[i-1][j-1]) and (page[i][j] > page[i-1][j]) and (page[i][j] > page[i-1][j+1]) and (page[i][j] > page[i][j-1]) and (page[i][j] > page[i][j+1]) and (page[i][j] > page[i+1][j-1]) and (page[i][j] > page[i+1][j]) and (page[i][j] > page[i+1][j+1]): #rewrite using chains
                lmax += [(i,index_pairs[j][0],index_pairs[j][1])]
    if (lmax == []):
        return "empty"
    return lmax

#takes the maximum of the local maxima in space
def globalmax(page):
    gmax = [(0.,0.,0.)]
    lmax = localmax(page)
    if lmax != "empty":
        for i in range(len(lmax)):
            if (lmax[i][2] > gmax[0][2]):
                gmax = [lmax[i]]
            elif (lmax[i][2] == gmax[0][2]):
                gmax += [lmax[i]]
                if gmax[0] == 0.:
                    gmax.pop(0)
    else:
        return "empty"
    if len(gmax)>1:
        gmax.pop(1)
    return gmax

#import membrane file for printing onto the plot
    cell_membrane = np.loadtxt('./data/shape-'+load.f_shape+'/membrane_files/membrane-'+load.f_param1+'-'+load.f_param2+'-'+load.f_param3+'-'+load.f_param4+'-'+load.f_param5+'.dat')

#returns a list of displacement vectors between points that are global maxima in space (in one data file) and
#local maxima in time (in all the data files)
def maxvectorplot(protein):
    positions = []
    displacements = []
    nonempty = []
    for t in range(protein.tsteps):
        if globalmax(protein.dataset[t]) != "empty":
            nonempty += globalmax(protein.dataset[t])
    for i in range(1,len(nonempty)-1):
        if (nonempty[i][2]>nonempty[i-1][2]) and (nonempty[i][2]>nonempty[i+1][2]):
            positions += [[nonempty[i][0],nonempty[i][1]]]
    for i in range(1,len(positions)):
        displacements += [[positions[i-1][0],positions[i-1][1],positions[i][0]-positions[i-1][0],positions[i][1]-positions[i-1][1]]]
    return displacements

#generate the plots based on maxvectorplot() output, and saves them.
for p in [nflE, nflD, nATP, nE, nADP, Nd]:
    X,Y,U,V = zip(*maxvectorplot(p))
    plt.figure()
    plt.axes().set_aspect('equal', 'datalim')
    plt.ax = plt.gca()
    plt.ax.quiver(X,Y,U,V,scale_units='xy',angles='xy',scale=1)
    plt.xlim((0,p.datashape[0]))
    plt.ylim((0,p.datashape[1]))
    plt.scatter(unzip_membrane(cell_membrane)[0],unzip_membrane(cell_membrane)[1])
    plt.savefig('./data/shape-'+f_shape+'/plots/'+load.hires_str+load.m_str+'extrema-'+str(p.protein)+'-'+f_shape+'-'+f_param1+'-'+f_param2+'-'+f_param3+'-'+f_param4+'-'+f_param5+'.pdf')
    print "extrema plot generated: " + str(p.protein)

