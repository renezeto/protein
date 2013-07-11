from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import sys
import file_loader as load
#
f_shape = sys.argv[1]
f_param1 = sys.argv[2]
f_param2 = sys.argv[3]
f_param3 = sys.argv[4]
f_param4 = sys.argv[5]
f_param5 = sys.argv[6]

cell_membrane = np.loadtxt('./data/shape-'+f_shape+'/membrane_files/membrane-'+f_param1+'-'+f_param2+'-'+f_param3+'-'+f_param4+'-'+f_param5+'.dat')

nflE = load.data(protein="nflE")
nflD = load.data(protein="nflD")
nATP = load.data(protein="nATP")
nE = load.data(protein="nE")
nADP = load.data(protein="nADP")
Nd = load.data(protein="Nd")

def unzip_membrane(memdat):
    x = []
    y = []
    for i in range(len(memdat)):
        for j in range(len(memdat[i])):
            if memdat[i][j] != 0:
                x += [i]
                y += [j]
    return x,y

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

for p in [nflE, nflD, nATP, nE, nADP, Nd]:
    X,Y,U,V = zip(*maxvectorplot(p))
    plt.figure()
    plt.axes().set_aspect('equal', 'datalim')
    plt.ax = plt.gca()
    plt.ax.quiver(X,Y,U,V,scale_units='xy',angles='xy',scale=1)
    plt.xlim((0,p.datashape[0]))
    plt.ylim((0,p.datashape[1]))
    plt.scatter(unzip_membrane(cell_membrane)[0],unzip_membrane(cell_membrane)[1])
    plt.savefig('./data/shape-'+f_shape+'/plots/extrema-'+str(p.protein)+'-'+f_shape+'-'+f_param1+'-'+f_param2+'-'+f_param3+'-'+f_param4+'-'+f_param5+'.pdf')
    print "extrema plot generated: " + str(p.protein)
