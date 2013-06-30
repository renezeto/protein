from __future__ import division
from numpy import *
from matplotlib import *
from pylab import *
import sys
import glob
from file_loader import *

cell_membrane = loadtxt('./data/shape-'+f_shape+'/membrane_files/membrane-'+f_param1+'-'+f_param2+'-'+f_param3+'-'+f_param4+'-'+f_param5+'.dat')

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
            if (page[i][j] > page[i-1][j-1]) and (page[i][j] > page[i-1][j]) and (page[i][j] > page[i-1][j+1]) and (page[i][j] > page[i][j-1]) and (page[i][j] > page[i][j+1]) and (page[i][j] > page[i+1][j-1]) and (page[i][j] > page[i+1][j]) and (page[i][j] > page[i+1][j+1]):
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

def maxvectorplot(dataset):
    positions = []
    displacements = []
    nonempty = []
    for t in range(t_steps):
        if globalmax(dataset[t]) != "empty":
            nonempty += globalmax(dataset[t])
    for i in range(1,len(nonempty)-1):
        if (nonempty[i][2]>nonempty[i-1][2]) and (nonempty[i][2]>nonempty[i+1][2]):
            positions += [[nonempty[i][0],nonempty[i][1]]]
    for i in range(1,len(positions)):
        displacements += [[positions[i-1][0],positions[i-1][1],positions[i][0]-positions[i-1][0],positions[i][1]-positions[i-1][1]]]
    X,Y,U,V = zip(*displacements)
    figure()
    ax = plt.gca()
    ax.quiver(X,Y,U,V,scale_units='xy',angles='xy',scale=1)
    xlim((0,data_shape[0]))
    ylim((0,data_shape[1]))
    scatter(unzip_membrane(cell_membrane)[0],unzip_membrane(cell_membrane)[1])
    savefig('./data/shape-'+f_shape+'/plots/extrema-natp-'+f_shape+'-'+f_param1+'-'+f_param2+'-'+f_param3+'-'+f_param4+'-'+f_param5+'.pdf')
    #this only does it for NATP right now. will need to automate for other protein types later.
    print "extrema plot generated."
    return 0

maxvectorplot(data_natp_set)
