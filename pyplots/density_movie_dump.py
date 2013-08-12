from __future__ import division
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import sys
import time
import file_loader as load

#create data objects (see file_loader.py)
NflE = load.data(protein="NflE")
NflD = load.data(protein="NflD")
nATP = load.data(protein="nATP")
nE = load.data(protein="nE")
nADP = load.data(protein="nADP")
Nd = load.data(protein="Nd")
Nde = load.data(protein="Nde")

f_shape = sys.argv[1]
f_param1 = sys.argv[2]
f_param2 = sys.argv[3]
f_param3 = sys.argv[4]
f_param4 = sys.argv[5]
f_param5 = sys.argv[6]

#computes the maximum of a two dimensional array - REPLACE with method
def maxnum(page):
    Z = [0. for i in range(page.shape[0])]
    for i in range(page.shape[0]):
        Z[i] = max(page[i])
    maxval = max(Z)
    return maxval

#computes the minimum of a two dimensional array - REPLACE with method
def minnum(page):
    Z = [0. for i in range(page.shape[0])]
    for i in range(page.shape[0]):
        Z[i] = min(page[i])
        if Z[i] == 0:
            Z[i] = 50000000 #eh
    minval = min(Z)
    return minval

#computes the global maximum over a set of two dimensional arrays (stored as files)
def timemax(protein):
    Z = [0. for i in range(protein.tsteps)]
    for i in (range(protein.tsteps-1)):
        Z[i+1] = maxnum(protein.dataset[i+1])
    maxval = max(Z)
    return maxval

#computes the global minimum over a set of two dimensional arrays (stored as files)
def timemin(protein):
    Z = [0. for i in range(protein.tsteps)]
    for i in (range(protein.tsteps-1)):
        Z[i+1] = minnum(protein.dataset[i+1])
    minval = min(Z)
    return minval

#function to carry out the animation generation
def contourplt(protein):
    #get extrema values for color bar (global extrema in time)
    maxval = timemax(protein)
    minval = timemin(protein)
    plt.figure(1)

    #shell command to clean up any previous .png's, just in case (perhaps a process was cancelled midway)
    os.system("rm -f ./data/shape-"+f_shape+"/plots/"+load.hires_str+load.m_str+str(protein.protein)+"-tmp_*" \
    +"-"+str(protein.protein)+"-"+f_shape+"-"+f_param1+"-"+f_param2 \
    +"-"+f_param3+"-"+f_param4+"-"+f_param5+".png")

    Z, Y = np.meshgrid(np.arange(0,protein.datashape[1],1), np.arange(0,protein.datashape[0],1))
    #generate a sequence of .png's for each file (printed time step). these will be used to create a gif.
    for k in range(len(protein.dataset)): #fig.dpi method
        page = protein.dataset[k]
        plt.clf()
        plt.axes().set_aspect('equal', 'datalim')
        CS = plt.contourf(Z, Y, page, cmap=plt.cm.jet,origin='lower',levels=np.arange(minval,maxval+10,5))
        cbar = plt.colorbar(CS)
        plt.clim(minval,maxval)
        plt.title(str(protein.protein)+" volume density at time: "+str(k*5)+" sec")
        plt.xlabel("Z position")
        plt.ylabel("Y position")
        plt.savefig('./data/shape-'+f_shape+'/plots/'+load.hires_str+load.m_str+str(protein.protein)+'-tmp_%03d'%k+'-'+str(protein.protein)+ \
                        '-'+f_shape+'-'+f_param1+'-'+f_param2+'-'+f_param3+'-'+f_param4+'-'+f_param5+'.png',dpi=50)

    #shell command to convert all of the recently generated .png's to a single .gif using convert utility
    os.system("convert -delay 8 ./data/shape-"+f_shape+"/plots/"+load.hires_str+load.m_str+str(protein.protein)+"-tmp_*" \
    +"-"+str(protein.protein)+"-"+f_shape+"-"+f_param1+"-"+f_param2 \
    +"-"+f_param3+"-"+f_param4+"-"+f_param5+".png ./data/shape-"+f_shape \
    +"/plots/"+load.hires_str+load.m_str+"density_movie-"+str(protein.protein)+"-"+f_shape \
    +"-"+f_param1+"-"+f_param2+"-"+f_param3+"-"+f_param4+"-"+f_param5+".gif")

    #shell command to clean up the .png's
    os.system("rm -f ./data/shape-"+f_shape+"/plots/"+load.hires_str+load.m_str+str(protein.protein)+"-tmp_*" \
    +"-"+str(protein.protein)+"-"+f_shape+"-"+f_param1+"-"+f_param2 \
    +"-"+f_param3+"-"+f_param4+"-"+f_param5+".png")
    return 0

contourplt(Nde)
print "Generating flourescent Min D tagging plot:"
contourplt(NflD)
print "Done! Generating flourescent Min E tagging plot:"
contourplt(NflE)
print "Generating nATP plot:"
contourplt(nATP)
print "Done! Generating nE plot:"
contourplt(nE)
print "Done! Generating nADP plot:"
contourplt(nADP)
print "Done! Generating nD plot:"
contourplt(Nd)
print "Finished generating plots."
