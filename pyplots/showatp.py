from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import time
import file_loader as load

f_shape = sys.argv[1]
f_param1 = sys.argv[2]
f_param2 = sys.argv[3]
f_param3 = sys.argv[4]
f_param4 = sys.argv[5]
f_param5 = sys.argv[6]

nflE = load.data(protein="nflE")
nflD = load.data(protein="nflD")
nATP = load.data(protein="nATP")
nE = load.data(protein="nE")
nADP = load.data(protein="nADP")
Nd = load.data(protein="Nd")

def maxnum(page):
    Z = [0. for i in range(page.shape[0])]
    for i in range(page.shape[0]):
        Z[i] = max(page[i])
    maxval = max(Z)
    return maxval

def minnum(page):
    Z = [0. for i in range(page.shape[0])]
    for i in range(page.shape[0]):
        Z[i] = min(page[i])
        if Z[i] == 0:
            Z[i] = 50000000 #eh
    minval = min(Z)
    return minval

def timemax(protein):
    Z = [0. for i in range(protein.tsteps)]
    for i in (range(protein.tsteps-1)):
        Z[i+1] = maxnum(protein.dataset[i+1])
    maxval = max(Z)
    return maxval

def timemin(protein):
    Z = [0. for i in range(protein.tsteps)]
    for i in (range(protein.tsteps-1)):
        Z[i+1] = minnum(protein.dataset[i+1])
    minval = min(Z)
    return minval

def contourplt(protein):
    maxval = timemax(protein)
    minval = timemin(protein)
    plt.figure(1)
    os.system("rm -f ./data/shape-"+f_shape+"/plots/"+str(protein.protein)+"-tmp_*.png")
    Z, Y = np.meshgrid(np.arange(0,protein.datashape[1],1), np.arange(0,protein.datashape[0],1))
    for k in range(len(protein.dataset)): #fig.dpi method
        page = protein.dataset[k]
        plt.clf()
        plt.axes().set_aspect('equal', 'datalim')
        CS = plt.contourf(Z, Y, page, cmap=plt.cm.jet,origin='lower',levels=np.arange(minval,maxval+10,5))
        cbar = plt.colorbar(CS)
        plt.clim(minval,maxval)
        plt.title(str(protein.protein)+" volume density at time: "+str(k/2)+" sec")
        plt.xlabel("Z position")
        plt.ylabel("Y position")
        if k<10:
            plt.savefig('./data/shape-'+f_shape+'/plots/'+str(protein.protein)+'-tmp_00'+str(k)+'-'+str(protein.protein)+ \
                            '-'+f_shape+'-'+f_param1+'-'+f_param2+'-'+f_param3+'-'+f_param4+'-'+f_param5+'.png',dpi=50)
        elif k<100:
            plt.savefig('./data/shape-'+f_shape+'/plots/'+str(protein.protein)+'-tmp_0'+str(k)+'-'+str(protein.protein)+ \
                            '-'+f_shape+'-'+f_param1+'-'+f_param2+'-'+f_param3+'-'+f_param4+'-'+f_param5+'.png',dpi=50)
        else:
            plt.savefig('./data/shape-'+f_shape+'/plots/'+str(protein.protein)+'-tmp_'+str(k)+'-'+str(protein.protein)+ \
                            '-'+f_shape+'-'+f_param1+'-'+f_param2+'-'+f_param3+'-'+f_param4+'-'+f_param5+'.png',dpi=50)
    os.system("convert -delay 8 ./data/shape-"+f_shape+"/plots/-tmp_*" \
    +"-"+str(protein.protein)+"-"+f_shape+"-"+f_param1+"-"+f_param2 \
    +"-"+f_param3+"-"+f_param4+"-"+f_param5+".png ./data/shape-"+f_shape \
    +"/plots/density_movie-"+str(protein.protein)+"-"+f_shape \
    +"-"+f_param1+"-"+f_param2+"-"+f_param3+"-"+f_param4+"-"+f_param5+".gif")
    os.system("rm -f ./data/shape-"+f_shape+"/plots/"+str(protein.protein)+"-tmp_*.png")
    return 0

print "Generating flourescent Min D tagging plot:"
contourplt(nflD)
print "Done! Generating flourescent Min E tagging plot:"
contourplt(nflE)
print "Generating nATP plot:"
contourplt(nATP)
print "Done! Generating nE plot:"
contourplt(nE)
print "Done! Generating nADP plot:"
contourplt(nADP)
print "Done! Generating nD plot:"
contourplt(Nd)
print "Finished generating plots."
