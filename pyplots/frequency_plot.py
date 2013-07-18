from __future__ import division
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab
import os
import sys
import time
import file_loader as load


dx = load.dx
#wip

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

f_shape = sys.argv[1]
f_param1 = sys.argv[2]
f_param2 = sys.argv[3]
f_param3 = sys.argv[4]
f_param4 = sys.argv[5]
f_param5 = sys.argv[6]

par1 = str(int(10*(float(f_param1))))
par2 = str(int(10*(float(f_param2))))
par3 = str(int(10*(float(f_param3))))
par4 = str(int(10*(float(f_param4))))
par5 = str(int(10*(float(f_param5))))


nflE = load.data(protein="nflE")
nflD = load.data(protein="nflD")
nATP = load.data(protein="nATP")
nE = load.data(protein="nE")
nADP = load.data(protein="nADP")
Nd = load.data(protein="Nd")
Nde = load.data(protein="Nde")

mid = False
if "-mid" in sys.argv:
    mid = True
    print "The -mid flag was set so we'll add the middle point data to the plot"

norm = False
if "-norm" in sys.argv:
    norm = True
    print "We are normalizing the plots so the highest density for each protein will be one"


proteins = [nflD, nATP, nE, nADP, Nd, Nde]#normally have nflD and nflE too
proteins_strings = ["nATP", "nE", "nADP", "Nd", "Nde"] #normally "nflD" and "nflE"
if mid:
    proteins_strings = ["nflD", "middle-nflD", "nATP", "middle-nATP","nE", "middle-nE","nADP", "middle-nADP"] #normally "nflD" and "nflE"



for p in proteins:
    for i in p.dataset:
        if (globalmax(i) != "empty"):
            p.pointInfo = (globalmax(i)[0][0], globalmax(i)[0][1])
            break
    print p.pointInfo

mid_z = 2 + round((float(f_param1)-dx)/dx)
mid_y = 2 + round((float(f_param2)-dx)/dx)


for p in proteins:
    p.densityInfo = []
    p.densityInfoTwo = []
    for file in p.dataset:
        p.densityInfo += [ file[p.pointInfo[0]][p.pointInfo[1]] ]
        p.densityInfoTwo += [ file[mid_y][mid_z] ]
    time = np.arange(0,len(p.densityInfo)/2,.5)
    plt.figure()
    plt.plot(time,p.densityInfo)
    plt.plot(time,p.densityInfoTwo)
    plt.title('Max point = '+str(p.pointInfo[0]*dx)+','+str(p.pointInfo[0]*dx)
              +'   Middle point = '+str(mid_y*dx)+','+str(mid_z*dx))

    plt.savefig('./data/shape-'+f_shape+'/plots/'+load.hires_str+load.m_str+'frequency_plot'+str(p.protein)+'-'+f_shape+'-'+par1+'-'+par2+'-'+par3+'-'+par4+'-'+par5+'.pdf')


plt.figure()
for p in proteins:
    pointInfoWall = Nd.pointInfo
    if p!=Nd and p!=Nde:
        pointInfoWall = (int(Nd.pointInfo[0])-1,int(Nd.pointInfo[1]))
    p.densityInfo = []
    p.densityInfoTwo = []
    for file in p.dataset:
        p.densityInfo += [ file[pointInfoWall[0]][pointInfoWall[1]] ]
        p.densityInfoTwo += [ file[mid_y][mid_z] ]
    if norm:
        for i in np.arange(len(p.densityInfo)):
            p.densityInfo[i]=p.densityInfo[i]/max(p.densityInfo)
    time = np.arange(0,len(p.densityInfo)/2,.5)
    if mid:
        if p!=Nd and p!=Nde:
            plt.plot(time,p.densityInfo)
            plt.plot(time,p.densityInfoTwo)
    else:
        plt.plot(time,p.densityInfo)
if norm:
    plt.ylim(0,1.5)
plt.title('Max point = '+str(pointInfoWall[0]*dx)+','+str(pointInfoWall[0]*dx))
plt.legend(proteins_strings,loc='best')
#pylab.legend(loc='best')
if norm:
    plt.savefig('./data/shape-'+f_shape+'/plots/'+load.hires_str+load.m_str+'frequency_plot_all_norm-'+f_shape+'-'+par1+'-'+par2+'-'+par3+'-'+par4+'-'+par5+'.pdf')
elif mid:
    plt.savefig('./data/shape-'+f_shape+'/plots/'+load.hires_str+load.m_str+'frequency_plot_all_middle-'+f_shape+'-'+par1+'-'+par2+'-'+par3+'-'+par4+'-'+par5+'.pdf')
else:
    plt.savefig('./data/shape-'+f_shape+'/plots/'+load.hires_str+load.m_str+'frequency_plot_all-'+f_shape+'-'+par1+'-'+par2+'-'+par3+'-'+par4+'-'+par5+'.pdf')
