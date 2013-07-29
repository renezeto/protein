from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import sys
import file_loader as load

f_shape = sys.argv[1]
f_param1 = sys.argv[2]
f_param2 = sys.argv[3]
f_param3 = sys.argv[4]
f_param4 = sys.argv[5]
f_param5 = sys.argv[6]

natp = load.data(protein="natp")
ne = load.data(protein="ne")
nadp = load.data(protein="nadp")
nd = load.data(protein="nd")

areaFname = "data/shape-" + f_shape +"/area_rating_two-" +f_param1+ "-" +f_param2 + "-" +f_param3+ "-" +f_param4 + "-" +f_param5+ ".dat"
areaFile = np.loadtxt(areaFname,dtype=float)


def splitdata(protein):
    areaFname = "data/shape-" + f_shape +"/area_rating_two-" +f_param1+ "-" +f_param2 + "-" +f_param3+ "-" +f_param4 + "-" +f_param5+ ".dat"
    areaFile = np.loadtxt(areaFname,dtype=float)
    
    proteinFname = "data/shape-" + f_shape +"/"+str(protein.protein)+"-avg_density-" +f_shape+ "-" +f_param1+ "-" +f_param2 + "-" +f_param3+ "-" +f_param4 + "-" +f_param5+ ".dat"
    proteinFile = np.loadtxt(proteinFname,dtype=float)
    
    avgProtein = []
    areaRating = []

    for i in range(len(areaFile)):
        if not (proteinFile[i][3]==0 and areaFile[i][3]==0):
            avgProtein += [proteinFile[i][3]]
            areaRating += [areaFile[i][3]]
    return areaRating, avgProtein

y = np.zeros(natp.datasize[0]/load.dx)
z = np.zeros(natp.datasize[1]/load.dx)

Y,Z = np.meshgrid(y,z) 
area_rating = np.zeros_like(Y)

print area_rating.shape
def unzip(areaFile):
    for i in range(areaFile.shape[0]):
        y = (round(areaFile[i][0])/load.dx)
        z = (round(areaFile[i][1])/load.dx)
        area_rating[y][z] = areaFile[i][3]
    return area_rating

plt.figure()
plt.contourf(y,z,unzip(areaFile))
plt.show()
