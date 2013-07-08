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

for p in [natp, ne, nadp, nd]:
    plt.figure()
    plt.title(str(p.protein))
    plt.ylabel("average protein density (proteins per micron^3)")
    plt.xlabel("area rating")
    plt.scatter(splitdata(p)[0], splitdata(p)[1])
    print './data/shape-'+f_shape+'/plots/avg_density-'+str(p.protein)+'-'+f_shape+'-'+f_param1+'-'+f_param2+'-'+f_param3+'-'+f_param4+'-'+f_param5+'.pdf'
    plt.savefig('./data/shape-'+f_shape+'/plots/avg_density-'+str(p.protein)+'-'+f_shape+'-'+f_param1+'-'+f_param2+'-'+f_param3+'-'+f_param4+'-'+f_param5+'.pdf')
