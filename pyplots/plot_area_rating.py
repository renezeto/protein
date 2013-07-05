from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import sys
import file_loader

f_shape = sys.argv[1]
f_param1 = sys.argv[2]
f_param2 = sys.argv[3]
f_param3 = sys.argv[4]
f_param4 = sys.argv[5]
f_param5 = sys.argv[6]

areaFname = "data/shape-" + f_shape +"/area_rating_two-" +f_param1+ "-" +f_param2 + "-" +f_param3+ "-" +f_param4 + "-" +f_param5+ ".dat"
areaFile = np.loadtxt(areaFname,dtype=float)

proteinFname = "data/shape-" + f_shape +"/avg_density-" +f_shape+ "-" +f_param1+ "-" +f_param2 + "-" +f_param3+ "-" +f_param4 + "-" +f_param5+ ".out" #change to .dat
proteinFile = np.loadtxt(proteinFname,dtype=float)

avgProtein = []
areaRating = []
for i in range(areaFile.shape[0]):
    areaRating += [areaFile[i][3]]
for i in range(avgProtein.shape[0]):
    avgProtein += [proteinFile[i][3]]
print areaRating
