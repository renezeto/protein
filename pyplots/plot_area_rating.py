from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import file_loader

f_shape = sys.argv[1]
f_param1 = sys.argv[2]
f_param2 = sys.argv[3]
f_param3 = sys.argv[4]
f_param4 = sys.argv[5]
f_param5 = sys.argv[6]

fname = "data/shape-" + f_shape +"area_rating_two-" +f_param1+ "-" +f_param2 + "-" +f_param3+ "-" +f_param4 + "-" +f_param5+ ".dat"
areaFile = np.loadtxt('fname',dtype=float)

x = np.array([],dtype='float')
y = np.array([],dtype='float')
z = np.array([],dtype='float')
areaRating=[]
for i in range(areaFile.shape[0]):
    x = append(x,(round(areaFile[i][0]/file_loader.dx)))
    y = append(y,(round(areaFile[i][1]/file_loader.dx)))
    z = append(z,(round(areaFile[i][2]/file_loader.dx)))
    areaRating[y][z] += [areaFile[i][3]]
areaRating = np.array(areaRating,dtype='float')
