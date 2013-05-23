import pylab, sys, math
from pylab import *
import matplotlib.pyplot as plt
import numpy as np
import sys
#print sys.argv[1]
#shape-randst/area_rating-0.14-6.00-8.00-98.00.dat

fname = "shape-" + str(sys.argv[1]) + "/area_rating-" + str(sys.argv[2]) + "-" + str(sys.argv[3]) + "-" + str(sys.argv[4]) + "-" + str(sys.argv[5]) + ".dat"
print fname
A = np.loadtxt(fname, dtype = float)
8/.05 +4
max_y = float(sys.argv[4])
max_z = float(sys.argv[3])
print "max_z" + str(max_z)
width_y = int(max_y/0.05+4.0)
width_z = int(max_z/0.05+4.0)
print len(A[:,0])
print width_y
print width_z
print "width_y*width_z = " + str(width_y*width_z)
Y = reshape(A[:,0],(width_y,width_z))
Z = reshape(A[:,1],(width_y,width_z))
area_rating = reshape(A[:,2],(width_y,width_z))
mylevel = np.arange(1.3,2.3,.001)
contour(Y,Z,area_rating,cmap=plt.cm.jet,levels=mylevel)
colorbar()
show()
