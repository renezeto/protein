import pylab, sys, math
from pylab import *
import matplotlib.pyplot as plt
import numpy as np
import sys


fname = "shape-" + str(sys.argv[1]) + "/area_rating_two-" + str(sys.argv[2]) + "-" + str(sys.argv[3]) + "-" + str(sys.argv[4]) + "-" + str(sys.argv[5]) + "-" + str(sys.argv[6]) + ".dat"
print fname
A = np.loadtxt(fname, dtype = float)
8/.05 +4
if sys.argv[1] == "randst" or sys.argv[1] == "triangle":
    max_y = float(sys.argv[4])
    max_z = float(sys.argv[3])
if sys.argv[1] == "p":
    max_y = 2*float(sys.argv[3])
    max_z = float(sys.argv[2])+2*float(sys.argv[3])
print "max_y " + str(max_y)
width_y = int(max_y/0.05+4.0)
width_z = int(max_z/0.05+4.0)
print len(A[:,0])
print width_y
print width_z
print "width_y*width_z = " + str(width_y*width_z)

y = np.zeros(width_y)
z = np.zeros(width_z)
for i in range(width_y):
    y[i] = i*.05
for i in range(width_z):
    z[i] = i*.05

Y,Z = meshgrid(y,z) #2d array of zeros
area_rating = np.zeros_like(Y)

def unzip(A):
    for i in range(A.shape[0]):
        y = (round(A[i][0]/.05))
        z = (round(A[i][1]/.05))
        area_rating[y][z] = A[i][2]
    return area_rating
max_rating = 0
min_rating = 5.0
for i in range(A.shape[0]):
    if (A[i][2] > max_rating):
        max_rating = A[i][2]
    if (A[i][2] < min_rating and A[i][2] != 0.0):
        min_rating = A[i][2]

area_rating = unzip(A)
mylevel = np.arange(0,0.1,.001)
contourf(Y,Z,area_rating,cmap=plt.cm.jet,levels=mylevel)
colorbar()
show()
