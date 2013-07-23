import pylab, sys, math
from pylab import *
import matplotlib.pyplot as plt
import numpy as np
import sys
import file_loader as load

dx = load.dx

f_shape = sys.argv[1]
f_param1 = sys.argv[2]
f_param2 = sys.argv[3]
f_param3 = sys.argv[4]
f_param4 = sys.argv[5]
f_param5 = sys.argv[6]

fname = "data/shape-" + str(sys.argv[1]) + "/area_rating-" + str(sys.argv[2]) + "-" + str(sys.argv[3]) + "-" + str(sys.argv[4]) + "-" + str(sys.argv[5]) + "-" + str(sys.argv[6]) + ".dat"
#print fname
A = np.loadtxt(fname, dtype = float)
8/.05 +4
#print A[:,0]
print "HERE"

width_y = 0
for i in range(len(A[:,0])):
    if A[i,1]==0:
        width_y +=1
width_z = len(A[:,0])/width_y

print width_y
print width_z
print len(A[:,0])
print "This place"

y = np.zeros(width_y)
z = np.zeros(width_z)
for i in range(width_y):
    y[i] = i*.05
for i in range(width_z):
    z[i] = i*.05

print "widths of Grid"
print width_y*dx, width_y
print width_z*dx, width_z
print "length of y"
print len(y)
Y,Z = meshgrid(y,z) #2d array of zeros
area_rating = np.zeros_like(Y)
count=0
print "HEEEEEEEEEER"
print len(Z)

for i in range(len(A[:,0])):
    #print A[i,3]
    count+=1
print count
print "Now here"

def unzip(A):
    print len(A[:,0]-1)
    print A[2,3]
    for i in range(len(A[:,0])-1):
        y = (round(A[i,1]/dx))
        z = (round(A[i,2]/dx))
        area_rating[y][z] = A[i,3]
        #print area_rating[y][z], i
    return area_rating
max_rating = 0
min_rating = 5.0
for i in range(len(A[:,0])):
    if (A[i,3] > max_rating):
        max_rating = A[i,3]
    if (A[i,3] < min_rating):# and A[i,3] != 0.0):
        min_rating = A[i,3]
print min_rating

area_rating = unzip(A)
mylevel = np.arange(500,800,1)
contourf(Y,Z,area_rating,cmap=plt.cm.jet,levels=mylevel)
colorbar()
savefig('./data/shape-'+f_shape+'/plots/show_area_rating-natp-'+f_shape+'-'+f_param1+'-'+f_param2+'-'+f_param3+'-'+f_param4+'-'+f_param5+'.pdf')
#this only does it for NATP right now. will need to automate for other protein types later.
print "show_area_rating plot generated."

#show()
