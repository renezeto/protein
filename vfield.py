from __future__ import division
from numpy import *
from matplotlib import *
from random import *

#this is designed for 3 dimensional data.
data = zeros((20,20,80),dtype='float') #loadtxt('filename.dat',delimiter=',')
length_x1 = 1 #import from filename string
length_x2 = 1 #import from filename string
length_x3 = 4 #import from filename string
dx = .05 #microns
n = data.shape[0]
m = data.shape[1]
o = data.shape[2] #n pages by m rows by o columns (3 dimensional array)

x1 = arange(0,length_x1,dx)
x2 = arange(0,length_x2,dx)
x3 = arange(0,length_x3,dx)

for i in range(n):
    for j in range(m):
        for k in range(o):
            data[i][j][k] = random() #warning: not continuous
data1 = [[[0.]*(o-2)]*(m-2)]*(n-2)
print data1[1][1][1]
data1[1][1][1]= 2
print data1[1][1][1],data1[2][1][1],data1[3][1][1]
'''
gradStore =  [[[0.]*(o-2)]*(m-2)]*(n-2) #zeros_like(data)
    
def gradient(data):
    for i in range(1,n-2):
        for j in range(1,m-2):
            for k in range(1,o-2):
                partial_x1 = (data[i][j][k+1] - data[i][j][k-1])/(2*dx)
                partial_x2 = (data[i][j+1][k] - data[i][j-1][k])/(2*dx)
                partial_x3 = (data[i+1][j][k] - data[i-1][j][k])/(2*dx)
                gradStore[i][j][k] = (partial_x1, partial_x2, partial_x3)
                if (j==1 and i ==1 and k==40):
                    print gradStore[i][j][k]
                    print "here"
    return 0
#print "here"
gradient(data)
for x in range(1,m-2):
    print gradStore[1][x][40]
    print "heretwo"
'''
#problems:
# - can't store a vector in an array position
# - partials will have boundary value problems. fix later with real data.

#to do:
# - animate: create a function to run gradient on each freshly loaded data
#file, and redraw.
