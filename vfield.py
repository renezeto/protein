from __future__ import division
from numpy import *
from matplotlib import *
from pylab import *
from random import *
from inspect import *

dx = .05 #microns
data = zeros((20,20),dtype='float') #loadtxt('filename.dat',delimiter=',')
data_shape = [data.shape[n] for n in range(len(data.shape))]
data_size = [data.shape[n]*dx for n in range(len(data.shape))]
axis = [arange(0,data_size[n],dx) for n in range(len(data.shape))]



#generate fake 2d data - works for 3d too.
for i in range(data_shape[0]):
	for j in range(data_shape[1]):
		data[i][j] = axis[1][i] + axis[0][j]


def gradient(data):
	if len(data.shape) == 2:
		gradStore =  [[(0.,0.) for y in range((data_shape[0]-2))] for z in range((data_shape[1]-2))]
		for i in range(2,data_shape[0]-2):
			for j in range(2,data_shape[1]-2):
				partial_x1 = (data[i][j+1] - data[i][j-1])/(2*dx)
				partial_x2 = (data[i+1][j] - data[i-1][j])/(2*dx)
				gradStore[i][j] = (partial_x1, partial_x2)
	else:
		gradStore =  [[[(0.,0.,0.) for x in range((data_shape[0]-2))] for y in range((data_shape[1]-2))] for z in range((data_shape[2]-2))]
		for i in range(2,data_shape[0]-2):
			for j in range(2,data_shape[1]-2):
				for k in range(2,data_shape[0]-2):
					partial_x1 = (data[i][j][k+1] - data[i][j][k-1])/(2*dx)
					partial_x2 = (data[i][j+1][k] - data[i][j-1][k])/(2*dx)
					partial_x3 = (data[i+1][j][k] - data[i-1][j][k])/(2*dx)
					gradStore[i][j][k] = (partial_x1, partial_x2, partial_x3)
	return gradStore

def mapping(data, page): #kind of like meshgrid(X,Y) for separating the partials.
	if len(data.shape) == 2:
		x2_component = zeros_like(data)
		x3_component = zeros_like(data)
		gradf = gradient(data)
		for i in range(data_shape[1]-2):
			for j in range(data_shape[0]-2):
				x2_component[i][j] = gradf[i][j][0]
				x3_component[i][j] = gradf[i][j][1]
	else:
		x2_component = zeros_like(data[page])
		x3_component = zeros_like(data[page])
		gradf = gradient(data)
		for i in range(data_shape[1]-2):
    			for j in range(data_shape[0]-2):
        			x2_component[i][j] = gradf[page][i][j][0]
				x3_component[i][j] = gradf[page][i][j][1]
	return x2_component, x3_component


X1, X2 = meshgrid(axis[0],axis[1])
quiver(X1,X2,mapping(data,1)[0],mapping(data,1)[1])
show()

