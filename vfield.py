from __future__ import division
from numpy import *
from matplotlib import *
from pylab import *
import sys
import glob
from matplotlib.widgets import Slider, RadioButtons


f_shape = sys.argv[1]
f_param1 = sys.argv[2]
f_param2 = sys.argv[3]
f_param3 = sys.argv[4]
f_param4 = sys.argv[5]
#f_density = sys.argv[6]
f_seed = sys.argv[6]

dat_filenames = []
for fn in glob.iglob('shape-'+f_shape+'/natp-'+f_shape+'-'+f_param1+'-'+f_param2+'-'+f_param3+'-'+f_param4+'-'+f_seed+'*.dat'):
        dat_filenames.append(fn)
t_steps = len(dat_filenames)
data_natp_set = array([np.loadtxt(dat_filenames[i]) for i in range(t_steps)])

print data_natp_set[0]

dx = .05 #microns
data = data_natp_set[0] #retrieves data format.
data_shape = [data.shape[n] for n in range(len(data.shape))]
data_size = [data.shape[n]*dx for n in range(len(data.shape))]
axis = [arange(0,data_size[n],dx) for n in range(len(data.shape))]


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
                                if (abs(x3_component[i][j])>1000 or abs(x2_component[i][j])>1000):
                                        x2_component[i][j]=0
                                        x3_component[i][j]=0
                                print x2_component[i][j]
                return x2_component, x3_component

        else:
		x2_component = zeros_like(data[page])
		x3_component = zeros_like(data[page])
		gradf = gradient(data)
		for i in range(data_shape[1]-2):
    			for j in range(data_shape[0]-2):
        			x2_component[i][j] = gradf[page][i][j][0]
				x3_component[i][j] = gradf[page][i][j][1]
        return x2_component, x3_component

t0 = 0
x0 = 0
Q = quiver(axis[0],axis[1],mapping(data_natp_set[0],x0)[0],mapping(data_natp_set[0],x0)[1])

#this is the time slider
t0_ax = axes([0.25, 0, 0.5, 0.03], axisbg='slategray')
t0_slider = Slider(t0_ax, 'time step', 0, t_steps, valinit = 0, valfmt='%0.0f')

def update_time(val):
        global t0
        t0 = round(t0_slider.val)
        Q.set_UVC(mapping(data_natp_set[t0],x0)[0],mapping(data_natp_set[t0],x0)[1])

t0_slider.on_changed(update_time)


#this would be used for going through 3d data along the x axis
depth_ax = axes([0.25, 0.03, 0.5, 0.03], axisbg='slategray')
depth_slider = Slider(depth_ax, 'depth', 0, data_shape[0], valinit = 0, valfmt='%0.1f')

def update_depth(val):
        global x0
        x0 = depth_slider.val
        vector_plot()

depth_slider.on_changed(update_depth)

show()

