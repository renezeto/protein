from __future__ import division
from matplotlib.pyplot import *
from numpy import *
from time import *
from random import *
import glob

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

def average_location(): #this function only works with a plane of data right now
    tsum = zeros_like(data_natp_set[0])
    for t in range(t_steps):
        tsum += data_natp_set[t]
    return tsum/t_steps

contourf(axis[0],axis[1],average_location(),500)
show()