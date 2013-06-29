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
f_param5 = sys.argv[6]

dat_filenames = []
for fn in glob.iglob('./data/shape-'+f_shape+'/m_natp-'+f_shape+'-'+f_param1+'-'+f_param2+'-'+f_param3+'-'+f_param4+'-'+f_param5+'*.dat'):
        dat_filenames.append(fn)
dat_filenames = sorted(dat_filenames)
dat_filenames.pop(0)

t_steps = len(dat_filenames)
data_natp_set = array([np.loadtxt(dat_filenames[i]) for i in range(t_steps)])
dx = .05 #microns

data = data_natp_set[0] #retrieves data format.
data_shape = [data.shape[n] for n in range(len(data.shape))]
data_size = [data.shape[n]*dx for n in range(len(data.shape))]
axis = [arange(0,data_size[n],dx) for n in range(len(data.shape))]

def average_location(dataset):
    tsum = zeros_like(dataset[0])
    for t in range(t_steps):
            tsum += dataset[t]
    return tsum/(t_steps)

#matshow(average_location())
figure()
contourf(axis[1],axis[0],average_location(data_natp_set),500)
colorbar()
savefig('./data/shape-'+f_shape+'/plots/time_map-natp-'+f_shape+'-'+f_param1+'-'+f_param2+'-'+f_param3+'-'+f_param4+'-'+f_param5+'.pdf')
#this only does it for NATP right now. will need to automate for other protein types later.
print "time_map plot generated."
