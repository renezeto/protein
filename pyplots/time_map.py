from __future__ import division
from matplotlib.pyplot import *
from numpy import *
from time import *
from random import *
import glob
from file_loader import *

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
print "time_map plot generated."
