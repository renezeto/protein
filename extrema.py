from __future__ import division
from numpy import *
from matplotlib import *
from pylab import *
import sys
import glob
from matplotlib.widgets import Slider, RadioButtons
'''
dat_filenames = []
for fn in glob.iglob('shape-'+f_shape+'/natp-'+f_shape+'-'+f_param1+'-'+f_param2+'-'+f_param3+'-'+f_param4+'-'+f_seed+'*.dat'):
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

#between tstep-1 and tstep, track the location of the maximum. first, write something to find the maximum.
'''

def maxloc(array):
    index_pairs = list(enumerate(array))
    maxinfo=index_pairs[0]
    for i in range(len(array)):
        if index_pairs[i][1] > maxinfo[1]:
            maxinfo = index_pairs[i]
    return maxinfo

a = [1,2,3,4,5,6]
b = [2,4,12,8,10,12]
c = [a,b]

def extrema(page):
    maxstore = [(0.,0.) for i in range(len(page))]
    for i in range(len(page)):
        maxstore[i] = (i,maxloc(page[i])[0])
    return maxstore
