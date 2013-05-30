from __future__ import division
from numpy import *
from matplotlib import *
from pylab import *
import sys
import glob

dat_filenames = []
for fn in glob.iglob('../data/shape-c/membrane_files/*.dat'):
    dat_filenames.append(fn)

print dat_filenames
a = loadtxt(dat_filenames[0])
print a
