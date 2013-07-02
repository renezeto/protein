import numpy as np
import sys
import glob

#test

f_shape = sys.argv[1]
f_param1 = sys.argv[2]
f_param2 = sys.argv[3]
f_param3 = sys.argv[4]
f_param4 = sys.argv[5]
f_param5 = sys.argv[6]

dx=.05

class data(object):
    def __init__(self,protein):
        self.protein = protein
        self.filenames = self.get_filenames(protein)
        self.tsteps = len(self.filenames)
        self.dataset = np.array([np.loadtxt(self.filenames[i]) for i in range(self.tsteps)])
        self.datashape = [self.dataset[0].shape[n] for n in range(len(self.dataset[0].shape))]
        self.datasize = [self.datashape[n]*dx for n in range(len(self.datashape))]
        self.axes = [np.arange(0,self.datasize[n],dx) for n in range(len(self.datashape))]

    @staticmethod
    def get_filenames(protein):
        dat_filenames = []
        for fn in glob.iglob('./data/shape-'+f_shape+'/m_'+protein+'-'+f_shape+'-'+f_param1+'-'+f_param2+'-'+f_param3+'-'+f_param4+'-'+f_param5+'*.dat'):
            dat_filenames.append(fn)
        dat_filenames = sorted(dat_filenames)
        if (dat_filenames == []):
            print "File loading error: filename list is empty."
            exit(1)
        else:
            dat_filenames.pop(0)
            return dat_filenames
