import numpy as np
import sys
import glob

f_shape = sys.argv[1]
f_param1 = sys.argv[2]
f_param2 = sys.argv[3]
f_param3 = sys.argv[4]
f_param4 = sys.argv[5]
f_param5 = sys.argv[6]

#get_filenames method sloppy but it works.

if "-hires" in sys.argv:
    dx=.05
else:
    dx=.15

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
        if "-hires" in sys.argv:
            if "-slice" in sys.argv:
                for fn in glob.iglob('./data/shape-'+f_shape+'/hires-'+protein+'-'+f_shape+'-'+f_param1+'-'+f_param2+'-'+f_param3+'-'+f_param4+'-'+f_param5+'*.dat'):
                    dat_filenames.append(fn)
            else:
                for fn in glob.iglob('./data/shape-'+f_shape+'/hires-m-'+protein+'-'+f_shape+'-'+f_param1+'-'+f_param2+'-'+f_param3+'-'+f_param4+'-'+f_param5+'*.dat'):
                    dat_filenames.append(fn)
        else:
            if "-slice" in sys.argv:
                for fn in glob.iglob('./data/shape-'+f_shape+'/'+protein+'-'+f_shape+'-'+f_param1+'-'+f_param2+'-'+f_param3+'-'+f_param4+'-'+f_param5+'*.dat'):
                    dat_filenames.append(fn)
            else:
                for fn in glob.iglob('./data/shape-'+f_shape+'/m-'+protein+'-'+f_shape+'-'+f_param1+'-'+f_param2+'-'+f_param3+'-'+f_param4+'-'+f_param5+'*.dat'):
                    dat_filenames.append(fn)
        dat_filenames = sorted(dat_filenames)
        i=0
        while i<round(len(dat_filenames)/10):
            dat_filenames.pop(0)
            i+=1
        if (dat_filenames == []):
            print "File loading error: filename list is empty."
            exit(1)
        else:
            return dat_filenames
