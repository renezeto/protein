import numpy as np
import sys
import glob

f_shape = sys.argv[1]
f_param1 = sys.argv[2]
f_param2 = sys.argv[3]
f_param3 = sys.argv[4]
f_param4 = sys.argv[5]
f_param5 = sys.argv[6]

class data(object):
    def __init__(self,protein):
        self.protein = protein
        self.filenames = get_filenames(protein)
        self.tsteps = 

    def get_filenames(protein):
        dat_filenames = []
        for fn in glob.iglob('./data/shape-'+f_shape+'/m_natp-'+f_shape+'-'+f_param1+'-'+f_param2+'-'+f_param3+'-'+f_param4+'-'+f_param5+'*.dat'):
            dat_filenames.append(fn)
        dat_filenames = sorted(dat_filenames)
        return 
