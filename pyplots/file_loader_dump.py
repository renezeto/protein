import numpy as np
import sys
import glob

#get the shape an physical parameters of the cell (last argument is density lopsidedness)
f_shape = sys.argv[1]
f_param1 = sys.argv[2]
f_param2 = sys.argv[3]
f_param3 = sys.argv[4]
f_param4 = sys.argv[5]
f_param5 = sys.argv[6]

dx=.15

#for the specific file reading, these are the protein strings
proteinList = ["D_nATP","D_nADP","D_E_NDE","E_nE","D_ND"]

#initialize as empty lists incase arguments not present in sys.argv
hires_str=""
slice_str=""

#change them to print the correct file names if arguments present
if "-hires" in sys.argv:
    dx=.05
    hires_str="hires-"
if "-slice" in sys.argv:
    slice_str="slice-"

#data object, to create in another plot: proteinname=data(protein="proteinname")
class data(object):
    def __init__(self,protein):
        self.protein = protein
        self.filenames = self.get_filenames(protein)
        self.tsteps = len(self.filenames)
        self.dataset = np.array([np.loadtxt(file) for file in self.filenames])
        self.datashape = self.dataset[0].shape
        self.axes = [[i*dx for i in range(self.datashape[1])],[j*dx for j in range(self.datashape[0])]]

#loads the files for creating the data object.
    @staticmethod
    def get_filenames(protein):
        dat_filenames = []
        for fn in glob.iglob("./data/shape-"+f_shape+"/"+hires_str+slice_str+protein+"-"+f_shape+"-"+f_param1+"-"+f_param2+"-"+f_param3+"-"+f_param4+"-"+f_param5+"*.dat"):
            dat_filenames.append(fn)
        dat_filenames = sorted(dat_filenames)
        i=0 #pop the first 10% of file names to let things equillibriate a bit
        while i<round(len(dat_filenames)/10):
            dat_filenames.pop(0)
            i+=1
        if (dat_filenames == []):
            print "File loading error: filename list is empty."
            print "./data/shape-"+f_shape+"/"+hires_str+slice_str+protein+"-"+f_shape+"-"+f_param1+"-"+f_param2+"-"+f_param3+"-"+f_param4+"-"+f_param5+"*.dat"
            exit(1)
        else:
            return dat_filenames

#function for easier plot name printing. probably should be renamed itself.
def print_string(plot_name,p):
    arg = [str(int(100*(float(i)))) for i in sys.argv[2:7]]
    filename = "./data/shape-%s/plots/%s%s%s-%s-%s-%s-%s-%s-%s-%s.pdf"%(f_shape, hires_str, slice_str, plot_name, p.protein, f_shape, arg[0], arg[1], arg[2], arg[3], arg[4])
    return filename
