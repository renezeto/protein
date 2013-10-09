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

dx=.05

#for the specific file reading, these are the protein strings
proteinList = ["D_nATP","D_nADP","D_E_NDE","E_nE","D_ND","NflD","NflE"]

#initialize as empty lists incase arguments not present in sys.argv
hires_str=""
slice_str=""
debug_str=""

#change them to print the correct file names if arguments present
if "-hires" in sys.argv:
    dx=.025
    hires_str="hires-"
if "-slice" in sys.argv:
    slice_str="slice-"
if "-debug" in sys.argv:
    debug_str="debug-"

expected_arg_number = 7
for arg in sys.argv:
    if "-" in arg:
        expected_arg_number+=1

if len(sys.argv)>expected_arg_number:
    f_param6 = sys.argv[7]
    f_param7 = sys.argv[8]

filename_tuple = (f_shape,hires_str,slice_str,f_shape,f_param1,f_param2,f_param3,f_param4,f_param5)

#data object, to create in another plot: proteinname=data(protein="proteinname")
class data(object):
    def __init__(self,protein):
        self.protein = protein
        self.filenames = self.get_filenames(protein,0,0)
        self.tsteps = len(self.filenames)
        self.dataset = np.array([np.loadtxt(file) for file in self.filenames])
        self.datashape = self.dataset[0].shape
        self.axes = [[i*dx for i in range(self.datashape[1])],[j*dx for j in range(self.datashape[0])]]
    # def __init__(self,protein,start_time,end_time):
    #     self.protein = protein
    #     self.filenames = self.get_filenames(protein,start_time,end_time)
    #     self.tsteps = len(self.filenames)
    #     self.dataset = np.array([np.loadtxt(file) for file in self.filenames])
    #     self.datashape = self.dataset[0].shape
    #     self.axes = [[i*dx for i in range(self.datashape[1])],[j*dx for j in range(self.datashape[0])]]

#loads the files for creating the data object.
    @staticmethod
    def get_filenames(protein,start_time,end_time):
        dat_filenames = []
        if end_time == 0:
            for fn in glob.iglob("./data/shape-%s/%s%s%s%s-%s-%s-%s-%s-%s-%s-*.dat"%(f_shape,debug_str,hires_str,slice_str,protein,f_shape,f_param1,f_param2,f_param3,f_param4,f_param5)):
                dat_filenames.append(fn)
        if end_time != 0:
            for i in np.arange(start_time,end_time,2.5):
                file_num = round(i/2.5)
                dat_filenames.append("./data/shape-%s/%s%s%s%s-%s-%s-%s-%s-%s-%s-%d.dat"%
                                     (f_shape,debug_str,hires_str,slice_str,protein,f_shape,f_param1,f_param2,f_param3,f_param4,f_param5,file_num))
        dat_filenames = sorted(dat_filenames)
        i=0 #pop the first 10% of file names to let things equillibriate a bit
        while i<round(len(dat_filenames)/10):
            dat_filenames.pop(0)
            i+=1
        if (dat_filenames == []):
            print "File loading error: filename list is empty."
            print "./data/shape-%s/%s%s%s%s-%s-%s-%s-%s-%s-%s-*.dat"%(f_shape,debug_str,hires_str,slice_str,protein,f_shape,f_param1,f_param2,f_param3,f_param4,f_param5)
            exit(1)
        else:
            return dat_filenames

#function for easier plot name printing. probably should be renamed itself.
def print_string(plot_name,p):
    arg = [str(int(100*(float(i)))) for i in sys.argv[2:7]]
    filename = "./data/shape-%s/plots/%s%s%s%s-%s-%s-%s-%s-%s-%s-%s.pdf"%(f_shape, debug_str, hires_str, slice_str, plot_name, p, f_shape, arg[0], arg[1], arg[2], arg[3], arg[4])
    return filename

#"./data/shape-"+f_shape+"/"debug_str+hires_str+slice_str+protein+"-"+f_shape+"-"+f_param1+"-"+f_param2+"-"+f_param3+"-"+f_param4+"-"+f_param5+"*.dat"
