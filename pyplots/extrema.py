from __future__ import division
from numpy import *
from matplotlib import *
from pylab import *
import sys
import glob

f_shape = sys.argv[1]
f_param1 = sys.argv[2]
f_param2 = sys.argv[3]
f_param3 = sys.argv[4]
f_param4 = sys.argv[5]
f_param5 = sys.argv[6]

dat_filenames = []
for fn in glob.iglob('./data/shape-'+f_shape+'/natp-'+f_shape+'-'+f_param1+'-'+f_param2+'-'+f_param3+'-'+f_param4+'-'+f_param5+'*.dat'):
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


cell_membrane = loadtxt('./data/shape-'+f_shape+'/membrane_files/membrane-'+f_param1+'-'+f_param2+'-'+f_param3+'-'+f_param4+'-'+f_param5+'.dat')

def unzip_membrane(memdat):
    x = []
    y = []
    for i in range(len(memdat)):
        for j in range(len(memdat[i])):
            if memdat[i][j] != 0:
                x += [i]
                y += [j]
    return x,y

def localmax(page):
    lmax = []
    for i in range(1,len(page)-1):
        index_pairs = list(enumerate(page[i]))
        for j in range(1,len(page[0])-1):
            if (page[i][j] > page[i-1][j-1]) and (page[i][j] > page[i-1][j]) and (page[i][j] > page[i-1][j+1]) and (page[i][j] > page[i][j-1]) and (page[i][j] > page[i][j+1]) and (page[i][j] > page[i+1][j-1]) and (page[i][j] > page[i+1][j]) and (page[i][j] > page[i+1][j+1]):
                lmax += [(i,index_pairs[j][0],index_pairs[j][1])]
    if (lmax == []):
        return "empty"
    return lmax

def globalmax(page):
    gmax = [(0.,0.,0.)]
    lmax = localmax(page)
    if lmax != "empty":
        for i in range(len(lmax)):
            if (lmax[i][2] > gmax[0][2]):
                gmax = [lmax[i]]
            elif (lmax[i][2] == gmax[0][2]):
                gmax += [lmax[i]]
                if gmax[0] == 0.:
                    gmax.pop(0) #need to check this 
    else:
        return "empty"
    if len(gmax)>1:
        gmax.pop(1)
    return gmax

def maxvectorplot(dataset):
    positions = []
    displacements = []
    nonempty = []
    for t in range(t_steps):
        if globalmax(dataset[t]) != "empty":
            nonempty += globalmax(dataset[t])
    for i in range(1,len(nonempty)-1):
        if (nonempty[i][2]>nonempty[i-1][2]) and (nonempty[i][2]>nonempty[i+1][2]):
            positions += [[nonempty[i][0],nonempty[i][1]]]
    for i in range(1,len(positions)):
        displacements += [[positions[i-1][0],positions[i-1][1],positions[i][0]-positions[i-1][0],positions[i][1]-positions[i-1][1]]] #[tail_x, tail_y, head_x, head_y]
    testv = array([[0,0,1,2], [1,2,2,3], [2,3,4,5], [4,5,6,8]])
    X,Y,U,V = zip(*displacements)
    plt.figure()
    ax = plt.gca()
    ax.quiver(X,Y,U,V)
    plt.xlim((0,data_shape[0]))
    plt.ylim((0,data_shape[1]))
    scatter(unzip_membrane(cell_membrane)[0],unzip_membrane(cell_membrane)[1])
    show()
    return 0

maxvectorplot(data_natp_set)
