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

#cell_membrane = loadtxt('./data/shape-'+f_shape+'/membrane_files/membrane-'+f_param1+'-'+f_param2+'-'+f_param3+'-'+f_param4+'-'+f_param5+'.dat')

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

# X, Y = meshgrid(arange(0,6,dx),arange(0,6,dx))
# Z = [0. for i in range(20)]
# one = ones_like(X)
# for t in range(20):
#     Z[t] = exp(-(X-one*t*(6/20))**2 - (Y-one*t*(6/20))**2)
# #z is a moving guassian with constant (1.0) amplitude.

# print "Debugging extrema.py, with our friend, the moving guassian distribution."
# # print "global maxima:"
# # for i in range(20):
# #     print globalmax(Z[i])

# # print "local maxima:"
# # for i in range(20):
# #     print localmax(Z[i])

# print "now lets add some time dependence to the amplitude so we have local maxima in time as well:"
# for t in range(20):
#     Z[t] = exp(-(X-one*t*(6/20))**2 - (Y-one*t*(6/20))**2)*sin(one*t*(6/20)) #periodic amplitude

# print "global maxima:"
# for i in range(20):
#     print globalmax(Z[i])

# print "local maxima in time:"
# #prune the "empty"'s:
# A = []
# for t in range(20):
#     if globalmax(Z[t]) != "empty":
#         A += globalmax(Z[t])

# print A

# for t in range(1,len(A)-1):
#     if (A[t][2]>A[t-1][2]) and (A[t][2]>A[t+1][2]):
#         print A[t]

# #hooray! working algorithm.

def maxtracker(dataset):
    locs = zeros_like(dataset[0])
    for t in range(1,t_steps-1):
        if globalmax(dataset[t]) != "empty":
            for i in range(len(globalmax(dataset[t]))):
                if (globalmax(dataset[t])[i] > globalmax(dataset[t+1])[i]) and (globalmax(dataset[t])[i] > globalmax(dataset[t-1])[i]):
                    x, y, m = globalmax(dataset[t])[i] #[i] accounting for possibility of multiple maxima
                    locs[x][y] = m
    return locs

def maxvectorplot(dataset):
    positions = []
    displacements = []
    nonempty = []
    for t in range(t_steps):
        if globalmax(dataset[t]) != "empty":
            nonempty += globalmax(dataset[t])
    print nonempty
    for i in range(1,len(nonempty)-1):
        if (nonempty[i][2]>nonempty[i-1][2]) and (nonempty[i][2]>nonempty[i+1][2]):
            print "got here"
            positions += [[nonempty[i][0],nonempty[i][1]]]
    print positions
    for i in range(1,len(positions)):
        displacements += [[positions[i-1][0],positions[i-1][1],positions[i][0]-positions[i-1][0],positions[i][1]-positions[i-1][1]]] #[tail_x, tail_y, head_x, head_y]
    print displacements
    X,Y,U,V = zip(*displacements)
    plt.figure()
    ax = plt.gca()
    ax.quiver(X,Y,U,V,angles='xy',scale=1,linewidths=.05)
#    scatter(unzip_membrane(cell_membrane)[0],unzip_membrane(cell_membrane)[1], marker='.')
    plt.xlim((0,data_shape[0]))
    plt.ylim((0,data_shape[1]))
    print displacements
    show()
    return 0

maxvectorplot(data_natp_set)
