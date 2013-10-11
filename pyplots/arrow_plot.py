from __future__ import division
import matplotlib
import sys
if "show" not in sys.argv:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import file_loader as load
from matplotlib.font_manager import FontProperties

#membrane.dat or arrow.dat printed transposed, gotta align them:
cell_membrane_initial = np.loadtxt("./data/shape-%s/membrane_files/%s%s%ssections-%s-%s-%s-%s-%s-%s.dat"%
                                        (load.f_shape,load.debug_str,load.hires_str,load.slice_str,
                                         load.f_shape,load.f_param1,load.f_param2,load.f_param3,load.f_param4,load.f_param5))

#http://stackoverflow.com/questions/8486294/how-to-add-an-extra-column-to-an-numpy-array
#this is algorithm want
Ny = len(cell_membrane_initial[:,0])
Nz = len(cell_membrane_initial[0,:])

temp_array = np.zeros((Ny+1,Nz+1))
temp_array[:-1,:-1] = cell_membrane_initial
cell_membrane = temp_array


row_of_zeros = np.zeros(len(cell_membrane[:,0]))
np.vstack((cell_membrane[:,0],row_of_zeros))
np.vstack((row_of_zeros,cell_membrane[:,0]))

row_of_zeros = np.zeros(len(cell_membrane[0,:]))
np.vstack((cell_membrane[0,:],row_of_zeros))
np.vstack((row_of_zeros,cell_membrane[0,:]))


cell_membrane
if load.f_param4 != "94.00":
    np.transpose(cell_membrane)



for protein in load.proteinList:
    tails = np.loadtxt('./data/shape-%s/%s%s%sarrow-plot-%s-%s-%s-%s-%s-%s-%s.dat'%(load.f_shape,load.debug_str,load.hires_str,load.slice_str,protein,load.f_shape,load.f_param1,load.f_param2,load.f_param3,load.f_param4,load.f_param5))
    tails = [list(element) for element in tails]
    cut = len(tails)-10
    tails = [tail for tail in tails[cut:][:]]
    # z_max = max(tails[:][2])
    # z_min = min(tails[:][2])
    # print "z_max = ",z_max
    # print "z_min = ",z_min
    # i = 0
    # print "tails len = ",len(tails)
    # start_iter = 0
    # end_iter = len(tails)-1
    # for i in xrange(len(tails)):
    #     print i, tails[i]
    #     if tails[i][2] > .5*z_max+.5*z_min:
    #         start_iter = i
    #         break
    # start_time = tails[start_iter][0]
    # for i in xrange(start_iter,len(tails)):
    #     if tails[i][2] < .5*z_max+.5*z_min:
    #         middle_iter = i
    #         break
    # for i in xrange(middle_iter,len(tails)):
    #     if tails[i][2] > .5*z_max+.5*z_min:
    #         end_iter = i
    #         break
    # end_time = tails[end_iter][0]

    # tails = [[tails[i][1],tails[i][2]] for i in range(0,len(tails)-1) if (tails[i][0]>start_time and tails[i][0]<=end_time)]
    # if len(tails)<2:
    #     print "you need a bigger time period!"
    #     exit(1)
    if len(tails)==0:
        print "The tails in this file has no tails!!"
        exit(0)
    heads = [[load.dx*(tails[1][1]-tails[0][1]),load.dx*(tails[1][0]-tails[0][0])]]
    heads += [[load.dx*(tails[i+1][1]-tails[i][1]), load.dx*(tails[i+1][0]-tails[i][0])] for i in range(1,len(tails)-1)]
    tails = [[load.dx*tails[i][1],load.dx*tails[i][0]] for i in range(len(tails))]
    plt.figure()
    plt.axes().set_aspect('equal', 'datalim')
    plt.ax = plt.gca()
    font=FontProperties()
    font.set_family('serif')
    for i in range(len(tails)-1):
        radial_length = np.sqrt((tails[i][0]-load.dx*Ny/2.0)**2+(tails[i][1]-load.dx*Nz/2.0)**2)
        dir_z = (tails[i][0]-load.dx*Nz/2.0)/radial_length
        dir_y = (tails[i][1]-load.dx*Ny/2.0)/radial_length
        number_zpos = tails[i][0]+.18*dir_z
        number_ypos = tails[i][1]+.18*dir_y
        plt.ax.annotate('',xy=(tails[i+1][0],tails[i+1][1]),xytext=(tails[i][0],tails[i][1]),
                        fontsize=7,
                        arrowprops=dict(facecolor='red',shrink=0.01, width=.3, headwidth=5.))
        plt.ax.text(number_zpos,number_ypos,"%d"%(i+1),fontsize=30,fontproperties=font,color='red')
    cell_membrane[cell_membrane>0] = 1
    cell_x = np.linspace(0,load.dx*len(cell_membrane[0,:]),len(cell_membrane[0,:]))
    cell_y = np.linspace(0,load.dx*len(cell_membrane[:,0]),len(cell_membrane[:,0]))
    plt.contour(cell_x,cell_y,cell_membrane, linewidths=2,levels=[.99])
    #plt.ax.quiver(X,Y,U,V,scale_units='xy',angles='xy',scale=1)
    plt.ax.get_yaxis().set_visible(False)
    plt.ax.get_xaxis().set_visible(False)
    plt.xlim((0,load.dx*cell_membrane.shape[1]))
    plt.ylim((0,load.dx*cell_membrane.shape[0]))
    plt.xlabel("Z grid position")
    plt.ylabel("Y grid position")
    plt.title("Local temporal maxima, global spatial maxima view of %s"%(protein))
    plt.savefig(load.print_string("arrow-plot",protein))

plt.show()
