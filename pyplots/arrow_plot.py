from __future__ import division
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import file_loader as load
import sys

#membrane.dat or arrow.dat printed transposed, gotta align them:
cell_membrane = np.transpose(np.loadtxt("./data/shape-%s/membrane_files/%s%s%smembrane-%s-%s-%s-%s-%s-%s.dat"%(load.f_shape,load.debug_str,load.hires_str,load.slice_str,load.f_shape,load.f_param1,load.f_param2,load.f_param3,load.f_param4,load.f_param5)))

for protein in load.proteinList:
    tails = np.loadtxt('./data/shape-%s/%s%s%sarrow-plot-%s-%s-%s-%s-%s-%s-%s.dat'%(load.f_shape,load.debug_str,load.hires_str,load.slice_str,protein,load.f_shape,load.f_param1,load.f_param2,load.f_param3,load.f_param4,load.f_param5))
    tails = [list(element) for element in tails]
    heads = [[tails[1][0]-tails[0][0],tails[1][1]-tails[0][1]]]
    heads += [[tails[i+1][0]-tails[i][0], tails[i+1][1]-tails[i][1]] for i in range(1,len(tails)-1)]
    displacements = [tail_element+head_element for (tail_element,head_element) in zip(tails,heads)] #in format: [point_x, point_y, vector_component_x, vector_component_y]
    X,Y,U,V = zip(*displacements)
    plt.figure()
    plt.axes().set_aspect('equal', 'datalim')
    plt.ax = plt.gca()
    plt.contourf(cell_membrane, levels=[0.5, 1])
    plt.ax.quiver(X,Y,U,V,scale_units='xy',angles='xy',scale=1)
    plt.xlim((0,cell_membrane.shape[1]))
    plt.ylim((0,cell_membrane.shape[0]))
    plt.xlabel("Z grid position")
    plt.ylabel("Y grid position")
    plt.title("Local temporal maxima, global spatial maxima view of %s"%(protein))
    plt.savefig(load.print_string("arrow-plot",protein))
