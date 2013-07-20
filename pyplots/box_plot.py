from __future__ import division
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab
import os
import sys
import file_loader as load

#bug: very memory intensive, especially with fill_between

#only need these two data types for now
NflE = load.data(protein="NflE")
NflD = load.data(protein="NflD")

proteins = [NflE, NflD]

for p in proteins:
    #pick a relatively arbitrary box divider
    divider_left = round(p.datashape[1]/3)
    divider_right = round(p.datashape[1]*2/3)
    p.proteins_left = []
    p.proteins_mid = []
    p.proteins_right = []

    #collect the total number of proteins in each box, at each time step
    for file in p.dataset:
        left_sum = 0
        mid_sum = 0
        right_sum = 0
        for i in range(p.datashape[0]-1):
            for j in range(p.datashape[1]-1):
                if j <= divider_left:
                    left_sum += file[i][j]
                elif divider_left < j <= divider_right:
                    mid_sum += file[i][j]
                else:
                    right_sum += file[i][j]
        p.proteins_left += [left_sum]
        p.proteins_mid += [mid_sum]
        p.proteins_right += [right_sum]

        #plot them and save the figures.
        plt.figure()
        plt.xlabel("Time")
        plt.ylabel("%s count"%p.protein)

        #the list comprehensions used below are just one-liner ways to add the elements in two lists
        lowest_line = p.proteins_left
        middle_line = [i+j for (i,j) in zip(lowest_line,p.proteins_mid)]
        top_line = [i+j for (i,j) in zip(middle_line,p.proteins_right)]

        time_axis = list(np.arange(0,len(p.proteins_left)*.5,.5))

        # plt.fill_between(time_axis,0,lowest_line,alpha=0.5)
        # plt.fill_between(time_axis,lowest_line,middle_line,alpha=0.5)
        # plt.fill_between(time_axis,middle_line,top_line,alpha=0.5)

        plt.plot(time_axis,lowest_line)
        plt.plot(time_axis,middle_line)
        plt.plot(time_axis,top_line)
        plt.savefig('./data/shape-'+load.f_shape+'/plots/'+load.hires_str+load.m_str+'box-plot-'+str(p.protein)+'-'+load.f_shape+'-'+load.f_param1+'-'+load.f_param2+'-'+load.f_param3+'-'+load.f_param4+'-'+load.f_param5+'.pdf')
