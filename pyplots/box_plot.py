from __future__ import division
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
import pylab
import file_loader as load

#works, but also a work in progress.

def oscillation_period(density_list):
    peaklist = []
    index_position_pairs = list(enumerate(density_list))
    for i in index_position_pairs:
        if (abs(i[1]-max(density_list)) < .01*max(density_list)):
            peaklist += [i[0]]
    return (peaklist[0], peaklist[1])

NflE = load.data(protein="NflE")
NflD = load.data(protein="NflD")
nATP = load.data(protein="nATP")
nE = load.data(protein="nE")
nADP = load.data(protein="nADP")
Nd = load.data(protein="Nd")
Nde = load.data(protein="Nde")

proteins = [NflE, NflD]
sub_proteins_D = [nATP, nADP, Nd, Nde]
sub_proteins_E = [nE, Nde]

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

    #gather the individual proteins to stick on the same graph as NflD or NflE.
    if p==NflE:
        for individual_proteins in sub_proteins_E:
            #the following 3 lines convert densities to counts (i.e. nATP to NATP) so the graphs make sense
            dV = 1
            if individual_proteins == nE:
                dV = load.dx**3
            individual_proteins.proteins_left = []
            individual_proteins.proteins_mid = []
            individual_proteins.proteins_right = []

            for file in individual_proteins.dataset:
                left_sum = 0
                mid_sum = 0
                right_sum = 0
                for i in range(individual_proteins.datashape[0]-1):
                    for j in range(individual_proteins.datashape[1]-1):
                        if j <= divider_left:
                            left_sum += file[i][j]
                        elif divider_left < j <= divider_right:
                            mid_sum += file[i][j]
                        else:
                            right_sum += file[i][j]
                individual_proteins.proteins_left += [left_sum*dV]
                individual_proteins.proteins_mid += [mid_sum*dV]
                individual_proteins.proteins_right += [right_sum*dV]
        NflE.lines = [nE.proteins_left, Nde.proteins_left, \
                          nE.proteins_mid, Nde.proteins_mid, \
                          nE.proteins_right, Nde.proteins_right]

    if p==NflD:
        for individual_proteins in sub_proteins_D:
            dV = 1
            if (individual_proteins == nATP) or (individual_proteins == nADP):
                dV = load.dx**3
            individual_proteins.proteins_left = []
            individual_proteins.proteins_mid = []
            individual_proteins.proteins_right = []

            for file in individual_proteins.dataset:
                left_sum = 0
                mid_sum = 0
                right_sum = 0
                for i in range(individual_proteins.datashape[0]-1):
                    for j in range(individual_proteins.datashape[1]-1):
                        if j <= divider_left:
                            left_sum += file[i][j]
                        elif divider_left < j <= divider_right:
                            mid_sum += file[i][j]
                        else:
                            right_sum += file[i][j]
                individual_proteins.proteins_left += [left_sum*dV]
                individual_proteins.proteins_mid += [mid_sum*dV]
                individual_proteins.proteins_right += [right_sum*dV]
        NflD.lines = [nATP.proteins_left, nADP.proteins_left, Nd.proteins_left, Nde.proteins_left, \
                          nATP.proteins_mid, nADP.proteins_mid, Nd.proteins_mid, Nde.proteins_mid, \
                          nATP.proteins_right, nADP.proteins_right, Nd.proteins_right, Nde.proteins_right]

    #calculate the period
    (start, end) = (0, -1)
    if ("-all" in sys.argv) or ("-short" in sys.argv):
        (start, end) = oscillation_period(p.proteins_mid)

    #plot them and save the figures.
    plt.figure()
    plt.xlabel("Time")
    plt.ylabel("%s count"%p.protein)
    plt.title("%s count vs time"%p.protein)

    plot_lines = [p.lines[0]]
    for i in range(1,len(p.lines)):
        plot_lines += [[j+k for (j,k) in zip(plot_lines[i-1], p.lines[i])]]

    #the list comprehensions used below are just one-liner ways to add the elements in two lists
    lowest_line = p.proteins_left
    middle_line = [i+j for (i,j) in zip(lowest_line,p.proteins_mid)]
    top_line = [i+j for (i,j) in zip(middle_line,p.proteins_right)]

    time_axis = list(np.arange(0,len(p.proteins_left)*5,5))

    plt.fill_between(time_axis[start:end],0,lowest_line[start:end],alpha=0.5)
    plt.fill_between(time_axis[start:end],lowest_line[start:end],middle_line[start:end],alpha=0.5,facecolor="green")
    plt.fill_between(time_axis[start:end],middle_line[start:end],top_line[start:end],alpha=0.5,facecolor="red")

    plt.plot(time_axis[start:end],lowest_line[start:end],color="blue",linewidth=0.5)
    plt.plot(time_axis[start:end],middle_line[start:end],color="green",linewidth=0.5)
    plt.plot(time_axis[start:end],top_line[start:end],color="red",linewidth=0.5)

    all_string = ""
    if "-all" in sys.argv:
        for line in plot_lines:
            plt.plot(time_axis[start:end],line[start:end],color="black",linewidth=0.5)
        all_string = "-all"

    short_string = ""
    if "-short" in sys.argv:
        short_string = "-short"

    plt.xlim(time_axis[start],max(time_axis[start:end]))
    plt.ylim(0,1.4*max(top_line[start:end]))
    #plt.legend(["left density", "middle density", "right density"],loc="best")
    plt.savefig(load.print_string("box-plot"+short_string+all_string,p))
