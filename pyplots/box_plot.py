from __future__ import division
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
import pylab
import file_loader as load

#works, but also a work in progress.
#flags: -all: includes all the proteins of type D/E on the D/E graph, and truncates the graph to 1 period for readability
#       -short: truncates the graph to 1 period for readability.

#returns the indicies of the first two peaks.
def oscillation_period(density_list):
    peaklist = []
    index_position_pairs = list(enumerate(density_list))
    for i in index_position_pairs:
        if (abs(i[1]-max(density_list)) < .01*max(density_list)):
            peaklist += [i[0]]
    return (peaklist[0], peaklist[1])

#load protein data.
NflE = load.data(protein="NflE")
NflD = load.data(protein="NflD")
nATP = load.data(protein="nATP")
nE = load.data(protein="nE")
nADP = load.data(protein="nADP")
Nd = load.data(protein="Nd")
Nde = load.data(protein="Nde")

#proteins is the main protein list. this program makes one graph for each element in proteins
proteins = [NflE, NflD]

#sub proteins are the individual, optionally-plotted proteins.
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

    #gather the individual proteins to stick on the same graph as NflE.
    if p==NflE:
        for individual_proteins in sub_proteins_E:
            #density to numbe rconversion for nE.
            dV = 1
            if individual_proteins == nE:
                dV = load.dx**3
            #use the same method as before to store proteins in each bo.x
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

    #gather the individual proteins to stick on the same graph as NflD.
    if p==NflD:
        for individual_proteins in sub_proteins_D:
            #density to number conversion for nATP and nADP.
            dV = 1
            if (individual_proteins == nATP) or (individual_proteins == nADP):
                dV = load.dx**3
            #use the same method as before to store proteins in each box.
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
        #store them in one big list for plotting purposes later. the order is important, since we will be
        #adding each element in this list to the one before it, starting with the first one (for the stacked plot)
        NflD.lines = [nATP.proteins_left, nADP.proteins_left, Nd.proteins_left, Nde.proteins_left, \
                          nATP.proteins_mid, nADP.proteins_mid, Nd.proteins_mid, Nde.proteins_mid, \
                          nATP.proteins_right, nADP.proteins_right, Nd.proteins_right, Nde.proteins_right]

    #calculate the indicies of the first period or so
    #if -all and -short are both not present in the system arguments,
    #plot all of the data.
    (start, end) = (0, -1)
    if ("-all" in sys.argv) or ("-short" in sys.argv):
        (start, end) = oscillation_period(p.proteins_mid)

    #plot them and save the figures.
    plt.figure()
    plt.xlabel("Time")
    plt.ylabel("%s count"%p.protein)
    plt.title("%s count vs time"%p.protein)


    #plot_lines is the list where we will store each line to be plotted. the first element is set to
    #NflD/NflE.lines[0], which was defined in the above loop. the list comprehension just adds the ith
    #list to the (i-1)th list.
    plot_lines = [p.lines[0]]
    for i in range(1,len(p.lines)):
        plot_lines += [[j+k for (j,k) in zip(plot_lines[i-1], p.lines[i])]]

    #use the same method to get NflD/E stacked on top of each other. lowest = left, mid = mid, top = right
    lowest_line = p.proteins_left
    middle_line = [i+j for (i,j) in zip(lowest_line,p.proteins_mid)]
    top_line = [i+j for (i,j) in zip(middle_line,p.proteins_right)]

    #create a time axis
    time_axis = list(np.arange(0,len(p.proteins_left)*5,5))

    #matplotlib method for filling in color between two lines. left = blue, mid = green, right = red.
    if "-all" not in sys.argv:
        plt.fill_between(time_axis[start:end],0,lowest_line[start:end],alpha=0.5)
        plt.fill_between(time_axis[start:end],lowest_line[start:end],middle_line[start:end],alpha=0.5,facecolor="green")
        plt.fill_between(time_axis[start:end],middle_line[start:end],top_line[start:end],alpha=0.5,facecolor="red")

    if ("-all" in sys.argv) and (p==NflD):
        plt.fill_between(time_axis[start:end],[0]*(end-start),plot_lines[0][start:end],alpha=0.25,facecolor="blue")
        plt.fill_between(time_axis[start:end],plot_lines[0][start:end],plot_lines[1][start:end],alpha=0.5,facecolor="blue")
        plt.fill_between(time_axis[start:end],plot_lines[1][start:end],plot_lines[2][start:end],alpha=0.75,facecolor="blue")
        plt.fill_between(time_axis[start:end],plot_lines[2][start:end],plot_lines[3][start:end],alpha=1.0,facecolor="blue")
        plt.fill_between(time_axis[start:end],plot_lines[3][start:end],plot_lines[4][start:end],alpha=0.25,facecolor="red")
        plt.fill_between(time_axis[start:end],plot_lines[4][start:end],plot_lines[5][start:end],alpha=0.5,facecolor="red")
        plt.fill_between(time_axis[start:end],plot_lines[5][start:end],plot_lines[6][start:end],alpha=0.75,facecolor="red")
        plt.fill_between(time_axis[start:end],plot_lines[6][start:end],plot_lines[7][start:end],alpha=1.0,facecolor="red")
        plt.fill_between(time_axis[start:end],plot_lines[7][start:end],plot_lines[8][start:end],alpha=0.25,facecolor="green")
        plt.fill_between(time_axis[start:end],plot_lines[8][start:end],plot_lines[9][start:end],alpha=0.5,facecolor="green")
        plt.fill_between(time_axis[start:end],plot_lines[9][start:end],plot_lines[10][start:end],alpha=0.75,facecolor="green")
        plt.fill_between(time_axis[start:end],plot_lines[10][start:end],plot_lines[11][start:end],alpha=1.0,facecolor="green")
        plt.plot(time_axis[start:end],plot_lines[0][start:end],alpha=0.25,color="blue",linewidth=0.5)
        plt.plot(time_axis[start:end],plot_lines[0][start:end],alpha=0.5,color="blue",linewidth=0.5)
        plt.plot(time_axis[start:end],plot_lines[0][start:end],alpha=0.75,color="blue",linewidth=0.5)
        plt.plot(time_axis[start:end],plot_lines[0][start:end],alpha=1.0,color="blue",linewidth=0.5)
        plt.plot(time_axis[start:end],plot_lines[0][start:end],alpha=0.25,color="red",linewidth=0.5)
        plt.plot(time_axis[start:end],plot_lines[0][start:end],alpha=0.5,color="red",linewidth=0.5)
        plt.plot(time_axis[start:end],plot_lines[0][start:end],alpha=0.75,color="red",linewidth=0.5)
        plt.plot(time_axis[start:end],plot_lines[0][start:end],alpha=1.0,color="red",linewidth=0.5)
        plt.plot(time_axis[start:end],plot_lines[0][start:end],alpha=0.25,color="green",linewidth=0.5)
        plt.plot(time_axis[start:end],plot_lines[0][start:end],alpha=0.5,color="green",linewidth=0.5)
        plt.plot(time_axis[start:end],plot_lines[0][start:end],alpha=0.75,color="green",linewidth=0.5)
        plt.plot(time_axis[start:end],plot_lines[0][start:end],alpha=1.0,color="green",linewidth=0.5)
        # plt.legend(["nATP proteins left", "nADP proteins left", "Nd proteins left", "Nde proteins left", \
        #                   "nATP proteins mid", "nADP proteins mid", "Nd proteins mid", "Nde proteins mid", \
        #                   "nATP proteins right", "nADP proteins right", "Nd proteins right", "Nde proteins right"],loc="best")

    if ("-all" in sys.argv) and (p==NflE):
        plt.fill_between(time_axis[start:end],[0]*(end-start),plot_lines[0][start:end],alpha=0.3,facecolor="blue")
        plt.fill_between(time_axis[start:end],plot_lines[0][start:end],plot_lines[1][start:end],alpha=1.0,facecolor="blue")
        plt.fill_between(time_axis[start:end],plot_lines[1][start:end],plot_lines[2][start:end],alpha=0.3,facecolor="red")
        plt.fill_between(time_axis[start:end],plot_lines[2][start:end],plot_lines[3][start:end],alpha=1.0,facecolor="red")
        plt.fill_between(time_axis[start:end],plot_lines[3][start:end],plot_lines[4][start:end],alpha=0.3,facecolor="green")
        plt.fill_between(time_axis[start:end],plot_lines[4][start:end],plot_lines[5][start:end],alpha=1.0,facecolor="green")
        plt.plot(time_axis[start:end],plot_lines[0][start:end],alpha=0.3,color="blue",linewidth=0.5)
        plt.plot(time_axis[start:end],plot_lines[1][start:end],alpha=1.0,color="blue",linewidth=0.5)
        plt.plot(time_axis[start:end],plot_lines[2][start:end],alpha=0.3,color="red",linewidth=0.5)
        plt.plot(time_axis[start:end],plot_lines[3][start:end],alpha=1.0,color="red",linewidth=0.5)
        plt.plot(time_axis[start:end],plot_lines[4][start:end],alpha=0.3,color="green",linewidth=0.5)
        plt.plot(time_axis[start:end],plot_lines[5][start:end],alpha=1.0,color="green",linewidth=0.5)
        plt.legend(["nE proteins left", "Nde proteins left", "nE proteins mid", "Nde proteins mid","nE proteins right", "Nde proteins right"])

    plt.plot(time_axis[start:end],lowest_line[start:end],color="blue",linewidth=0.5)
    plt.plot(time_axis[start:end],middle_line[start:end],color="green",linewidth=0.5)
    plt.plot(time_axis[start:end],top_line[start:end],color="red",linewidth=0.5)

    #if -all is present, include all the lines, not just the usual 3, and flag it in the file name
    all_string = ""
    if "-all" in sys.argv:
        for line in plot_lines:
            plt.plot(time_axis[start:end],line[start:end],color="black",linewidth=0.5)
        all_string = "-all"

    #if short is present, flag it in the file name
    short_string = ""
    if "-short" in sys.argv:
        short_string = "-short"

    plt.xlim(time_axis[start],max(time_axis[start:end]))
    plt.ylim(0,1.4*max(top_line[start:end]))
    if "-all" not in sys.argv:
        plt.legend(["left density", "middle density", "right density"],loc="best")
    plt.savefig(load.print_string("box-plot"+short_string+all_string,p))

    plt.figure()
    plt.xlabel("Time")
    plt.ylabel("%s count"%p.protein)
    plt.title("%s count vs time at end caps of cell"%p.protein)

    if ("-all" in sys.argv) and (p==NflD):
        plt.fill_between(time_axis[start:end],[0]*(end-start),plot_lines[0][start:end],alpha=1.0,facecolor="blue")
        plt.fill_between(time_axis[start:end],plot_lines[0][start:end],plot_lines[1][start:end],alpha=1.0,facecolor="red")
        plt.fill_between(time_axis[start:end],plot_lines[1][start:end],plot_lines[2][start:end],alpha=1.0,facecolor="green")
        plt.fill_between(time_axis[start:end],plot_lines[2][start:end],plot_lines[3][start:end],alpha=1.0,facecolor="gold")
        plt.plot(time_axis[start:end],plot_lines[0][start:end],alpha=1.0,color="blue",linewidth=0.5)
        plt.plot(time_axis[start:end],plot_lines[0][start:end],alpha=1.0,color="red",linewidth=0.5)
        plt.plot(time_axis[start:end],plot_lines[0][start:end],alpha=1.0,color="green",linewidth=0.5)
        plt.plot(time_axis[start:end],plot_lines[0][start:end],alpha=1.0,color="gold",linewidth=0.5)
        plt.legend(["nATP proteins", "nADP proteins", "Nd proteins", "Nde proteins"],loc="best")

    if ("-all" in sys.argv) and (p==NflE):
        plt.fill_between(time_axis[start:end],[0]*(end-start),plot_lines[0][start:end],alpha=1.0,facecolor="blue")
        plt.fill_between(time_axis[start:end],plot_lines[0][start:end],plot_lines[1][start:end],alpha=1.0,facecolor="blue")
        plt.plot(time_axis[start:end],plot_lines[0][start:end],alpha=1.0,color="blue",linewidth=0.5)
        plt.plot(time_axis[start:end],plot_lines[1][start:end],alpha=1.0,color="blue",linewidth=0.5)
        plt.legend(["nE proteins", "Nde proteins"],loc="best")

    all_string = ""
    if "-all" in sys.argv:
        for i in range(4):
            plt.plot(time_axis[start:end],plot_lines[i][start:end],color="black",linewidth=0.5)
        all_string = "-all"

    #if short is present, flag it in the file name
    short_string = ""
    if "-short" in sys.argv:
        short_string = "-short"

    plt.xlim(time_axis[start],max(time_axis[start:end]))
    plt.ylim(0,1.4*max(top_line[start:end]))
    plt.savefig(load.print_string("endcap-plot"+short_string+all_string,p))
