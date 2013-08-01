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

hires_str = ""
if "-hires" in sys.argv:
    dx=.05
    hires_str="hires-"


f_shape = sys.argv[1]
f_param1 = sys.argv[2]
f_param2 = sys.argv[3]
f_param3 = sys.argv[4]
f_param4 = sys.argv[5]
f_param5 = sys.argv[6]

div_file =  "./data/shape-"+f_shape+"/"+hires_str+"m-divisions-"+f_shape+"-"+f_param1+"-"+f_param2+"-"+f_param3+"-"+f_param4+"-"+f_param5+".dat"

compare_file = "./data/shape-"+f_shape+"/"+hires_str+"m-compare-"+f_shape+"-"+f_param1+"-"+f_param2+"-"+f_param3+"-"+f_param4+"-"+f_param5+".dat"

divider_file = "./data/shape-"+f_shape+"/"+hires_str+"m-cell-dividers-"+f_shape+"-"+f_param1+"-"+f_param2+"-"+f_param3+"-"+f_param4+"-"+f_param5+".dat"

print div_file
div_data = np.loadtxt(div_file)
compare_data = np.loadtxt(compare_file)
#divider_data = np.loadtxt(divider_file)


#create a time axis
time_axis = np.arange(0,len(div_data[:,0]))

indices_from_cpp_D = [1,2,3,7,8,9,10,14,15,16,17]

b=0
plot_lines_D = [div_data[:,0]]
for a in indices_from_cpp_D:
    plot_lines_D += [[j+k for (j,k) in zip(plot_lines_D[b], div_data[:,a])]]
    print a
    b+=1


plt.figure()
plt.xlabel("Time")
plt.ylabel("Ave Number of MinD Proteins Per Area")
plt.title("Ave Number of MinD Proteins Per Area Vs Time")

plt.fill_between(time_axis,[0]*len(time_axis),plot_lines_D[0],alpha=0.25,facecolor="blue")
plt.fill_between(time_axis,plot_lines_D[0],plot_lines_D[1],alpha=0.5,facecolor="blue")
plt.fill_between(time_axis,plot_lines_D[1],plot_lines_D[2],alpha=0.75,facecolor="blue")
plt.fill_between(time_axis,plot_lines_D[2],plot_lines_D[3],alpha=1.0,facecolor="blue")

plt.fill_between(time_axis,plot_lines_D[3],plot_lines_D[4],alpha=0.25,facecolor="red")
plt.fill_between(time_axis,plot_lines_D[4],plot_lines_D[5],alpha=0.5,facecolor="red")
plt.fill_between(time_axis,plot_lines_D[5],plot_lines_D[6],alpha=0.75,facecolor="red")
plt.fill_between(time_axis,plot_lines_D[6],plot_lines_D[7],alpha=1.0,facecolor="red")

plt.fill_between(time_axis,plot_lines_D[7],plot_lines_D[8],alpha=0.25,facecolor="green")
plt.fill_between(time_axis,plot_lines_D[8],plot_lines_D[9],alpha=0.5,facecolor="green")
plt.fill_between(time_axis,plot_lines_D[9],plot_lines_D[10],alpha=0.75,facecolor="green")
plt.fill_between(time_axis,plot_lines_D[10],plot_lines_D[11],alpha=1.0,facecolor="green")

plt.plot(time_axis,plot_lines_D[0],alpha=0.25,color="blue",linewidth=0.5)
plt.plot(time_axis,plot_lines_D[1],alpha=0.5,color="blue",linewidth=0.5)
plt.plot(time_axis,plot_lines_D[2],alpha=0.75,color="blue",linewidth=0.5)
plt.plot(time_axis,plot_lines_D[3],alpha=1.0,color="blue",linewidth=0.5)

plt.plot(time_axis,plot_lines_D[4],alpha=0.25,color="red",linewidth=0.5)
plt.plot(time_axis,plot_lines_D[5],alpha=0.5,color="red",linewidth=0.5)
plt.plot(time_axis,plot_lines_D[6],alpha=0.75,color="red",linewidth=0.5)
plt.plot(time_axis,plot_lines_D[7],alpha=1.0,color="red",linewidth=0.5)

plt.plot(time_axis,plot_lines_D[8],alpha=0.25,color="green",linewidth=0.5)
plt.plot(time_axis,plot_lines_D[9],alpha=0.5,color="green",linewidth=0.5)
plt.plot(time_axis,plot_lines_D[10],alpha=0.75,color="green",linewidth=0.5)
plt.plot(time_axis,plot_lines_D[11],alpha=1.0,color="green",linewidth=0.5)


plt.savefig("./data/shape-"+f_shape+"/plots/m"+hires_str+"-divisions-MinD-"+f_shape+"-"+f_param1+"-"+f_param2+"-"+f_param3+"-"+f_param4+"-"+f_param5+".pdf")

indices_from_cpp_E = [4,10,11,17,18]

b=0
plot_lines_E = [div_data[:,3]]
for a in indices_from_cpp_E:
    plot_lines_E += [[j+k for (j,k) in zip(plot_lines_E[b], div_data[:,a])]]
    print a
    b+=1

plt.figure()
plt.xlabel("Time")
plt.ylabel("Ave Number of MinE Proteins Per Area")
plt.title("Ave Number of MinE Proteins Per Area Vs Time")

plt.fill_between(time_axis,[0]*len(time_axis),plot_lines_E[0],alpha=0.25,facecolor="blue")
plt.fill_between(time_axis,plot_lines_E[0],plot_lines_E[1],alpha=0.5,facecolor="blue")

plt.fill_between(time_axis,plot_lines_E[1],plot_lines_E[2],alpha=0.25,facecolor="red")
plt.fill_between(time_axis,plot_lines_E[2],plot_lines_E[3],alpha=0.5,facecolor="red")

plt.fill_between(time_axis,plot_lines_E[3],plot_lines_E[4],alpha=0.25,facecolor="green")
plt.fill_between(time_axis,plot_lines_E[4],plot_lines_E[5],alpha=0.5,facecolor="green")

plt.plot(time_axis,plot_lines_E[0],alpha=0.25,color="blue",linewidth=0.5)
plt.plot(time_axis,plot_lines_E[1],alpha=0.5,color="blue",linewidth=0.5)

plt.plot(time_axis,plot_lines_E[2],alpha=0.25,color="red",linewidth=0.5)
plt.plot(time_axis,plot_lines_E[3],alpha=0.5,color="red",linewidth=0.5)

plt.plot(time_axis,plot_lines_E[4],alpha=0.25,color="green",linewidth=0.5)
plt.plot(time_axis,plot_lines_E[5],alpha=0.5,color="green",linewidth=0.5)



plt.savefig("./data/shape-"+f_shape+"/plots/m"+hires_str+"-divisions-MinE-"+f_shape+"-"+f_param1+"-"+f_param2+"-"+f_param3+"-"+f_param4+"-"+f_param5+".pdf")
#     if ("-all" in sys.argv) and (p==NflD):
#         plt.fill_between(time_axis[start:end],[0]*(end-start),plot_lines[0][start:end],alpha=1.0,facecolor="blue")
#         plt.fill_between(time_axis[start:end],plot_lines[0][start:end],plot_lines[1][start:end],alpha=1.0,facecolor="red")
#         plt.fill_between(time_axis[start:end],plot_lines[1][start:end],plot_lines[2][start:end],alpha=1.0,facecolor="green")
#         plt.fill_between(time_axis[start:end],plot_lines[2][start:end],plot_lines[3][start:end],alpha=1.0,facecolor="gold")
#         plt.plot(time_axis[start:end],plot_lines[0][start:end],alpha=1.0,color="blue",linewidth=0.5)
#         plt.plot(time_axis[start:end],plot_lines[0][start:end],alpha=1.0,color="red",linewidth=0.5)
#         plt.plot(time_axis[start:end],plot_lines[0][start:end],alpha=1.0,color="green",linewidth=0.5)
#         plt.plot(time_axis[start:end],plot_lines[0][start:end],alpha=1.0,color="gold",linewidth=0.5)
#         plt.legend(["nADP proteins", "nATP proteins", "Nd proteins", "Nde proteins"],loc="best")

#     if ("-all" in sys.argv) and (p==NflE):
#         plt.fill_between(time_axis[start:end],[0]*(end-start),plot_lines[0][start:end],alpha=1.0,facecolor="blue")
#         plt.fill_between(time_axis[start:end],plot_lines[0][start:end],plot_lines[1][start:end],alpha=1.0,facecolor="blue")
#         plt.plot(time_axis[start:end],plot_lines[0][start:end],alpha=1.0,color="blue",linewidth=0.5)
#         plt.plot(time_axis[start:end],plot_lines[1][start:end],alpha=1.0,color="blue",linewidth=0.5)
#         plt.legend(["nE proteins", "Nde proteins"],loc="best")

#     all_string = ""
#     if "-all" in sys.argv:
#         for i in range(4):
#             plt.plot(time_axis[start:end],plot_lines[i][start:end],color="black",linewidth=0.5)
#         all_string = "-all"

#     #if short is present, flag it in the file name
#     short_string = ""
#     if "-short" in sys.argv:
#         short_string = "-short"
#     if "-all" in sys.argv:
#         plt.xlim(time_axis[start],max(time_axis[start:end]))
#         plt.ylim(0,1.4*max(top_line[start:end]))
#         plt.savefig(load.print_string("endcap-plot"+short_string+all_string,p))

print "Got here"





# sprintf(outfilenameDivisions, "

# #load protein data.
# NflE = load.data(protein="NflE")
# NflD = load.data(protein="NflD")
# nATP = load.data(protein="nATP")
# nE = load.data(protein="nE")
# nADP = load.data(protein="nADP")
# Nd = load.data(protein="Nd")
# Nde = load.data(protein="Nde")

# #proteins is the main protein list. this program makes one graph for each element in proteins
# proteins = [NflE, NflD]

# #sub proteins are the individual, optionally-plotted proteins.
# sub_proteins_D = [nADP, nATP, Nd, Nde]
# sub_proteins_E = [nE, Nde]

# for p in proteins:
#     #pick a relatively arbitrary box divider
#     divider_left = 2 + (round(p.datashape[0]/2)-2)
#     divider_right = p.datashape[1] - (round(p.datashape[0]/2)-2)-1
#     print "Here!"
#     print p.datashape[1]
#     print p.datashape[0]
#     print divider_left
#     print divider_right
#     p.proteins_left = []
#     p.proteins_mid = []
#     p.proteins_right = []

#     #collect the total number of proteins in each box, at each time step
#     for file in p.dataset:
#         left_sum = 0
#         mid_sum = 0
#         right_sum = 0
#         for i in range(p.datashape[0]-1):
#             for j in range(p.datashape[1]-1):
#                 if j <= divider_left:
#                     left_sum += file[i][j]
#                 elif divider_left < j <= divider_right:
#                     mid_sum += file[i][j]
#                 else:
#                     right_sum += file[i][j]
#         p.proteins_left += [left_sum]
#         p.proteins_mid += [mid_sum]
#         p.proteins_right += [right_sum]

#     #gather the individual proteins to stick on the same graph as NflE.
#     if p==NflE:
#         for individual_proteins in sub_proteins_E:
#             #density to numbe rconversion for nE.
#             dV = 1
#             if individual_proteins == nE:
#                 dV = load.dx**3
#             #use the same method as before to store proteins in each bo.x
#             individual_proteins.proteins_left = []
#             individual_proteins.proteins_mid = []
#             individual_proteins.proteins_right = []
#             for file in individual_proteins.dataset:
#                 left_sum = 0
#                 mid_sum = 0
#                 right_sum = 0
#                 for i in range(individual_proteins.datashape[0]-1):
#                     for j in range(individual_proteins.datashape[1]-1):
#                         if j <= divider_left:
#                             left_sum += file[i][j]
#                         elif divider_left < j <= divider_right:
#                             mid_sum += file[i][j]
#                         else:
#                             right_sum += file[i][j]
#                 individual_proteins.proteins_left += [left_sum*dV]
#                 individual_proteins.proteins_mid += [mid_sum*dV]
#                 individual_proteins.proteins_right += [right_sum*dV]
#         NflE.lines = [nE.proteins_left, Nde.proteins_left, \
#                           nE.proteins_mid, Nde.proteins_mid, \
#                           nE.proteins_right, Nde.proteins_right]

#     #gather the individual proteins to stick on the same graph as NflD.
#     if p==NflD:
#         for individual_proteins in sub_proteins_D:
#             #density to number conversion for nATP and nADP.
#             dV = 1
#             if (individual_proteins == nATP) or (individual_proteins == nADP):
#                 dV = load.dx**3
#             #use the same method as before to store proteins in each box.
#             individual_proteins.proteins_left = []
#             individual_proteins.proteins_mid = []
#             individual_proteins.proteins_right = []
#             for file in individual_proteins.dataset:
#                 left_sum = 0
#                 mid_sum = 0
#                 right_sum = 0
#                 for i in range(individual_proteins.datashape[0]-1):
#                     for j in range(individual_proteins.datashape[1]-1):
#                         if j <= divider_left:
#                             left_sum += file[i][j]
#                         elif divider_left < j <= divider_right:
#                             mid_sum += file[i][j]
#                         else:
#                             right_sum += file[i][j]
#                 individual_proteins.proteins_left += [left_sum*dV]
#                 individual_proteins.proteins_mid += [mid_sum*dV]
#                 individual_proteins.proteins_right += [right_sum*dV]
#         #store them in one big list for plotting purposes later. the order is important, since we will be
#         #adding each element in this list to the one before it, starting with the first one (for the stacked plot)
#         NflD.lines = [nADP.proteins_left, nATP.proteins_left, Nd.proteins_left, Nde.proteins_left, \
#                           nADP.proteins_mid, nATP.proteins_mid, Nd.proteins_mid, Nde.proteins_mid, \
#                           nADP.proteins_right, nATP.proteins_right, Nd.proteins_right, Nde.proteins_right]

#     #calculate the indicies of the first period or so
#     #if -all and -short are both not present in the system arguments,
#     #plot all of the data.
#     (start, end) = (0, -1)
#     if ("-all" in sys.argv) or ("-short" in sys.argv):
#         (start, end) = oscillation_period(proteins[0].proteins_mid)
#         print oscillation_period(proteins[0].proteins_mid)

#     #plot them and save the figures.
#     plt.figure()
#     plt.xlabel("Time")
#     plt.ylabel("%s count"%p.protein)
#     plt.title("%s count vs time"%p.protein)


#     #plot_lines is the list where we will store each line to be plotted. the first element is set to
#     #NflD/NflE.lines[0], which was defined in the above loop. the list comprehension just adds the ith
#     #list to the (i-1)th list.
#     plot_lines = [p.lines[0]]
#     for i in range(1,len(p.lines)):
#         plot_lines += [[j+k for (j,k) in zip(plot_lines[i-1], p.lines[i])]]

#     #use the same method to get NflD/E stacked on top of each other. lowest = left, mid = mid, top = right
#     lowest_line = p.proteins_left
#     middle_line = [i+j for (i,j) in zip(lowest_line,p.proteins_mid)]
#     top_line = [i+j for (i,j) in zip(middle_line,p.proteins_right)]

#     #create a time axis
#     time_axis = list(np.arange(0,len(p.proteins_left)*5,5))

#     #matplotlib method for filling in color between two lines. left = blue, mid = green, right = red.
#     if "-all" not in sys.argv:
#         plt.plot(time_axis[start:end],lowest_line[start:end],color="blue",linewidth=0.5)
#         plt.plot(time_axis[start:end],middle_line[start:end],color="green",linewidth=0.5)
#         plt.plot(time_axis[start:end],top_line[start:end],color="red",linewidth=0.5)
#         plt.fill_between(time_axis[start:end],0,lowest_line[start:end],alpha=0.5)
#         plt.fill_between(time_axis[start:end],lowest_line[start:end],middle_line[start:end],alpha=0.5,facecolor="green")
#         plt.fill_between(time_axis[start:end],middle_line[start:end],top_line[start:end],alpha=0.5,facecolor="red")

#     if ("-all" in sys.argv) and (p==NflD):
#         plt.fill_between(time_axis[start:end],[0]*(end-start),plot_lines[0][start:end],alpha=0.25,facecolor="blue")
#         plt.fill_between(time_axis[start:end],plot_lines[0][start:end],plot_lines[1][start:end],alpha=0.5,facecolor="blue")
#         plt.fill_between(time_axis[start:end],plot_lines[1][start:end],plot_lines[2][start:end],alpha=0.75,facecolor="blue")
#         plt.fill_between(time_axis[start:end],plot_lines[2][start:end],plot_lines[3][start:end],alpha=1.0,facecolor="blue")

#         plt.fill_between(time_axis[start:end],plot_lines[3][start:end],plot_lines[4][start:end],alpha=0.25,facecolor="red")
#         plt.fill_between(time_axis[start:end],plot_lines[4][start:end],plot_lines[5][start:end],alpha=0.5,facecolor="red")
#         plt.fill_between(time_axis[start:end],plot_lines[5][start:end],plot_lines[6][start:end],alpha=0.75,facecolor="red")
#         plt.fill_between(time_axis[start:end],plot_lines[6][start:end],plot_lines[7][start:end],alpha=1.0,facecolor="red")

#         plt.fill_between(time_axis[start:end],plot_lines[7][start:end],plot_lines[8][start:end],alpha=0.25,facecolor="green")
#         plt.fill_between(time_axis[start:end],plot_lines[8][start:end],plot_lines[9][start:end],alpha=0.5,facecolor="green")
#         plt.fill_between(time_axis[start:end],plot_lines[9][start:end],plot_lines[10][start:end],alpha=0.75,facecolor="green")
#         plt.fill_between(time_axis[start:end],plot_lines[10][start:end],plot_lines[11][start:end],alpha=1.0,facecolor="green")

#         plt.plot(time_axis[start:end],plot_lines[0][start:end],alpha=0.25,color="blue",linewidth=0.5)
#         plt.plot(time_axis[start:end],plot_lines[1][start:end],alpha=0.5,color="blue",linewidth=0.5)
#         plt.plot(time_axis[start:end],plot_lines[2][start:end],alpha=0.75,color="blue",linewidth=0.5)
#         plt.plot(time_axis[start:end],plot_lines[3][start:end],alpha=1.0,color="blue",linewidth=0.5)

#         plt.plot(time_axis[start:end],plot_lines[4][start:end],alpha=0.25,color="red",linewidth=0.5)
#         plt.plot(time_axis[start:end],plot_lines[5][start:end],alpha=0.5,color="red",linewidth=0.5)
#         plt.plot(time_axis[start:end],plot_lines[6][start:end],alpha=0.75,color="red",linewidth=0.5)
#         plt.plot(time_axis[start:end],plot_lines[7][start:end],alpha=1.0,color="red",linewidth=0.5)

#         plt.plot(time_axis[start:end],plot_lines[8][start:end],alpha=0.25,color="green",linewidth=0.5)
#         plt.plot(time_axis[start:end],plot_lines[9][start:end],alpha=0.5,color="green",linewidth=0.5)
#         plt.plot(time_axis[start:end],plot_lines[10][start:end],alpha=0.75,color="green",linewidth=0.5)
#         plt.plot(time_axis[start:end],plot_lines[11][start:end],alpha=1.0,color="green",linewidth=0.5)

#         plt.legend(["nADP proteins left", "nATP proteins left", "Nd proteins left", "Nde proteins left", \
#                           "nADP proteins mid", "nATP proteins mid", "Nd proteins mid", "Nde proteins mid", \
#                           "nADP proteins right", "nATP proteins right", "Nd proteins right", "Nde proteins right"],loc=9,ncol=2)


#     if ("-all" in sys.argv) and (p==NflE):
#         plt.fill_between(time_axis[start:end],[0]*(end-start),plot_lines[0][start:end],alpha=0.3,facecolor="blue")
#         plt.fill_between(time_axis[start:end],plot_lines[0][start:end],plot_lines[1][start:end],alpha=1.0,facecolor="blue")

#         plt.fill_between(time_axis[start:end],plot_lines[1][start:end],plot_lines[2][start:end],alpha=0.3,facecolor="red")
#         plt.fill_between(time_axis[start:end],plot_lines[2][start:end],plot_lines[3][start:end],alpha=1.0,facecolor="red")

#         plt.fill_between(time_axis[start:end],plot_lines[3][start:end],plot_lines[4][start:end],alpha=0.3,facecolor="green")
#         plt.fill_between(time_axis[start:end],plot_lines[4][start:end],plot_lines[5][start:end],alpha=1.0,facecolor="green")

#         plt.plot(time_axis[start:end],plot_lines[0][start:end],alpha=0.3,color="blue",linewidth=0.5)
#         plt.plot(time_axis[start:end],plot_lines[1][start:end],alpha=1.0,color="blue",linewidth=0.5)

#         plt.plot(time_axis[start:end],plot_lines[2][start:end],alpha=0.3,color="red",linewidth=0.5)
#         plt.plot(time_axis[start:end],plot_lines[3][start:end],alpha=1.0,color="red",linewidth=0.5)

#         plt.plot(time_axis[start:end],plot_lines[4][start:end],alpha=0.3,color="green",linewidth=0.5)
#         plt.plot(time_axis[start:end],plot_lines[5][start:end],alpha=1.0,color="green",linewidth=0.5)

#         plt.legend(["nE proteins left", "Nde proteins left", "nE proteins mid", "Nde proteins mid","nE proteins right", "Nde proteins right"],loc=9,ncol=2)

#     #if -all is present, include all the lines, not just the usual 3, and flag it in the file name
#     all_string = ""
#     if "-all" in sys.argv:
#         all_string = "-all"

#     #if short is present, flag it in the file name
#     short_string = ""
#     if "-short" in sys.argv:
#         short_string = "-short"

#     plt.xlim(time_axis[start],max(time_axis[start:end]))
#     plt.ylim(0,1.7*max(top_line[start:end]))
#     if "-all" not in sys.argv:
#         plt.legend(["left density", "middle density", "right density"],loc="best")
#     plt.savefig(load.print_string("box-plot"+short_string+all_string,p))

#     plt.figure()
#     plt.xlabel("Time")
#     plt.ylabel("%s count"%p.protein)
#     plt.title("%s count vs time at end caps of cell"%p.protein)

#     if ("-all" in sys.argv) and (p==NflD):
#         plt.fill_between(time_axis[start:end],[0]*(end-start),plot_lines[0][start:end],alpha=1.0,facecolor="blue")
#         plt.fill_between(time_axis[start:end],plot_lines[0][start:end],plot_lines[1][start:end],alpha=1.0,facecolor="red")
#         plt.fill_between(time_axis[start:end],plot_lines[1][start:end],plot_lines[2][start:end],alpha=1.0,facecolor="green")
#         plt.fill_between(time_axis[start:end],plot_lines[2][start:end],plot_lines[3][start:end],alpha=1.0,facecolor="gold")
#         plt.plot(time_axis[start:end],plot_lines[0][start:end],alpha=1.0,color="blue",linewidth=0.5)
#         plt.plot(time_axis[start:end],plot_lines[0][start:end],alpha=1.0,color="red",linewidth=0.5)
#         plt.plot(time_axis[start:end],plot_lines[0][start:end],alpha=1.0,color="green",linewidth=0.5)
#         plt.plot(time_axis[start:end],plot_lines[0][start:end],alpha=1.0,color="gold",linewidth=0.5)
#         plt.legend(["nADP proteins", "nATP proteins", "Nd proteins", "Nde proteins"],loc="best")

#     if ("-all" in sys.argv) and (p==NflE):
#         plt.fill_between(time_axis[start:end],[0]*(end-start),plot_lines[0][start:end],alpha=1.0,facecolor="blue")
#         plt.fill_between(time_axis[start:end],plot_lines[0][start:end],plot_lines[1][start:end],alpha=1.0,facecolor="blue")
#         plt.plot(time_axis[start:end],plot_lines[0][start:end],alpha=1.0,color="blue",linewidth=0.5)
#         plt.plot(time_axis[start:end],plot_lines[1][start:end],alpha=1.0,color="blue",linewidth=0.5)
#         plt.legend(["nE proteins", "Nde proteins"],loc="best")

#     all_string = ""
#     if "-all" in sys.argv:
#         for i in range(4):
#             plt.plot(time_axis[start:end],plot_lines[i][start:end],color="black",linewidth=0.5)
#         all_string = "-all"

#     #if short is present, flag it in the file name
#     short_string = ""
#     if "-short" in sys.argv:
#         short_string = "-short"
#     if "-all" in sys.argv:
#         plt.xlim(time_axis[start],max(time_axis[start:end]))
#         plt.ylim(0,1.4*max(top_line[start:end]))
#         plt.savefig(load.print_string("endcap-plot"+short_string+all_string,p))

