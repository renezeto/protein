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

dx = 0.15
hires_str = ""
if "-hires" in sys.argv:
    dx=.05
    hires_str="hires-"

difD = 2.5 # (um)^2 s^- 1
time_step = .1*dx*dx/difD

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
time_axis = np.arange(0,time_step*len(div_data[:,0]),time_step)

#These are what is printed out by the c++ file and some index numbers that state the column they're in in the dat file
# 0 = NATPtotal_left, 1 = NADPtotal_left, 2 = Ndtotal_left, 3 = Ndetotal_left, 4 = NEtotal_left, 5 = NflDtotal_left, 6 = NflEtotal_left,
# 7 = NATPtotal_mid,  8 = NADPtotal_mid, 9 = Ndtotal_mid, 10 = Ndetotal_mid, 11 = NEtotal_mid, 12 = NflDtotal_mid, 13 = NflEtotal_mid,
# 14 = NATPtotal_right, 15 = NADPtotal_right, 16 = Ndtotal_right, 17 = Ndetotal_right, 18 = NEtotal_right, 19 = NflDtotal_right, 20 = NflEtotal_right);

indices_from_cpp_D = [1,2,3,7,8,9,10,17,16,15,14]

b=0
plot_lines_D = [div_data[:,0]]
for a in indices_from_cpp_D:
    plot_lines_D += [[j+k for (j,k) in zip(plot_lines_D[b], div_data[:,a])]]
    print a
    b+=1


plt.figure()
plt.xlabel("Time")
plt.ylabel("Number of MinD Proteins in Section")
plt.title("Number of MinD Proteins in Section Vs Time")

plt.fill_between(time_axis,[0]*len(time_axis),plot_lines_D[0],alpha=0.25,facecolor="blue")
plt.fill_between(time_axis,plot_lines_D[0],plot_lines_D[1],alpha=0.5,facecolor="blue")
plt.fill_between(time_axis,plot_lines_D[1],plot_lines_D[2],alpha=0.75,facecolor="blue")
plt.fill_between(time_axis,plot_lines_D[2],plot_lines_D[3],alpha=1.0,facecolor="blue")

plt.fill_between(time_axis,plot_lines_D[3],plot_lines_D[4],alpha=0.25,facecolor="red")
plt.fill_between(time_axis,plot_lines_D[4],plot_lines_D[5],alpha=0.5,facecolor="red")
plt.fill_between(time_axis,plot_lines_D[5],plot_lines_D[6],alpha=0.75,facecolor="red")
plt.fill_between(time_axis,plot_lines_D[6],plot_lines_D[7],alpha=1.0,facecolor="red")

plt.fill_between(time_axis,plot_lines_D[7],plot_lines_D[8],alpha=1.0,facecolor="green")
plt.fill_between(time_axis,plot_lines_D[8],plot_lines_D[9],alpha=0.75,facecolor="green")
plt.fill_between(time_axis,plot_lines_D[9],plot_lines_D[10],alpha=0.5,facecolor="green")
plt.fill_between(time_axis,plot_lines_D[10],plot_lines_D[11],alpha=0.25,facecolor="green")

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

plt.ylim(0,1.4*plot_lines_D[11][100])
plt.legend(["nADP proteins left", "nATP proteins left", "Nd proteins left", "Nde proteins left", \
                          "nADP proteins mid", "nATP proteins mid", "Nd proteins mid", "Nde proteins mid", \
                          "nADP proteins right", "nATP proteins right", "Nd proteins right", "Nde proteins right"],loc=9,ncol=2,prop={'size':10})


plt.savefig("./data/shape-"+f_shape+"/plots/m"+hires_str+"-divisions-MinD-"+f_shape+"-"+f_param1+"-"+f_param2+"-"+f_param3+"-"+f_param4+"-"+f_param5+".pdf")

#These are what is printed out by the c++ file and some index numbers that state the column they're in in the dat file
# 0 = NATPtotal_left, 1 = NADPtotal_left, 2 = Ndtotal_left, 3 = Ndetotal_left, 4 = NEtotal_left, 5 = NflDtotal_left, 6 = NflEtotal_left,
# 7 = NATPtotal_mid,  8 = NADPtotal_mid, 9 = Ndtotal_mid, 10 = Ndetotal_mid, 11 = NEtotal_mid, 12 = NflDtotal_mid, 13 = NflEtotal_mid,
# 14 = NATPtotal_right, 15 = NADPtotal_right, 16 = Ndtotal_right, 17 = Ndetotal_right, 18 = NEtotal_right, 19 = NflDtotal_right, 20 = NflEtotal_right);

indices_from_cpp_E = [4,10,11,18,17]

b=0
plot_lines_E = [div_data[:,3]]
for a in indices_from_cpp_E:
    plot_lines_E += [[j+k for (j,k) in zip(plot_lines_E[b], div_data[:,a])]]
    print a
    b+=1

plt.figure()
plt.xlabel("Time")
plt.ylabel("Number of MinE Proteins in Section")
plt.title("Number of MinE Proteins in Section Vs Time")

plt.fill_between(time_axis,[0]*len(time_axis),plot_lines_E[0],alpha=0.5,facecolor="blue")
plt.fill_between(time_axis,plot_lines_E[0],plot_lines_E[1],alpha=1.0,facecolor="blue")

plt.fill_between(time_axis,plot_lines_E[1],plot_lines_E[2],alpha=0.5,facecolor="red")
plt.fill_between(time_axis,plot_lines_E[2],plot_lines_E[3],alpha=1.0,facecolor="red")

plt.fill_between(time_axis,plot_lines_E[3],plot_lines_E[4],alpha=1.0,facecolor="green")
plt.fill_between(time_axis,plot_lines_E[4],plot_lines_E[5],alpha=0.5,facecolor="green")

plt.plot(time_axis,plot_lines_E[0],alpha=0.25,color="blue",linewidth=0.5)
plt.plot(time_axis,plot_lines_E[1],alpha=0.5,color="blue",linewidth=0.5)

plt.plot(time_axis,plot_lines_E[2],alpha=0.25,color="red",linewidth=0.5)
plt.plot(time_axis,plot_lines_E[3],alpha=0.5,color="red",linewidth=0.5)

plt.plot(time_axis,plot_lines_E[4],alpha=0.25,color="green",linewidth=0.5)
plt.plot(time_axis,plot_lines_E[5],alpha=0.5,color="green",linewidth=0.5)

plt.ylim(0,1.1*plot_lines_D[5][100])
plt.legend(["NE proteins left", "Nde proteins left", "NE proteins middle", "Nde proteins middle", "NE proteins right", "Nde proteins right"],loc=9,ncol=2,prop={'size':14})


plt.savefig("./data/shape-"+f_shape+"/plots/m"+hires_str+"-divisions-MinE-"+f_shape+"-"+f_param1+"-"+f_param2+"-"+f_param3+"-"+f_param4+"-"+f_param5+".pdf")


time_axis = np.arange(0,time_step*len(compare_data[:,0]),time_step)
plt.figure()
plt.xlabel("Time")
plt.ylabel("Ave Number of MinD Proteins Per Area")

plt.plot(time_axis,compare_data[:,0],linewidth=0.5)
plt.plot(time_axis,compare_data[:,1],linewidth=0.5)
plt.title("Ave Number of MinD Proteins Per Area in Middle and on Right Vs Time")
plt.legend(["Right section", "Middle section"],loc='best',prop={'size':10})

plt.savefig("./data/shape-"+f_shape+"/plots/m"+hires_str+"-compare-MinD-"+f_shape+"-"+f_param1+"-"+f_param2+"-"+f_param3+"-"+f_param4+"-"+f_param5+".pdf")

print "DONE!!!! :) "
