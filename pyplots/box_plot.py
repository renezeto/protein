from __future__ import division
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
import pylab
import file_loader_dump as load

## WIP!! 

#this is the ordering of proteins in the box plot list
# sprintf(proteinList[0]->name,"nATP");
# sprintf(proteinList[1]->name,"nADP");
# sprintf(proteinList[2]->name,"nE");
# sprintf(proteinList[3]->name,"ND");
# sprintf(proteinList[4]->name,"NDE");
#starting with nATP, it's 3 lines (one for each box), break, then nADP, etc
#each column (moving right) represents a time step

#box_data = np.loadtxt("../data/shape-%s/box-plot--%s-%s-%s-%s-%s.dat"%(load.f_shape,load.f_param1,load.f_param2,load.f_param3,load.f_param4,load.f_param5))

def returnData(protein,boxName):
    with open("../data/shape-%s/box-plot--%s-%s-%s-%s-%s.dat"%(load.f_shape,load.f_param1,load.f_param2,load.f_param3,load.f_param4,load.f_param5),"r") as boxData:
        proteinsOverTime = [line for line in f if (protein in line) and (boxName in line)]
    return proteinsOvertime

#*args is a tuple of the function input separated by commas
def stackData(plotList):
    stackedPlotList = []
    for proteinData in plotList:
        splitString = proteinData.split('-')
        (protein, boxName) = (splitString[0], splitString[1])
        stackedPlotList.append(returnData(protein, boxName))
    stackedPlotList = ["list comprehension from old box plot wip"]
        

def main():
    returnData("nATP","left")
    return 0

if __name__ == '__main__':
    main()
