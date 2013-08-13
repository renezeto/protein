from __future__ import division
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
import pylab
import file_loader_dump as load

## WIP!!

#reads a special kind of data file printed by protein_microscopy.cpp
#format is:
# --------
# proteinType1 boxName1 n(t=0) n(t=1) n(t=2) ...
# proteinType1 boxName2 n(t=0) n(t=1) n(t=2) ...
# proteinType1 boxName3 n(t=0) n(t=1) n(t=2) ...
#
# proteinType2 boxName1 n(t=0) n(t=1) n(t=2) ...
# proteinType2 boxName2 n(t=0) n(t=1) n(t=2) ...
# proteinType2 boxName3 n(t=0) n(t=1) n(t=2) ...
# --------
#where n(t) is the number of proteins of one type at time t

#opens the file and grabs a particular line matching proteinType and boxName. returns list of protein counts at each time.
def returnData(proteinType,boxName):

    #open the data file, grab the line with the correct protein type and box partition, load it as a [string] (so we can use list comprehensions)
    with open("./data/shape-%s/box-plot--%s-%s-%s-%s-%s.dat"%(load.f_shape,load.f_param1,load.f_param2,load.f_param3,load.f_param4,load.f_param5),"r") as boxData:
        proteinsOverTime = [line for line in boxData if (proteinType in line) and (boxName in line)]

    #format the string so that it is a list of numbers (split on tab, pop off keywords and newlines, convert str -> float)
    proteinsOverTime = proteinsOverTime[0].split('\t')
    proteinsOverTime = proteinsOverTime[2:-1]
    proteinsOverTime = [float(i) for i in proteinsOverTime]

    return proteinsOverTime

#takes input format: ["proteinType1-boxNum1","proteinType1-boxnum2",proteinType2-boxnum1"...]. will return a list of lists
#in the stacking order specified by the input (first entry is at the bottom).
def stackData(plotList):

    #parse the input
    tempList = []
    for proteinData in plotList:
        splitString = proteinData.split('-')
        (protein, boxName) = (splitString[0], splitString[1])
        tempList += [returnData(protein, boxName)]

    #"stack" the lists
    stackedPlotList = [tempList[0]]
    for i in range(1,len(tempList)):
        stackedPlotList += [[j+k for (j,k) in zip(tempList[i-1],tempList[i])]]

    return stackedPlotList

def main():
    return 0

if __name__ == '__main__':
    main()
