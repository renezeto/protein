from __future__ import division
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
import pylab
import file_loader as load

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
def returnData(boxName,proteinType):

    #open the data file, grab the line with the correct protein type and box partition, load it as a [string] (so we can use list comprehensions)
    with open("./data/shape-%s/%s%s%sbox-plot--%s-%s-%s-%s-%s-%s.dat"%(load.f_shape,load.debug_str,load.hires_str,load.slice_str,load.f_shape,load.f_param1,load.f_param2,load.f_param3,load.f_param4,load.f_param5),"r") as boxData:
        proteinsOverTime = [line for line in boxData if (proteinType in line) and (boxName in line)]

    #format the string so that it is a list of numbers (split on tab, pop off keywords and newlines, convert str -> float)
    proteinsOverTime = proteinsOverTime[0].split('\t')
    proteinsOverTime = proteinsOverTime[2:-1]
    proteinsOverTime = [float(i) for i in proteinsOverTime]
    output = np.array(proteinsOverTime)
    return (output)


def main():

    with open("./data/shape-%s/%s%s%sbox-plot--%s-%s-%s-%s-%s-%s.dat"%(load.f_shape,load.debug_str,load.hires_str,load.slice_str,load.f_shape,load.f_param1,load.f_param2,load.f_param3,load.f_param4,load.f_param5),"r") as boxData:
        fileLines = boxData.readlines()

    #get number of boxes and protein types. little hokey but it works. in boxData.readlines(), there is exactly one '\n' newline string
    #for each protein type block. therefor, the number of protein types is equal to the number of times "\n" appears by itself in the list.
    numProteinTypes = len([line for line in fileLines if line=="\n"])
    numNewLines = numProteinTypes

    #it follows that the total number of lines in the data file, minus the number of blank lines in the data file, is equal to
    #the number of protein types * the number of box types. divide by number of protein types to get number of box types.
    numBoxes = (len(fileLines) - numNewLines)/numProteinTypes

    #grab the names of the proteins used, and the names of the boxes
    proteinTypeList = []
    boxList = []
    for line in fileLines:
        if (line != "\n"):
            proteinTypeList += [line.split("\t")[0]]
            boxList += [line.split("\t")[1]]

    #prune duplicates
    proteinTypeList = list(set(proteinTypeList))
    boxList = list(set(boxList))
 
    #generate list of proteinType and box combinations to feed into returnData
    plotNameList_D = []
    plotNameList_E = []
    numProteinTypes_D = 0
    numProteinTypes_E = 0
    numPlots_D = 0
    numPlots_E = 0
    for box in boxList:
        for proteinType in proteinTypeList:
            if "D_" in proteinType and "n" not in proteinType:
                plotNameList_D += ["%s-%s"%(box,proteinType)]
                numPlots_D += 1
            if "E_" in proteinType and "n" not in proteinType:
                plotNameList_E += ["%s-%s"%(box,proteinType)]
                numPlots_E += 1

    for proteinType in proteinTypeList:
        if "D_" in proteinType and "n" not in proteinType:
            numProteinTypes_D += 1
        if "E_" in proteinType and "n" not in proteinType:
            numProteinTypes_E +=1

    plotCurveList_D = []
    plotCurveList_E = []
    for proteinData in plotNameList_D:
        splitString = proteinData.split('-')
        (protein, boxName) = (splitString[0], splitString[1])
        plotCurveList_D += [returnData(boxName,protein)]
    for proteinData in plotNameList_E:
        splitString = proteinData.split('-')
        (protein, boxName) = (splitString[0], splitString[1])
        plotCurveList_E += [returnData(boxName,protein)]

    print len(plotCurveList_D[:][0])
#pass plotNameList through stackData to generate the list of line data to be plotted
    #plotCurveList_D = stackData(plotNameList_D)
    #plotCurveList_E = stackData(plotNameList_E)

    #get a time axis for the plot from the length of one of the data sets we have
    timeAxis = range(len(plotCurveList_D[0]))
    #print len(plotCurveList_D[:,0])

    #begin messy code (to deal with matplotlib) - don't judge me
    (start, end) = (int(5.2*int(len(timeAxis)/10)),int(6*int(len(timeAxis)/10)))

    #plot scales. colors limited for now.
    colorScale = ["b","g","r","c","m","y"]
    alphaScale_D = [n/numProteinTypes_D for n in range(1,numProteinTypes_D+1)]
    alphaScale_E = [n/numProteinTypes_E for n in range(1,numProteinTypes_E+1)]


                 #generate the plot
    plt.figure()
    j=0
    k=0
    for i in range(numPlots_D):
        if i%(numProteinTypes_D)==0:
            j+=1
            k=0
        print i
        print k
        plt.plot(timeAxis[start:end],
                 plotCurveList_D[i][start:end],
                 color=colorScale[j],alpha=alphaScale_D[k])
        k+=1
    plt.xlim(start,end+.40*(end-start))
    #plt.ylim(0,whatever)
    plt.title("Average MinE concentration on walls over time")
    plt.xlabel("Time (s)")
    plt.ylabel("proteins per uint Area")
    plt.legend(plotNameList_D,loc="best",prop={'size':10})
    plt.savefig(load.print_string("ave-plot_D",""))

    plt.figure()
    j=0
    k=0
    for i in range(numPlots_E):
        if i%(numProteinTypes_E)==0:
            j+=1
            k=0
        print i
        print len(plotCurveList_E[0])
        plt.plot(timeAxis[start:end],
                 plotCurveList_E[i][start:end],
                 color=colorScale[j])
        k+=1
    plt.xlim(start,end+.40*(end-start))
    #plt.ylim(0,500)
    plt.title("Average MinD concentration on walls over time")
    plt.xlabel("Time (s)")
    plt.ylabel("Proteins per unit Area")
    plt.legend(plotNameList_E,loc="best",prop={'size':10})
    plt.savefig(load.print_string("ave-plot_E",""))
    return 0

if __name__ == '__main__':
    main()


