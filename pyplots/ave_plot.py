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

def find_period(f):
    """
      Find the period of a function that is described by the input
      array f, and return indices for a start and end range for one
      period.  If we cannot find the period accurately, just return
      the entire range.
    """
    # first we look at the fft to get a guess at the period (probably
    # not *too* accurate or too bad).
    fk = np.fft.fft(f)
    fk[0] = 0
    kmax = 1
    fkmax = np.abs(fk[:int(len(fk)/2)]).max()
    for i in xrange(1,int(len(fk)/2)):
        if np.abs(fk[i]) == fkmax:
            kmax = i
            break
    #print 'kmax is', kmax
    period_estimate = len(f)/kmax
    #plt.plot(np.abs(fk))
    #plt.figure()
    if kmax < 5:
        return (0, len(f))
    # now we locate the final minimum of the function.
    lastmin = len(f)-2
    while f[lastmin] > f[lastmin+1] or f[lastmin] > f[lastmin-1]:
        lastmin -= 1
    # and last (but not least), we locate the second-to-last
    # (penultimate) minimum, which should have a very similar value to
    # the final minimum.
    penultimate_min = lastmin - int(period_estimate*.7)
    while f[penultimate_min] > f[penultimate_min+1] or f[penultimate_min] > f[penultimate_min-1] or np.abs(f[penultimate_min]/f[lastmin]-1) > 0.01:
        penultimate_min -= 1
    #return (0, len(f) - 1)
    if penultimate_min < 0:
        return (0, len(f))
    return (penultimate_min, lastmin)


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

    #generate list of proteinType and box combinations to feed into ckData
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

    for proteinType in proteinTypeList:
        if "D_" in proteinType and "n" not in proteinType:
            numProteinTypes_D += 1

    plotCurveList_D = []
    for proteinData in plotNameList_D:
        splitString = proteinData.split('-')
        (protein, boxName) = (splitString[0], splitString[1])
        plotCurveList_D += [returnData(boxName,protein)]

    #get a time axis for the plot from the length of one of the data sets we have
    #timeAxis = range(len(plotCurveList_D[0]))
    difD = 2.5 # (um)^2 s^- 1
    time_step = .1*load.dx*load.dx/difD #sec
    print_denominator = 1000 #This is from the c++ I wanted to format things the same here.
    box_time_step = time_step*print_denominator
    timeAxis = np.linspace(0,box_time_step*len(plotCurveList_D[0]),len(plotCurveList_D[0]))


    #begin messy code (to deal with matplotlib) - don't judge me
    #(start, end) = (int(5.3*int(len(timeAxis)/10)),int(5.7*int(len(timeAxis)/10)))
    (start, end) = find_period(plotCurveList_D[3])

    #integral of proteins over time
    for i in range(numPlots_D):
        integral = 0
        for j in range(len(plotCurveList_D[0])):
            integral += plotCurveList_D[i][j]

    #plot scales. colors limited for now.
    colorScale = ["b","g","r","c","m","y"]
    alphaScale_D = [n/numProteinTypes_D for n in range(1,numProteinTypes_D+1)]

    # yLimit = max(plotCurveList_D[0][start:end])

                 #generate the plot
    plt.figure()
    j=0
    k=0
    for i in range(numPlots_D):
        if i%(numProteinTypes_D)==0:
            j+=1
            k=0
        plt.plot(timeAxis[start:end],
                 plotCurveList_D[i][start:end],
                 color=colorScale[j],alpha=alphaScale_D[k])
        k+=1
    plt.xlim(timeAxis[start],timeAxis[end-1])
    #plt.ylim(0,10000)
    plt.title("Min D protein counts over time")
    plt.xlabel("Time (s)")
    plt.ylabel("Fraction of proteins")
    plt.legend(plotNameList_D,loc="best",prop={'size':10})
    plt.savefig(load.print_string("ave-plot_D",""))

    return 0

if __name__ == '__main__':
    main()


