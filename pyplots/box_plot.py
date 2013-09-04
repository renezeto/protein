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
    return proteinsOverTime

#takes input format: ["proteinType1-boxNum1","proteinType1-boxnum2",proteinType2-boxnum1"...]. will return a list of lists
#in the stacking order specified by the input (first entry is at the bottom).
def stackData(plotList):
    #parse the input
    tempList = []
    for proteinData in plotList:
        splitString = proteinData.split('-')
        (protein, boxName) = (splitString[0], splitString[1])
        tempList += [returnData(boxName,protein)]

    #"stack" the lists
    stackedPlotList = [tempList[0]]
    for i in range(1,len(tempList)):
        stackedPlotList += [[j+k for (j,k) in zip(stackedPlotList[i-1],tempList[i])]]

    output = np.array(stackedPlotList)
    return output/output[len(output[:,0])-1, 0] # normalize output as a fraction of total

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

    #generate list of proteinType and box combinations to feed into stackData
    plotNameList_D = []
    plotNameList_E = []
    numProteinTypes_D = 0
    numProteinTypes_E = 0
    for box in boxList:
        for proteinType in proteinTypeList:
            if "D_" in proteinType:
                plotNameList_D += ["%s-%s"%(box,proteinType)]
            if "E_" in proteinType:
                plotNameList_E += ["%s-%s"%(box,proteinType)]

    #pass plotNameList through stackData to generate the list of line data to be plotted
    plotCurveList_D = stackData(plotNameList_D)
    plotCurveList_E = stackData(plotNameList_E)

    #get a time axis for the plot from the length of one of the data sets we have
    timeAxis = range(len(plotCurveList_D[0]))

    #begin messy code (to deal with matplotlib) - don't judge me

    start_time_as_frac_of_ten = float(load.f_param6)
    end_time_as_frac_of_ten = float(load.f_param7)
    # print load.f_shape
    # print start_time_as_frac_of_ten
    # print end_time_as_frac_of_ten
    (start, end) = (int(start_time_as_frac_of_ten*len(timeAxis)/10),int(end_time_as_frac_of_ten*len(timeAxis)/10))

    #get num on each plot
    for proteinType in proteinTypeList:
        if "D_" in proteinType:
            numProteinTypes_D += 1
        if "E_" in proteinType:
            numProteinTypes_E +=1

    #plot scales. colors limited for now.
    colorScale = ["b","g","r","c","m","y"]
    alphaScale_D = [n/numProteinTypes for n in range(1,numProteinTypes_D+1)]
    alphaScale_E = [n/numProteinTypes for n in range(1,numProteinTypes_E+1)]

    #generate the plot
    f, (bax,sectionax) = plt.subplots(1, 2)

    # first plot the section data...
    sectiondata = np.loadtxt("data/shape-%s/membrane_files/%s%s%ssections-%s-%s-%s-%s-%s-%s.dat"
                             % (load.f_shape,load.debug_str,load.hires_str,load.slice_str,load.f_shape,
                                load.f_param1,load.f_param2,load.f_param3,load.f_param4,load.f_param5))
    def plot_sections(sectionax, sectiondata):
        dx = 0.1 # FIXME
        x = np.arange(sectiondata.shape[1]*1.0)*dx
        y = np.arange(sectiondata.shape[0]*1.0)*dx
        X,Y = np.meshgrid(x,y)
        inmembrane = np.zeros_like(sectiondata)
        inmembrane[sectiondata>0] = 1.0
        xmax = X[sectiondata>0].max()
        xmin = X[sectiondata>0].min()
        ymax = Y[sectiondata>0].max()
        ymin = Y[sectiondata>0].min()
        ymean = (Y*inmembrane).sum()/inmembrane.sum()
        xmean = (X*inmembrane).sum()/inmembrane.sum()
        yweighted = (Y*sectiondata).sum()/sectiondata.sum()
        xweighted = (X*sectiondata).sum()/sectiondata.sum()
        levels = [0.5, 1.5, 2.5, 3.5, 4.5]
        mycolors = ["w","g","r","c","m","y"]
        for i in xrange(min(4, len(boxList))):
            if boxList[i] == 'Right':
                mycolors[1] = colorScale[i]
            if boxList[i] == 'Mid':
                mycolors[2] = colorScale[i]
            if boxList[i] == 'Left':
                mycolors[3] = colorScale[i]
            if boxList[i] == 'UpperLeft':
                mycolors[1] = colorScale[i]
            if boxList[i] == 'UpperRight':
                mycolors[2] = colorScale[i]
            if boxList[i] == 'LowerLeft':
                mycolors[3] = colorScale[i]
            if boxList[i] == 'LowerRight':
                mycolors[4] = colorScale[i]
        mycolors = colorScale[1:]
        # here we rotate and flip so that the order of sections will
        # match the box plot as well as we can manage.
        if yweighted - ymean > abs(xweighted - xmean):
            sectionax.contourf(X, Y, sectiondata, levels=levels, colors=mycolors)
            sectionax.set_xlim(xmin, xmax)
            sectionax.set_ylim(ymin, ymax)
        elif ymean - yweighted > abs(xweighted - xmean):
            sectionax.contourf(X, ymax - (Y - ymin), sectiondata, levels=levels, colors=mycolors)
            sectionax.set_xlim(xmin, xmax)
            sectionax.set_ylim(ymin, ymax)
        elif xweighted > xmean:
            sectionax.contourf(Y, X, sectiondata, levels=levels, colors=mycolors)
            sectionax.set_xlim(ymin, ymax)
            sectionax.set_ylim(xmin, xmax)
        else:
            sectionax.contourf(Y, xmax - (X - xmin), sectiondata, levels=levels, colors=mycolors)
            sectionax.set_xlim(ymin, ymax)
            sectionax.set_ylim(xmin, xmax)
        sectionax.set_aspect('equal')
    plot_sections(sectionax, sectiondata)

    j=0
    k=0
    for i in range(len(plotCurveList_D[:,0])):
        if i%(numProteinTypes_D)==0:
            j+=1
            k=0
        if i==0:
            bax.plot(timeAxis[start:end],
                       plotCurveList_D[i, start:end],
                       color=colorScale[j],alpha=alphaScale_D[k])
            bax.fill_between(timeAxis[start:end],
                             [0 for x in range(len(timeAxis))[start:end]],
                             plotCurveList_D[i, start:end],
                             alpha=alphaScale_D[k],facecolor=colorScale[j])
        elif i!=0:
            bax.plot(timeAxis[start:end],
                     plotCurveList_D[i, start:end],
                     color=colorScale[j],alpha=alphaScale_D[k])
            bax.fill_between(timeAxis[start:end],
                             plotCurveList_D[i-1, start:end],
                             plotCurveList_D[i, start:end],
                             alpha=alphaScale_D[k],facecolor=colorScale[j])
        #print "i is ",i," || k is", k," || j is",j
        k+=1
    bax.set_xlim(start,end+.40*(end-start))
    bax.set_ylim(0, 1)
    bax.set_title("Min D protein counts over time")
    bax.set_xlabel("Time (s)")
    bax.set_ylabel("Fraction of proteins")
    bax.legend(plotNameList_D,loc="best",prop={'size':10})


    plt.savefig(load.print_string("box-plot_D",""))

    f, (bax,sectionax) = plt.subplots(1, 2)

    # First plot the section data...
    plot_sections(sectionax, sectiondata)

    j=0
    k=0
    for i in range(len(plotCurveList_E)):
        if i%(numProteinTypes_E)==0:
            j+=1
            k=0
        if i==0:
            bax.plot(timeAxis[start:end],plotCurveList_E[i][start:end],color=colorScale[j],alpha=alphaScale_E[k])
            bax.fill_between(timeAxis[start:end],[0 for x in range(len(timeAxis))[start:end]],plotCurveList_E[i][start:end],alpha=alphaScale_E[k],facecolor=colorScale[j])
        elif i!=0:
            bax.plot(timeAxis[start:end],plotCurveList_E[i][start:end],color=colorScale[j],alpha=alphaScale_E[k])
            bax.fill_between(timeAxis[start:end],plotCurveList_E[i-1][start:end],plotCurveList_E[i][start:end],alpha=alphaScale_E[k],facecolor=colorScale[j])
        #print "i is ",i," || k is", k," || j is",j
        k+=1
    bax.set_xlim(start,end+.40*(end-start))
    bax.set_ylim(0, 1)
    bax.set_title("Min E protein counts over time")
    bax.set_xlabel("Time (s)")
    bax.set_ylabel("Fraction of proteins")
    bax.legend(plotNameList_E,loc="best",prop={'size':10})
    plt.savefig(load.print_string("box-plot_E",""))

    plt.plot([1,2,3], [4,5,6])
    return 0

if __name__ == '__main__':
    main()
