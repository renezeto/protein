from __future__ import division
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
import pylab
import file_loader as load
import matplotlib.patheffects
from matplotlib.font_manager import FontProperties


from mpl_toolkits.axes_grid.anchored_artists import AnchoredSizeBar

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

#opens the file and grabs a particular line matching proteinType and
#boxName. returns list of protein counts at each time.
def returnData(boxName,proteinType):

    #open the data file, grab the line with the correct protein type and box partition, load it as a [string] (so we can use list comprehensions)
    with open("./data/shape-%s/%s%s%sbox-plot--%s-%s-%s-%s-%s-%s.dat"%
              (load.f_shape,load.debug_str,load.hires_str,load.slice_str,
               load.f_shape,load.f_param1,load.f_param2,load.f_param3,
               load.f_param4,load.f_param5),"r") as boxData:
        proteinsOverTime = [line for line in boxData if (proteinType in line) and (boxName in line)]


    #format the string so that it is a list of numbers (split on tab, pop off keywords and newlines, convert str -> float)
    proteinsOverTime = proteinsOverTime[0].split('\t')
    proteinsOverTime = proteinsOverTime[2:-1]
    proteinsOverTime = [float(i) for i in proteinsOverTime]
    return proteinsOverTime

#takes input format:
#["proteinType1-boxNum1","proteinType1-boxnum2",proteinType2-boxnum1"...]. will
#return a list of lists in the stacking order specified by the input
#(first entry is at the bottom).
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

    print ""
    print "plotNameList before ", plotNameList_D, "\n"
    new_plotNameList_D = [0]*len(plotNameList_D)
    P_Ord = [3,0,2,1,7,4,6,5,11,8,10,9]
    if load.f_param4 == '97.00':
        P_Ord = [3,0,2,1,11,8,10,9,15,12,14,13,7,4,6,5]
    if load.f_param4 == '96.00':
        P_Ord = [15,12,14,13,3,0,2,1,7,4,6,5,11,8,10,9]
    for i in range(len(P_Ord)):
        new_plotNameList_D[i] = plotNameList_D[P_Ord[i]]
    for i in range(len(plotNameList_D)):
        plotNameList_D[i] = new_plotNameList_D[i]
    print "plotNameList after ",plotNameList_D,"\n"
    plotProteinLabels = ['MinD:ATP (cyto)','MinD:ATP (mem)','MinE:MinD:ATP','MinD:ADP (cyto)']



    #pass plotNameList through stackData to generate the list of line data to be plotted
    plotCurveList_D = stackData(plotNameList_D)
    plotCurveList_E = stackData(plotNameList_E)

    #get a time axis for the plot from the length of one of the data sets we have
    difD = 2.5 # (um)^2 s^- 1
    time_step = .1*load.dx*load.dx/difD #sec
    print_denominator = 1000 #This is from the c++ I wanted to format things the same here.
    box_time_step = time_step*print_denominator
    timeAxis = np.linspace(0,box_time_step*len(plotCurveList_D[0]),len(plotCurveList_D[0]))

    #begin messy code (to deal with matplotlib) - don't judge me

    start_time_as_frac_of_ten = float(load.f_param6)
    end_time_as_frac_of_ten = float(load.f_param7)
    tot_time = len(plotCurveList_D[0])*box_time_step
    start = int(tot_time*start_time_as_frac_of_ten/10.0/box_time_step)
    end = int(tot_time*end_time_as_frac_of_ten/10.0/box_time_step)
    (start, end) = find_period(plotCurveList_D[3])
    if load.f_param4 == '97.00' or (load.f_shape == 'triangle' and load.f_param3 == '6.01'):
        if load.f_param4 == '97.00':
            start_time_as_frac_of_ten = 0
            end_time_as_frac_of_ten = 2.3
        if load.f_shape == 'triangle' and load.f_param3 == '6.01':
            start_time_as_frac_of_ten = 5.00
            end_time_as_frac_of_ten = 9.00
        start = int(tot_time*start_time_as_frac_of_ten/10.0/box_time_step)
        end = int(tot_time*end_time_as_frac_of_ten/10.0/box_time_step)


    #print set(plotCurveList_D[1]).union(set(plotCurveList_D[2]))

    arrow_start = start+int(0.05*(end-start))
    arrow_x = [arrow_start]*len(plotCurveList_D)
    for i in range(arrow_start,arrow_start + int((end-start)/100)):
            if ((plotCurveList_D[0][arrow_x[0]]-0) < (plotCurveList_D[0][i]-0)):
                arrow_x[0] = i
    for type in range(1,len(plotCurveList_D)):
        for i in range(arrow_start,arrow_start + int((end-start)/100)):
            if ((plotCurveList_D[type][arrow_x[type]]-plotCurveList_D[type-1][arrow_x[type]]) < (plotCurveList_D[type][i]-plotCurveList_D[type-1][i])):
                arrow_x[type] = i

    # arrow_x = [100]*len(plotCurveList_D)
    # print len(plotCurveList_D)
    # for type in range(len(plotCurveList_D)):
    #     for i in range(start,end):
    #         if ((plotCurveList_D[type][i]-plotCurveList_D[type][i-1])>.1):
    #             arrow_x[type] = i
    #             break
    #             #print type,i
    # print arrow_x


    # arrow_start, arrow_end = find_period(plotCurveList_D[3])
    # arrow_x = start + int(.25*(arrow_end - arrow_start))
    # print arrow_x

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
    #f, (bax,sectionax) = plt.subplots(1, 2)
    bax = plt.subplot2grid((2,5), (0,0), colspan=4, rowspan=2)
    sectionax = plt.subplot2grid((2,5), (0,4), colspan=1,rowspan=2)

    # first plot the section data...
    sectiondata = np.loadtxt("data/shape-%s/membrane_files/%s%s%ssections-%s-%s-%s-%s-%s-%s.dat"
                             % (load.f_shape,load.debug_str,load.hires_str,load.slice_str,load.f_shape,
                                load.f_param1,load.f_param2,load.f_param3,load.f_param4,load.f_param5))

    def plot_sections(sectionax, sectiondata):
        dx = load.dx
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
        mycolors = ["w","g","r","m","c","y"]
        for i in xrange(min(4, len(boxList))):
            if boxList[i] == 'Right':
                mycolors[1] = colorScale[i]
            if boxList[i] == 'Mid':
                mycolors[2] = colorScale[i]
            if boxList[i] == 'Left':
                mycolors[3] = colorScale[i]
        mycolors = colorScale[1:]
        if load.f_param4 == '97.00':
            mycolors = ['g','r','m','c']
        if load.f_param4 == '96.00':
            #rightup = 2, rightdown = 1, leftup = 4, leftdown = 3
            mycolors = ['g','r','c','m']
        print mycolors
        # here we rotate so that the order of sections will match the
        # box plot.
        xdir, ydir = xweighted - xmean, yweighted - ymean
        xdir, ydir = xdir/np.sqrt(xdir**2+ydir**2), ydir/np.sqrt(xdir**2+ydir**2)
        extraxspace = .5
        extrayspace = 0
        Yrotated = X*xdir + Y*ydir
        Xrotated = X*ydir - Y*xdir
        sectionax.contourf(Xrotated, Yrotated, sectiondata, levels=levels, colors=mycolors)
        xmin = Xrotated[sectiondata>0].min()
        xmax = Xrotated[sectiondata>0].max()
        ymin = Yrotated[sectiondata>0].min()
        ymax = Yrotated[sectiondata>0].max()
        sectionax.set_xlim(xmin-extraxspace, xmax)
        sectionax.set_ylim(ymin-extrayspace, ymax)
        sectionax.set_aspect('equal')
        sectionax.set_frame_on(False)
        sectionax.axes.get_xaxis().set_visible(False)
        sectionax.axes.get_yaxis().set_visible(False)
        sectionax.add_artist(AnchoredSizeBar(
                sectionax.transData,
                1.00, # length of the bar in the data reference
                "1$\mu$", # label of the bar
                #bbox_to_anchor=(0.,0.,1.,1.),
                loc=8, # 'best', # location (lower right)
                pad=-(ymax-ymin)/2.0 + 0.4, borderpad=0.25, sep=3,
                frameon=False
                ))
    plot_sections(sectionax, sectiondata)
    section_names = ['Left Section','Center Section','Right Section']
    if load.f_param4 == '97.00':
        section_names = ['Lower Section','Middle Left Section','Middle Right Section','Upper Section']
#        section_names = ['rightup','mid','left','rightdown']
    if load.f_param4 == '96.00':
        section_names = ['Lower Left Section','Lower Right Section','Upper Left Section','Upper Right Section']
#        section_names = ['rightdown','rightup','leftdown','leftup']
    font=FontProperties()
    font.set_family('serif')
    text_adjust = -.2*box_time_step*(end-start)
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
            y_text_label = i*.8/len(plotCurveList_D[:,0]) + .1*np.floor(i/numProteinTypes_D)
            if load.f_param4 == '97.00' or load.f_param4 == '96.00':
                y_text_label = i*.8/len(plotCurveList_D[:,0]) + .07*np.floor(i/numProteinTypes_D)
            y_label = (plotCurveList_D[i, arrow_x[i]])/2.0
            bax.annotate('%s'%plotProteinLabels[i],xy=(arrow_x[i]*box_time_step,y_label),xytext=(start*box_time_step+text_adjust,y_text_label),
                         fontsize=7,
                         fontproperties=font,
                         arrowprops=dict(facecolor='black',shrink=0.05, width=.3, headwidth=5.))
            bax.fill_between(timeAxis[start:end],
                             [0 for x in range(len(timeAxis))[start:end]],
                             plotCurveList_D[i, start:end],
                             alpha=alphaScale_D[k],facecolor=colorScale[j])
        elif i!=0:
            bax.plot(timeAxis[start:end],
                     plotCurveList_D[i,start:end],
                     color=colorScale[j],alpha=alphaScale_D[k])
            y_text_label = i*.8/len(plotCurveList_D[:,0]) + .1*np.floor(i/numProteinTypes_D)
            y_label = (plotCurveList_D[i, arrow_x[i]] + plotCurveList_D[i-1, arrow_x[i]])/2.0
            if load.f_param4 == '97.00' or load.f_param4 == '96.00':
                y_text_label = i*.8/len(plotCurveList_D[:,0]) + .07*np.floor(i/numProteinTypes_D)
            bax.annotate('%s'%plotProteinLabels[i%numProteinTypes_D],xy=(arrow_x[i]*box_time_step,y_label),xytext=(start*box_time_step+text_adjust,y_text_label),
                         fontsize=7,
                         fontproperties=font,
                         arrowprops=dict(facecolor='black',shrink=0.05, width=.3, headwidth=5.))
            if (i+1)%(numProteinTypes_D)==0:
                bax.text(-0.2,y_text_label+.04,section_names[int(np.floor(i/numProteinTypes_D))],transform=bax.transAxes,fontsize=9,fontproperties=font,)
            bax.fill_between(timeAxis[start:end],
                             plotCurveList_D[i-1, start:end],
                             plotCurveList_D[i, start:end],
                             alpha=alphaScale_D[k],facecolor=colorScale[j])
            k+=1
    bax.set_xlim(box_time_step*start,box_time_step*end)
    bax.get_yaxis().set_visible(False)
    bax.set_ylim(0, 1)
    bax.set_title("MinD protein counts over time")
    bax.set_xlabel("Time (s)")
    rax = bax.twinx()
    rax.set_ylabel('Fraction of proteins in each stage and section',labelpad=-15)
    rax.yaxis.set_ticklabels([0,"","","","",1.0])
#bax.set_ylabel("Fraction of proteins")


# 'A', xy=(Az, Ax), xytext=(1.2,-3.5),
#          path_effects=texteff,
#          arrowprops=dict(shrink=0.01, width=1,
#                         headwidth=hw, path_effects=arroweff))

    #bax.legend(plotNameList_D,bbox_to_anchor=(0.3,-0.05,1.0,1.0),loc=4,prop={'size':8}).draw_frame(False)


    plt.savefig(load.print_string("box-plot_D",""))
    plt.figure()

    #f, (bax,sectionax) = plt.subplots(1, 2)
    bax = plt.subplot2grid((2,5), (0,0), colspan=4, rowspan=2)
    sectionax = plt.subplot2grid((2,5), (0,4), colspan=1,rowspan=2)

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
    bax.set_xlim(box_time_step*start,box_time_step*end)
    bax.set_ylim(0, 1)
    bax.set_title("MinE protein counts over time")
    bax.set_xlabel("Time (s)")
    bax.set_ylabel("Fraction of proteins")
    bax.legend(plotNameList_E,bbox_to_anchor=(0.3,-0.05,1.0,1.0),loc="lower right",prop={'size':8}).draw_frame(False)
    plt.savefig(load.print_string("box-plot_E",""))

    plt.show()
    return 0

if __name__ == '__main__':
    main()
