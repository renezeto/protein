from __future__ import division
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import sys
import time
import file_loader as load
import Image

f_shape = sys.argv[1]
f_param1 = sys.argv[2]
f_param2 = sys.argv[3]
f_param3 = sys.argv[4]
f_param4 = sys.argv[5]
f_param5 = sys.argv[6]
start_time = float(sys.argv[7])
end_time = float(sys.argv[8])


dump_time_step = 2.5
#computes the maximum of a two dimensional array - REPLACE with method
def maxnum(page):
    Z = [0. for i in range(page.shape[0])]
    for i in range(page.shape[0]):
        Z[i] = max(page[i])
    maxval = max(Z)
    return maxval

#computes the minimum of a two dimensional array - REPLACE with method
def minnum(page):
    Z = [0. for i in range(page.shape[0])]
    for i in range(page.shape[0]):
        for j in range(len(page[i])):
            if page[i][j] == 0:
             page[i][j] = 50000000 #eh
        Z[i] = min(page[i])
    minval = min(Z)
    return minval

#computes the global maximum over a set of two dimensional arrays (stored as files)
def timemax(protein):
    Z = [0. for i in range(protein.tsteps)]
    for i in (range(int(end_time/dump_time_step)-int(start_time/dump_time_step))):
        Z[i] = maxnum(protein.dataset[i])
    maxval = max(Z)
    return maxval

#computes the global minimum over a set of two dimensional arrays (stored as files)
def timemin(protein):
    Z = [0. for i in range(protein.tsteps)]
    for i in (range(int(end_time/dump_time_step)-int(start_time/dump_time_step))):
        Z[i] = minnum(protein.dataset[i])
    minval = min(Z)
    return minval

def build(proteins,proteinList):
    total_length = int((end_time/dump_time_step)+1)- int(start_time/dump_time_step)
    spread = Image.new("RGB", (120+300*total_length, 300*len(proteins)), "white")
    for i in range(len(proteins)):
        if (proteinList[i]=="ND" or proteinList[i]=="NDE"):
            maxval = timemax(proteins[i])
        else:
            maxval = timemax(proteins[i])
        minval = timemin(proteins[i])
        plt.figure()
        Z, Y = np.meshgrid(np.arange(0,proteins[i].datashape[1],1), np.arange(0,proteins[i].datashape[0],1))
    #generate a sequence of .png's for each file (printed time step). these will be used to create a gif.
        for k in range(int(end_time/dump_time_step)-int(start_time/dump_time_step)):
            print "1st prints"
            print k
            print maxval
            print minval
            sys.stdout.flush()
            page = proteins[i].dataset[k]
            plt.clf()
            plt.axes().set_aspect('equal', 'datalim')
            CS = plt.contourf(Z, Y, page, cmap=plt.cm.jet,origin='lower',levels=np.arange(minval,maxval,.1))
        #plt.axis('off')
            plt.axes().get_xaxis().set_ticks([])
            plt.axes().get_yaxis().set_ticks([])
            plt.xlabel('%s  %gs' %  (proteinList[i], float(k*dump_time_step)) )
            if (k == int(start_time/dump_time_step)):
                plt.ylabel('%s' % proteinList[i])
            plt.savefig('./data/shape-'+f_shape+'/pngs/'+load.debug_str+load.hires_str+load.slice_str+'tmp_%03d'%k+'-'+str(proteins[i].protein)+ \
                            '-'+f_shape+'-'+f_param1+'-'+f_param2+'-'+f_param3+'-'+f_param4+'-'+f_param5+'.png',dpi=50)
            orig = Image.open("data/shape-"+f_shape+"/pngs/"+load.debug_str+load.hires_str+load.slice_str+'tmp_%03d'%k \
                               +"-"+str(proteins[i].protein)+"-"+f_shape+"-"+f_param1+"-"+f_param2 \
                               +"-"+f_param3+"-"+f_param4+"-"+f_param5+".png")
            init_box = orig.getbbox()
            cut = 60
            print "2nd prints"
            print cut
            print init_box[2]
            print init_box[3]
            sys.stdout.flush()
            box = (cut,0,init_box[2]-cut,init_box[3])
            if (k == int(start_time/dump_time_step)):
                box = (0,0,init_box[2]-cut,init_box[3])
            cropped = orig.crop(box)
            del orig
            if (k == int(start_time/dump_time_step)):
                spread.paste(cropped, (0,box[3]*i))
            else:
                spread.paste(cropped, (cut+(init_box[2]-2*cut)*(k-int(start_time/dump_time_step)),box[3]*i))
    spread.save(load.print_string("image-plot",""))
    #spread.save("data/shape-"+f_shape +"/plots/"+load.debug_str+load.hires_str+load.slice_str+"image-" \
    #                        +f_shape+"-"+f_param1+"-"+f_param2 \
    #                        +"-"+f_param3+"-"+f_param4+"-"+f_param5+".pdf")

proteinList = ["nADP","nATP","ND","nE","NDE"]#["nE","nATP","nADP","ND"]#,"NDE"]
proteins = [0]*len(proteinList)

for i in range(len(proteinList)):
    proteins[i] = load.data(proteinList[i],start_time,end_time)

build(proteins,proteinList)
#may need to use --mem-per-cpu=4000 after the srun for this

    # NflE = load.data(protein="NflE")
    #NflD = load.data(protein="NflD")
    # nATP = load.data(protein="nATP")
    # nE = load.data(protein="nE")
    # nADP = load.data(protein="nADP")
    # ND = load.data(protein="ND")
#NDE = load.data(protein="NDE")


