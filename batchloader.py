from __future__ import division
import numpy as np
import sys
import os
##
batch_pill = []
i=1
while i<20:
    a = abs(i/2)
    b = abs(5-i/2)
    if ([a, b] not in batch_pill) and ([b, a] not in batch_pill):
        batch_pill+=[[a,b]]
    i+=1

for job in batch_pill:
    os.system("srun ./protein_microscopy p "+str(job[0])+" "+str(job[1])+" 0.00 0.00 15.00 > /dev/null")
