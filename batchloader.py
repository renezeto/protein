from __future__ import division
import numpy as np
import sys
import os
import subprocess
import time

#set up simulations to be run:
batch_pill_simulations = []
i=1
while i<20:
    a = abs(i/2)
    b = abs(5-i/2)
    if ([a, b] not in batch_pill_simulations) and ([b, a] not in batch_pill_simulations):
        batch_pill_simulations+=[[a,b]]
    i+=1

#simulation subprocess management:
if '-sim' in sys.argv:
    processes = set()
    max_processes = 7
    for job in batch_pill_simulations:
        processes.add(subprocess.Popen(['./protein_microscopy','p',str(job[0]),str(job[1]),'0.00','0.00','15.00']))
        if len(processes) >= max_processes:
            os.wait()
            processes.difference_update(p for p in processes if p.poll() is not None)

if '-plot' in sys.argv:
    processes = set()
    max_processes = 14
    for job in batch_pill_simulations:
        processes.add(subprocess.Popen(['python','pyplots/extrema.py','p',str(job[0]),str(job[1]),'0.00','0.00','15.00']))
        processes.add(subprocess.Popen(['python','pyplots/time_map.py','p',str(job[0]),str(job[1]),'0.00','0.00','15.00']))
        if len(processes) >= max_processes:
            os.wait()
            processes.difference_update(p for p in processes if p.poll() is not None)
