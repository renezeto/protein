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

batch_randst_simulations = [[1.00, 1.00, 1.00, 96.00], \
                                [1.00, 1.00, 1.00, 97.00], \
                                [1.00, 1.00, 1.00, 98.00], \
                                [1.00, 1.00, 1.00, 99.00]]

#simulation subprocess management:
if '-sim' in sys.argv:
    if 'p' in sys.argv:
        processes = set()
        max_processes = 7
        for job in batch_pill_simulations:
            processes.add(subprocess.Popen(['./protein_microscopy','p',str(job[0]),str(job[1]),'0.00','0.00','15.00']))
            if len(processes) >= max_processes:
                os.wait()
                processes.difference_update(p for p in processes if p.poll() is not None)
    if 'randst' in sys.argv:
        processes = set()
        max_processes = 7
        for job in batch_randst_simulations:
            processes.add(subprocess.Popen(['./protein_microscopy','randst',str(job[0]),str(job[1]),str(job[2]),str(job[3]),'15.00']))
            if len(processes) >= max_processes:
                os.wait()
                processes.difference_update(p for p in processes if p.poll() is not None)

#plots the recent simulations:
if '-plot' in sys.argv:
    processes = set()
    max_processes = 7
    if 'p' in sys.argv:
        for job in batch_pill_simulations:
            processes.add(subprocess.Popen(['python','pyplots/extrema.py','p',str(job[0])+'0',str(job[1])+'0','0.00','0.00','15.00']))
            processes.add(subprocess.Popen(['python','pyplots/time_map.py','p',str(job[0])+'0',str(job[1])+'0','0.00','0.00','15.00']))
            if len(processes) >= max_processes:
                os.wait()
                processes.difference_update(p for p in processes if p.poll() is not None)
    if 'randst' in sys.argv:
        for job in batch_pill_simulations:
            processes.add(subprocess.Popen(['python','pyplots/extrema.py','p',str(job[0])+'0',str(job[1])+'0','0.00','0.00','15.00']))
            processes.add(subprocess.Popen(['python','pyplots/time_map.py','randst',str(job[0]),str(job[1]),str(job[2]),str(job[3]),'15.00']))
            if len(processes) >= max_processes:
                os.wait()
                processes.difference_update(p for p in processes if p.poll() is not None)
