from __future__ import division
import numpy as np
import sys
import os
import subprocess
import time

#set up simulations to be run:
batch_randst_simulations = []
i=1

batch_pill_simulations = [
    [10.00, 0.50], \
        [4.00, 0.50], \
        [4.00, 2.00], \
        [4.00, 3.00] ]

batch_randst_simulations = [
    [1.00,6.00,6.00,99.00], \
    [1.00,8.00,6.00,98.00], \
    [1.00,6.00,8.00,98.00], \
    [1.00,6.00,6.00,97.00], \
    [1.00,6.00,6.00,96.00]
    ]

#simulation subprocess management:
if '-sim' in sys.argv:
    processes = set()
    max_processes = 100
    if 'p' in sys.argv:
        for job in batch_pill_simulations:
            processes.add(subprocess.Popen(['srun','./protein_microscopy','p','%.2f'%job[0],'%.2f'%job[1],'0.00','0.00','15.00']))
            if len(processes) >= max_processes:
                os.wait()
                processes.difference_update(p for p in processes if p.poll() is not None)
    if 'randst' in sys.argv:
        for job in batch_randst_simulations:
            processes.add(subprocess.Popen(['srun','./protein_microscopy','randst','%.2f'%job[0],'%.2f'%job[1],'%.2f'%job[2],'%.2f'%job[3],'15.00']))
            if len(processes) >= max_processes:
                os.wait()
                processes.difference_update(p for p in processes if p.poll() is not None)

#plots the recent simulations:
if '-plot' in sys.argv:
    processes = set()
    max_processes = 100
    if 'p' in sys.argv:
        for job in batch_pill_simulations:
            processes.add(subprocess.Popen(['srun','python','pyplots/extrema.py','p','%.2f'%job[0],'%.2f'%job[1],'0.00','0.00','15.00']))
            processes.add(subprocess.Popen(['srun','python','pyplots/time_map.py','p','%.2f'%job[0],'%.2f'%job[1],'0.00','0.00','15.00']))
            processes.add(subprocess.Popen(['srun','python','pyplots/showatp.py','p','%.2f'%job[0],'%.2f'%job[1],'0.00','0.00','15.00']))
            if len(processes) >= max_processes:
                os.wait()
                processes.difference_update(p for p in processes if p.poll() is not None)
    if 'randst' in sys.argv:
        for job in batch_randst_simulations:
            processes.add(subprocess.Popen(['srun','python','pyplots/extrema.py','randst','%.2f'%job[0],'%.2f'%job[1],'%.2f'%job[2],'%.2f'%job[3],'15.00']))
            processes.add(subprocess.Popen(['srun','python','pyplots/time_map.py','randst','%.2f'%job[0],'%.2f'%job[1],'%.2f'%job[2],'%.2f'%job[3],'15.00']))
            processes.add(subprocess.Popen(['srun','python','pyplots/showatp.py','randst','%.2f'%job[0],'%.2f'%job[1],'%.2f'%job[2],'%.2f'%job[3],'15.00']))
            if len(processes) >= max_processes:
                os.wait()
                processes.difference_update(p for p in processes if p.poll() is not None)
print batch_randst_simulations
print processes

if '-area' in sys.argv:
    processes = set()
    max_processes = 7
    if 'p' in sys.argv:
        for job in batch_pill_simulations:
            processes.add(subprocess.Popen(['./protein_microscopy','p','%.2f'%job[0],'%.2f'%job[1],'0.00','0.00','15.00','-area']))
            if len(processes) >= max_processes:
                os.wait()
                processes.difference_update(p for p in processes if p.poll() is not None)
    if 'randst' in sys.argv:
        for job in batch_randst_simulations:
            processes.add(subprocess.Popen(['./protein_microscopy','randst','%.2f'%job[0],'%.2f'%job[1],'%.2f'%job[2],'%.2f'%job[3],'15.00','-area']))
            if len(processes) >= max_processes:
                os.wait()
                processes.difference_update(p for p in processes if p.poll() is not None)
