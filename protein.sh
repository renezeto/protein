#!/bin/bash
#for d in 99 98 97 96; do
#srun -p debian ./protein randst 0.5 6.0 6.0 $d 
#done

srun ./protein randst 0.5 6.0 6.0 99 &
srun ./protein randst 0.5 6.0 6.0 98 &
srun ./protein randst 0.5 6.0 6.0 97 &
srun ./protein randst 0.5 6.0 6.0 96 &


