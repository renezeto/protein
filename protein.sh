#!/bin/bash
#for d in 99 98 97 96; do
#srun -p debian ./protein randst 0.5 6.0 6.0 $d
#done

srun ./protein randst 0.50 6.0 6.0 99 &
srun ./protein randst 0.50 6.0 8.0 98 &
srun ./protein randst 0.50 6.0 6.0 97 &
srun ./protein randst 0.50 6.0 6.0 96 &
srun ./protein p 8.0 0.5 0.0 0.0 &
srun ./protein p 8.0 0.25 0.0 0.0 &
srun ./protein p 10.0 0.5 0.0 0.0 &
srun ./protein p 6.0 0.5 0.0 0.0 &
srun ./protein p 6.0 0.25 0.0 0.0 &
srun ./protein p 2.0 0.5 0.0 0.0 &
srun ./protein p 1.0 0.5 0.0 0.0 &
srun ./protein triangle 0.50 6.0 6.0 0.0 &
srun ./protein triangle 0.50 6.0 6.0 0.0 &





