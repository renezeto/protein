#!/bin/bash
srun -p debian ./protein randst 0.5 3.5 3.5 12 &
srun -p debian ./protein randst 0.5 3.5 3.5 13 &
srun -p debian ./protein randst 0.5 3.5 3.5 14 &
srun -p debian ./protein randst 0.5 3.5 3.5 15 &
srun -p debian ./protein randst 0.5 3.5 3.5 16 &
srun -p debian ./protein randst 0.5 3.5 3.5 17 &
srun -p debian ./protein randst 0.5 3.5 3.5 18 &
srun -p debian ./protein randst 0.5 3.5 3.5 19 &




