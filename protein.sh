#!/bin/bash
srun -p debian ./protein randst 0.5 3.5 3.5 2 &
srun -p debian ./protein randst 0.5 3.5 3.5 3 &
srun -p debian ./protein randst 0.5 3.5 3.5 4 &
srun -p debian ./protein randst 0.5 3.5 3.5 5 &
srun -p debian ./protein randst 0.5 3.5 3.5 6 &
srun -p debian ./protein randst 0.5 3.5 3.5 7 &
srun -p debian ./protein randst 0.5 3.5 3.5 8 &
srun -p debian ./protein randst 0.5 3.5 3.5 9 &




