protein
=======

protein_microscopy
--------------------------------
The command line argument syntax for running a simulation is as follows:

    ./protein_microscopy shape A B C D density_lopsidedness -optional_flags

Where the optional flags are:

    -hires: use a finer grid resolution (default: .05 microns^3, hires: .025 microns^3)
    -slice: print data for the middle yz plane (x = Nx/2) rather than data summed over x
    -dump: print out the state of the cell every few thousand time steps

The possible shapes are:

    p: pill
    randst: one of 4 custom shapes, specified by argument D
    triangle: triangle
    st: stadium shape. think a flattened pill.

The parameters A, B, C, and D have different meanings for each shape. They are:

batch submit
--------------------------------
To add a sim or plot job, edit 'jobs' and add a new line with the simulation parameters in the same format as if you were running protein_microscopy. Example:

   p 4.00 0.50 0.00 0.00 15.00 -slice -hires
   p 4.00 0.50 0.00 0.00 15.00 -hires
   randst 0.50 6.00 6.00 96.00 15.00
   triangle 0.50 1.00 1.00 1.00 20.00 -dump

The command to then run all of these jobs is:

   ./batch sim plot srun

Which will then simulate each of the jobs, wait for them to finish, and then plot them (using slurm). Both the sim and plot arguments are optional (i.e., you can sim and not plot, or plot and not sim). It can filter jobs if you like, according to shape types and plot types:

   ./batch sim plot shapes="p randst triangle" plots="box_plot time_map" srun
   ./batch sim plot shapes="st triangle" plots="arrow_plot" srun

It will only run the jobs on slurm if the srun command is present (after ./batch).

plotting programs
--------------------------------
The three main plotting programs used are located in pyplots. They are box_plot.py, arrow_plot.py, and time_map.py. They take command line arguments in exactly the same way as running protein_microscopy. For example:

   python pyplots/box_plot.py p 4.00 0.50 0.00 0.00 15.00 -slice -hires
   python pyplots/box_plot.py p 4.00 0.50 0.00 0.00 15.00 -hires

   python pyplots/time_map.py p 4.00 0.50 0.00 0.00 15.00 -slice -hires
   python pyplots/time_map.py p 4.00 0.50 0.00 0.00 15.00 -slice

Each plot works this way.
