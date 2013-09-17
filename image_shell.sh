#!/bin/bash

srun python pyplots/image_plot.py randst 0.25 8.00 6.00 97.00 15.00 701.00 861.00
srun python pyplots/image_plot.py randst 0.25 6.00 6.00 96.00 15.00 701.00 861.00
srun python pyplots/image_plot.py randst 0.25 6.00 8.00 98.00 15.00 701.00 861.00
srun python pyplots/image_plot.py triangle 0.25 4.00 4.00 4.00 15.00 701.00 861.00
srun python pyplots/image_plot.py triangle 0.25 5.00 5.00 3.00 15.00 701.00 861.00
srun python pyplots/image_plot.py triangle 0.25 3.60 4.80 6.00 15.00 701.00 861.00
