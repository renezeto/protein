#!/bin/sh
python pyplots/time_map.py $1 $2 $3 $4 $5 $6
python pyplots/extrema.py $1 $2 $3 $4 $5 $6
python pyplots/show_area_rating.py $1 $2 $3 $4 $5 $6