CXXFLAGS = -g -Wall -Werror

.SUFFIXES: .tex .dvi .ps .bib .bbl .pdf .fig .eps .aux .jpg .png .svg \
		.gp .mf .2602gf .pl .xgr

all: sim #paper/paper.pdf ${ALL_FIGURES}

test: test.cpp

clean: rm -f protein_microscopy paper/paper.pdf

sim: protein_microscopy

#ALL_FIGURES = \
	data/shape-p/plots/debug-box-plot_D--p-150-60-0-0-1500.pdf

#paper/paper.pdf: paper/paper.tex #${ALL_FIGURES}
#	cd paper && pdflatex paper.tex && bibtex paper && pdflatex paper.tex && pdflatex paper.tex


#${ALL_FIGURES}: batch jobs $(wildcard data/shape-*/*.dat)
#	./batch plot include="-paper"

#time maps
data/shape-randst/plots/time-map-compare-randst-10-60-60-990-150.pdf: pyplots/time_map.py
	python pyplots/time_map.py randst 1.00 6.00 6.00 99.00 15.00
data/shape-randst/plots/time-map-compare-randst-10-80-60-980-150.pdf: pyplots/time_map.py
	python pyplots/time_map.py randst 1.00 8.00 6.00 98.00 15.00
data/shape-randst/plots/time-map-compare-randst-10-60-60-970-150.pdf: pyplots/time_map.py
	python pyplots/time_map.py randst 1.00 6.00 6.00 97.00 15.00
data/shape-randst/plots/time-map-compare-randst-10-60-60-960-150.pdf: pyplots/time_map.py
	python pyplots/time_map.py randst 1.00 6.00 6.00 96.00 15.00
data/shape-p/plots/time-map-compare-p-100-5-0-0-150.pdf: pyplots/time_map.py
	python pyplots/time_map.py p 10.00 0.50 0.00 0.00 15.00
data/shape-p/plots/time-map-compare-p-40-5-0-0-150.pdf: pyplots/time_map.py
	python pyplots/time_map.py p 4.00 0.50 0.00 0.00 15.00
data/shape-p/plots/time-map-compare-p-40-20-0-0-150.pdf: pyplots/time_map.py
	python pyplots/time_map.py p 4.00 2.00 0.00 0.00 15.00
data/shape-p/plots/time-map-compare-p-40-30-0-0-150.pdf: pyplots/time_map.py
	python pyplots/time_map.py p 4.00 3.00 0.00 0.00 15.00

#arrow plots
#box plots
#frequency plots
