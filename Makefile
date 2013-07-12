#CXXFLAGS = -g -O3
CXXFLAGS = -g -Wall -Werror

.SUFFIXES: .tex .dvi .ps .bib .bbl .pdf .fig .eps .aux .jpg .png .svg \
		.gp .mf .2602gf .pl .xgr

#OBJECTS = grid.h showDensity.

all : protein_microscopy paper/paper.pdf

#FIGURES=$(patsubst %.pdf)

clean:
	rm -f *.o protein

# plots:
# 	sh runplots.sh randst 0.50 6.00 8.00 98.00 15.00
# 	sh runplots.sh p 1.5 1.0 1.1 0.0 15.0
# 	sh runplots.sh p 1.5 1.0 0.0 0.0 15.0
# ne-randst-0.50-6.00-8.00-98.00-15.00-296.dat

#protein: test.o #foo.o

ALL_FIGURES = data/shape-randst/plots/time_map_compare-randst-10-60-60-990-150.pdf \

paper/paper.pdf: paper/paper.tex $(ALL_FIGURES)
	sh paper/protein_paper.sh

data/shape-randst/plots/time_map_compare-randst-10-60-60-990-150.pdf: pyplots/time_map.py
	python pyplots/time_map.py randst 1.00 6.00 6.00 99.00 15.00

figs/dadz-3-3.pdf: figs/plot-da_dz.py
	python figs/plot-da_dz.py .3 3.005


protein: protein.cpp protein.h
	g++ -g protein.cpp protein.h -o protein -I /usr/include/eigen2 -fpermissive
test.o: test.h

#foo.o: protein.h

.PHONY: clean all

#CXXFLAGS=-g -O3

#OBJECTS=Grid.o foo.o

#all: grid

#clean:
#	rm -f *.o grid

#grid: $(OBJECTS)
#	g++ -o grid $(OBJECTS)


#Grid.o: foo.h
#foo.o: foo.h

#.PHONY: clean all

