#CXXFLAGS = -g -O3
CXXFLAGS = -g -Wall -Werror

.SUFFIXES: .tex .dvi .ps .bib .bbl .pdf .fig .eps .aux .jpg .png .svg \
		.gp .mf .2602gf .pl .xgr

#OBJECTS = grid.h showDensity.

all : protein_microscopy paper/paper.pdf

#FIGURES=$(patsubst %.pdf)

clean:
	rm -f *.o protein

plots:
	sh runplots.sh randst 0.50 6.00 8.00 98.00 15.00
	sh runplots.sh p 1.5 1.0 1.1 0.0 15.0
	sh runplots.sh p 1.5 1.0 0.0 0.0 15.0
	ne-randst-0.50-6.00-8.00-98.00-15.00-296.dat

protein: test.o #foo.o

# ALL_FIGURES = data/shape-randst/plots/time-map-compare-randst-10-60-60-990-150.pdf \
# 	data/shape-randst/plots/time-map-compare-randst-10-80-60-980-150.pdf \
# 	data/shape-randst/plots/time-map-compare-randst-10-60-60-970-150.pdf \
# 	data/shape-randst/plots/time-map-compare-randst-10-60-60-960-150.pdf \
# 	data/shape-p/plots/time-map-compare-p-100-5-0-0-150.pdf \
# 	data/shape-p/plots/time-map-compare-p-40-5-0-0-150.pdf \
# 	data/shape-p/plots/time-map-compare-p-40-20-0-0-150.pdf \
# 	data/shape-p/plots/time-map-compare-p-40-30-0-0-150.pdf \


paper/paper.pdf: paper/paper.tex $(ALL_FIGURES)
	cd paper && pdflatex paper.tex && bibtex paper && pdflatex paper.tex && pdflatex paper.tex

# data/shape-randst/plots/time-map-compare-randst-10-60-60-990-150.pdf: pyplots/time_map.py
# 	python pyplots/time_map.py randst 1.00 6.00 6.00 99.00 15.00
# data/shape-randst/plots/time-map-compare-randst-10-80-60-980-150.pdf: pyplots/time_map.py
# 	python pyplots/time_map.py randst 1.00 8.00 6.00 98.00 15.00
# data/shape-randst/plots/time-map-compare-randst-10-60-60-970-150.pdf: pyplots/time_map.py
# 	python pyplots/time_map.py randst 1.00 6.00 6.00 97.00 15.00
# data/shape-randst/plots/time-map-compare-randst-10-60-60-960-150.pdf: pyplots/time_map.py
# 	python pyplots/time_map.py randst 1.00 6.00 6.00 96.00 15.00
# data/shape-p/plots/time-map-compare-p-100-5-0-0-150.pdf: pyplots/time_map.py
# 	python pyplots/time_map.py p 10.00 0.50 0.00 0.00 15.00
# data/shape-p/plots/time-map-compare-p-40-5-0-0-150.pdf: pyplots/time_map.py
# 	python pyplots/time_map.py p 4.00 0.50 0.00 0.00 15.00
# data/shape-p/plots/time-map-compare-p-40-20-0-0-150.pdf: pyplots/time_map.py
# 	python pyplots/time_map.py p 4.00 2.00 0.00 0.00 15.00
# data/shape-p/plots/time-map-compare-p-40-30-0-0-150.pdf: pyplots/time_map.py
# 	python pyplots/time_map.py p 4.00 3.00 0.00 0.00 15.00


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

