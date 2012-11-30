#!/bin/bash
pdflatex protein.tex
bibtex protein
pdflatex protein.tex
pdflatex protein.tex
evince protein.pdf
