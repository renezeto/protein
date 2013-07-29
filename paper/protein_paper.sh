#!/bin/bash
pdflatex paper/paper.tex
bibtex paper/paper
pdflatex paper/paper.tex
pdflatex paper/paper.tex

