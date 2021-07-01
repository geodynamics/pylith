#!/bin/bash
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2017 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------
#

if [ $# == 1 ]; then
    if [ $1 == "clean" ]; then
        latexmk -C
    elif [ $1 == "cover" ]; then
        pdflatex coveronly.tex && convert coveronly.pdf -background white -flatten -resize 250 -quality 95 cover/cover_small.jpg
    fi
else
    latexmk -pdf -pdflatex="pdflatex -interaction=nonstopmode" -use-make userguide.tex
fi
									       
# End of file
