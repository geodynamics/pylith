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
# Copyright (c) 2010 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

if [ $1 == "clean" ]; then
    latexmk -C
else
    latexmk -pdf -pdflatex="pdflatex -interaction=nonstopmode" -use-make userguide.tex
fi
									       
# End of file
