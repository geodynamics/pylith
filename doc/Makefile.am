# -*- Makefile -*-
#
# ----------------------------------------------------------------------
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# <LicenseText>
#
# ----------------------------------------------------------------------
#

manual.pdf: manual.xml
	docbook2pdf manual.xml

view:
	xpdf manual.pdf&

clean:
	-rm manual.pdf

.PHONY: view clean

# version
# $Id$

# End of file 
