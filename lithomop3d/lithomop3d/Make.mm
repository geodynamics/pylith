# -*- Makefile -*-
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                             Charles A. Williams
#                       Rensselaer Polytechnic Institute
#                        (C) 2003  All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

PROJECT = lithomop3d

#--------------------------------------------------------------------------
#

all: export

release: tidy
	cvs release .

update: clean
	cvs update .

#--------------------------------------------------------------------------
#
# export

EXPORT_PYTHON_MODULES = \
	Application.py \
	Lithomop3d_run.py \
	Lithomop3d_scan.py \
	Lithomop3d_setup.py \
	__init__.py


export:: export-python-modules


# version
# $Id: Make.mm,v 1.1 2004/04/14 21:22:47 willic3 Exp $

# Generated automatically by MakeMill on Tue Mar  2 17:05:23 2004

# End of file 
