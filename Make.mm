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
PROJECT = pylith
PACKAGE = pylith

RECURSE_DIRS = \
	pylith3d \

# pylith2d \
# pylithop \
# pylithax \
# pylithsph \

#--------------------------------------------------------------------------
#

all:
	BLD_ACTION="all" $(MM) recurse

PROJ_CLEAN =
clean::
	BLD_ACTION="clean" $(MM) recurse

tidy::
	BLD_ACTION="tidy" $(MM) recurse


# version
# $Id: Make.mm,v 1.1 2004/04/14 21:07:29 willic3 Exp $

# Generated automatically by MakeMill on Tue Mar  2 16:26:04 2004

# End of file 
