# -*- Makefile -*-
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                             Charles A. Williams
#                       Rensselaer Polytechnic Institute
#                        (C) 2005  All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

test: clean
	lithomop3dapp.py -typos=relaxed -lm3dscan.fileRoot=bm1b -log_summary -petsc_solver 1 -pc_type ilu -ksp_monitor -ksp_view
	# lithomop3dapp.py -typos=relaxed -lm3dscan.fileRoot=bm1b -log_summary -petsc_solver 1 -pc_type ilu -ksp_monitor -ksp_view -start_in_debugger

clean::
	$(RM_F) bm1b.ascii bm1b.plot bm1b*.inp


# version
# $Id: Make.mm,v 1.1 2005/04/01 22:55:27 willic3 Exp $

# Generated automatically by MakeMill on Tue Mar  2 17:05:23 2004

# End of file 
