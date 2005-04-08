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

test: clean
	lithomop3dapp.py -typos=relaxed -lm3dscan.fileRoot=bm1a -log_summary -petsc_solver 1 -pc_type ilu -ksp_monitor -ksp_view -ksp_rtol 1e-9

	# lithomop3dapp.py -typos=relaxed -lm3dscan.fileRoot=bm1a -log_summary -petsc_solver 1 -pc_type ilu -ksp_monitor -ksp_view -ksp_rtol 1e-9 -start_in_debugger

clean::
	$(RM_F) bm1a.ascii bm1a.plot bm1a*.inp


# version
# $Id: Make.mm,v 1.6 2005/04/08 00:47:49 willic3 Exp $

# Generated automatically by MakeMill on Tue Mar  2 17:05:23 2004

# End of file 
