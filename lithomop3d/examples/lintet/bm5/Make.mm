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
	lithomop3dapp.py -typos=relaxed -lm3dscan.fileRoot=bm5 -log_summary -petsc_solver 1 -pc_type ilu -ksp_monitor -ksp_view
	# lithomop3dapp.py -typos=relaxed -lm3dscan.fileRoot=bm5 -log_summary -petsc_solver 1 -pc_type ilu -ksp_monitor -ksp_view -start_in_debugger noxterm

clean::
	$(RM_F) bm5.ascii bm5.plot bm5*.inp


# version
# $Id: Make.mm,v 1.1 2005/04/05 23:14:54 willic3 Exp $

# Generated automatically by MakeMill on Tue Mar  2 17:05:23 2004

# End of file 
