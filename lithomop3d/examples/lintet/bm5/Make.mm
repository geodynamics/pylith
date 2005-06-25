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
	lithomop3dapp.py -typos=relaxed -lm3dscan.fileRoot=bm5 -lm3dscan.ucdOutput=binary -log_summary -pc_type ilu -ksp_monitor -ksp_view -ksp_rtol 1e-9

# test: clean
	# lithomop3dapp.py -typos=relaxed -lm3dscan.fileRoot=bm5 -lm3dscan.ucdOutput=ascii -log_summary -pc_type ilu -ksp_monitor -ksp_view -ksp_rtol 1e-9 -start_in_debugger

clean::
	$(RM_F) bm5.ascii bm5.plot bm5*.inp


# version
# $Id: Make.mm,v 1.3 2005/06/24 20:26:18 willic3 Exp $

# Generated automatically by MakeMill on Tue Mar  2 17:05:23 2004

# End of file 
