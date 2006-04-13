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
	pylith3dapp.py -pl3dscan.fileRoot=../bm1a/bm1a -pl3dscan.bcInputFile=bm2a.bc -pl3dscan.keywordEqualsValueFile=bm2a.keyval -pl3dscan.timeStepInputFile=bm2a.time -pl3dscan.asciiOutputFile=bm2a.ascii -pl3dscan.ucdOutputRoot=bm2a

clean::
	$(RM_F) bm2a.ascii bm2a.plot bm2a*.inp


# version
# $Id: Make.mm,v 1.1 2005/01/06 16:33:56 willic3 Exp $

# Generated automatically by MakeMill on Tue Mar  2 17:05:23 2004

# End of file 
