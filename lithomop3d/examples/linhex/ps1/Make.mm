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
	lithomop3dapp.py --lm3dscan.fileRoot=elas-test1/ps1 --lm3dscan.coordinateInputFile=ps1.coord --lm3dscan.bcInputFile=ps1.bc --lm3dscan.stateVariableInputFile=ps1.statevar --lm3dscan.connectivityInputFile=ps1.connect --lm3dscan.fullOutputInputFile=ps1.fuldat
	lithomop3dapp.py --lm3dscan.fileRoot=elas-test2/ps1 --lm3dscan.coordinateInputFile=ps1.coord --lm3dscan.bcInputFile=ps1.bc --lm3dscan.stateVariableInputFile=ps1.statevar --lm3dscan.connectivityInputFile=ps1.connect --lm3dscan.fullOutputInputFile=ps1.fuldat
	lithomop3dapp.py --lm3dscan.fileRoot=elas-test3/ps1 --lm3dscan.coordinateInputFile=ps1.coord --lm3dscan.bcInputFile=ps1.bc --lm3dscan.stateVariableInputFile=ps1.statevar --lm3dscan.connectivityInputFile=ps1.connect --lm3dscan.fullOutputInputFile=ps1.fuldat
	lithomop3dapp.py --lm3dscan.fileRoot=elas-test4/ps1 --lm3dscan.coordinateInputFile=ps1.coord --lm3dscan.bcInputFile=ps1.bc --lm3dscan.stateVariableInputFile=ps1.statevar --lm3dscan.connectivityInputFile=ps1.connect --lm3dscan.fullOutputInputFile=ps1.fuldat
	lithomop3dapp.py --lm3dscan.fileRoot=elas-test5/ps1 --lm3dscan.coordinateInputFile=ps1.coord --lm3dscan.bcInputFile=ps1.bc --lm3dscan.stateVariableInputFile=ps1.statevar --lm3dscan.connectivityInputFile=ps1.connect --lm3dscan.fullOutputInputFile=ps1.fuldat
	lithomop3dapp.py --lm3dscan.fileRoot=elas-test6/ps1 --lm3dscan.coordinateInputFile=ps1.coord --lm3dscan.bcInputFile=ps1.bc --lm3dscan.stateVariableInputFile=ps1.statevar --lm3dscan.connectivityInputFile=ps1.connect --lm3dscan.fullOutputInputFile=ps1.fuldat
	lithomop3dapp.py --lm3dscan.fileRoot=elas-test7/ps1 --lm3dscan.coordinateInputFile=ps1.coord --lm3dscan.bcInputFile=ps1.bc --lm3dscan.stateVariableInputFile=ps1.statevar --lm3dscan.connectivityInputFile=ps1.connect --lm3dscan.fullOutputInputFile=ps1.fuldat
	lithomop3dapp.py --lm3dscan.fileRoot=elas-test8/ps1 --lm3dscan.coordinateInputFile=ps1.coord --lm3dscan.bcInputFile=ps1.bc --lm3dscan.stateVariableInputFile=ps1.statevar --lm3dscan.connectivityInputFile=ps1.connect --lm3dscan.fullOutputInputFile=ps1.fuldat
	lithomop3dapp.py --lm3dscan.fileRoot=elas-test9/ps1 --lm3dscan.coordinateInputFile=ps1.coord --lm3dscan.bcInputFile=ps1.bc --lm3dscan.stateVariableInputFile=ps1.statevar --lm3dscan.connectivityInputFile=ps1.connect --lm3dscan.fullOutputInputFile=ps1.fuldat
	lithomop3dapp.py --lm3dscan.fileRoot=elas-test10/ps1 --lm3dscan.coordinateInputFile=ps1.coord --lm3dscan.bcInputFile=ps1.bc --lm3dscan.stateVariableInputFile=ps1.statevar --lm3dscan.connectivityInputFile=ps1.connect --lm3dscan.fullOutputInputFile=ps1.fuldat
	lithomop3dapp.py --lm3dscan.fileRoot=elas-test11/ps1 --lm3dscan.coordinateInputFile=ps1.coord --lm3dscan.bcInputFile=ps1.bc --lm3dscan.stateVariableInputFile=ps1.statevar --lm3dscan.connectivityInputFile=ps1.connect --lm3dscan.fullOutputInputFile=ps1.fuldat
	lithomop3dapp.py --lm3dscan.fileRoot=elas-test12/ps1 --lm3dscan.coordinateInputFile=ps1.coord --lm3dscan.bcInputFile=ps1.bc --lm3dscan.stateVariableInputFile=ps1.statevar --lm3dscan.connectivityInputFile=ps1.connect --lm3dscan.fullOutputInputFile=ps1.fuldat
	lithomop3dapp.py --lm3dscan.fileRoot=vis-test1/ps1 --lm3dscan.coordinateInputFile=ps1.coord --lm3dscan.bcInputFile=ps1.bc --lm3dscan.stateVariableInputFile=ps1.statevar --lm3dscan.connectivityInputFile=ps1.connect --lm3dscan.fullOutputInputFile=ps1.fuldat
	lithomop3dapp.py --lm3dscan.fileRoot=vis-test2/ps1 --lm3dscan.coordinateInputFile=ps1.coord --lm3dscan.bcInputFile=ps1.bc --lm3dscan.stateVariableInputFile=ps1.statevar --lm3dscan.connectivityInputFile=ps1.connect --lm3dscan.fullOutputInputFile=ps1.fuldat
	lithomop3dapp.py --lm3dscan.fileRoot=vis-test3/ps1 --lm3dscan.coordinateInputFile=ps1.coord --lm3dscan.bcInputFile=ps1.bc --lm3dscan.stateVariableInputFile=ps1.statevar --lm3dscan.connectivityInputFile=ps1.connect --lm3dscan.fullOutputInputFile=ps1.fuldat
	lithomop3dapp.py --lm3dscan.fileRoot=vis-test4/ps1 --lm3dscan.coordinateInputFile=ps1.coord --lm3dscan.bcInputFile=ps1.bc --lm3dscan.stateVariableInputFile=ps1.statevar --lm3dscan.connectivityInputFile=ps1.connect --lm3dscan.fullOutputInputFile=ps1.fuldat
	lithomop3dapp.py --lm3dscan.fileRoot=vis-test5/ps1 --lm3dscan.coordinateInputFile=ps1.coord --lm3dscan.bcInputFile=ps1.bc --lm3dscan.stateVariableInputFile=ps1.statevar --lm3dscan.connectivityInputFile=ps1.connect --lm3dscan.fullOutputInputFile=ps1.fuldat
	lithomop3dapp.py --lm3dscan.fileRoot=vis-test6/ps1 --lm3dscan.coordinateInputFile=ps1.coord --lm3dscan.bcInputFile=ps1.bc --lm3dscan.stateVariableInputFile=ps1.statevar --lm3dscan.connectivityInputFile=ps1.connect --lm3dscan.fullOutputInputFile=ps1.fuldat
	lithomop3dapp.py --lm3dscan.fileRoot=vis-test7/ps1 --lm3dscan.coordinateInputFile=ps1.coord --lm3dscan.bcInputFile=ps1.bc --lm3dscan.stateVariableInputFile=ps1.statevar --lm3dscan.connectivityInputFile=ps1.connect --lm3dscan.fullOutputInputFile=ps1.fuldat
	lithomop3dapp.py --lm3dscan.fileRoot=vis-test8/ps1 --lm3dscan.coordinateInputFile=ps1.coord --lm3dscan.bcInputFile=ps1.bc --lm3dscan.stateVariableInputFile=ps1.statevar --lm3dscan.connectivityInputFile=ps1.connect --lm3dscan.fullOutputInputFile=ps1.fuldat
	lithomop3dapp.py --lm3dscan.fileRoot=vis-test9/ps1 --lm3dscan.coordinateInputFile=ps1.coord --lm3dscan.bcInputFile=ps1.bc --lm3dscan.stateVariableInputFile=ps1.statevar --lm3dscan.connectivityInputFile=ps1.connect --lm3dscan.fullOutputInputFile=ps1.fuldat
	lithomop3dapp.py --lm3dscan.fileRoot=vis-test10/ps1 --lm3dscan.coordinateInputFile=ps1.coord --lm3dscan.bcInputFile=ps1.bc --lm3dscan.stateVariableInputFile=ps1.statevar --lm3dscan.connectivityInputFile=ps1.connect --lm3dscan.fullOutputInputFile=ps1.fuldat
	lithomop3dapp.py --lm3dscan.fileRoot=vis-test11/ps1 --lm3dscan.coordinateInputFile=ps1.coord --lm3dscan.bcInputFile=ps1.bc --lm3dscan.stateVariableInputFile=ps1.statevar --lm3dscan.connectivityInputFile=ps1.connect --lm3dscan.fullOutputInputFile=ps1.fuldat
	lithomop3dapp.py --lm3dscan.fileRoot=vis-test12/ps1 --lm3dscan.coordinateInputFile=ps1.coord --lm3dscan.bcInputFile=ps1.bc --lm3dscan.stateVariableInputFile=ps1.statevar --lm3dscan.connectivityInputFile=ps1.connect --lm3dscan.fullOutputInputFile=ps1.fuldat

clean::
	$(RM_F) elas-test*/ps1.ascii elas-test*/ps1.plot elas-test*/ps1*.inp \
	vis-test*/ps1.ascii vis-test*/ps1.plot vis-test*/ps1*.inp


# version
# $Id: Make.mm,v 1.1 2005/02/24 00:52:48 willic3 Exp $

# Generated automatically by MakeMill on Tue Mar  2 17:05:23 2004

# End of file 
