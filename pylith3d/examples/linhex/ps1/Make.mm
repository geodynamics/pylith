# -*- Makefile -*-
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  PyLith by Charles A. Williams
#  Copyright (c) 2003-2006 Rensselaer Polytechnic Institute
#
#  Permission is hereby granted, free of charge, to any person obtaining
#  a copy of this software and associated documentation files (the
#  "Software"), to deal in the Software without restriction, including
#  without limitation the rights to use, copy, modify, merge, publish,
#  distribute, sublicense, and/or sell copies of the Software, and to
#  permit persons to whom the Software is furnished to do so, subject to
#  the following conditions:
#
#  The above copyright notice and this permission notice shall be
#  included in all copies or substantial portions of the Software.
#
#  THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
#  EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
#  MERCHANTABILITY,    FITNESS    FOR    A   PARTICULAR    PURPOSE    AND
#  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
#  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
#  OF CONTRACT, TORT OR OTHERWISE,  ARISING FROM, OUT OF OR IN CONNECTION
#  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

test: clean
	pylith3dapp.py --pl3dscan.fileRoot=elas-test1/ps1 --pl3dscan.coordinateInputFile=ps1.coord --pl3dscan.bcInputFile=ps1.bc --pl3dscan.stateVariableInputFile=ps1.statevar --pl3dscan.connectivityInputFile=ps1.connect --pl3dscan.fullOutputInputFile=ps1.fuldat
	pylith3dapp.py --pl3dscan.fileRoot=elas-test2/ps1 --pl3dscan.coordinateInputFile=ps1.coord --pl3dscan.bcInputFile=ps1.bc --pl3dscan.stateVariableInputFile=ps1.statevar --pl3dscan.connectivityInputFile=ps1.connect --pl3dscan.fullOutputInputFile=ps1.fuldat
	pylith3dapp.py --pl3dscan.fileRoot=elas-test3/ps1 --pl3dscan.coordinateInputFile=ps1.coord --pl3dscan.bcInputFile=ps1.bc --pl3dscan.stateVariableInputFile=ps1.statevar --pl3dscan.connectivityInputFile=ps1.connect --pl3dscan.fullOutputInputFile=ps1.fuldat
	pylith3dapp.py --pl3dscan.fileRoot=elas-test4/ps1 --pl3dscan.coordinateInputFile=ps1.coord --pl3dscan.bcInputFile=ps1.bc --pl3dscan.stateVariableInputFile=ps1.statevar --pl3dscan.connectivityInputFile=ps1.connect --pl3dscan.fullOutputInputFile=ps1.fuldat
	pylith3dapp.py --pl3dscan.fileRoot=elas-test5/ps1 --pl3dscan.coordinateInputFile=ps1.coord --pl3dscan.bcInputFile=ps1.bc --pl3dscan.stateVariableInputFile=ps1.statevar --pl3dscan.connectivityInputFile=ps1.connect --pl3dscan.fullOutputInputFile=ps1.fuldat
	pylith3dapp.py --pl3dscan.fileRoot=elas-test6/ps1 --pl3dscan.coordinateInputFile=ps1.coord --pl3dscan.bcInputFile=ps1.bc --pl3dscan.stateVariableInputFile=ps1.statevar --pl3dscan.connectivityInputFile=ps1.connect --pl3dscan.fullOutputInputFile=ps1.fuldat
	pylith3dapp.py --pl3dscan.fileRoot=elas-test7/ps1 --pl3dscan.coordinateInputFile=ps1.coord --pl3dscan.bcInputFile=ps1.bc --pl3dscan.stateVariableInputFile=ps1.statevar --pl3dscan.connectivityInputFile=ps1.connect --pl3dscan.fullOutputInputFile=ps1.fuldat
	pylith3dapp.py --pl3dscan.fileRoot=elas-test8/ps1 --pl3dscan.coordinateInputFile=ps1.coord --pl3dscan.bcInputFile=ps1.bc --pl3dscan.stateVariableInputFile=ps1.statevar --pl3dscan.connectivityInputFile=ps1.connect --pl3dscan.fullOutputInputFile=ps1.fuldat
	pylith3dapp.py --pl3dscan.fileRoot=elas-test9/ps1 --pl3dscan.coordinateInputFile=ps1.coord --pl3dscan.bcInputFile=ps1.bc --pl3dscan.stateVariableInputFile=ps1.statevar --pl3dscan.connectivityInputFile=ps1.connect --pl3dscan.fullOutputInputFile=ps1.fuldat
	pylith3dapp.py --pl3dscan.fileRoot=elas-test10/ps1 --pl3dscan.coordinateInputFile=ps1.coord --pl3dscan.bcInputFile=ps1.bc --pl3dscan.stateVariableInputFile=ps1.statevar --pl3dscan.connectivityInputFile=ps1.connect --pl3dscan.fullOutputInputFile=ps1.fuldat
	pylith3dapp.py --pl3dscan.fileRoot=elas-test11/ps1 --pl3dscan.coordinateInputFile=ps1.coord --pl3dscan.bcInputFile=ps1.bc --pl3dscan.stateVariableInputFile=ps1.statevar --pl3dscan.connectivityInputFile=ps1.connect --pl3dscan.fullOutputInputFile=ps1.fuldat
	pylith3dapp.py --pl3dscan.fileRoot=elas-test12/ps1 --pl3dscan.coordinateInputFile=ps1.coord --pl3dscan.bcInputFile=ps1.bc --pl3dscan.stateVariableInputFile=ps1.statevar --pl3dscan.connectivityInputFile=ps1.connect --pl3dscan.fullOutputInputFile=ps1.fuldat
	pylith3dapp.py --pl3dscan.fileRoot=vis-test1/ps1 --pl3dscan.coordinateInputFile=ps1.coord --pl3dscan.bcInputFile=ps1.bc --pl3dscan.stateVariableInputFile=ps1.statevar --pl3dscan.connectivityInputFile=ps1.connect --pl3dscan.fullOutputInputFile=ps1.fuldat
	pylith3dapp.py --pl3dscan.fileRoot=vis-test2/ps1 --pl3dscan.coordinateInputFile=ps1.coord --pl3dscan.bcInputFile=ps1.bc --pl3dscan.stateVariableInputFile=ps1.statevar --pl3dscan.connectivityInputFile=ps1.connect --pl3dscan.fullOutputInputFile=ps1.fuldat
	pylith3dapp.py --pl3dscan.fileRoot=vis-test3/ps1 --pl3dscan.coordinateInputFile=ps1.coord --pl3dscan.bcInputFile=ps1.bc --pl3dscan.stateVariableInputFile=ps1.statevar --pl3dscan.connectivityInputFile=ps1.connect --pl3dscan.fullOutputInputFile=ps1.fuldat
	pylith3dapp.py --pl3dscan.fileRoot=vis-test4/ps1 --pl3dscan.coordinateInputFile=ps1.coord --pl3dscan.bcInputFile=ps1.bc --pl3dscan.stateVariableInputFile=ps1.statevar --pl3dscan.connectivityInputFile=ps1.connect --pl3dscan.fullOutputInputFile=ps1.fuldat
	pylith3dapp.py --pl3dscan.fileRoot=vis-test5/ps1 --pl3dscan.coordinateInputFile=ps1.coord --pl3dscan.bcInputFile=ps1.bc --pl3dscan.stateVariableInputFile=ps1.statevar --pl3dscan.connectivityInputFile=ps1.connect --pl3dscan.fullOutputInputFile=ps1.fuldat
	pylith3dapp.py --pl3dscan.fileRoot=vis-test6/ps1 --pl3dscan.coordinateInputFile=ps1.coord --pl3dscan.bcInputFile=ps1.bc --pl3dscan.stateVariableInputFile=ps1.statevar --pl3dscan.connectivityInputFile=ps1.connect --pl3dscan.fullOutputInputFile=ps1.fuldat
	pylith3dapp.py --pl3dscan.fileRoot=vis-test7/ps1 --pl3dscan.coordinateInputFile=ps1.coord --pl3dscan.bcInputFile=ps1.bc --pl3dscan.stateVariableInputFile=ps1.statevar --pl3dscan.connectivityInputFile=ps1.connect --pl3dscan.fullOutputInputFile=ps1.fuldat
	pylith3dapp.py --pl3dscan.fileRoot=vis-test8/ps1 --pl3dscan.coordinateInputFile=ps1.coord --pl3dscan.bcInputFile=ps1.bc --pl3dscan.stateVariableInputFile=ps1.statevar --pl3dscan.connectivityInputFile=ps1.connect --pl3dscan.fullOutputInputFile=ps1.fuldat
	pylith3dapp.py --pl3dscan.fileRoot=vis-test9/ps1 --pl3dscan.coordinateInputFile=ps1.coord --pl3dscan.bcInputFile=ps1.bc --pl3dscan.stateVariableInputFile=ps1.statevar --pl3dscan.connectivityInputFile=ps1.connect --pl3dscan.fullOutputInputFile=ps1.fuldat
	pylith3dapp.py --pl3dscan.fileRoot=vis-test10/ps1 --pl3dscan.coordinateInputFile=ps1.coord --pl3dscan.bcInputFile=ps1.bc --pl3dscan.stateVariableInputFile=ps1.statevar --pl3dscan.connectivityInputFile=ps1.connect --pl3dscan.fullOutputInputFile=ps1.fuldat
	pylith3dapp.py --pl3dscan.fileRoot=vis-test11/ps1 --pl3dscan.coordinateInputFile=ps1.coord --pl3dscan.bcInputFile=ps1.bc --pl3dscan.stateVariableInputFile=ps1.statevar --pl3dscan.connectivityInputFile=ps1.connect --pl3dscan.fullOutputInputFile=ps1.fuldat
	pylith3dapp.py --pl3dscan.fileRoot=vis-test12/ps1 --pl3dscan.coordinateInputFile=ps1.coord --pl3dscan.bcInputFile=ps1.bc --pl3dscan.stateVariableInputFile=ps1.statevar --pl3dscan.connectivityInputFile=ps1.connect --pl3dscan.fullOutputInputFile=ps1.fuldat

clean::
	$(RM_F) elas-test*/ps1.ascii elas-test*/ps1.plot elas-test*/ps1*.inp \
	vis-test*/ps1.ascii vis-test*/ps1.plot vis-test*/ps1*.inp


# version
# $Id: Make.mm,v 1.1 2005/02/24 00:52:48 willic3 Exp $

# Generated automatically by MakeMill on Tue Mar  2 17:05:23 2004

# End of file 
