#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  PyLith by Charles A. Williams, Brad Aagaard, and Matt Knepley
#
#  Copyright (c) 2004-2006 Rensselaer Polytechnic Institute
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
# The function of this code is to call the elastic and time-dependent solution
# drivers.  To do this, a number of previously-defined parameters need to be
# bundled into lists.  This portion of the code requires access to all of the
# information previously defined in Pylith3d_scan.py and Pylith3d_setup.py.
#


from pyre.components.Component import Component


class Pylith3d_run(Component):

    def initialize(self, scanner, setup):

        pl3dscan = scanner
        pl3dsetup = setup

	print ""
        print "Hello from pl3drun.initialize (begin)!"
        print "Importing information from other modules:"

        # The only parameters required from Pylith3d_scan are those in the
        # inventory.  All others have been imported into Pylith3d_setup, and
        # possibly modified there.  Get all required info from the inventory.
        
        self.analysisType = pl3dscan.inventory.analysisType

        # Import all necessary pointers, etc. from Pylith3d_setup.
	self.memorySize = pl3dsetup.memorySize
	self.intSize = pl3dsetup.intSize
	self.doubleSize = pl3dsetup.doubleSize
        
        self.numberTimeStepGroups = pl3dsetup.numberTimeStepGroups

        self.pointerToAlnz = pl3dsetup.pointerToAlnz
        self.pointerToPcg = pl3dsetup.pointerToPcg
        self.pointerToZcg = pl3dsetup.pointerToZcg
        self.pointerToJa = pl3dsetup.pointerToJa

        self.pointerToB = pl3dsetup.pointerToB
        self.pointerToBtot = pl3dsetup.pointerToBtot
        self.pointerToBres = pl3dsetup.pointerToBres
        self.pointerToPvec = pl3dsetup.pointerToPvec
        self.pointerToGvec1 = pl3dsetup.pointerToGvec1
        self.pointerToGvec2 = pl3dsetup.pointerToGvec2
        self.pointerToListArrayGrav = pl3dsetup.pointerToListArrayGrav

        self.pointerToX = pl3dsetup.pointerToX
        self.pointerToD = pl3dsetup.pointerToD
        self.pointerToDeld = pl3dsetup.pointerToDeld
        self.pointerToDprev = pl3dsetup.pointerToDprev
        self.pointerToDcur = pl3dsetup.pointerToDcur
        self.pointerToId = pl3dsetup.pointerToId
        self.pointerToIwink = pl3dsetup.pointerToIwink
        self.pointerToWink = pl3dsetup.pointerToWink
        self.pointerToListArrayNsysdat = pl3dsetup.pointerToListArrayNsysdat

        self.pointerToIbond = pl3dsetup.pointerToIbond
        self.pointerToBond = pl3dsetup.pointerToBond

        self.pointerToDx = pl3dsetup.pointerToDx
        self.pointerToDeldx = pl3dsetup.pointerToDeldx
        self.pointerToDxcur = pl3dsetup.pointerToDxcur
        self.pointerToDiforc = pl3dsetup.pointerToDiforc
        self.pointerToIdx = pl3dsetup.pointerToIdx
        self.pointerToIwinkx = pl3dsetup.pointerToIwinkx
        self.pointerToWinkx = pl3dsetup.pointerToWinkx
        self.pointerToIdslp = pl3dsetup.pointerToIdslp
        self.pointerToIpslp = pl3dsetup.pointerToIpslp
        self.pointerToIdhist = pl3dsetup.pointerToIdhist

        self.pointerToFault = pl3dsetup.pointerToFault
        self.pointerToNfault = pl3dsetup.pointerToNfault
        self.pointerToDfault = pl3dsetup.pointerToDfault
        self.pointerToTfault = pl3dsetup.pointerToTfault

        self.pointerToS = pl3dsetup.pointerToS
        self.pointerToStemp = pl3dsetup.pointerToStemp

        self.pointerToState = pl3dsetup.pointerToState
        self.pointerToDstate = pl3dsetup.pointerToDstate
        self.pointerToDmat = pl3dsetup.pointerToDmat
        self.pointerToIen = pl3dsetup.pointerToIen
        self.pointerToLm = pl3dsetup.pointerToLm
        self.pointerToLmx = pl3dsetup.pointerToLmx
        self.pointerToLmf = pl3dsetup.pointerToLmf
        self.pointerToInfiel = pl3dsetup.pointerToInfiel
        self.pointerToListArrayIddmat = pl3dsetup.pointerToListArrayIddmat
        self.pointerToListArrayNpar = pl3dsetup.pointerToListArrayNpar

        self.pointerToIelno = pl3dsetup.pointerToIelno
        self.pointerToIside = pl3dsetup.pointerToIside
        self.pointerToIhistry = pl3dsetup.pointerToIhistry
        self.pointerToPres = pl3dsetup.pointerToPres
        self.pointerToPdir = pl3dsetup.pointerToPdir

        self.pointerToListArrayPropertyList = pl3dsetup.pointerToListArrayPropertyList
        self.pointerToMhist = pl3dsetup.pointerToMhist
        self.pointerToMaterialInfo = pl3dsetup.pointerToMaterialInfo
        self.pointerToMaterialModelInfo = pl3dsetup.pointerToMaterialModelInfo
        self.pointerToMaterialModelStates = pl3dsetup.pointerToMaterialModelStates

        self.pointerToGauss = pl3dsetup.pointerToGauss
        self.pointerToSh = pl3dsetup.pointerToSh
        self.pointerToShj = pl3dsetup.pointerToShj
        self.pointerToElementTypeInfo = pl3dsetup.pointerToElementTypeInfo

        self.pointerToHistry = pl3dsetup.pointerToHistry
        self.pointerToListArrayRtimdat = pl3dsetup.pointerToListArrayRtimdat
        self.pointerToListArrayNvisdat = pl3dsetup.pointerToListArrayNvisdat
        self.pointerToMaxstp = pl3dsetup.pointerToMaxstp
        self.pointerToDelt = pl3dsetup.pointerToDelt
        self.pointerToAlfa = pl3dsetup.pointerToAlfa
        self.pointerToMaxit = pl3dsetup.pointerToMaxit
        self.pointerToNtdinit = pl3dsetup.pointerToNtdinit
        self.pointerToLgdef = pl3dsetup.pointerToLgdef
        self.pointerToUtol = pl3dsetup.pointerToUtol
        self.pointerToFtol = pl3dsetup.pointerToFtol
        self.pointerToEtol = pl3dsetup.pointerToEtol
        self.pointerToItmax = pl3dsetup.pointerToItmax

        self.pointerToListArrayRgiter = pl3dsetup.pointerToListArrayRgiter
        self.pointerToListArrayRmin = pl3dsetup.pointerToListArrayRmin
        self.pointerToListArrayRmult = pl3dsetup.pointerToListArrayRmult
        self.pointerToListArrayNsiter = pl3dsetup.pointerToListArrayNsiter

        self.pointerToSkew = pl3dsetup.pointerToSkew

        self.pointerToIprint = pl3dsetup.pointerToIprint
        self.pointerToListArrayNcodat = pl3dsetup.pointerToListArrayNcodat
        self.pointerToListArrayNunits = pl3dsetup.pointerToListArrayNunits
        self.pointerToListArrayNprint = pl3dsetup.pointerToListArrayNprint
        self.pointerToIstatout = pl3dsetup.pointerToIstatout

        self.asciiOutputFile = pl3dsetup.asciiOutputFile
        self.plotOutputFile = pl3dsetup.plotOutputFile

	print ""
        print "Hello from pl3drun.initialize (end)!"

        return

    def run(self):
        import pylith3d
        
        # First define all of the lists that maintain variable values.  The
        # variables in these lists are altered during the running of the code
        # and should not be accessed directly except as a member of the list.
        # They should not have been defined previously.

	print ""
        print "Hello from pl3drun.run (begin)!"
        print "Beginning problem solution:"

        # ntimdat array
        self.currentTimeStep = 0
        self.currentIterationsBetweenReform = 0
        self.currentStepsBetweenReform = 0
        self.currentLargeDeformationFlag = 0
        self.currentBbarFlag = 0
        self.currentMaximumIterations = 0
        self.currentNumberTotalIterations = 0
        self.currentNumberReforms = 0
        self.currentNumberTotalPcgIterations = 0
        self.reformFlagInt = 0
        self.listNtimdat = [
            self.currentTimeStep,
            self.currentIterationsBetweenReform,
            self.currentStepsBetweenReform,
            self.currentLargeDeformationFlag,
            self.currentMaximumIterations,
            self.currentNumberTotalIterations,
            self.currentNumberReforms,
            self.currentNumberTotalPcgIterations,
            self.reformFlagInt]
        self.pointerToListArrayNtimdat = pylith3d.intListToArray(
            self.listNtimdat)
	self.memorySize += 9*self.intSize

        # gcurr array
        self.currentDisplacementNorm = 0.0
        self.currentForceNorm = 0.0
        self.currentEnergyNorm = 0.0
        self.listGcurr = [
            self.currentDisplacementNorm,
            self.currentForceNorm,
            self.currentEnergyNorm]
        self.pointerToListArrayGcurr = pylith3d.doubleListToArray(
            self.listGcurr)
	self.memorySize += 3*self.doubleSize

        # gi array
        self.initialDisplacementNorm = 0.0
        self.initialForceNorm = 0.0
        self.initialEnergyNorm = 0.0
        self.listGi = [
            self.initialDisplacementNorm,
            self.initialForceNorm,
            self.initialEnergyNorm]
        self.pointerToListArrayGi = pylith3d.doubleListToArray(
            self.listGi)
	self.memorySize += 3*self.doubleSize
        
        # gprev array
        self.previousDisplacementNorm = 0.0
        self.previousForceNorm = 0.0
        self.previousEnergyNorm = 0.0
        self.listGprev = [
            self.previousDisplacementNorm,
            self.previousForceNorm,
            self.previousEnergyNorm]
        self.pointerToListArrayGprev = pylith3d.doubleListToArray(
            self.listGprev)
	self.memorySize += 3*self.doubleSize
        
        # gtol array
        self.currentDisplacementTolerance = 0.0
        self.currentForceTolerance = 0.0
        self.currentEnergyTolerance = 0.0
        self.listGtol = [
            self.currentDisplacementTolerance,
            self.currentForceTolerance,
            self.currentEnergyTolerance]
        self.pointerToListArrayGtol = pylith3d.doubleListToArray(
            self.listGtol)
	self.memorySize += 3*self.doubleSize

        # Output approximate memory usage
        self.memorySizeMB =0.0
        self.memorySizeMB=self.memorySize/(1024.0*1024.0)

	print ""
	print "Approximate memory allocation for f77 arrays (MB): %g" % self.memorySizeMB
	# print "Just before pylith3d.elastc:"

        # Perform elastic solution, if requested.

        if self.analysisType == "elasticSolution" or self.analysisType == "fullSolution":
            pylith3d.elastc(
                self.pointerToAlnz,
                self.pointerToPcg,
                self.pointerToZcg,
                self.pointerToJa,
                self.pointerToB,
                self.pointerToBtot,
                self.pointerToBres,
                self.pointerToPvec,
                self.pointerToGvec1,
                self.pointerToGvec2,
                self.pointerToListArrayGrav,
                self.pointerToX,
                self.pointerToD,
                self.pointerToDeld,
                self.pointerToDprev,
                self.pointerToDcur,
                self.pointerToId,
                self.pointerToIwink,
                self.pointerToWink,
                self.pointerToListArrayNsysdat,
                self.pointerToIbond,
                self.pointerToBond,
                self.pointerToDx,
                self.pointerToDeldx,
                self.pointerToDxcur,
                self.pointerToDiforc,
                self.pointerToIdx,
                self.pointerToIwinkx,
                self.pointerToWinkx,
                self.pointerToIdslp,
                self.pointerToIpslp,
                self.pointerToIdhist,
                self.pointerToFault,
                self.pointerToNfault,
                self.pointerToDfault,
                self.pointerToTfault,
                self.pointerToS,
                self.pointerToStemp,
                self.pointerToState,
                self.pointerToDstate,
                self.pointerToDmat,
                self.pointerToIen,
                self.pointerToLm,
                self.pointerToLmx,
                self.pointerToLmf,
                self.pointerToInfiel,
                self.pointerToListArrayIddmat,
                self.pointerToListArrayNpar,
                self.pointerToIelno,
                self.pointerToIside,
                self.pointerToIhistry,
                self.pointerToPres,
                self.pointerToPdir,
                self.pointerToListArrayPropertyList,
                self.pointerToMhist,
                self.pointerToMaterialInfo,
                self.pointerToMaterialModelInfo,
                self.pointerToMaterialModelStates,
                self.pointerToGauss,
                self.pointerToSh,
                self.pointerToShj,
                self.pointerToElementTypeInfo,
                self.pointerToHistry,
                self.pointerToListArrayRtimdat,
                self.pointerToListArrayNtimdat,
                self.pointerToListArrayNvisdat,
                self.pointerToMaxstp,
                self.pointerToDelt,
                self.pointerToAlfa,
                self.pointerToMaxit,
                self.pointerToNtdinit,
                self.pointerToLgdef,
                self.pointerToUtol,
                self.pointerToFtol,
                self.pointerToEtol,
                self.pointerToItmax,
                self.pointerToListArrayRgiter,
                self.pointerToListArrayGcurr,
                self.pointerToListArrayGi,
                self.pointerToListArrayGprev,
                self.pointerToListArrayGtol,
                self.pointerToListArrayRmin,
                self.pointerToListArrayRmult,
                self.pointerToListArrayNsiter,
                self.pointerToSkew,
                self.pointerToListArrayNcodat,
                self.pointerToListArrayNunits,
                self.pointerToListArrayNprint,
                self.pointerToIstatout,
                self.asciiOutputFile,
                self.plotOutputFile)

        # Perform time-dependent solution, if requested.

        if self.analysisType == "fullSolution" and self.numberTimeStepGroups > 1:
            pylith3d.viscos(
                self.pointerToAlnz,
                self.pointerToPcg,
                self.pointerToZcg,
                self.pointerToJa,
                self.pointerToB,
                self.pointerToBtot,
                self.pointerToBres,
                self.pointerToPvec,
                self.pointerToGvec1,
                self.pointerToGvec2,
                self.pointerToListArrayGrav,
                self.pointerToX,
                self.pointerToD,
                self.pointerToDeld,
                self.pointerToDprev,
                self.pointerToDcur,
                self.pointerToId,
                self.pointerToIwink,
                self.pointerToWink,
                self.pointerToListArrayNsysdat,
                self.pointerToIbond,
                self.pointerToBond,
                self.pointerToDx,
                self.pointerToDeldx,
                self.pointerToDxcur,
                self.pointerToDiforc,
                self.pointerToIdx,
                self.pointerToIwinkx,
                self.pointerToWinkx,
                self.pointerToIdslp,
                self.pointerToIpslp,
                self.pointerToIdhist,
                self.pointerToFault,
                self.pointerToNfault,
                self.pointerToDfault,
                self.pointerToTfault,
                self.pointerToS,
                self.pointerToStemp,
                self.pointerToState,
                self.pointerToDstate,
                self.pointerToDmat,
                self.pointerToIen,
                self.pointerToLm,
                self.pointerToLmx,
                self.pointerToLmf,
                self.pointerToInfiel,
                self.pointerToListArrayIddmat,
                self.pointerToListArrayNpar,
                self.pointerToIelno,
                self.pointerToIside,
                self.pointerToIhistry,
                self.pointerToPres,
                self.pointerToPdir,
                self.pointerToListArrayPropertyList,
                self.pointerToMhist,
                self.pointerToMaterialInfo,
                self.pointerToMaterialModelInfo,
                self.pointerToMaterialModelStates,
                self.pointerToGauss,
                self.pointerToSh,
                self.pointerToShj,
                self.pointerToElementTypeInfo,
                self.pointerToHistry,
                self.pointerToListArrayRtimdat,
                self.pointerToListArrayNtimdat,
                self.pointerToListArrayNvisdat,
                self.pointerToMaxstp,
                self.pointerToDelt,
                self.pointerToAlfa,
                self.pointerToMaxit,
                self.pointerToNtdinit,
                self.pointerToLgdef,
                self.pointerToUtol,
                self.pointerToFtol,
                self.pointerToEtol,
                self.pointerToItmax,
                self.pointerToListArrayRgiter,
                self.pointerToListArrayGcurr,
                self.pointerToListArrayGi,
                self.pointerToListArrayGprev,
                self.pointerToListArrayGtol,
                self.pointerToListArrayRmin,
                self.pointerToListArrayRmult,
                self.pointerToListArrayNsiter,
                self.pointerToSkew,
                self.pointerToIprint,
                self.pointerToListArrayNcodat,
                self.pointerToListArrayNunits,
                self.pointerToListArrayNprint,
                self.pointerToIstatout,
                self.asciiOutputFile,
                self.plotOutputFile)
                          
	print ""
        print "Hello from pl3drun.run (end)!"
        return


    def __init__(self):
        Component.__init__(self, "pl3drun", "solver")

	print ""
        print "Hello from pl3drun.__init__!"

        return


# version
# $Id: Pylith3d_run.py,v 1.1 2004/09/23 17:37:55 willic3 Exp $

# End of file 
