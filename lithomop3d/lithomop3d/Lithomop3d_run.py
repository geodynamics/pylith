#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                             Charles A. Williams
#                       Rensselaer Polytechnic Institute
#                        (C) 2004  All Rights Reserved
#
#  Copyright 2004 Rensselaer Polytechnic Institute.
#  All worldwide rights reserved.  A license to use, copy, modify and
#  distribute this software for non-commercial research purposes only
#  is hereby granted, provided that this copyright notice and
#  accompanying disclaimer is not modified or removed from the software.
#
#  DISCLAIMER:  The software is distributed "AS IS" without any express
#  or implied warranty, including but not limited to, any implied
#  warranties of merchantability or fitness for a particular purpose
#  or any warranty of non-infringement of any current or pending patent
#  rights.  The authors of the software make no representations about
#  the suitability of this software for any particular purpose.  The
#  entire risk as to the quality and performance of the software is with
#  the user.  Should the software prove defective, the user assumes the
#  cost of all necessary servicing, repair or correction.  In
#  particular, neither Rensselaer Polytechnic Institute, nor the authors
#  of the software are liable for any indirect, special, consequential,
#  or incidental damages related to the software, to the maximum extent
#  the law permits.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# The function of this code is to call the elastic and time-dependent solution
# drivers.  To do this, a number of previously-defined parameters need to be
# bundled into lists.  This portion of the code requires access to all of the
# information previously defined in Lithomop3d_scan.py and Lithomop3d_setup.py.
#


from pyre.components.Component import Component


class Lithomop3d_run(Component):

    def initialize(self, scanner, setup):

        lm3dscan = scanner
        lm3dsetup = setup

	print ""
        print "Hello from lm3drun.initialize (begin)!"
        print "Importing information from other modules:"

        # The only parameters required from Lithomop3d_scan are those in the
        # inventory.  All others have been imported into Lithomop3d_setup, and
        # possibly modified there.  Get all required info from the inventory.

        # PETSc logging
        self.elasticStage = lm3dsetup.elasticStage
        self.viscousStage = lm3dsetup.viscousStage
        self.iterateEvent = lm3dsetup.iterateEvent
        self.A = lm3dsetup.A
        
        self.analysisType = lm3dscan.inventory.analysisType

        # Import all necessary pointers, etc. from Lithomop3d_setup.
	self.memorySize = lm3dsetup.memorySize
	self.intSize = lm3dsetup.intSize
	self.doubleSize = lm3dsetup.doubleSize
        
        self.numberTimeStepGroups = lm3dsetup.numberTimeStepGroups

        self.pointerToAlnz = lm3dsetup.pointerToAlnz
        self.pointerToPcg = lm3dsetup.pointerToPcg
        self.pointerToZcg = lm3dsetup.pointerToZcg
        self.pointerToDprev = lm3dsetup.pointerToDprev
        self.pointerToJa = lm3dsetup.pointerToJa

        # self.pointerToB = lm3dsetup.pointerToB
        # self.pointerToBtot = lm3dsetup.pointerToBtot
        # self.pointerToBres = lm3dsetup.pointerToBres
        # self.pointerToPvec = lm3dsetup.pointerToPvec
        # self.pointerToGvec1 = lm3dsetup.pointerToGvec1
        # self.pointerToGvec2 = lm3dsetup.pointerToGvec2
        self.pointerToBextern = lm3dsetup.pointerToBextern
        self.pointerToBtraction = lm3dsetup.pointerToBtraction
        self.pointerToBgravity = lm3dsetup.pointerToBgravity
        self.pointerToBconcForce = lm3dsetup.pointerToBconcForce
        # self.pointerToBprestress = lm3dsetup.pointerToBprestress
        self.pointerToBintern = lm3dsetup.pointerToBintern
        self.pointerToBresid = lm3dsetup.pointerToBresid
        self.pointerToBwork = lm3dsetup.pointerToBwork
        self.pointerToDispVec = lm3dsetup.pointerToDispVec
        self.externFlag = lm3dsetup.externFlag
        self.tractionFlag = lm3dsetup.tractionFlag
        self.gravityFlag = lm3dsetup.gravityFlag
        self.concForceFlag = lm3dsetup.concForceFlag
        self.prestressFlag = lm3dsetup.prestressFlag
        self.usePreviousDisplacementFlag = lm3dsetup.usePreviousDisplacementFlag
        self.pointerToListArrayGrav = lm3dsetup.pointerToListArrayGrav

        self.pointerToX = lm3dsetup.pointerToX
        self.pointerToD = lm3dsetup.pointerToD
        self.pointerToDeld = lm3dsetup.pointerToDeld
        self.pointerToDcur = lm3dsetup.pointerToDcur
        self.pointerToId = lm3dsetup.pointerToId
        self.pointerToIwink = lm3dsetup.pointerToIwink
        self.pointerToWink = lm3dsetup.pointerToWink
        self.pointerToListArrayNsysdat = lm3dsetup.pointerToListArrayNsysdat

        self.pointerToIbond = lm3dsetup.pointerToIbond
        self.pointerToBond = lm3dsetup.pointerToBond

        self.pointerToDx = lm3dsetup.pointerToDx
        self.pointerToDeldx = lm3dsetup.pointerToDeldx
        self.pointerToDxcur = lm3dsetup.pointerToDxcur
        self.pointerToDiforc = lm3dsetup.pointerToDiforc
        self.pointerToIdx = lm3dsetup.pointerToIdx
        self.pointerToIwinkx = lm3dsetup.pointerToIwinkx
        self.pointerToWinkx = lm3dsetup.pointerToWinkx
        self.pointerToIdslp = lm3dsetup.pointerToIdslp
        self.pointerToIpslp = lm3dsetup.pointerToIpslp
        self.pointerToIdhist = lm3dsetup.pointerToIdhist

        self.pointerToFault = lm3dsetup.pointerToFault
        self.pointerToNfault = lm3dsetup.pointerToNfault
        self.pointerToDfault = lm3dsetup.pointerToDfault
        self.pointerToTfault = lm3dsetup.pointerToTfault

        self.pointerToS = lm3dsetup.pointerToS
        self.pointerToStemp = lm3dsetup.pointerToStemp

        self.pointerToState = lm3dsetup.pointerToState
        self.pointerToDstate = lm3dsetup.pointerToDstate
        self.pointerToState0 = lm3dsetup.pointerToState0
        self.pointerToDmat = lm3dsetup.pointerToDmat
        self.pointerToIen = lm3dsetup.pointerToIen
        self.pointerToLm = lm3dsetup.pointerToLm
        self.pointerToLmx = lm3dsetup.pointerToLmx
        self.pointerToLmf = lm3dsetup.pointerToLmf
        self.pointerToInfiel = lm3dsetup.pointerToInfiel
        self.pointerToListArrayIddmat = lm3dsetup.pointerToListArrayIddmat
        self.pointerToListArrayNpar = lm3dsetup.pointerToListArrayNpar

        self.prestressAutoComputeInt = lm3dsetup.prestressAutoComputeInt

        self.pointerToIelno = lm3dsetup.pointerToIelno
        self.pointerToIside = lm3dsetup.pointerToIside
        self.pointerToIhistry = lm3dsetup.pointerToIhistry
        self.pointerToPres = lm3dsetup.pointerToPres
        self.pointerToPdir = lm3dsetup.pointerToPdir

        self.pointerToListArrayPropertyList = lm3dsetup.pointerToListArrayPropertyList
        self.pointerToMhist = lm3dsetup.pointerToMhist
        self.pointerToMaterialInfo = lm3dsetup.pointerToMaterialInfo
        self.pointerToMaterialModelInfo = lm3dsetup.pointerToMaterialModelInfo
        self.pointerToMaterialModelStates = lm3dsetup.pointerToMaterialModelStates

        self.pointerToGauss = lm3dsetup.pointerToGauss
        self.pointerToSh = lm3dsetup.pointerToSh
        self.pointerToShj = lm3dsetup.pointerToShj
        self.pointerToElementTypeInfo = lm3dsetup.pointerToElementTypeInfo

        self.pointerToHistry = lm3dsetup.pointerToHistry
        self.pointerToListArrayRtimdat = lm3dsetup.pointerToListArrayRtimdat
        self.pointerToListArrayNvisdat = lm3dsetup.pointerToListArrayNvisdat
        self.pointerToMaxstp = lm3dsetup.pointerToMaxstp
        self.pointerToDelt = lm3dsetup.pointerToDelt
        self.pointerToAlfa = lm3dsetup.pointerToAlfa
        self.pointerToMaxit = lm3dsetup.pointerToMaxit
        self.pointerToNtdinit = lm3dsetup.pointerToNtdinit
        self.pointerToLgdef = lm3dsetup.pointerToLgdef
        self.pointerToUtol = lm3dsetup.pointerToUtol
        self.pointerToFtol = lm3dsetup.pointerToFtol
        self.pointerToEtol = lm3dsetup.pointerToEtol
        self.pointerToItmax = lm3dsetup.pointerToItmax

        self.pointerToListArrayRgiter = lm3dsetup.pointerToListArrayRgiter
        self.pointerToListArrayRmin = lm3dsetup.pointerToListArrayRmin
        self.pointerToListArrayRmult = lm3dsetup.pointerToListArrayRmult
        self.pointerToListArrayNsiter = lm3dsetup.pointerToListArrayNsiter

        self.pointerToSkew = lm3dsetup.pointerToSkew

        self.pointerToIprint = lm3dsetup.pointerToIprint
        self.pointerToListArrayNcodat = lm3dsetup.pointerToListArrayNcodat
        self.pointerToListArrayNunits = lm3dsetup.pointerToListArrayNunits
        self.pointerToListArrayNprint = lm3dsetup.pointerToListArrayNprint
        self.pointerToIstatout = lm3dsetup.pointerToIstatout

        self.asciiOutputFile = lm3dsetup.asciiOutputFile
        self.plotOutputFile = lm3dsetup.plotOutputFile
        self.ucdOutputRoot = lm3dsetup.ucdOutputRoot

	print ""
        print "Hello from lm3drun.initialize (end)!"

        return

    def run(self):
        import lithomop3d
        
        # First define all of the lists that maintain variable values.  The
        # variables in these lists are altered during the running of the code
        # and should not be accessed directly except as a member of the list.
        # They should not have been defined previously.

	print ""
        print "Hello from lm3drun.run (begin)!"
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
        self.pointerToListArrayNtimdat = lithomop3d.intListToArray(
            self.listNtimdat)
	# print "After listNtimdat"
	self.memorySize += 9*self.intSize

        # gcurr array
        self.currentDisplacementNorm = 0.0
        self.currentForceNorm = 0.0
        self.currentEnergyNorm = 0.0
        self.listGcurr = [
            self.currentDisplacementNorm,
            self.currentForceNorm,
            self.currentEnergyNorm]
        self.pointerToListArrayGcurr = lithomop3d.doubleListToArray(
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
        self.pointerToListArrayGi = lithomop3d.doubleListToArray(
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
        self.pointerToListArrayGprev = lithomop3d.doubleListToArray(
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
        self.pointerToListArrayGtol = lithomop3d.doubleListToArray(
            self.listGtol)
	self.memorySize += 3*self.doubleSize

        # Nforce array
        self.listNforce = [
            self.externFlag,
            self.tractionFlag,
            self.gravityFlag,
            self.concForceFlag,
            self.prestressFlag,
            self.usePreviousDisplacementFlag]
        self.pointerToListArrayNforce = lithomop3d.intListToArray(
            self.listNforce)
	self.memorySize += 6*self.intSize

        # Output approximate memory usage
        self.memorySizeMB =0.0
        self.memorySizeMB=self.memorySize/(1024.0*1024.0)

	print ""
	print "Approximate memory allocation for f77 arrays (MB): %g" % self.memorySizeMB
	# print "Just before lithomop3d.autoprestr:"

        # Compute gravitational prestresses, if requested.
        if self.analysisType == "elasticSolution" or self.analysisType == "fullSolution":
            if self.prestressAutoComputeInt == 1:
                lithomop3d.autoprestr(
                    self.A,
                    self.pointerToAlnz,
                    self.pointerToPcg,
                    self.pointerToZcg,
                    self.pointerToDprev,
                    self.pointerToJa,
                    self.pointerToBextern,
                    self.pointerToBtraction,
                    self.pointerToBgravity,
                    self.pointerToBconcForce,
                    # self.pointerToBprestress,
                    self.pointerToBintern,
                    self.pointerToBresid,
                    self.pointerToBwork,
                    self.pointerToDispVec,
                    self.pointerToListArrayNforce,
                    self.pointerToListArrayGrav,
                    self.pointerToX,
                    self.pointerToD,
                    self.pointerToDeld,
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
                    self.pointerToState0,
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
                    self.plotOutputFile,
                    self.ucdOutputRoot,
                    self.iterateEvent)

            # print "Just before lithomop3d.elastc:"

            # Perform elastic solution, if requested.
            
            lithomop3d.elastc(
                self.A,
                self.pointerToAlnz,
                self.pointerToPcg,
                self.pointerToZcg,
                self.pointerToDprev,
                self.pointerToJa,
                self.pointerToBextern,
                self.pointerToBtraction,
                self.pointerToBgravity,
                self.pointerToBconcForce,
                # self.pointerToBprestress,
                self.pointerToBintern,
                self.pointerToBresid,
                self.pointerToBwork,
                self.pointerToDispVec,
                self.pointerToListArrayNforce,
                self.pointerToListArrayGrav,
                self.pointerToX,
                self.pointerToD,
                self.pointerToDeld,
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
                self.pointerToState0,
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
                self.plotOutputFile,
                self.ucdOutputRoot,
                self.elasticStage,
                self.iterateEvent)

        # Perform time-dependent solution, if requested.

        if self.analysisType == "fullSolution" and self.numberTimeStepGroups > 1:
            lithomop3d.viscos(
                self.A,
                self.pointerToAlnz,
                self.pointerToPcg,
                self.pointerToZcg,
                self.pointerToDprev,
                self.pointerToJa,
                self.pointerToBextern,
                self.pointerToBtraction,
                self.pointerToBgravity,
                self.pointerToBconcForce,
                # self.pointerToBprestress,
                self.pointerToBintern,
                self.pointerToBresid,
                self.pointerToBwork,
                self.pointerToDispVec,
                self.pointerToListArrayNforce,
                self.pointerToListArrayGrav,
                self.pointerToX,
                self.pointerToD,
                self.pointerToDeld,
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
                self.pointerToState0,
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
                self.plotOutputFile,
                self.ucdOutputRoot,
                self.viscousStage,
                self.iterateEvent)
        lithomop3d.destroyPETScMat(self.A)
        lithomop3d.PetscFinalize()
	print ""
        print "Hello from lm3drun.run (end)!"
        return


    def __init__(self):
        Component.__init__(self, "lm3drun", "solver")

	print ""
        print "Hello from lm3drun.__init__!"

        return


# version
# $Id: Lithomop3d_run.py,v 1.12 2005/03/10 01:10:37 knepley Exp $

# End of file 
