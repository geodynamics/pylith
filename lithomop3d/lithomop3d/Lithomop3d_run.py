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

        print "Hello from lm3drun.initialize (begin)!"

        # The only parameters required from Lithomop3d_scan are those in the
        # inventory.  All others have been imported into Lithomop3d_setup, and
        # possibly modified there.  Get all required info from the inventory.
        
        self.analysisType = lm3dscan.inventory.analysisType

        # Import all necessary pointers, etc. from Lithomop3d_setup.
        
        self.numberTimeStepGroups = lm3dsetup.numberTimeStepGroups

        self.pointerToAlnz = lm3dsetup.pointerToAlnz
        self.pointerToPcg = lm3dsetup.pointerToPcg
        self.pointerToZcg = lm3dsetup.pointerToZcg
        self.pointerToJa = lm3dsetup.pointerToJa
        self.pointerToB = lm3dsetup.pointerToB
        self.pointerToBtot = lm3dsetup.pointerToBtot
        self.pointerToBres = lm3dsetup.pointerToBres
        self.pointerToPvec = lm3dsetup.pointerToPvec
        self.pointerToGvec1 = lm3dsetup.pointerToGvec1
        self.pointerToGvec2 = lm3dsetup.pointerToGvec2
        self.pointerToX = lm3dsetup.pointerToX
        self.pointerToD = lm3dsetup.pointerToD
        self.pointerToDx = lm3dsetup.pointerToDx
        self.pointerToDeld = lm3dsetup.pointerToDeld
        self.pointerToDeldx = lm3dsetup.pointerToDeldx
        self.pointerToDprev = lm3dsetup.pointerToDprev
        self.pointerToDcur = lm3dsetup.pointerToDcur
        self.pointerToDxcur = lm3dsetup.pointerToDxcur
        self.pointerToId = lm3dsetup.pointerToId
        self.pointerToIdx = lm3dsetup.pointerToIdx
        self.pointerToSkew = lm3dsetup.pointerToSkew
        self.pointerToHistry = lm3dsetup.pointerToHistry
        self.pointerToIen = lm3dsetup.pointerToIen
        self.pointerToInfin = lm3dsetup.pointerToInfin
        self.pointerToMat = lm3dsetup.pointerToMat
        self.pointerToLm = lm3dsetup.pointerToLm
        self.pointerToLmx = lm3dsetup.pointerToLmx
        self.pointerToLmf = lm3dsetup.pointerToLmf
        self.pointerToProp = lm3dsetup.pointerToProp
        self.pointerToListArrayGauss = lm3dsetup.pointerToListArrayGauss
        self.pointerToIbond = lm3dsetup.pointerToIbond
        self.pointerToBond = lm3dsetup.pointerToBond
        self.pointerToDmat = lm3dsetup.pointerToDmat
        self.pointerToStn = lm3dsetup.pointerToStn
        self.pointerToScur = lm3dsetup.pointerToScur
        self.pointerToSt0 = lm3dsetup.pointerToSt0
        self.pointerToEps = lm3dsetup.pointerToEps
        self.pointerToDeps = lm3dsetup.pointerToDeps
        self.pointerToBeta = lm3dsetup.pointerToBeta
        self.pointerToDbeta = lm3dsetup.pointerToDbeta
        self.pointerToBetb = lm3dsetup.pointerToBetb
        self.pointerToDbetb = lm3dsetup.pointerToDbetb
        self.pointerToListArrayIddmat = lm3dsetup.pointerToListArrayIddmat
        self.pointerToIelno = lm3dsetup.pointerToIelno
        self.pointerToIside = lm3dsetup.pointerToIside
        self.pointerToIhistry = lm3dsetup.pointerToIhistry
        self.pointerToPres = lm3dsetup.pointerToPres
        self.pointerToPdir = lm3dsetup.pointerToPdir
        self.pointerToMaxstp = lm3dsetup.pointerToMaxstp
        self.pointerToDelt = lm3dsetup.pointerToDelt
        self.pointerToAlfa = lm3dsetup.pointerToAlfa
        self.pointerToMaxit = lm3dsetup.pointerToMaxit
        self.pointerToMaxitc = lm3dsetup.pointerToMaxitc
        self.pointerToLgdef = lm3dsetup.pointerToLgdef
        self.pointerToIbbar = lm3dsetup.pointerToIbbar
        self.pointerToUtol = lm3dsetup.pointerToUtol
        self.pointerToFtol = lm3dsetup.pointerToFtol
        self.pointerToEtol = lm3dsetup.pointerToEtol
        self.pointerToItmax = lm3dsetup.pointerToItmax
        self.pointerToIprint = lm3dsetup.pointerToIprint
        self.pointerToFault = lm3dsetup.pointerToFault
        self.pointerToNfault = lm3dsetup.pointerToNfault
        self.pointerToDfault = lm3dsetup.pointerToDfault
        self.pointerToTfault = lm3dsetup.pointerToTfault
        self.pointerToIdftn = lm3dsetup.pointerToIdftn
        self.pointerToIdslp = lm3dsetup.pointerToIdslp
        self.pointerToIpslp = lm3dsetup.pointerToIpslp
        self.pointerToDiforc = lm3dsetup.pointerToDiforc
        self.pointerToIdhist = lm3dsetup.pointerToIdhist
        self.pointerToIwink = lm3dsetup.pointerToIwink
        self.pointerToWink = lm3dsetup.pointerToWink
        self.pointerToIwinkx = lm3dsetup.pointerToIwinkx
        self.pointerToWinkx = lm3dsetup.pointerToWinkx
        self.pointerToS = lm3dsetup.pointerToS
        self.pointerToStemp = lm3dsetup.pointerToStemp
        self.pointerToListArrayGrav = lm3dsetup.pointerToListArrayGrav
        self.pointerToListArrayNcodat = lm3dsetup.pointerToListArrayNcodat
        self.pointerToListArrayNconsts = lm3dsetup.pointerToListArrayNconsts
        self.pointerToListArrayNdimens = lm3dsetup.pointerToListArrayNdimens
        self.pointerToListArrayNpar = lm3dsetup.pointerToListArrayNpar
        self.pointerToListArrayNprint = lm3dsetup.pointerToListArrayNprint
        self.pointerToListArrayNsiter = lm3dsetup.pointerToListArrayNsiter
        self.pointerToListArrayNsysdat = lm3dsetup.pointerToListArrayNsysdat
        self.pointerToListArrayNunits = lm3dsetup.pointerToListArrayNunits
        self.pointerToListArrayNvisdat = lm3dsetup.pointerToListArrayNvisdat
        self.pointerToListArrayRconsts = lm3dsetup.pointerToListArrayRconsts
        self.pointerToListArrayRgiter = lm3dsetup.pointerToListArrayRgiter
        self.pointerToListArrayRmin = lm3dsetup.pointerToListArrayRmin
        self.pointerToListArrayRmult = lm3dsetup.pointerToListArrayRmult
        self.pointerToListArrayRtimdat = lm3dsetup.pointerToListArrayRtimdat
        self.asciiOutputFile = lm3dsetup.asciiOutputFile
        self.plotOutputFile = lm3dsetup.plotOutputFile
        self.viscousFlagInt = lm3dsetup.viscousFlagInt
        self.plasticFlagInt = lm3dsetup.plasticFlagInt
        self.materialHistoryFlagInt = lm3dsetup.materialHistoryFlagInt

        print "Hello from lm3drun.initialize (end)!"

        return

    def run(self):
        import lithomop3d
        
        # First define all of the lists that maintain variable values.  The
        # variables in these lists are altered during the running of the code
        # and should not be accessed directly except as a member of the list.
        # They should not have been defined previously.

        print "Hello from lm3drun.run (begin)!"

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

        # ntimdat array
        self.currentTimeStep = 0
        self.currentStepsBetweenReform = 0
        self.currentIterationsBetweenReform = 0
        self.currentLargeDeformationFlag = 0
        self.currentBbarFlag = 0
        self.currentMaximumIterations = 0
        self.currentNumberTotalIterations = 0
        self.currentNumberReforms = 0
        self.currentNumberTotalPcgIterations = 0
        self.reformFlagInt = 0
        self.listNtimdat = [
            self.currentTimeStep,
            self.currentStepsBetweenReform,
            self.currentIterationsBetweenReform,
            self.currentLargeDeformationFlag,
            self.currentBbarFlag,
            self.currentMaximumIterations,
            self.currentNumberTotalIterations,
            self.currentNumberReforms,
            self.currentNumberTotalPcgIterations,
            self.reformFlagInt,
            self.viscousFlagInt,
            self.plasticFlagInt,
            self.materialHistoryFlagInt]
        self.pointerToListArrayNtimdat = lithomop3d.intListToArray(
            self.listNtimdat)


        # Perform elastic solution, if requested.

        if self.analysisType == "elasticSolution" or self.analysisType == "fullSolution":
            lithomop3d.elastc(
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
                self.pointerToX,
                self.pointerToD,
                self.pointerToDx,
                self.pointerToDeld,
                self.pointerToDeldx,
                self.pointerToDprev,
                self.pointerToDcur,
                self.pointerToDxcur,
                self.pointerToId,
                self.pointerToIdx,
                self.pointerToSkew,
                self.pointerToHistry,
                self.pointerToIen,
                self.pointerToInfin,
                self.pointerToMat,
                self.pointerToLm,
                self.pointerToLmx,
                self.pointerToLmf,
                self.pointerToProp,
                self.pointerToListArrayGauss,
                self.pointerToIbond,
                self.pointerToBond,
                self.pointerToDmat,
                self.pointerToStn,
                self.pointerToScur,
                self.pointerToSt0,
                self.pointerToEps,
                self.pointerToDeps,
                self.pointerToBeta,
                self.pointerToDbeta,
                self.pointerToBetb,
                self.pointerToDbetb,
                self.pointerToListArrayIddmat,
                self.pointerToIelno,
                self.pointerToIside,
                self.pointerToIhistry,
                self.pointerToPres,
                self.pointerToPdir,
                self.pointerToMaxstp,
                self.pointerToDelt,
                self.pointerToAlfa,
                self.pointerToMaxit,
                self.pointerToMaxitc,
                self.pointerToLgdef,
                self.pointerToIbbar,
                self.pointerToUtol,
                self.pointerToFtol,
                self.pointerToEtol,
                self.pointerToItmax,
                self.pointerToFault,
                self.pointerToNfault,
                self.pointerToDfault,
                self.pointerToTfault,
                self.pointerToIdftn,
                self.pointerToIdslp,
                self.pointerToIpslp,
                self.pointerToDiforc,
                self.pointerToIdhist,
                self.pointerToIwink,
                self.pointerToWink,
                self.pointerToIwinkx,
                self.pointerToWinkx,
                self.pointerToS,
                self.pointerToStemp,
                self.pointerToListArrayGcurr,
                self.pointerToListArrayGi,
                self.pointerToListArrayGprev,
                self.pointerToListArrayGrav,
                self.pointerToListArrayGtol,
                self.pointerToListArrayNcodat,
                self.pointerToListArrayNconsts,
                self.pointerToListArrayNdimens,
                self.pointerToListArrayNpar,
                self.pointerToListArrayNprint,
                self.pointerToListArrayNsiter,
                self.pointerToListArrayNsysdat,
                self.pointerToListArrayNtimdat,
                self.pointerToListArrayNunits,
                self.pointerToListArrayNvisdat,
                self.pointerToListArrayRconsts,
                self.pointerToListArrayRgiter,
                self.pointerToListArrayRmin,
                self.pointerToListArrayRmult,
                self.pointerToListArrayRtimdat,
                self.asciiOutputFile,
                self.plotOutputFile)

        # Perform time-dependent solution, if requested.

        if self.analysisType == "fullSolution" and self.numberTimeStepGroups > 1:
            lithomop3d.viscos(
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
                self.pointerToX,
                self.pointerToD,
                self.pointerToDx,
                self.pointerToDeld,
                self.pointerToDeldx,
                self.pointerToDprev,
                self.pointerToDcur,
                self.pointerToDxcur,
                self.pointerToId,
                self.pointerToIdx,
                self.pointerToSkew,
                self.pointerToHistry,
                self.pointerToIen,
                self.pointerToInfin,
                self.pointerToMat,
                self.pointerToLm,
                self.pointerToLmx,
                self.pointerToLmf,
                self.pointerToProp,
                self.pointerToListArrayGauss,
                self.pointerToIbond,
                self.pointerToBond,
                self.pointerToDmat,
                self.pointerToStn,
                self.pointerToScur,
                self.pointerToSt0,
                self.pointerToEps,
                self.pointerToDeps,
                self.pointerToBeta,
                self.pointerToDbeta,
                self.pointerToBetb,
                self.pointerToDbetb,
                self.pointerToListArrayIddmat,
                self.pointerToIelno,
                self.pointerToIside,
                self.pointerToIhistry,
                self.pointerToPres,
                self.pointerToPdir,
                self.pointerToMaxstp,
                self.pointerToDelt,
                self.pointerToAlfa,
                self.pointerToMaxit,
                self.pointerToMaxitc,
                self.pointerToLgdef,
                self.pointerToIbbar,
                self.pointerToUtol,
                self.pointerToFtol,
                self.pointerToEtol,
                self.pointerToItmax,
                self.pointerToIprint,
                self.pointerToFault,
                self.pointerToNfault,
                self.pointerToDfault,
                self.pointerToTfault,
                self.pointerToIdftn,
                self.pointerToIdslp,
                self.pointerToIpslp,
                self.pointerToDiforc,
                self.pointerToIdhist,
                self.pointerToIwink,
                self.pointerToWink,
                self.pointerToIwinkx,
                self.pointerToWinkx,
                self.pointerToS,
                self.pointerToStemp,
                self.pointerToListArrayGcurr,
                self.pointerToListArrayGi,
                self.pointerToListArrayGprev,
                self.pointerToListArrayGrav,
                self.pointerToListArrayGtol,
                self.pointerToListArrayNcodat,
                self.pointerToListArrayNconsts,
                self.pointerToListArrayNdimens,
                self.pointerToListArrayNpar,
                self.pointerToListArrayNprint,
                self.pointerToListArrayNsiter,
                self.pointerToListArrayNsysdat,
                self.pointerToListArrayNtimdat,
                self.pointerToListArrayNunits,
                self.pointerToListArrayNvisdat,
                self.pointerToListArrayRconsts,
                self.pointerToListArrayRgiter,
                self.pointerToListArrayRmin,
                self.pointerToListArrayRmult,
                self.pointerToListArrayRtimdat,
                self.asciiOutputFile,
                self.plotOutputFile)
                          
        print "Hello from lm3drun.run (end)!"
        return


    def __init__(self):
        Component.__init__(self, "lm3drun", "solver")

        print "Hello from lm3drun.__init__!"

        return


# version
# $Id: Lithomop3d_run.py,v 1.1 2004/04/14 21:22:47 willic3 Exp $

# End of file 
