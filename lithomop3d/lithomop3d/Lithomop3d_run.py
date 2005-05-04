#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                             Charles A. Williams
#                       Rensselaer Polytechnic Institute
#                        (C) 2005  All Rights Reserved
#
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
        self.autoprestrStage = lm3dsetup.autoprestrStage
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

        self.pointerToBextern = lm3dsetup.pointerToBextern
        self.pointerToBtraction = lm3dsetup.pointerToBtraction
        self.pointerToBgravity = lm3dsetup.pointerToBgravity
        self.pointerToBconcForce = lm3dsetup.pointerToBconcForce
        self.pointerToBintern = lm3dsetup.pointerToBintern
        self.pointerToBresid = lm3dsetup.pointerToBresid
        self.pointerToBwink = lm3dsetup.pointerToBwink
        self.pointerToBwinkx = lm3dsetup.pointerToBwinkx
        self.pointerToDispVec = lm3dsetup.pointerToDispVec
        self.pointerToDprev = lm3dsetup.pointerToDprev
        self.pointerToListArrayNforce = lm3dsetup.pointerToListArrayNforce
        self.pointerToListArrayGrav = lm3dsetup.pointerToListArrayGrav

        self.pointerToX = lm3dsetup.pointerToX
        self.pointerToD = lm3dsetup.pointerToD
        self.pointerToDeld = lm3dsetup.pointerToDeld
        self.pointerToDcur = lm3dsetup.pointerToDcur
        self.pointerToId = lm3dsetup.pointerToId
        self.pointerToIwink = lm3dsetup.pointerToIwink
        self.pointerToWink = lm3dsetup.pointerToWink
        self.pointerToListArrayNsysdat = lm3dsetup.pointerToListArrayNsysdat
        self.pointerToListArrayIddmat = lm3dsetup.pointerToListArrayIddmat

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
        self.pointerToIens = lm3dsetup.pointerToIens
        self.pointerToLm = lm3dsetup.pointerToLm
        self.pointerToLmx = lm3dsetup.pointerToLmx
        self.pointerToLmf = lm3dsetup.pointerToLmf
        self.pointerToIvfamily = lm3dsetup.pointerToIvfamily
        self.pointerToListArrayNpar = lm3dsetup.pointerToListArrayNpar

        self.prestressAutoComputeInt = lm3dsetup.prestressAutoComputeInt

        self.pointerToIelno = lm3dsetup.pointerToIelno
        self.pointerToIside = lm3dsetup.pointerToIside
        self.pointerToIhistry = lm3dsetup.pointerToIhistry
        self.pointerToPres = lm3dsetup.pointerToPres
        self.pointerToPdir = lm3dsetup.pointerToPdir

        self.pointerToListArrayPropertyList = lm3dsetup.pointerToListArrayPropertyList
        self.pointerToMaterialModelInfo = lm3dsetup.pointerToMaterialModelInfo

        self.pointerToGauss = lm3dsetup.pointerToGauss
        self.pointerToSh = lm3dsetup.pointerToSh
        self.pointerToShj = lm3dsetup.pointerToShj
        self.pointerToListArrayElementTypeInfo = lm3dsetup.pointerToListArrayElementTypeInfo

        self.pointerToHistry = lm3dsetup.pointerToHistry
        self.pointerToListArrayRtimdat = lm3dsetup.pointerToListArrayRtimdat
        self.pointerToListArrayNtimdat = lm3dsetup.pointerToListArrayNtimdat
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

        self.pointerToSkew = lm3dsetup.pointerToSkew

        self.pointerToIprint = lm3dsetup.pointerToIprint
        self.pointerToListArrayNcodat = lm3dsetup.pointerToListArrayNcodat
        self.pointerToListArrayNunits = lm3dsetup.pointerToListArrayNunits
        self.pointerToListArrayNprint = lm3dsetup.pointerToListArrayNprint
        self.pointerToIstatout = lm3dsetup.pointerToIstatout
        self.pointerToNstatout = lm3dsetup.pointerToNstatout

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
                    self.pointerToBextern,
                    self.pointerToBtraction,
                    self.pointerToBgravity,
                    self.pointerToBconcForce,
                    self.pointerToBintern,
                    self.pointerToBresid,
                    self.pointerToBwink,
                    self.pointerToBwinkx,
                    self.pointerToDispVec,
                    self.pointerToDprev,
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
                    self.pointerToListArrayIddmat,
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
                    self.pointerToIens,
                    self.pointerToLm,
                    self.pointerToLmx,
                    self.pointerToLmf,
                    self.pointerToIvfamily,
                    self.pointerToListArrayNpar,
                    self.pointerToIelno,
                    self.pointerToIside,
                    self.pointerToIhistry,
                    self.pointerToPres,
                    self.pointerToPdir,
                    self.pointerToListArrayPropertyList,
                    self.pointerToMaterialModelInfo,
                    self.pointerToGauss,
                    self.pointerToSh,
                    self.pointerToShj,
                    self.pointerToListArrayElementTypeInfo,
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
                    self.pointerToSkew,
                    self.pointerToListArrayNcodat,
                    self.pointerToListArrayNunits,
                    self.pointerToListArrayNprint,
                    self.pointerToIstatout,
                    self.pointerToNstatout,
                    self.asciiOutputFile,
                    self.plotOutputFile,
                    self.ucdOutputRoot,
                    self.autoprestrStage,
                    self.iterateEvent)

            # Perform elastic solution, if requested.
            
            lithomop3d.elastc(
                self.A,
                self.pointerToBextern,
                self.pointerToBtraction,
                self.pointerToBgravity,
                self.pointerToBconcForce,
                self.pointerToBintern,
                self.pointerToBresid,
                self.pointerToBwink,
                self.pointerToBwinkx,
                self.pointerToDispVec,
                self.pointerToDprev,
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
                self.pointerToListArrayIddmat,
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
                self.pointerToIens,
                self.pointerToLm,
                self.pointerToLmx,
                self.pointerToLmf,
                self.pointerToIvfamily,
                self.pointerToListArrayNpar,
                self.pointerToIelno,
                self.pointerToIside,
                self.pointerToIhistry,
                self.pointerToPres,
                self.pointerToPdir,
                self.pointerToListArrayPropertyList,
                self.pointerToMaterialModelInfo,
                self.pointerToGauss,
                self.pointerToSh,
                self.pointerToShj,
                self.pointerToListArrayElementTypeInfo,
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
                self.pointerToSkew,
                self.pointerToListArrayNcodat,
                self.pointerToListArrayNunits,
                self.pointerToListArrayNprint,
                self.pointerToIstatout,
                self.pointerToNstatout,
                self.asciiOutputFile,
                self.plotOutputFile,
                self.ucdOutputRoot,
                self.elasticStage,
                self.iterateEvent)

        # Perform time-dependent solution, if requested.

        if self.analysisType == "fullSolution" and self.numberTimeStepGroups > 1:
            lithomop3d.viscos(
                self.A,
                self.pointerToBextern,
                self.pointerToBtraction,
                self.pointerToBgravity,
                self.pointerToBconcForce,
                self.pointerToBintern,
                self.pointerToBresid,
                self.pointerToBwink,
                self.pointerToBwinkx,
                self.pointerToDispVec,
                self.pointerToDprev,
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
                self.pointerToListArrayIddmat,
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
                self.pointerToIens,
                self.pointerToLm,
                self.pointerToLmx,
                self.pointerToLmf,
                self.pointerToIvfamily,
                self.pointerToListArrayNpar,
                self.pointerToIelno,
                self.pointerToIside,
                self.pointerToIhistry,
                self.pointerToPres,
                self.pointerToPdir,
                self.pointerToListArrayPropertyList,
                self.pointerToMaterialModelInfo,
                self.pointerToGauss,
                self.pointerToSh,
                self.pointerToShj,
                self.pointerToListArrayElementTypeInfo,
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
                self.pointerToSkew,
                self.pointerToIprint,
                self.pointerToListArrayNcodat,
                self.pointerToListArrayNunits,
                self.pointerToListArrayNprint,
                self.pointerToIstatout,
                self.pointerToNstatout,
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
# $Id: Lithomop3d_run.py,v 1.17 2005/05/03 18:47:35 willic3 Exp $

# End of file 
