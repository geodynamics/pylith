#!/usr/bin/env python
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

        # PETSc logging
        self.autoprestrStage = pl3dsetup.autoprestrStage
        self.elasticStage = pl3dsetup.elasticStage
        self.viscousStage = pl3dsetup.viscousStage
        self.iterateEvent = pl3dsetup.iterateEvent
        self.A = pl3dsetup.A
        self.rhs = pl3dsetup.rhs
        self.sol = pl3dsetup.sol
        
        self.analysisType = pl3dscan.inventory.analysisType

        # Import all necessary pointers, etc. from Pylith3d_setup.
        self.memorySize = pl3dsetup.memorySize
        self.intSize = pl3dsetup.intSize
        self.doubleSize = pl3dsetup.doubleSize
        
        self.numberTimeStepGroups = pl3dsetup.numberTimeStepGroups

        self.pointerToBextern = pl3dsetup.pointerToBextern
        self.pointerToBtraction = pl3dsetup.pointerToBtraction
        self.pointerToBgravity = pl3dsetup.pointerToBgravity
        self.pointerToBconcForce = pl3dsetup.pointerToBconcForce
        self.pointerToBintern = pl3dsetup.pointerToBintern
        self.pointerToBresid = pl3dsetup.pointerToBresid
        self.pointerToBwink = pl3dsetup.pointerToBwink
        self.pointerToBwinkx = pl3dsetup.pointerToBwinkx
        self.pointerToDispVec = pl3dsetup.pointerToDispVec
        self.pointerToDprev = pl3dsetup.pointerToDprev
        self.pointerToListArrayNforce = pl3dsetup.pointerToListArrayNforce
        self.pointerToListArrayGrav = pl3dsetup.pointerToListArrayGrav

        self.pointerToX = pl3dsetup.pointerToX
        self.pointerToD = pl3dsetup.pointerToD
        self.pointerToDeld = pl3dsetup.pointerToDeld
        self.pointerToDcur = pl3dsetup.pointerToDcur
        self.pointerToId = pl3dsetup.pointerToId
        self.pointerToIwink = pl3dsetup.pointerToIwink
        self.pointerToWink = pl3dsetup.pointerToWink
        self.pointerToListArrayNsysdat = pl3dsetup.pointerToListArrayNsysdat
        self.pointerToListArrayIddmat = pl3dsetup.pointerToListArrayIddmat

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
        self.pointerToState0 = pl3dsetup.pointerToState0
        self.pointerToDmat = pl3dsetup.pointerToDmat
        self.pointerToIens = pl3dsetup.pointerToIens
        self.pointerToLm = pl3dsetup.pointerToLm
        self.pointerToLmx = pl3dsetup.pointerToLmx
        self.pointerToLmf = pl3dsetup.pointerToLmf
        self.pointerToIvfamily = pl3dsetup.pointerToIvfamily
        self.pointerToListArrayNpar = pl3dsetup.pointerToListArrayNpar

        self.prestressAutoComputeInt = pl3dsetup.prestressAutoComputeInt

        self.pointerToIelno = pl3dsetup.pointerToIelno
        self.pointerToIside = pl3dsetup.pointerToIside
        self.pointerToIhistry = pl3dsetup.pointerToIhistry
        self.pointerToPres = pl3dsetup.pointerToPres
        self.pointerToPdir = pl3dsetup.pointerToPdir

        self.pointerToListArrayPropertyList = pl3dsetup.pointerToListArrayPropertyList
        self.pointerToMaterialModelInfo = pl3dsetup.pointerToMaterialModelInfo

        self.pointerToGauss = pl3dsetup.pointerToGauss
        self.pointerToSh = pl3dsetup.pointerToSh
        self.pointerToShj = pl3dsetup.pointerToShj
        self.pointerToListArrayElementTypeInfo = pl3dsetup.pointerToListArrayElementTypeInfo

        self.pointerToHistry = pl3dsetup.pointerToHistry
        self.pointerToListArrayRtimdat = pl3dsetup.pointerToListArrayRtimdat
        self.pointerToListArrayNtimdat = pl3dsetup.pointerToListArrayNtimdat
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

        self.pointerToSkew = pl3dsetup.pointerToSkew

        self.pointerToIprint = pl3dsetup.pointerToIprint
        self.pointerToListArrayNcodat = pl3dsetup.pointerToListArrayNcodat
        self.pointerToListArrayNunits = pl3dsetup.pointerToListArrayNunits
        self.pointerToListArrayNprint = pl3dsetup.pointerToListArrayNprint
        self.pointerToIstatout = pl3dsetup.pointerToIstatout
        self.pointerToNstatout = pl3dsetup.pointerToNstatout

        self.asciiOutputFile = pl3dsetup.asciiOutputFile
        self.plotOutputFile = pl3dsetup.plotOutputFile
        self.ucdOutputRoot = pl3dsetup.ucdOutputRoot

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

        # Output approximate memory usage
        self.memorySizeMB =0.0
        self.memorySizeMB=self.memorySize/(1024.0*1024.0)

        print ""
        print "Approximate memory allocation for f77 arrays (MB): %g" % self.memorySizeMB
        # print "Just before pylith3d.autoprestr:"

        # Compute gravitational prestresses, if requested.
        if self.analysisType == "elasticSolution" or self.analysisType == "fullSolution":
            if self.prestressAutoComputeInt == 1:
                pylith3d.autoprestr(
                    self.A,self.rhs,self.sol,
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
            
            pylith3d.elastc(
                self.A,self.rhs,self.sol,
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
            pylith3d.viscos(
                self.A,self.rhs,self.sol,
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
        pylith3d.destroyPETScMat(self.A,self.rhs,self.sol)
        pylith3d.PetscFinalize()
        print ""
        print "Hello from pl3drun.run (end)!"
        return


    def __init__(self):
        Component.__init__(self, "pl3drun", "solver")

        print ""
        print "Hello from pl3drun.__init__!"

        return


# version
# $Id: Pylith3d_run.py,v 1.17 2005/05/03 18:47:35 willic3 Exp $

# End of file 
