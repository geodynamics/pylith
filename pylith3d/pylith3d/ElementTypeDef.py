#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                              Charles A. Williams
#                        Rensselaer Polytechnic Institute
#                      (C) 2005  All Rights Reserved
#
#  Copyright 2005 Rensselaer Polytechnic Institute.
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

class ElementTypeDef:

    def getdef(self,
               elementType,
               quadratureOrderInt,
               numberSpaceDimensions,
               numberDegreesFreedom):

	import pylith3d

        print "Hello from ElementTypeDef.getdef!"

        self.numberVolumeElementNodes = self.elementNodes[elementType - 1]

        if quadratureOrderInt == 1:
            self.numberVolumeElementGaussPoints = \
                                          self.elementFullGauss[elementType -1]
        elif quadratureOrderInt == 2:
            self.numberVolumeElementGaussPoints =  \
                                          self.elementReducedGauss[elementType -1]
        elif quadratureOrderInt == 3:
            self.numberVolumeElementGaussPoints = \
                                           self.elementFullGauss[elementType -1]

        self.numberVolumeElementEquations = \
                                    numberDegreesFreedom * \
                                    self.numberVolumeElementNodes

        self.numberVolumeElementCoordinates = \
                                      numberSpaceDimensions * \
                                      self.numberVolumeElementNodes

        self.elementTypeInfo = [self.numberVolumeElementNodes,
                                self.numberVolumeElementGaussPoints,
                                self.numberVolumeElementEquations,
                                self.numberVolumeElementCoordinates]

        self.pointerToSh = pylith3d.allocateDouble(
            (numberSpaceDimensions+1)*
            self.numberVolumeElementNodes*
            self.numberVolumeElementGaussPoints)
            
        self.pointerToShj = pylith3d.allocateDouble(
            (numberSpaceDimensions+1)*
            self.numberVolumeElementNodes*
            self.numberVolumeElementGaussPoints)
            
        self.pointerToGauss = pylith3d.allocateDouble(
            (numberSpaceDimensions+1)*
            self.numberVolumeElementGaussPoints)

        pylith3d.preshape(
            self.pointerToSh,
            self.pointerToShj,
            self.pointerToGauss,
            quadratureOrderInt,
            elementType,
            self.numberVolumeElementNodes,
            self.numberVolumeElementGaussPoints)
        
        return


    def __init__(self):
	print ""
        print "Hello from ElementTypeDef.__init__!"
        self.numberVolumeElementNodes = 0
        self.numberVolumeElementGaussPoints = 0
        self.numberVolumeElementEquations = 0
        self.numberVolumeElementCoordinates = 0
        self.elementTypeInfo = [0, 0, 0, 0]
        self.pointerToListArrayElementTypeInfo = None
        self.pointerToSh = None
        self.pointerToShj = None
        self.pointerToGauss = None
        self.elementNodes = [ 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, \
                              8, 8, 8, 8, 8, 8, 8, 8, 8, 8, \
                              8, 8, 8, 8, 8, 8, 8, 7, 6, 5, \
                              4,20,20,20,20,20,20,20,20,20, \
                             20,20,20,20,20,20,20,20,20,20, \
                             20,20,20,20,20,20,20,20,18,15, \
                             13,10]
        self.elementFullGauss = [ 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, \
                                  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, \
                                  8, 8, 8, 8, 8, 8, 8, 8, 2, 5, \
                                  1,27,27,27,27,27,27,27,27,27, \
                                 27,27,27,27,27,27,27,27,27,27, \
                                 27,27,27,27,27,27,27,27,27, 9, \
                                 13, 4]
        self.elementReducedGauss = [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, \
                                     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, \
                                     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, \
                                     1, 8, 8, 8, 8, 8, 8, 8, 8, 8, \
                                     8, 8, 8, 8, 8, 8, 8, 8, 8, 8, \
                                     8, 8, 8, 8, 8, 8, 8, 8, 8, 2, \
                                     5, 1]
        return

# version
# $Id: ElementTypeDef.py,v 1.2 2005/04/01 23:40:46 willic3 Exp $

# End of file 
