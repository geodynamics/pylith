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

        print "Hello from ElementTypeDef.getdef!"

        self.numberElementNodes = self.elementNodes[elementType - 1]


        if quadratureOrderInt == 1:
            self.numberElementGaussPoints = \
                                          self.elementFullGauss[elementType -1]
        elif quadratureOrderInt == 2:
            self.numberElementGaussPoints =  \
                                          self.elementReducedGauss[elementType -1]
        elif quadratureOrderInt == 3:
            self.numberElementGaussPoints = \
                                           self.elementFullGauss[elementType -1]


        self.numberElementEquations = \
                                    numberDegreesFreedom * \
                                    self.numberElementNodes
        self.numberElementCoordinates = \
                                      numberSpaceDimensions * \
                                      self.numberElementNodes

        self.elementTypeInfo = [self.numberElementNodes,
                                self.numberElementGaussPoints,
                                self.numberElementEquations,
                                self.numberElementCoordinates]

        self.pointerToSh = lithomop3d.allocateDouble(
            (numberSpaceDimensions+1)*
            numberElementNodes*
            numberElementGaussPoints)
            
        self.pointerToShj = lithomop3d.allocateDouble(
            (numberSpaceDimensions+1)*
            numberElementNodes*
            numberElementGaussPoints)
            
        self.pointerToGauss = lithomop3d.allocateDouble(
            (numberSpaceDimensions+1)*
            numberElementGaussPoints)

        lithomop3d.preshape(
            self.pointerToSh,
            self.pointerToShj,
            self.pointerToGauss,
            quadratureOrderInt,
            elementType,
            self.numberElementNodes,
            self.numberElementGaussPoints)
        
        return


    def __init__(self):
	print ""
        print "Hello from ElementTypeDef.__init__!"
        self.numberElementNodes = 0
        self.numberElementGaussPoints = 0
        self.numberElementEquations = 0
        self.numberElementCoordinates = 0
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
# $Id: ElementTypeDef.py,v 1.1 2005/03/22 02:26:26 willic3 Exp $

# End of file 
