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
