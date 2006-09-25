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

class ElementTypeDef:

    def getdef(self,
               elementType,
               quadratureOrderInt,
               numberSpaceDimensions,
               numberDegreesFreedom):

	import pylith3d

        print "Hello from ElementTypeDef.getdef!"

        self.numberVolumeElementNodes = self.elementNodes[elementType - 1]
        self.numberSurfaceElementNodes = self.elementNodes2d[elementType - 1]

        if quadratureOrderInt == 1:
            self.numberVolumeElementGaussPoints = \
                                          self.elementFullGauss[elementType -1]
            self.numberSurfaceElementGaussPoints = \
                                          self.elementFullGauss2d[elementType -1]
        elif quadratureOrderInt == 2:
            self.numberVolumeElementGaussPoints =  \
                                          self.elementReducedGauss[elementType -1]
            self.numberSurfaceElementGaussPoints =  \
                                          self.elementReducedGauss2d[elementType -1]
        elif quadratureOrderInt == 3:
            self.numberVolumeElementGaussPoints = \
                                           self.elementFullGauss[elementType -1]
            self.numberSurfaceElementGaussPoints = \
                                           self.elementFullGauss2d[elementType -1]

        self.numberVolumeElementEquations = \
                                    numberDegreesFreedom * \
                                    self.numberVolumeElementNodes

        self.numberSurfaceElementEquations = \
                                    numberDegreesFreedom * \
                                    self.numberSurfaceElementNodes

        self.numberSurfaceElementCoordinates = \
                                      numberSpaceDimensions * \
                                      self.numberSurfaceElementNodes

        self.elementTypeInfo = [self.numberVolumeElementNodes,
                                self.numberVolumeElementGaussPoints,
                                self.numberVolumeElementEquations,
                                self.numberVolumeElementCoordinates]

        self.elementTypeInfo2d = [self.numberSurfaceElementNodes,
                                  self.numberSurfaceElementGaussPoints,
                                  self.numberSurfaceElementEquations,
                                  self.numberSurfaceElementCoordinates]

        self.pointerToSh = pylith3d.allocateDouble(
            (numberSpaceDimensions+1)*
            self.numberVolumeElementNodes*
            self.numberVolumeElementGaussPoints)

        self.pointerToSh2d = pylith3d.allocateDouble(
            numberSpaceDimensions*
            self.numberSurfaceElementNodes*
            self.numberSurfaceElementGaussPoints)
            
        self.pointerToShj = pylith3d.allocateDouble(
            (numberSpaceDimensions+1)*
            self.numberVolumeElementNodes*
            self.numberVolumeElementGaussPoints)
            
        self.pointerToGauss = pylith3d.allocateDouble(
            (numberSpaceDimensions+1)*
            self.numberVolumeElementGaussPoints)
            
        self.pointerToGauss2d = pylith3d.allocateDouble(
            numberSpaceDimensions*
            self.numberSurfaceElementGaussPoints)

        pylith3d.preshape(
            self.pointerToSh,
            self.pointerToShj,
            self.pointerToGauss,
            quadratureOrderInt,
            elementType,
            self.numberVolumeElementNodes,
            self.numberVolumeElementGaussPoints)

        pylith3d.preshape2d(
            self.pointerToSh2d,
            self.pointerToGauss2d,
            quadratureOrderInt,
            elementType,
            self.numberSurfaceElementNodes,
            self.numberSurfaceElementGaussPoints)
        
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
        self.numberSurfaceElementNodes = 0
        self.numberSurfaceElementGaussPoints = 0
        self.numberSurfaceElementEquations = 0
        self.numberSurfaceElementCoordinates = 0
        self.elementTypeInfo2d = [0, 0, 0, 0]
        self.pointerToListArrayElementTypeInfo2d = None
        self.pointerToSh2d = None
        self.pointerToGauss2d = None
        self.elementNodes2d = [ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, \
                                4, 4, 4, 4, 4, 4, 4, 4, 4, 4, \
                                4, 4, 4, 4, 4, 4, 4, 4, 4, 3, \
                                3, 8, 8, 8, 8, 8, 8, 8, 8, 8, \
                                8, 8, 8, 8, 8, 8, 8, 8, 8, 8, \
                                8, 8, 8, 8, 8, 8, 8, 8, 8, 8, \
                                6, 6]
        self.elementFullGauss2d = [ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, \
                                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, \
                                    4, 4, 4, 4, 4, 4, 4, 4, 4, 1, \
                                    1, 9, 9, 9, 9, 9, 9, 9, 9, 9, \
                                    9, 9, 9, 9, 9, 9, 9, 9, 9, 9, \
                                    9, 9, 9, 9, 9, 9, 9, 9, 9, 8, \
                                    3, 3]
        self.elementReducedGauss2d = [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, \
                                       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, \
                                       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, \
                                       1, 4, 4, 4, 4, 4, 4, 4, 4, 4, \
                                       4, 4, 4, 4, 4, 4, 4, 4, 4, 4, \
                                       4, 4, 4, 4, 4, 4, 4, 4, 4, 4, \
                                       1, 1]
        return

# version
# $Id: ElementTypeDef.py,v 1.2 2005/04/01 23:40:46 willic3 Exp $

# End of file 
