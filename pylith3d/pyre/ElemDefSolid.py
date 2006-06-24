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

# This could probably be set up as a Pyre component, but it doesn't seem
# necessary right now.

class ElemDefSolid:
  """Definitions for solid elements used by pylith3d."""

  def getdef(self, elemTypeInt, quadratureOrderInt):
      """Gets dimension and quadrature info, depending on the global
      quadrature order setting."""

    import pylith3d as pl3d

    print "Hello from ElemDefSolid.getdef!"

    numElemNodes = self.elemNodes[elemTypeInt - 1]
    self.elemInfoSolid['numElemNodes'] = numElemNodes
    self.elemInfoSolid['numElemEquations'] = self.numDegreesFreedom * numElemNodes
    self.elemInfoSolid['numElemCoords'] = self.numSpaceDims * numElemNodes

    # Full integration order
    if quadratureOrderInt == 1:
      numElemGaussPoints = self.elemFullGauss[elementType -1]
      self.elemInfoSolid['fptrBmatrix'] = pl3d.bmatrixn
      self.elemInfoSolid['fptrGetShape'] = pl3d.getshapn

    # Reduced integration order
    elif quadratureOrderInt == 2:
      numElemGaussPoints = self.elemReducedGauss[elementType -1]
      self.elemInfoSolid['fptrBmatrix'] = pl3d.bmatrixn
      self.elemInfoSolid['fptrGetShape'] = pl3d.getshapn

    # Selective (B-bar) integration order
    elif quadratureOrderInt == 3:
      numElemGaussPoints = self.elemFullGauss[elementType -1]
      self.elemInfoSolid['fptrBmatrix'] = pl3d.bmatrixb
      self.elemInfoSolid['fptrGetShape'] = pl3d.getshapb

    self.elemInfoSolid['numElemGaussPoints'] = numElemGaussPoints

    ptrShape = pl3d.allocateDouble(
      (self.numSpaceDims+1)*
      numElemNodes*
      numElemGaussPoints)
            
    ptrShapej = pl3d.allocateDouble(
      (self.numSpaceDims+1)*
      numElemNodes*
      numElemGaussPoints)
            
    ptrGauss = pl3d.allocateDouble(
      (self.numSpaceDims+1)*
      numElemGaussPoints)

    pl3d.preshape(
      ptrShape,
      ptrShapej,
      ptrGauss,
      quadratureOrderInt,
      elemTypeInt,
      numElemNodes,
      numElemGaussPoints)

    self.elemInfoSolid['ptrShape'] = ptrShape
    self.elemInfoSolid['ptrShapej'] = ptrShapej
    self.elemInfoSolid['ptrGauss'] = ptrGauss
        
    return self.elemInfoSolid


    def __init__(self):
      """Constructor."""
      print ""
      print "Hello from ElemDefSolid.__init__!"

      # I am leaving these hard-wired for now, since solid elements will always
      # have 3 DOF and 3 spatial dimensions, but it might be better to leave
      # these as variables for more generality.  Ideally, we would get this
      # info from global definitions.
      self.numDegreesFreedom = 3
      self.numSpaceDims = 3

      self.elemInfoSolid = {'numElemNodes': 0,
                            'numElemEquations': 0,
                            'numElemCoords': 0,
                            'numElemGaussPoints': 0,
                            'ptrShape': None,
                            'ptrShapej': None,
                            'ptrGauss': None,
                            'fptrBmatrix': None,
                            'fptrGetShape': None}

      self.elemNodes = [ 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, \
                         8, 8, 8, 8, 8, 8, 8, 8, 8, 8, \
                         8, 8, 8, 8, 8, 8, 8, 7, 6, 5, \
                         4,20,20,20,20,20,20,20,20,20, \
                         20,20,20,20,20,20,20,20,20,20, \
                         20,20,20,20,20,20,20,20,18,15, \
                         13,10]
      self.elemFullGauss = [ 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, \
                             8, 8, 8, 8, 8, 8, 8, 8, 8, 8, \
                             8, 8, 8, 8, 8, 8, 8, 8, 2, 5, \
                             1,27,27,27,27,27,27,27,27,27, \
                             27,27,27,27,27,27,27,27,27,27, \
                             27,27,27,27,27,27,27,27,27, 9, \
                             13, 4]
      self.elemReducedGauss = [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, \
                                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, \
                                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, \
                                1, 8, 8, 8, 8, 8, 8, 8, 8, 8, \
                                8, 8, 8, 8, 8, 8, 8, 8, 8, 8, \
                                8, 8, 8, 8, 8, 8, 8, 8, 8, 2, \
                                5, 1]
      return

# version
# $Id: ElemDefSolid.py,v 1.2 2005/04/01 23:40:46 willic3 Exp $

# End of file 
