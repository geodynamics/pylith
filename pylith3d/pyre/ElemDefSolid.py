#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                              Charles A. Williams
#                        Rensselaer Polytechnic Institute
#                      (C) 2006  All Rights Reserved
#
#  Copyright 2006 Rensselaer Polytechnic Institute.
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

# This could probably be set up as a Pyre component, but it doesn't seem
# necessary right now.

class ElemIntegrationSolid:
  """Integration information for solid elements used by pylith3d."""

  def getintegration(self, elemType, quadratureOrder):
      """Gets dimension and quadrature info, depending on the element type
      and the global quadrature order setting."""

    import pylith3d as pl3d

    print "Hello from ElemIntegrationSolid.getintegration!"

    numElemNodes = self.elemNodes[elemType - 1]
    self.elemIntegrationSolid['numElemNodes'] = numElemNodes
    self.elemIntegrationSolid['numElemEquations'] = self.numDegreesFreedom * numElemNodes
    self.elemIntegrationSolid['numElemCoords'] = self.numSpaceDims * numElemNodes

    # Full integration order
    if quadratureOrder == "Full":
      numElemGaussPoints = self.elemFullGauss[elementType -1]
      self.elemIntegrationSolid['fptrBmatrix'] = pl3d.bmatrixn
      self.elemIntegrationSolid['fptrGetShape'] = pl3d.getshapn

    # Reduced integration order
    elif quadratureOrder == "Reduced":
      numElemGaussPoints = self.elemReducedGauss[elementType -1]
      self.elemIntegrationSolid['fptrBmatrix'] = pl3d.bmatrixn
      self.elemIntegrationSolid['fptrGetShape'] = pl3d.getshapn

    # Selective (B-bar) integration order
    elif quadratureOrder == "Selective":
      numElemGaussPoints = self.elemFullGauss[elementType -1]
      self.elemIntegrationSolid['fptrBmatrix'] = pl3d.bmatrixb
      self.elemIntegrationSolid['fptrGetShape'] = pl3d.getshapb

    self.elemIntegrationSolid['numElemGaussPoints'] = numElemGaussPoints

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
      quadratureOrder,
      elemType,
      numElemNodes,
      numElemGaussPoints)

    self.elemIntegrationSolid['ptrShape'] = ptrShape
    self.elemIntegrationSolid['ptrShapej'] = ptrShapej
    self.elemIntegrationSolid['ptrGauss'] = ptrGauss
        
    return self.elemIntegrationSolid


    def __init__(self):
      """Constructor."""
      print ""
      print "Hello from ElemIntegrationSolid.__init__!"

      # I am leaving these hard-wired for now, since solid elements will always
      # have 3 DOF and 3 spatial dimensions, but it might be better to leave
      # these as variables for more generality.  Ideally, we would get this
      # info from global definitions.
      self.numDegreesFreedom = 3
      self.numSpaceDims = 3

      self.elemIntegrationSolid = {'numElemNodes': 0,
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
# $Id: ElemIntegrationSolid.py,v 1.2 2005/04/01 23:40:46 willic3 Exp $

# End of file 
