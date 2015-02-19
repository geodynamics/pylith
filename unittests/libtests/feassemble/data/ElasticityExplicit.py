#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file unittests/libtests/feassemble/data/ElasticityExplicit.py

## @brief Python application for generating C++ data files for testing
## C++ ElasticityExplicit object.

from pyre.components.Component import Component

import numpy

# ----------------------------------------------------------------------

# ElasticityExplicit class
class ElasticityExplicit(Component):
  """
  Python application for generating C++ data files for testing C++
  ElasticityExplicit object.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="elasticityexplicit"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="formulation")
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def calculateResidual(self, integrator):
    """
    Calculate contribution to residual of operator for integrator.

    {r} = (1/dt**2)[M](-{u(t+dt)} + 2 {u(t)} - {u(t-dt)}) -
          [K]{u(t)}
    """
    K = integrator._calculateStiffnessMat()    
    M = integrator._calculateMassMat()

    vel = (integrator.fieldT + integrator.fieldTIncr - integrator.fieldTmdt) / (2.0*integrator.dt)
    acc = (integrator.fieldTIncr - integrator.fieldT + integrator.fieldTmdt) / (integrator.dt**2) 
    dispAdj = integrator.fieldT + integrator.dt*integrator.normViscosity*vel
    residual = -numpy.dot(M, acc) - numpy.dot(K, dispAdj)
    return residual.flatten()


  def calculateResidual(self, integrator):
    """
    Calculate contribution to residual of operator for integrator.

    {r} = (1/dt**2)[M](-{u(t+dt)} + 2 {u(t)} - {u(t-dt)}) -
          [K]{u(t)}
    """
    K = integrator._calculateStiffnessMat()    
    M = integrator._calculateMassMat()
    Ml = self._lumpMatrix(M, integrator.numBasis, integrator.spaceDim)
    
    vel = (integrator.fieldT + integrator.fieldTIncr - integrator.fieldTmdt) / (2.0*integrator.dt)
    acc = (integrator.fieldTIncr - integrator.fieldT + integrator.fieldTmdt) / (integrator.dt**2) 
    acc = acc.flatten()
    dispAdj = integrator.fieldT + integrator.dt*integrator.normViscosity*vel
    residual = -Ml*acc - numpy.dot(K, dispAdj).flatten()
    return residual.flatten()


  def calculateJacobian(self, integrator):
    """
    Calculate contribution to Jacobian matrix of operator for integrator.

    [A] = (1/dt**2)[M]
    """
    M = integrator._calculateMassMat()
    Ml = self._lumpMatrix(M, integrator.numBasis, integrator.spaceDim)

    jacobian = 1.0/integrator.dt**2 * Ml
    return jacobian.flatten()


  def _lumpMatrix(self, matrix, numBasis, spaceDim):
    """
    Lump matrix.
    """
    (nrows,ncols) = matrix.shape
    assert(numBasis * spaceDim == nrows)
    assert(nrows == ncols)
    vector = numpy.zeros( (ncols), dtype=numpy.float64)

    for iBasis in xrange(numBasis):
      for iDim in xrange(spaceDim):
        v = 0.0
        for jBasis in xrange(numBasis):
          v += matrix[iBasis*spaceDim+iDim,jBasis*spaceDim+iDim]
        vector[iBasis*spaceDim+iDim] = v
    return vector


# FACTORY //////////////////////////////////////////////////////////////
def formulation():
  return ElasticityExplicit()


# End of file 
