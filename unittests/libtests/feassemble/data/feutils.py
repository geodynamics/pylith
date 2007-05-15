#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# <LicenseText>
#
# ----------------------------------------------------------------------
#

## @file unittests/libtests/feassemble/data/feutils.py

## @brief Python routines for doing simple finite-element related
## calculations.

import numpy

# ----------------------------------------------------------------------
def calculateJacobian(quadrature, vertices):
  """
  Calculate jacobian, its determinant, and its inverse at quadrature
  points for a given cell.

  @param quadrature Quadrature information
  @param vertices Coordinates of cell's vertices
  """
  jacobian = numpy.zeros( (quadrature.numQuadPts,
                           quadrature.spaceDim, quadrature.cellDim),
                          dtype=numpy.float64)
  jacobianInv = numpy.zeros( (quadrature.numQuadPts,
                              quadrature.spaceDim, quadrature.cellDim),
                             dtype=numpy.float64)
  jacobianDet = numpy.zeros( (quadrature.numQuadPts,), dtype=numpy.float64)
    
  iQuad = 0
  for q in quadrature.quadPtsRef:
    # Jacobian at quadrature points
    deriv = quadrature.basisDeriv[iQuad]
    j = numpy.dot(vertices.transpose(), deriv)
    jacobian[iQuad] = j

    # Determinant of Jacobian and Jacobian inverse at quadrature points
    if quadrature.spaceDim == quadrature.cellDim:
      jacobianDet[iQuad] = numpy.linalg.det(j)
      jacobianInv[iQuad] = numpy.linalg.inv(j)
    else:
      det = numpy.linalg.det(numpy.dot(j.transpose(), j))**0.5
      jacobianDet[iQuad] = det

      if 1 == quadrature.cellDim:
        jacobianInv[iQuad] = 1.0 / j
      elif 2 == quadrature.cellDim:
        minJacobian = 1.0e-06
        jj10 = j[[0,1],:]
        jj21 = j[[1,2],:]
        jj20 = j[[0,2],:]
        det10 = numpy.linalg.det(jj10)
        det21 = numpy.linalg.det(jj21)
        det20 = numpy.linalg.det(jj20)
        if abs(det10) > minJacobian:
          ij10 = numpy.linalg.inv(jj10)
          if abs(det21) > minJacobian:
            ij21 = numpy.linalg.inv(jj21)
            jacobianInv[iQuad] = numpy.array([ [ij10[0,0], ij10[1,0]],
                                               [ij10[0,1], ij10[1,1]],
                                               [ij21[0,1], ij21[1,1]] ],
                                             dtype=numpy.float64)
          elif abs(det20) > minJacobian:
            ij20 = numpy.linalg.inv(jj20)
            jacobianInv[iQuad] = numpy.array([ [ij10[0,0], ij10[1,0]],
                                               [ij10[0,1], ij10[1,1]],
                                               [ij20[0,1], ij20[1,1]] ],
                                             dtype=numpy.float64)
          else:
            jacobianInv[iQuad] = numpy.array([ [ij10[0,0], ij10[0,1]],
                                               [ij10[1,0], ij10[1,1]],
                                               [      0.0,       0.0] ],
                                             dtype=numpy.float64)
        elif abs(det20) > minJacobian:
          ij20 = numpy.linalg.inv(jj20)
          if abs(det21) > minJacobian:
            ij21 = numpy.linalg.inv(jj21)
            jacobianInv[iQuad] = numpy.array([ [ij20[0,0], ij20[1,0]],
                                               [ij21[0,0], ij21[1,0]],
                                               [ij20[0,1], ij20[1,1]] ],
                                             dtype=numpy.float64)
          else:
            jacobianInv[iQuad] = numpy.array([ [ij20[0,0], ij20[1,0]],
                                               [      0.0,       0.0],
                                               [ij20[0,1], ij20[1,1]] ],
                                             dtype=numpy.float64)
        elif abs(det21) > minJacobian:
          ij21 = numpy.linalg.inv(jj21)
          jacobianInv[iQuad] = numpy.array([ [      0.0,       0.0],
                                             [ij21[0,0], ij21[1,0]],
                                             [ij21[0,1], ij21[1,1]] ],
                                           dtype=numpy.float64)
        else:
          raise ValueError("Could not find inverse of Jacobian.")
      else:
        raise ValueError("Could not find inverse of Jacobian.")
    iQuad += 1
  return (jacobian, jacobianInv, jacobianDet)
    

# ----------------------------------------------------------------------
def assembleMat(globalMat, cellMat, cell, fiberDim):
  """
  Assemble cell matrix into global matrix.
  """
  (nrows, ncols) = cellMat.shape
  for iR in xrange(nrows/fiberDim):
    ibeginL = iR * fiberDim
    iendL = ibeginL + fiberDim
    ibeginG = cell[iR] * fiberDim
    iendG = ibeginG + fiberDim
    for iC in xrange(ncols/fiberDim):
      jbeginL = iC * fiberDim
      jendL = jbeginL + fiberDim
      jbeginG = cell[iC] * fiberDim
      jendG = jbeginG + fiberDim
      globalMat[ibeginG:iendG, jbeginG:jendG] += cellMat[ibeginL:iendL,
                                                         jbeginL:jendL]
  return

  
# ----------------------------------------------------------------------
def assembleVec(globalVec, cellVec, cell, fiberDim):
  """
  Assemble cell vector into global vector.
  """
  (nrows,) = cellVec.shape
  for iR in xrange(nrows/fiberDim):
    ibeginL = iR * fiberDim
    iendL = ibeginL + fiberDim
    ibeginG = cell[iR] * fiberDim
    iendG = ibeginG + fiberDim
    globalVec[ibeginG:iendG] += cellVec[ibeginL:iendL]
  return

  
# End of file 
