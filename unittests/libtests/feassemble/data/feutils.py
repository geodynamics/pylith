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
  (basis, basisDerivRef) = quadrature.calculateBasis()

  numQuadPts = quadrature.numQuadPts
  cellDim = quadrature.cellDim
  spaceDim = quadrature.spaceDim
  numBasis = quadrature.numBasis
  jacobian = numpy.zeros( (numQuadPts, spaceDim, cellDim),
                          dtype=numpy.float64)
  jacobianInv = numpy.zeros( (numQuadPts, cellDim, spaceDim),
                             dtype=numpy.float64)
  jacobianDet = numpy.zeros( (numQuadPts,), dtype=numpy.float64)
  basisDeriv = numpy.zeros( (numQuadPts, numBasis, spaceDim),
                            dtype=numpy.float64)

  iQuad = 0
  for q in quadrature.quadPtsRef:
    # Jacobian at quadrature points
    deriv = basisDerivRef[iQuad]
    j = numpy.dot(vertices.transpose(), deriv)
    jacobian[iQuad] = j

    # Determinant of Jacobian and Jacobian inverse at quadrature points
    if spaceDim == cellDim:
      jacobianDet[iQuad] = numpy.linalg.det(j)
      jacobianInv[iQuad] = numpy.linalg.inv(j)
    else:
      det = numpy.linalg.det(numpy.dot(j.transpose(), j))**0.5
      jacobianDet[iQuad] = det

      if 1 == quadrature.cellDim:
        jacobianInv[iQuad] = 1.0 / j.transpose()
      elif 2 == quadrature.cellDim:
        minJacobian = 1.0e-06
        jj01 = j[[0,1],:]
        jj12 = j[[1,2],:]
        jj02 = j[[0,2],:]
        det01 = numpy.linalg.det(jj01)
        det12 = numpy.linalg.det(jj12)
        det02 = numpy.linalg.det(jj02)
        if abs(det01) > minJacobian:
          ij01 = numpy.linalg.inv(jj01)
          if abs(det12) > minJacobian:
            ij12 = numpy.linalg.inv(jj12)
            jacobianInv[iQuad] = numpy.array([ [ij01[0,0], ij01[0,1], ij12[0,1]],
                                               [ij01[1,0], ij01[1,1], ij12[1,1]] ],
                                             dtype=numpy.float64)
          elif abs(det02) > minJacobian:
            ij02 = numpy.linalg.inv(jj02)
            jacobianInv[iQuad] = numpy.array([ [ij01[0,0], ij01[0,0], ij02[0,1]],
                                               [ij01[1,0], ij11[1,1], ij02[1,1]] ],
                                             dtype=numpy.float64)
          else:
            jacobianInv[iQuad] = numpy.array([ [ij01[0,0], ij01[0,1], 0.0],
                                               [ij01[1,0], ij01[1,1], 0.0] ],
                                             dtype=numpy.float64)
        elif abs(det02) > minJacobian:
          ij02 = numpy.linalg.inv(jj02)
          if abs(det12) > minJacobian:
            ij12 = numpy.linalg.inv(jj12)
            jacobianInv[iQuad] = numpy.array([ [ij02[0,0], ij12[0,0], ij02[0,1]],
                                               [ij02[1,0], ij12[1,0], ij02[1,1]] ],
                                             dtype=numpy.float64)
          else:
            jacobianInv[iQuad] = numpy.array([ [ij02[0,0], 0.0, ij02[0,1]],
                                               [ij02[1,0], 0.0, ij02[1,1]] ],
                                             dtype=numpy.float64)
        elif abs(det12) > minJacobian:
          ij12 = numpy.linalg.inv(jj12)
          jacobianInv[iQuad] = numpy.array([ [0.0, ij12[0,0], ij12[0,1]],
                                             [0.0, ij12[1,0], ij12[1,1]] ],
                                           dtype=numpy.float64)
        else:
          raise ValueError("Could not find inverse of Jacobian.")
      else:
        raise ValueError("Could not find inverse of Jacobian.")

    # Compute derivatives of basis functions with respect to global
    # coordinates using derivatives of basis functions with respect
    # to local cell coordinates and inverse of Jacobian matrix
    for iBasis in xrange(numBasis):
      matBasis = numpy.array(basisDerivRef[iQuad,iBasis], dtype=numpy.float64)
      jInv = jacobianInv[iQuad]
      basisDeriv[iQuad,iBasis,:] = numpy.dot(matBasis, jInv)
    iQuad += 1
  return (jacobian, jacobianInv, jacobianDet, basisDeriv)
    

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
