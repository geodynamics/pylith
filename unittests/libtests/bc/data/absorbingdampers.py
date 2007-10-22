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

## @file unittests/libtests/bc/data/absorbingdampers.py

## @brief Python routines for calculating values for test of absorbing
## dampers.

import math

# ----------------------------------------------------------------------
def calcTri3():
  """
  Calculate damping constants, residual, and Jacobian values for
  absorbing dampers on mesh with tri3 cells.
  """

  density = 2500.0
  vp = 5000.0
  vs = 3000.0
  velocityX = 0.3
  velocityY = 0.2
  edgeLen = 0.5**0.5
  dt = 0.25
  normal = [0.5**0.5, -0.5**0.5]
  N0 = 0.5
  N1 = 0.5

  constNormal = density*vp
  constTangential = density*vs
  dampingConsts = [abs(constNormal*normal[0] - constTangential*normal[1]),
                   abs(constNormal*normal[1] + constTangential*normal[0])]
  residualX = -dampingConsts[0]*velocityX*edgeLen/(2.0*dt)
  residualY = -dampingConsts[1]*velocityY*edgeLen/(2.0*dt)
  residual = [residualX, residualY,
              residualX, residualY]

  j00 = 2*edgeLen*N0**2 / (2.0*dt)
  j01 = 2*edgeLen*N0*N1 / (2.0*dt)
  j10 = j01
  j11 = 2*edgeLen*N1**2 / (2.0*dt)
  jacobian = [dampingConsts[0]*j00, dampingConsts[1]*j00,
              dampingConsts[0]*j01, dampingConsts[1]*j01,
              dampingConsts[0]*j10, dampingConsts[1]*j10,
              dampingConsts[0]*j11, dampingConsts[1]*j11]

  print "Absorbing boundary for tri3 mesh"
  print "damping constants:"
  for v in dampingConsts:
      print "  %16.8e" % v
  print "values for residual:"
  for v in residual:
      print "  %16.8e" % v
  print "values for jacobian:"
  for j in jacobian:
      print "  %16.8e" % j

  
# ----------------------------------------------------------------------
calcTri3()

  
# End of file 
