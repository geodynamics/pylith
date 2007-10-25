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

  dt = 0.25
  density = 2500.0
  vp = 5000.0
  vs = 3000.0
  edgeLen = 2.0**0.5
  N0 = 0.5
  N1 = 0.5
  velocityX = 0.5*(N0+N1)*(N0*(1.5-1.1)+N1*(2.1-1.3)) / (2.0*dt)
  velocityY = 0.5*(N0+N1)*(N0*(2.4-1.8)+N1*(2.4-2.2)) / (2.0*dt)
  normal = [0.5**0.5, -0.5**0.5]

  constNormal = density*vp
  constTangential = density*vs
  dampingConsts = [abs(constNormal*normal[0] - constTangential*normal[1]),
                   abs(constNormal*normal[1] + constTangential*normal[0])]
  residualX = -dampingConsts[0]*velocityX*edgeLen
  residualY = -dampingConsts[1]*velocityY*edgeLen
  residual = [residualX, residualY,
              residualX, residualY]

  j00 = edgeLen*N0**2 / (2.0*dt)
  j01 = edgeLen*N0*N1 / (2.0*dt)
  j10 = j01
  j11 = edgeLen*N1**2 / (2.0*dt)
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
def calcQuad4():
  """
  Calculate damping constants, residual, and Jacobian values for
  absorbing dampers on mesh with quad4 cells.
  """

  density = 2500.0
  vp = 5000.0
  vs = 3000.0
  dt = 0.25
  N0 = 0.5
  N1 = 0.5
  velocity1X = 0.5*(N0+N1)*(N0*(1.5-1.1) + N1*(1.2-1.0)) / (2.0*dt)
  velocity1Y = 0.5*(N0+N1)*(N0*(2.4-1.8) + N1*(1.6-2.4)) / (2.0*dt)
  normal1 = [-1.0, 0.0]

  velocity2X = (N0*(2.4-1.4) + N1*(2.7-1.5)) / (2.0*dt)
  velocity2Y = (N0*(2.4-2.4) + N1*(3.4-1.6)) / (2.0*dt)
  normal2 = [1.0, 0.0]

  edgeLen = 2.0

  constNormal = density*vp
  constTangential = density*vs
  dampingConsts = [abs(constNormal*normal1[0] - constTangential*normal1[1]),
                   abs(constNormal*normal1[1] + constTangential*normal1[0]),
                   abs(constNormal*normal2[0] - constTangential*normal2[1]),
                   abs(constNormal*normal2[1] + constTangential*normal2[0])]
  residual1X = -dampingConsts[0]*velocity1X*edgeLen
  residual1Y = -dampingConsts[1]*velocity1Y*edgeLen
  residual2X = -dampingConsts[2]*velocity2X*edgeLen
  residual2Y = -dampingConsts[3]*velocity2Y*edgeLen
  residual = [residual1X, residual1Y,
              residual1X, residual1Y,
              residual2X, residual2Y,
              residual2X, residual2Y]

  j00 = edgeLen*N0**2 / (2.0*dt)
  j01 = edgeLen*N0*N1 / (2.0*dt)
  j10 = j01
  j11 = edgeLen*N1**2 / (2.0*dt)
  jacobian = [dampingConsts[0]*j00, dampingConsts[1]*j00,
              dampingConsts[0]*j01, dampingConsts[1]*j01,
              dampingConsts[0]*j10, dampingConsts[1]*j10,
              dampingConsts[0]*j11, dampingConsts[1]*j11,
              dampingConsts[2]*j00, dampingConsts[3]*j00,
              dampingConsts[2]*j01, dampingConsts[3]*j01,
              dampingConsts[2]*j10, dampingConsts[3]*j10,
              dampingConsts[2]*j11, dampingConsts[3]*j11]
  print "Absorbing boundary for quad4mesh"
  print "damping constants:"
  for v in dampingConsts:
      print "  %16.8e" % v
  print "velocity:"
  print "  velocity1: ",velocity1X,"  ",velocity1Y
  print "  velocity2: ",velocity2X,"  ",velocity2Y
  print "values for residual:"
  for v in residual:
      print "  %16.8e" % v
  print "values for jacobian:"
  for j in jacobian:
      print "  %16.8e" % j

  
# ----------------------------------------------------------------------
calcTri3()
calcQuad4()

  
# End of file 
