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
  velocity = 1.2
  edgeLen = 0.5**0.5
  dt = 0.25
  normal = [0.5**0.5, -0.5**0.5]

  dampingConsts = [density*vp, density*vs]
  residualNormal = dampingConsts[0]*velocity*edgeLen/(2.0*dt)
  residualTangential = dampingConsts[1]*velocity*edgeLen/(2.0*dt)
  residual = [residualNormal*normal[0] - residualTangential*normal[1],
              residualNormal*normal[1] + residualTangential*normal[0],
              residualNormal*normal[0] - residualTangential*normal[1],
              residualNormal*normal[1] + residualTangential*normal[0]]

  jacobian = []

  print "Absorbing boundary for tri3 mesh"
  print "damping constants:"
  for v in dampingConsts:
      print "  %16.8e" % v
  print "residual:"
  for v in residual:
      print "  %16.8e" % v

  
# ----------------------------------------------------------------------
calcTri3()

  
# End of file 
