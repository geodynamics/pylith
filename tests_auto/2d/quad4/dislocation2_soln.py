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

## @file tests/2d/quad4/dislocation2_soln.py
##
## @brief Analytical solution to dislocation2 (two faults, rigid motion).

import numpy

# Physical properties
p_density = 2500.0
p_vs = 3000.0
p_vp = 5291.502622129181

p_mu = p_density*p_vs**2
p_lambda = p_density*p_vp**2 - 2*p_mu

# Uniform stress field (plane strain)
sxx = 0.0
sxy = 0.0
syy = 0.0
szz = p_lambda/(2*p_lambda+2*p_mu)*(sxx+syy)

# Uniform strain field
exx = 1.0/(2*p_mu) * (sxx - p_lambda/(3*p_lambda+2*p_mu) * (sxx+syy+szz))
eyy = 1.0/(2*p_mu) * (syy - p_lambda/(3*p_lambda+2*p_mu) * (sxx+syy+szz))
ezz = 1.0/(2*p_mu) * (szz - p_lambda/(3*p_lambda+2*p_mu) * (sxx+syy+szz))

exy = 1.0/(2*p_mu) * (sxy)

#print exx,eyy,exy,ezz
#print -exx*p_lambda/(p_lambda+2*p_mu)

# ----------------------------------------------------------------------
class AnalyticalSoln(object):
  """
  Analytical solution to axial/shear displacement problem.
  """

  def __init__(self):
    return


  def displacement(self, locs, nlocsO):
    """
    Compute displacement field at locations.
    """
    (nlocs, dim) = locs.shape

    disp = numpy.zeros( (nlocs, 3), dtype=numpy.float64)
    maskP = locs[:,0] >= 0.0
    maskP[nlocsO:nlocs] = False
    maskN = numpy.bitwise_and(locs[:,0] <= 0.0, ~maskP)
    disp[:,1] = \
        maskN*(+1.0) + \
        maskP*(-1.0)
    return disp


  def strain(self, locs):
    """
    Compute strain field at locations.
    """
    (nlocs, dim) = locs.shape
    strain = numpy.zeros( (nlocs, 3), dtype=numpy.float64)
    strain[:,0] = exx
    strain[:,1] = eyy
    strain[:,2] = exy
    return strain
  

  def stress(self, locs):
    """
    Compute stress field at locations.
    """
    (nlocs, dim) = locs.shape
    stress = numpy.zeros( (nlocs, 3), dtype=numpy.float64)
    stress[:,0] = sxx
    stress[:,1] = syy
    stress[:,2] = sxy
    return stress


# End of file 
