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
# Copyright (c) 2010-2014 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file tests/3dnew/hex8/rigidbody_soln.py
##
## @brief Analytical solution to rigid body motion problem.

import numpy

# Physical properties
p_density = 2500.0
p_vs = 3000.0
p_vp = 5291.502622129181

p_mu = p_density*p_vs**2
p_lambda = p_density*p_vp**2 - 2*p_mu

# Uniform stress field (plane strain)
sxx = 0.0
syy = 0.0
szz = 0.0
sxy = 0.0
syz = 0.0
sxz = 0.0

# Uniform strain field
exx = 0.0
eyy = 0.0
ezz = 0.0
exy = 0.0
eyz = 0.0
exz = 0.0

# ----------------------------------------------------------------------
class AnalyticalSoln(object):
  """
  Analytical solution to axial/shear displacement problem.
  """

  def __init__(self):
    return


  def displacement(self, locs):
    """
    Compute displacement field at locations.
    """
    u0 = 400.0
    v0 = -300.0
    w0 = 500.0
    from math import pi
    theta = 60.0/180.0*pi

    (npts, dim) = locs.shape
    disp = numpy.zeros( (npts, 3), dtype=numpy.float64)
    disp[:,0] = u0 -locs[:,0] + \
        u0 + numpy.cos(theta)*locs[:,0] + numpy.sin(theta)*locs[:,2]
    disp[:,1] = v0
    disp[:,2] = -locs[:,2] + \
        w0 - numpy.sin(theta)*locs[:,0] + numpy.cos(theta)*locs[:,2]
    return disp


  def strain(self, locs):
    """
    Compute strain field at locations.
    """
    (npts, dim) = locs.shape
    strain = numpy.zeros( (npts, 6), dtype=numpy.float64)
    strain[:,0] = exx
    strain[:,1] = eyy
    strain[:,2] = ezz
    strain[:,3] = exy
    strain[:,4] = eyz
    strain[:,5] = exz
    return strain
  

  def stress(self, locs):
    """
    Compute stress field at locations.
    """
    (npts, dim) = locs.shape
    stress = numpy.zeros( (npts, 6), dtype=numpy.float64)
    stress[:,0] = sxx
    stress[:,1] = syy
    stress[:,2] = szz
    stress[:,3] = sxy
    stress[:,4] = syz
    stress[:,5] = sxz
    return stress


# End of file 
