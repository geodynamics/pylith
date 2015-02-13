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

## @file tests/2d/quad4/rigidbody_soln.py
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
sxy = 0.0
syy = 0.0
szz = p_lambda/(2*p_lambda+2*p_mu)*(sxx+syy)

# Uniform strain field
exx = 0.0
eyy = 0.0
ezz = 0.0

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


  def displacement(self, locs):
    """
    Compute displacement field at locations.
    """ 
    u0 = 200.0
    v0 = 300.0
    from math import pi
    theta = -25.0/180.0*pi

    (npts, dim) = locs.shape
    disp = numpy.zeros( (1, npts, 2), dtype=numpy.float64)
    disp[0,:,0] = -locs[:,0] + \
        u0 + numpy.cos(theta)*locs[:,0] + numpy.sin(theta)*locs[:,1]
    disp[0,:,1] = -locs[:,1] + \
        v0 - numpy.sin(theta)*locs[:,0] + numpy.cos(theta)*locs[:,1]
    return disp


  def strain(self, locs):
    """
    Compute strain field at locations.
    """
    (npts, dim) = locs.shape
    strain = numpy.zeros( (1, npts, 3), dtype=numpy.float64)
    strain[0,:,0] = exx
    strain[0,:,1] = eyy
    strain[0,:,2] = exy
    return strain
  

  def stress(self, locs):
    """
    Compute stress field at locations.
    """
    (npts, dim) = locs.shape
    stress = numpy.zeros( (1, npts, 3), dtype=numpy.float64)
    stress[0,:,0] = sxx
    stress[0,:,1] = syy
    stress[0,:,2] = sxy
    return stress


# End of file 
