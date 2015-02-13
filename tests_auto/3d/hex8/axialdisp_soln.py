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

## @file tests/2d/tri3/axialdisp_soln.py
##
## @brief Analytical solution to axial displacement problem.

import numpy

# Physical properties
p_density = 2500.0
p_vs = 3000.0
p_vp = 5291.502622129181

p_mu = p_density*p_vs**2
p_lambda = p_density*p_vp**2 - 2*p_mu

# Uniform stress field (plane stress)
sxx = 1.0e+7
syy = 0.0
szz = 0.0
sxy = 0.0
syz = 0.0
sxz = 0.0

# Uniform strain field
exx = 1.0/(2*p_mu) * (sxx - p_lambda/(3*p_lambda+2*p_mu) * (sxx+syy+szz))
eyy = 1.0/(2*p_mu) * (syy - p_lambda/(3*p_lambda+2*p_mu) * (sxx+syy+szz))
ezz = 1.0/(2*p_mu) * (szz - p_lambda/(3*p_lambda+2*p_mu) * (sxx+syy+szz))

exy = 1.0/(2*p_mu) * (sxy)
eyz = 1.0/(2*p_mu) * (syz)
exz = 1.0/(2*p_mu) * (sxz)

#print exx,eyy,ezz,exy,eyz,exz

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
    (npts, dim) = locs.shape
    disp = numpy.zeros( (1, npts, 3), dtype=numpy.float64)
    disp[0,:,0] = exx*locs[:,0] + exy*locs[:,1] + exz*locs[:,2]
    disp[0,:,1] = eyy*locs[:,1] + exy*locs[:,0] + eyz*locs[:,2]
    disp[0,:,2] = ezz*locs[:,2] + exz*locs[:,0] + eyz*locs[:,1]
    return disp


  def strain(self, locs):
    """
    Compute strain field at locations.
    """
    (npts, dim) = locs.shape
    strain = numpy.zeros( (1, npts, 6), dtype=numpy.float64)
    strain[0,:,0] = exx
    strain[0,:,1] = eyy
    strain[0,:,2] = ezz
    strain[0,:,3] = exy
    strain[0,:,4] = eyz
    strain[0,:,5] = exz
    return strain
  

  def stress(self, locs):
    """
    Compute stress field at locations.
    """
    (npts, dim) = locs.shape
    stress = numpy.zeros( (1, npts, 6), dtype=numpy.float64)
    stress[0,:,0] = sxx
    stress[0,:,1] = syy
    stress[0,:,2] = szz
    stress[0,:,3] = sxy
    stress[0,:,4] = syz
    stress[0,:,5] = sxz
    return stress


# End of file 
