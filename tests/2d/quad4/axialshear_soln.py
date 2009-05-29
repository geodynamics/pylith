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

## @file tests/2d/quad4/axialsheardisp_soln.py
##
## @brief Analytical solution to axial/shear displacement problem.

import numpy

# Physical properties
p_density = 2500.0
p_vs = 3000.0
p_vp = 5291.502622129181

p_mu = p_density*p_vs**2
p_lambda = p_density*p_vp**2 - 2*p_mu

# Uniform stress field
sxx = 1.0e+7
sxy = 8.0e+6
syy = 0.0

# Uniform strain field
exx = 1.0/(2*p_mu) * (sxx - p_lambda/(3*p_lambda+2*p_mu) * (sxx+syy))
eyy = 1.0/(2*p_mu) * (syy - p_lambda/(3*p_lambda+2*p_mu) * (sxx+syy))
exy = 1.0/(2*p_mu) * (sxy)
      
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
    disp = numpy.zeros( (npts, 2), dtype=numpy.float64)
    disp[:,0] = exx*locs[:,0] + exy*locs[:,1]
    disp[:,1] = eyy*locs[:,1] + exy*locs[:,0]
    return disp


  def strain(self, locs):
    """
    Compute strain field at locations.
    """
    (npts, dim) = locs.shape
    strain = numpy.zeros( (npts, 3), dtype=numpy.float64)
    strain[:,0] = exx
    strain[:,1] = eyy
    strain[:,2] = exy
    return
  

  def stress(self, locs):
    """
    Compute stress field at locations.
    """
    (npts, dim) = locs.shape
    stress = numpy.zeros( (npts, 3), dtype=numpy.float64)
    stress[:,0] = self.sxx
    stress[:,1] = self.syy
    stress[:,2] = self.sxy
    return


# End of file 
