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
# Copyright (c) 2010-2017 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file tests/2d/quad4/compressrotate_soln.py
##
## @brief Analytical solution to compress and rigid body rotation.

## Axial compression in the x-direction + 90 degree CCW rigid body rotation.
##
## Compression
## Ux1 =  ex*(x-x0)
## Uy1 = ey*(y-y0)
## x1 = x + Ux1
## y1 = y + Uy1
##
## Ux2 = xr + (x1-xr)*cos(theta) + (y1-yr)*sin(theta) - x
## Uy2 = yr - (x1-xr)*sin(theta) + (y1-yr)*cos(theta) - y
##
## Dirichlet boundaries on -x, +x, and -y edges.

import numpy

# Physical properties
p_density = 2500.0
p_vs = 3000.0
p_vp = 5291.502622129181

p_mu = p_density*p_vs**2
p_lambda = p_density*p_vp**2 - 2*p_mu

# Uniform strain field
ex = -0.04
exx = ex + 0.5*ex*ex
eyy = -p_lambda/(p_lambda+2*p_mu)*exx
ey = -1+(1.0+2.0*eyy)**0.5

ezz = 0.0
exy = 0.0

# Uniform stress field (plane strain) in undeformed (original) configuration
sxx = p_lambda*(exx+eyy+ezz) + 2.0*p_mu*exx
syy = p_lambda*(exx+eyy+ezz) + 2.0*p_mu*eyy
szz = p_lambda*(exx+eyy+ezz) + 2.0*p_mu*ezz
sxy = 2.0*p_mu*exy

theta_d = -2.0
x0 = 0.0
y0 = -500.0
xr = -1000.0
yr = 0.0

#print ex,ey
#print exx,eyy,exy,ezz
#print sxx, syy, sxy

# ----------------------------------------------------------------------
class AnalyticalSoln(object):
  """
  Analytical solution to axial/shear displacement problem.
  """

  def __init__(self):
    from math import pi,cos,sin
    theta = theta_d / 180.0*pi
    self.R = numpy.array([[cos(theta), sin(theta)], [-sin(theta), cos(theta)]])
    self.U = numpy.array([[1.0+ex, 0], [0, 1.0+ey]])
    return

  def _transform(self, m1):
    import numpy.linalg
    X = numpy.dot(self.R, self.U)
    detX = numpy.linalg.det(X)
    m2 = numpy.dot(numpy.dot(X, m1), X.transpose()) / detX
    return m2


  def displacement(self, locs):
    """
    Compute displacement field at locations.
    """
    (npts, dim) = locs.shape
    disp = numpy.zeros( (1, npts, 2), dtype=numpy.float64)
    from math import pi,cos,sin
    x = locs[:,0]
    y = locs[:,1]
    ux1 = ex*(x-x0)
    uy1 = ey*(y-y0)
    x1 = x + ux1
    y1 = y + uy1
    theta = theta_d / 180.0*pi
    ux2 = xr + (x1-xr)*cos(theta) + (y1-yr)*sin(theta) - x
    uy2 = yr - (x1-xr)*sin(theta) + (y1-yr)*cos(theta) - y

    disp[0,:,0] = ux2
    disp[0,:,1] = uy2
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


  def cauchy_stress(self, locs):
    """
    Compute stress field at locations.
    """
    (npts, dim) = locs.shape

    stress1 = numpy.array([[sxx, sxy], [sxy, syy]])
    stress2 = self._transform(stress1)

    stress = numpy.zeros( (1, npts, 3), dtype=numpy.float64)
    stress[0,:,0] = stress2[0,0]
    stress[0,:,1] = stress2[1,1]
    stress[0,:,2] = stress2[0,1]
    return stress


# End of file 
