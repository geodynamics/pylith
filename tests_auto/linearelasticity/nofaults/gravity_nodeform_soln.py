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
# Copyright (c) 2010-2018 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file tests_auto/linearelasticity/nofaultsgravity_nodeform_soln.py
##
## @brief Analytical solution to gravity with initial stress and no displacement.

## 2-D uniform shear test.
##
##             --->
##          ----------
##          |        |
##        | |        | ^
##        v |        | |
##          |        |
##          ----------
##             <--
## 
## Dirichlet boundary conditions
##   Ux(-4000,y) = 0
##   Ux(+4000,y) = 0
##   Ux(x,-4000) = 0

import numpy

# Physical properties
p_density = 2500.0 # kg/m**3
p_vs = 3000.0 # m/s
p_vp = 5291.502622129181 # m/s

p_mu = p_density*p_vs**2
p_lambda = p_density*p_vp**2 - 2*p_mu

gacc = 9.80665 # m/s
ytop = +4000.0 # m

# Uniform strain field
exx = 0
eyy = 0
ezz = 0
exy = 0

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
    (npts, dim) = locs.shape
    disp = numpy.zeros( (1, npts, 2), dtype=numpy.float64)
    disp[0,:,0] = exx*locs[:,0] + exy*locs[:,1]
    disp[0,:,1] = eyy*locs[:,1] + exy*locs[:,0]
    return disp


  def strain(self, locs):
    """
    Compute strain field at locations.
    """
    (npts, dim) = locs.shape
    strain = numpy.zeros( (1, npts, 4), dtype=numpy.float64)
    strain[0,:,0] = exx
    strain[0,:,1] = eyy
    strain[0,:,2] = ezz
    strain[0,:,3] = exy
    return strain
  

  def stress(self, locs):
    """
    Compute stress field at locations.
    """
    syy = -p_density * gacc * (ytop-locs[:,1])
    sxx = syy
    szz = syy
    sxy = 0
    
    (npts, dim) = locs.shape
    stress = numpy.zeros( (1, npts, 4), dtype=numpy.float64)
    stress[0,:,0] = sxx
    stress[0,:,1] = syy
    stress[0,:,2] = szz
    stress[0,:,3] = sxy
    return stress


  def matprops(self, locs):
      """Get material properties at locations.
      """
      (npts, dim) = locs.shape
      
      props = {
          "density": p_density * numpy.ones(npts),
          "vp": p_vp * numpy.ones(npts),
          "vs": p_vs * numpy.ones(npts),
          }
      return props

  
# End of file 
