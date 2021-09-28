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

## @file tests/fullscale/linearelasticity/cylinder-2d/wedgetraction_soln.py

## @brief Analytical solution for elastic cylinder with pressure on outer boundary.

# Neumann boundary conditions.
# boundary_outer
#   Tr(r=b) = -100 MPa
#
# Dirichlet boundary conditions.
#   Ux(0,y) = 0
#   Uy(x,0) = 0


import math
import numpy

# ----------------------------------------------------------------------
# Solution parameters.
b = 20000.0 # Outer cylinder radius.
Pb = 1.0e8 # Pressure applied at r=b.

# PyLith solution parameters.
p_density = 2500.0
p_vs = 3464.1016
p_vp = 6000.0

p_mu = p_density * p_vs**2
p_lambda = p_density * p_vp**2 - 2 * p_mu

p_young = p_mu*(3.0*p_lambda + 2.0*p_mu)/(p_lambda + p_mu)
p_poisson = 0.5*p_lambda/(p_lambda + p_mu)

C = -0.5*Pb

# Uniform stress field.
sxx = 2.0*C
syy = 2.0*C
szz = p_poisson*2.0*sxx
sxy = 0.0

# Uniform strain field.
exx = 2.0*C*(1.0 + p_poisson)*(1.0 - 2.0*p_poisson)/p_young
eyy = exx
ezz = 0.0
exy = 0.0

# ----------------------------------------------------------------------
class AnalyticalSoln(object):
    """Analytical solution to pressurized cylinder problem.
    """
    SPACE_DIM = 2
    TENSOR_SIZE = 4

    def __init__(self):
        self.fields = {
            "displacement": self.displacement,
            "density": self.density,
            "shear_modulus": self.shear_modulus,
            "bulk_modulus": self.bulk_modulus,
            "cauchy_strain": self.strain,
            "cauchy_stress": self.stress,
            "initial_amplitude": {
                "bc_xneg": self.zero_vector,
                "bc_yneg": self.zero_vector,
                "bc_outer": self.bc_outer_traction
                }
        }
        return

    def getField(self, name, mesh_entity, pts):
        if name == "initial_amplitude":
            field = self.fields[name][mesh_entity](pts)
        else:
            field = self.fields[name](pts)
        return field
    
    def displacement(self, locs):
        """
        Compute analytical displacement solution for a given set of coordinates.
        """
        (npts, dim) = locs.shape
        r = numpy.linalg.norm(locs[:,0:2], axis=1)
        ur = (1.0 + p_poisson)*(2.0*(1.0 - 2.0*p_poisson)*C*r)/p_young

        disp = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        
        (disp[0, :, 0], disp[0, :, 1]) = self.cylToCartDisp(ur, locs)

        return disp

    def zero_vector(self, locs):
        (npts, dim) = locs.shape
        return numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)

    def density(self, locs):
        """Compute density field at locations.
        """
        (npts, dim) = locs.shape
        density = p_density * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return density

    def shear_modulus(self, locs):
        """Compute shear modulus field at locations.
        """
        (npts, dim) = locs.shape
        shear_modulus = p_mu * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return shear_modulus

    def bulk_modulus(self, locs):
        """Compute bulk modulus field at locations.
        """
        (npts, dim) = locs.shape
        bulk_modulus = (p_lambda + 2.0 / 3.0 * p_mu) * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return bulk_modulus

    def strain(self, locs):
        """Compute strain field at locations.
        """
        (npts, dim) = locs.shape
        strain = numpy.zeros((1, npts, self.TENSOR_SIZE), dtype=numpy.float64)
        strain[0,:, 0] = exx
        strain[0,:, 1] = eyy
        strain[0,:, 2] = ezz
        strain[0,:, 3] = exy
        return strain

    def stress(self, locs):
        """
        Compute analytical stress solution for a given set of coordinates.
        """
        (npts, dim) = locs.shape
        stress = numpy.zeros((1, npts, self.TENSOR_SIZE), dtype=numpy.float64)
        stress[0,:, 0] = sxx
        stress[0,:, 1] = syy
        stress[0,:, 2] = szz
        stress[0,:, 3] = sxy

        return (stress)

    def cylToCartDisp(self, ur, locs):
        """
        Convert displacement increments in cylindrical coordinates to Cartesian.
        Assumption is that tangential displacements are zero.
        """
        angs = numpy.arctan2(locs[:,1], locs[:,0])
        ca = numpy.cos(angs)
        sa = numpy.sin(angs)

        ux = ur*ca
        uy = ur*sa

        return (ux, uy)

    def bc_outer_traction(self, locs):
        """Compute initial traction at locations.
        """
        (npts, dim) = locs.shape
        traction = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        traction[0,:, 0] = 0.0
        traction[0,:, 1] = -Pb
        return traction


# End of file 
