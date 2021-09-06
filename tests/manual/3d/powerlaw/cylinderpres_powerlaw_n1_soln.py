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

## @file tests/manual/3d/powerlaw/cylinderpres_soln.py

## @brief Analytical steady-state solution for pressurized power-law viscoelastic cylinder.

# Neumann boundary conditions.
# Note that in the analytical solution, positive boundary tractions are considered to be compressive.
# boundary_inner
#   Tr(r=a) = 0 MPa
# boundary_outer
#   Tr(r=b) = -100 MPa
#
# Dirichlet boundary conditions.
#   Ux(0,y,z) = 0
#   Uy(x,0,z) = 0
#   Uz(x,y,-1000) = 0
#   Uz(x,y,1000) = 0


import math
import numpy

# ----------------------------------------------------------------------
# Solution parameters.
# NOTE:  A negative difference (Pb - Pa) will yield imaginary velocity values for fractional exponents.
a = 2000.0 # Inner cylinder radius.
b = 20000.0 # Outer cylinder radius.
Pa = 0.0 # Pressure applied at r=a.
Pb = 1.0e8 # Pressure applied at r=b.

# PyLith solution parameters.
p_density = 2500.0
p_vs = 3464.1016
p_vp = 6000.0
p_powerLawExponent = 1.0
p_powerLawReferenceStrainRate = 1.0e-6
p_powerLawReferenceStress = 5680367950354.795

p_mu = p_density * p_vs**2
p_lambda = p_density * p_vp**2 - 2 * p_mu

# Parameters converted to FLAC/analytical form.
AT = p_powerLawReferenceStrainRate/(p_powerLawReferenceStress**p_powerLawExponent)
A = 2.0*AT/(math.sqrt(3.0)**(3.0 - p_powerLawExponent))

# ----------------------------------------------------------------------
class AnalyticalSoln(object):
    """Analytical solution to pressurized power-law cylinder problem.
    """
    SPACE_DIM = 3
    TENSOR_SIZE = 6

    def __init__(self):
        self.fields = {
            "velocity": self.velocity,
            "density": self.density,
            "shear_modulus": self.shear_modulus,
            "bulk_modulus": self.bulk_modulus,
            "cauchy_stress": self.stress,
        }
        return
    
    def velocity(self, locs):
        """
        Compute analytical velocity solution for a given set of coordinates.
        """
        (npts, dim) = locs.shape
        r = numpy.linalg.norm(locs[:,0:2], axis=1)
        k1 = 2.0/p_powerLawExponent
        k3 = (3.0/4.0)**((p_powerLawExponent + 1.0)/2.0)
        urDot = -A*k3*((Pb - Pa)*(k1/((b/a)**k1 - 1.0)))**p_powerLawExponent*(b**2.0/r)

        vel = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        
        (vel[0, :, 0], vel[0, :, 1]) = self.cylToCartDisp(urDot, locs)

        return vel

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

    def stress(self, locs):
        """
        Compute analytical stress solution for a given set of coordinates.
        """
        npts = locs.shape[0]
        r = numpy.linalg.norm(locs[:,0:2], axis=1)
        k1 = 2.0/p_powerLawExponent
        k2 = 1.0/p_powerLawExponent
        k3 = (3.0/4.0)**((p_powerLawExponent + 1.0)/2.0)
        srr = -Pb + (Pb - Pa)*((b/r)**k1 - 1.0)/((b/a)**k1 - 1.0)
        stt = -Pb - (Pb - Pa)*((k1 - 1.0)*(b/r)**k1 + 1.0)/((b/a)**k1 - 1.0)
        szz = -Pb - (Pb - Pa)*((k2 - 1.0)*(b/r)**k1 + 1.0)/((b/a)**k1 - 1.0)

        (sxx, syy, sxy) = self.cylToCartStress(srr, stt, locs)
        s0 = numpy.zeros_like(sxx)
        stressVecCart = numpy.column_stack((sxx, syy, szz, sxy, s0, s0)).reshape(1,npts,6)
        stressVecCyl = numpy.column_stack((srr, stt, szz, s0, s0, s0)).reshape(1,npts,6)

        return (stressVecCart, stressVecCyl)

    def cylToCartStress(self, srr, stt, locs):
        """
        Convert stresses in cylindrical coordinates to Cartesian.
        Assumption is that shear stresses are zero.
        """
        angs = numpy.arctan2(locs[:,1], locs[:,0])
        ca = numpy.cos(angs)
        sa = numpy.sin(angs)

        sxx = srr*ca*ca + stt*sa*sa
        syy = srr*sa*sa + stt*ca*ca
        sxy = srr*ca*sa - stt*ca*sa

        return (sxx, syy, sxy)


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


# End of file 
