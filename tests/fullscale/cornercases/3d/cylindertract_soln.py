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
# @file tests/fullscale/cornercases/3d/cylindertract_soln.py
#
# @brief Analytical solution to thick-walled cylinder subjected to internal pressure.
#
# Neumann boundary conditions.
# boundary_inner
#   Tr(r=a) = -100 MPa
#
# Dirichlet boundary conditions.
#   Ux(0,y,z) = 0
#   Uy(x,0,z) = 0
#   Uz(x,y,-1000) = 0
#   Uz(x,y,1000) = 0

import math
import numpy
# import pdb

# Physical properties
p_density = 2500.0
p_vs = 3000.0
p_vp = 5291.502622129181

p_mu = p_density * p_vs**2
p_lambda = p_density * p_vp**2 - 2 * p_mu

p_youngs = p_mu*(3.0*p_lambda + 2.0*p_mu)/(p_lambda + p_mu)
p_poissons = 0.5*p_lambda/(p_lambda + p_mu)

p_fac1 = (1.0 + p_poissons)/p_youngs
p_fac2 = 2.0*(1.0 - 2.0*p_poissons)

# Cylinder parameters and BC.
a = 2000.0 # Cylinder inner radius.
b = 5000.0 # Cylinder outer radius.
Pa = 1.0e8 # Pressure applied to inner wall of cylinder.
A = -Pa/(1.0/a**2 - 1.0/b**2)
C = -0.5*A/b**2

# ----------------------------------------------------------------------
class AnalyticalSoln(object):
    """
    Analytical solution to axial extension problem.
    """
    SPACE_DIM = 3
    TENSOR_SIZE = 6

    def __init__(self):
        self.fields = {
            "displacement": self.displacement,
            "density": self.density,
            "shear_modulus": self.shear_modulus,
            "bulk_modulus": self.bulk_modulus,
            "cauchy_strain": self.strain,
            "cauchy_stress": self.stress,
            "initial_amplitude": self.traction,
        }
        return

    def getField(self, name, pts):
        return self.fields[name](pts)

    def displacement(self, locs):
        """
        Compute displacement field at locations.
        """
        # pdb.set_trace()
        (npts, dim) = locs.shape
        r = numpy.linalg.norm(locs[:,0:2], axis=1)
        ur = p_fac1*(-A/r + p_fac2*C*r)
        angs = numpy.arctan2(locs[:,1], locs[:,0])
        ca = numpy.cos(angs)
        sa = numpy.sin(angs)
        ux = ur*ca
        uy = ur*sa
        disp = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        disp[0,:,0] = ux
        disp[0,:,1] = uy
        disp[0,:,2] = numpy.zeros_like(ux)
        return disp

    def traction(self, locs):
        """
        Return applied traction on Neumann boundary.
        """
        # pdb.set_trace()
        (npts, dim) = locs.shape
        traction = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        traction[0,:,0] = 0.0
        traction[0,:,1] = 0.0
        traction[0,:,2] = -Pa
        return traction

    def density(self, locs):
        """
        Compute density field at locations.
        """
        (npts, dim) = locs.shape
        density = p_density * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return density

    def shear_modulus(self, locs):
        """
        Compute shear modulus field at locations.
        """
        (npts, dim) = locs.shape
        shear_modulus = p_mu * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return shear_modulus

    def bulk_modulus(self, locs):
        """
        Compute bulk modulus field at locations.
        """
        (npts, dim) = locs.shape
        bulk_modulus = (p_lambda + 2.0 / 3.0 * p_mu) * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return bulk_modulus

    def strain(self, locs):
        """
        Compute strain field at locations.
        """
        # pdb.set_trace()
        (npts, dim) = locs.shape
        strain = numpy.zeros((1, npts, self.TENSOR_SIZE), dtype=numpy.float64)
        r = numpy.linalg.norm(locs[:,0:2], axis=1)
        angs = numpy.arctan2(locs[:,1], locs[:,0])
        ca = numpy.cos(angs)
        sa = numpy.sin(angs)
        err = p_fac1*(A/r**2 + p_fac2*C)
        ett = p_fac1*(-A/r**2 + p_fac2*C)
        ezz = numpy.zeros_like(err)
        exx = err*ca*ca + ett*sa*sa
        eyy = err*sa*sa + ett*ca*ca
        strain[0, :, 0] = exx
        strain[0, :, 1] = eyy
        strain[0, :, 2] = ezz
        strain[0, :, 3] = ezz
        strain[0, :, 4] = ezz
        strain[0, :, 5] = ezz
        return strain

    def stress(self, locs):
        """
        Compute stress field at locations.
        """
        # pdb.set_trace()
        (npts, dim) = locs.shape
        stress = numpy.zeros((1, npts, self.TENSOR_SIZE), dtype=numpy.float64)
        r = numpy.linalg.norm(locs[:,0:2], axis=1)
        angs = numpy.arctan2(locs[:,1], locs[:,0])
        ca = numpy.cos(angs)
        sa = numpy.sin(angs)
        srr = A/r**2 + 2.0*C
        stt = -A/r**2 + 2.0*C
        szz = 4.0*p_poissons*C*numpy.ones_like(srr)
        sxx = srr*ca*ca + stt*sa*sa
        syy = srr*sa*sa + stt*ca*ca
        shear = numpy.zeros_like(sxx)
        stress[0, :, 0] = sxx
        stress[0, :, 1] = syy
        stress[0, :, 2] = szz
        stress[0, :, 3] = shear
        stress[0, :, 4] = shear
        stress[0, :, 5] = shear
        return stress


# End of file
