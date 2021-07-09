# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2021 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------
#
# @file tests/fullscale/linearelasticity/nofaults-2d/gravity_incompressible_soln.py
#
# @brief Analytical solution to gravitational body foces for incompressible linear elasticity.
#
#          p=0
#       ----------
#       |        |
# Ux=0  |        | Ux=0
#       |        |
#       |        |
#       ----------
#         Uy=0
#
# Dirichlet boundary conditions
# Ux(+-4000,0) = 0
# Uy(x,-4000) = 0
# p(x,+4000) = 0

import numpy


# Physical properties
p_density = 2500.0  # kg/m**3
p_vs = 3000.0  # m/s
p_vp = 1.0e+15  # m/s

p_mu = p_density * p_vs**2
p_lambda = p_density * p_vp**2 - 2 * p_mu

gacc = 9.80665  # m/s
ymax = 0.0  # m


# ----------------------------------------------------------------------
class AnalyticalSoln(object):
    """Analytical solution to gravitational body forces for incompssible linear elasticity.
    """
    SPACE_DIM = 3
    TENSOR_SIZE = 6

    def __init__(self):
        self.fields = {
            "displacement": self.zero_vector,
            "pressure": self.pressure,
            "density": self.density,
            "shear_modulus": self.shear_modulus,
            "bulk_modulus": self.bulk_modulus,
            "cauchy_strain": self.strain,
            "cauchy_stress": self.stress,
            "gravitational_acceleration": self.gacc,
            "initial_amplitude": {
                "bc_xneg": self.zero_vector,
                "bc_xpos": self.zero_vector,
                "bc_yneg": self.zero_vector,
                "bc_ypos": self.zero_vector,
                "bc_zneg": self.zero_vector,
                "bc_zpos": self.zero_scalar,
            },
        }
        return

    def getField(self, name, mesh_entity, pts):
        if name in "initial_amplitude":
            field = self.fields[name][mesh_entity](pts)
        else:
            field = self.fields[name](pts)
        return field

    def displacement(self, locs):
        """Compute displacement field at locations.
        """
        strain = self.strain(locs)

        (npts, dim) = locs.shape
        disp = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        return disp

    def pressure(self, locs):
        """Compute pressure field at locations.
        """
        stress = self.stress(locs)
        (npts, dim) = locs.shape
        pressure = -stress[:,:, 2].reshape(1, npts, 1)
        return pressure

    def zero_vector(self, locs):
        (npts, dim) = locs.shape
        return numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)

    def zero_scalar(self, locs):
        (npts, dim) = locs.shape
        return numpy.zeros((1, npts, 1), dtype=numpy.float64)

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
        return strain

    def stress(self, locs):
        """Compute stress field at locations.
        """
        szz = p_density * gacc * (locs[:, 2] - ymax)
        sxx = szz
        syy = szz
        sxy = 0.0
        syz = 0.0
        sxz = 0.0

        (npts, dim) = locs.shape
        stress = numpy.zeros((1, npts, self.TENSOR_SIZE), dtype=numpy.float64)
        stress[0,:, 0] = sxx
        stress[0,:, 1] = syy
        stress[0,:, 2] = szz
        stress[0,:, 3] = sxy
        stress[0,:, 4] = syz
        stress[0,:, 5] = sxz
        return stress

    def gacc(self, locs):
        """Compute gravitational acceleration at locations.
        """
        (npts, dim) = locs.shape
        gravacc = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        gravacc[0,:, 2] = -gacc
        return gravacc


# End of file
