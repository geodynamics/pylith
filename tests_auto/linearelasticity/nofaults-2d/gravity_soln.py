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
# @file tests_auto/linearelasticity/nofaults-2d/gravity_soln.py
#
# @brief Analytical solution to gravitational body foces (no initial stress).
#
# 2-D gravitational body forces tiquadrilateral cells.
#
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
# Uy(-4000,y) = 0

import numpy


# Physical properties
p_density = 2500.0
p_vs = 3000.0
p_vp = 5291.502622129181
gacc = 9.80665
ymax = +4000.0

p_mu = p_density * p_vs**2
p_lambda = p_density * p_vp**2 - 2 * p_mu


# ----------------------------------------------------------------------
class AnalyticalSoln(object):
    """
    Analytical solution to axial extension problem.
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
            "gravitational_acceleration": self.gacc,
            "initial_amplitude": self.zero_vector,
        }
        return

    def getField(self, name, pts):
        return self.fields[name](pts)

    def displacement(self, locs):
        """
        Compute displacement field at locations.
        """
        strain = self.strain(locs)

        (npts, dim) = locs.shape
        disp = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        disp[:, :, 0] = strain[:, :, 0] * locs[:, 0] + strain[:, :, 3] * locs[:, 1]
        disp[:, :, 1] = strain[:, :, 3] * locs[:, 0] + strain[:, :, 1] * (ymax - locs[:, 1])
        return disp

    def zero_vector(self, locs):
        (npts, dim) = locs.shape
        return numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)

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
        stress = self.stress(locs)
        sxx = stress[:, :, 0]
        syy = stress[:, :, 1]
        szz = stress[:, :, 2]
        sxy = stress[:, :, 3]

        (npts, dim) = locs.shape
        strain = numpy.zeros((1, npts, self.TENSOR_SIZE), dtype=numpy.float64)
        strain[:, :, 0] = 1.0 / (2 * p_mu) * (sxx - p_lambda / (3 * p_lambda + 2 * p_mu) * (sxx + syy + szz))
        strain[:, :, 1] = 1.0 / (2 * p_mu) * (syy - p_lambda / (3 * p_lambda + 2 * p_mu) * (sxx + syy + szz))
        strain[:, :, 2] = 1.0 / (2 * p_mu) * (szz - p_lambda / (3 * p_lambda + 2 * p_mu) * (sxx + syy + szz))
        strain[:, :, 3] = 1.0 / (2 * p_mu) * (sxy)
        return strain

    def stress(self, locs):
        """
        Compute stress field at locations.
        """
        syy = -p_density * gacc * (ymax - locs[:, 1])
        sxx = p_lambda / (p_lambda + 2 * p_mu) * syy
        sxy = 0.0
        szz = p_lambda / (2 * p_lambda + 2 * p_mu) * (sxx + syy)

        (npts, dim) = locs.shape
        stress = numpy.zeros((1, npts, self.TENSOR_SIZE), dtype=numpy.float64)
        stress[0, :, 0] = sxx
        stress[0, :, 1] = syy
        stress[0, :, 2] = szz
        stress[0, :, 3] = sxy
        return stress

    def gacc(self, locs):
        """Compute gravitational acceleration at locations.
        """
        (npts, dim) = locs.shape
        gravacc = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        gravacc[0, :, 1] = -gacc
        return gravacc


# End of file
