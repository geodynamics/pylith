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
# @file tests/fullscale/linearelasticity/nofaults-3d/gravity_soln.py
#
# @brief Analytical solution to gravitational body foces (no initial stress).
#
# Dirichlet boundary conditions on lateral sides and bottom
# boundary +-x: Ux(+-4000,0,z) = 0
# boundary +-y: Uy(x,-4000,z) = 0
# boundary -z: Uz(x,y,-9000) = 0

import numpy


# Physical properties
p_density = 2500.0  # kg/m**3
p_vs = 3000.0  # m/s
p_vp = 5291.5026  # m/s

p_mu = p_density * p_vs**2
p_lambda = p_density * p_vp**2 - 2 * p_mu

gacc = 9.80665  # m/s
zmax = 0.0
zmin = -9000.0  # m


# ----------------------------------------------------------------------
class AnalyticalSoln(object):
    """Analytical solution to gravitational body forces (no initial stress).
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
            "gravitational_acceleration": self.gacc,
            "initial_amplitude": self.zero_vector,
        }
        return

    def getField(self, name, mesh_entity, pts):
        return self.fields[name](pts)

    def displacement(self, locs):
        """Compute displacement field at locations.
        """
        strain = self.strain(locs)

        (npts, dim) = locs.shape
        disp = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        disp[:,:, 2] = p_density * gacc / (p_lambda + 2 * p_mu) * \
            (0.5 * (locs[:, 2]**2 - zmin**2) - zmax * (locs[:, 2] - zmin))
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
        exx = 0.0
        eyy = 0.0
        ezz = p_density * gacc * (locs[:, 2] - zmax) / (p_lambda + 2 * p_mu)
        exy = 0.0
        eyz = 0.0
        exz = 0.0

        (npts, dim) = locs.shape
        strain = numpy.zeros((1, npts, self.TENSOR_SIZE), dtype=numpy.float64)
        strain[:,:, 0] = exx
        strain[:,:, 1] = eyy
        strain[:,:, 2] = ezz
        strain[:,:, 3] = exy
        strain[:,:, 4] = eyz
        strain[:,:, 5] = exz
        return strain

    def stress(self, locs):
        """Compute stress field at locations.
        """
        szz = p_density * gacc * (locs[:, 2] - zmax)
        sxx = p_lambda / (p_lambda + 2 * p_mu) * szz
        syy = p_lambda / (p_lambda + 2 * p_mu) * szz
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
