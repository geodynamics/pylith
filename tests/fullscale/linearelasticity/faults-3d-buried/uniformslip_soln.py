# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2022 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------
#
# @brief Uniform slip on a buried fault (no analytical solution).
#
# We just check info fields.

import numpy

# Physical properties
p_density = 2500.0
p_vs = 3000.0
p_vp = 5291.502622129181

p_mu = p_density * p_vs**2
p_lambda = p_density * p_vp**2 - 2 * p_mu

# ----------------------------------------------------------------------
class AnalyticalSoln(object):
    """Analytical solution to axial/shear displacement problem.
    """
    SPACE_DIM = 3
    TENSOR_SIZE = 6

    def __init__(self):
        self.fields = {
            "density": self.density,
            "shear_modulus": self.shear_modulus,
            "bulk_modulus": self.bulk_modulus,
            "initial_amplitude": self.zero_displacement,
            "slip": self.slip,
        }
        return

    def getField(self, name, mesh_entity, pts):
        return self.fields[name](pts)

    def getMask(self, name, mesh_entity, pts):
        mask = None
        if name == "displacement":
            fdist = numpy.dot(pts-fault_center, fault_normal)
            index = numpy.argmin(numpy.abs(numpy.min(fdist)))
            mask = numpy.abs(fdist) < 1.0
        return mask
    
    def zero_displacement(self, locs):
        """Compute displacement field at locations.
        """
        (npts, dim) = locs.shape

        disp = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        return disp

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

    def slip(self, locs):
        """Compute slip field at locations.
        """
        (npts, dim) = locs.shape
        slip = numpy.zeros((1, npts, 3), dtype=numpy.float64)
        slip[:,:,1] = +2.0
        return slip


# End of file
