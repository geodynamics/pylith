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
"""Checks for output files from Green's function problem.
"""

import numpy

# Physical properties
p_density = 2500.0
p_vs = 3000.0
p_vp = 5291.502622129181

p_mu = p_density * p_vs**2
p_lambda = p_density * p_vp**2 - 2 * p_mu


# ----------------------------------------------------------------------
class AnalyticalSoln(object):
    """Analytical solution for material properties and initial conditions.
    """
    SPACE_DIM = 2
    TENSOR_SIZE = 4

    def __init__(self):
        self.fields = {
            "density": self.density,
            "shear_modulus": self.shear_modulus,
            "bulk_modulus": self.bulk_modulus,
            "initial_amplitude": self.initial_displacement,
        }
        return

    def getField(self, name, mesh_entity, pts):
        return self.fields[name](pts)

    def initial_displacement(self, locs):
        """Compute displacement field at locations.
        """
        (npts, dim) = locs.shape
        disp = numpy.zeros((1, npts, 2), dtype=numpy.float64)
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


# ----------------------------------------------------------------------
class SolnDims(object):
    """Sizes of output fields.
    """
    SPACE_DIM = 2
    TENSOR_SIZE = 4

    def __init__(self, nimpulses, basis_order):
        self.nimpulses = nimpulses
        self.basis_order = basis_order
        self.fields = {
            "displacement": self.displacement,
            "cauchy_stress": self.cauchy_stress,
            "cauchy_strain": self.cauchy_strain,
            "slip": self.slip,
        }
        return

    def getDims(self, name, mesh):
        basis_order = self.basis_order[name]
        npoints = {
            0: mesh.ncells,
            1: mesh.nvertices,
            2: mesh.nvertices + 0,
        }
        return self.fields[name](npoints[basis_order])

    def displacement(self, npoints):
        return (self.nimpulses, npoints, self.SPACE_DIM)

    def cauchy_stress(self, npoints):
        return (self.nimpulses, npoints, self.TENSOR_SIZE)

    def cauchy_strain(self, npoints):
        return (self.nimpulses, npoints, self.TENSOR_SIZE)

    def slip(self, npoints):
        return (self.nimpulses, npoints, self.SPACE_DIM)


# End of file
