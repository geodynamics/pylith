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
# @file tests/fullscale/linearelasticity/faults/axialdisp_soln.py
#
# @brief Analytical solution to two rigid blocks displacement problem.
#
# 2-D rigid motion of blocks using Dirichlet BC and fault slip.
#
#  -----
#  |   |-----
#  |   ||   |
#  |   ||   |
#  -----|   |
#       -----
#
# Ux(x,y) = +1.0 if x < 0
#           -1.0 if x > 0
#
# Dirichlet boundary conditions
#   Uy(-4000,y) = +1.0
#   Ux(-4000,y) =  0.0
#
#   Uy(+4000,y) = -1.0
#   Ux(+4000,y) = 0.0

import numpy

# Physical properties
p_density = 2500.0
p_vs = 3000.0
p_vp = 5291.502622129181

p_mu = p_density * p_vs**2
p_lambda = p_density * p_vp**2 - 2 * p_mu

# Uniform stress field (plane strain)
sxx = 0.0
sxy = 0.0
syy = 0.0
szz = p_lambda / (2 * p_lambda + 2 * p_mu) * (sxx + syy)

# Uniform strain field
exx = 1.0 / (2 * p_mu) * (sxx - p_lambda / (3 * p_lambda + 2 * p_mu) * (sxx + syy + szz))
eyy = 1.0 / (2 * p_mu) * (syy - p_lambda / (3 * p_lambda + 2 * p_mu) * (sxx + syy + szz))
ezz = 1.0 / (2 * p_mu) * (szz - p_lambda / (3 * p_lambda + 2 * p_mu) * (sxx + syy + szz))
exy = 1.0 / (2 * p_mu) * (sxy)

# Time steps
ts = 0.5  # year
nts = 3
tsteps = numpy.arange(0, ts * nts, ts) # yaer

# ----------------------------------------------------------------------
class AnalyticalSoln(object):
    """Analytical solution to axial/shear displacement problem.
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
            "initial_amplitude": self.zero_vector,
        }
        return

    def getField(self, name, pts):
        return self.fields[name](pts)

    def getMask(self, name, pts):
        mask = None
        if name == "displacement":
            mask = pts[:, 0] == 0.0
        return mask
    
    def zero_vector(self, locs):
        (npts, dim) = locs.shape
        return numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)

    def displacement(self, locs):
        """Compute displacement field at locations.
        """
        (npts, dim) = locs.shape
        ntpts = tsteps.shape[0]
        disp = numpy.zeros((ntpts, npts, 2), dtype=numpy.float64)

        maskN = locs[:, 0] < 0
        maskP = locs[:, 0] > 0

        t_track = 0
        for t in tsteps:
            if t < 0.0:
                disp[t_track, :, 1] = 0.0
            else:
                disp[t_track, :, 1] = (+1.0 * maskN - 1.0 * maskP) * t_track/(ntpts-1)

            t_track += 1

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

    def strain(self, locs):
        """Compute strain field at locations.
        """
        (npts, dim) = locs.shape
        strain = numpy.zeros((ntpts, npts, self.TENSOR_SIZE), dtype=numpy.float64)

        strain[:,:, 0] = exx
        strain[:,:, 1] = eyy
        strain[:,:, 2] = ezz
        strain[:,:, 3] = exy

        return strain

    def stress(self, locs):
        """Compute stress field at locations.
        """
        (npts, dim) = locs.shape
        stress = numpy.zeros((ntpts, npts, self.TENSOR_SIZE), dtype=numpy.float64)
    
        stress[:,:, 0] = sxx
        stress[:,:, 1] = syy
        stress[:,:, 2] = szz
        stress[:,:, 3] = sxy

        return stress


# End of file
