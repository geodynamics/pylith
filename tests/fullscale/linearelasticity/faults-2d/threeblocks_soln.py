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
# @brief Analytical solution to three rigid blocks displacement problem.
#
# 2-D rigid motion of blocks using Dirichlet BC and fault slip.
#
#       -----
#  -----|   |-----
#  |   ||   ||   |
#  |   ||   ||   |
#  |   |-----|   |
#  -----     -----
#
# Uy(x,y) = 0.0 if x < 2 or x > 0
#          +2.0 if x > 2 and x < 0
#
# Dirichlet boundary conditions
#   Ux(-4000,y) =  0.0
#   Uy(-4000,y) =  0.0
#
#   Ux(+4000,y) = 0.0
#   Uy(+4000,y) = 0.0

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
            "initial_amplitude": self.displacement,
            "slip": {
                "fault_xmid": self.slip_xmid,
                "fault_xneg": self.slip_xneg,
            },
            "traction_change": self.traction_change,
            "normal_dir": self.orientation_dir((+1, 0)),
            "strike_dir": self.orientation_dir((0, +1)),
        }

    def getField(self, name, mesh_entity, pts):
        if isinstance(self.fields[name], dict):
            field = self.fields[name][mesh_entity](pts)
        else:
            field = self.fields[name](pts)
        return field

    def getMask(self, name, mesh_entity, pts):
        mask = None
        if name == "displacement":
            mask = numpy.logical_or(pts[:, 0] == 0.0, pts[:,0] == -2.0e+3)
        return mask
    
    def displacement(self, locs):
        """Compute displacement field at locations.
        """
        (npts, dim) = locs.shape
        disp = numpy.zeros((1, npts, 2), dtype=numpy.float64)
        mask = numpy.logical_and(locs[:, 0] > -2.0e+3, locs[:,0] < 0.0)
        disp[0,:, 1] = +2.0 * mask
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
        strain = numpy.zeros((1, npts, self.TENSOR_SIZE), dtype=numpy.float64)
        strain[0,:, 0] = exx
        strain[0,:, 1] = eyy
        strain[0,:, 2] = ezz
        strain[0,:, 3] = exy
        return strain

    def stress(self, locs):
        """Compute stress field at locations.
        """
        (npts, dim) = locs.shape
        stress = numpy.zeros((1, npts, self.TENSOR_SIZE), dtype=numpy.float64)
        stress[0,:, 0] = sxx
        stress[0,:, 1] = syy
        stress[0,:, 2] = szz
        stress[0,:, 3] = sxy
        return stress

    def slip_xmid(self, locs):
        """Compute slip field on fault xmid.
        """
        (npts, dim) = locs.shape
        slip = numpy.zeros((1, npts, 2), dtype=numpy.float64)
        slip[0,:, 1] = -2.0
        return slip

    def slip_xneg(self, locs):
        """Compute slip field on fault xneg.
        """
        (npts, dim) = locs.shape
        slip = numpy.zeros((1, npts, 2), dtype=numpy.float64)
        slip[0,:, 1] = +2.0
        return slip

    def traction_change(self, locs):
        """Compute change in traction on faults.
        """
        (npts, dim) = locs.shape
        traction = numpy.zeros((1, npts, 2), dtype=numpy.float64)
        return traction

    def orientation_dir(self, vector):
        def fn_dir(locs):
            (npts, dim) = locs.shape
            values = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
            for d in range(self.SPACE_DIM):
                values[:,:,d] = vector[d]
            return values
        return fn_dir


# End of file
