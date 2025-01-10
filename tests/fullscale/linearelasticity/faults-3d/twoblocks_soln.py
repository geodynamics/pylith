# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================
# @brief Analytical solution to three rigid blocks displacement problem.
#
# 2D rigid motion of blocks using Dirichlet BC and fault slip.
#
#       -----
#  -----|   |
#  |   ||   |
#  |   ||   |
#  |   |-----
#  -----     
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

fault_dip = 0.4*numpy.pi
fault_normal = numpy.array([numpy.sin(fault_dip), 0.0, numpy.cos(fault_dip)])
fault_center = numpy.array([[0.0, 0.0, -4000]])

# ----------------------------------------------------------------------
class AnalyticalSoln(object):
    """Analytical solution to axial/shear displacement problem.
    """
    SPACE_DIM = 3
    TENSOR_SIZE = 6
    FAULT_DIP_ANGLE = 0.4*numpy.pi

    def __init__(self):
        sindip = numpy.sin(self.FAULT_DIP_ANGLE)
        cosdip = numpy.cos(self.FAULT_DIP_ANGLE)
        self.fields = {
            "displacement": self.displacement,
            "density": self.density,
            "shear_modulus": self.shear_modulus,
            "bulk_modulus": self.bulk_modulus,
            "cauchy_strain": self.strain,
            "cauchy_stress": self.stress,
            "initial_amplitude": self.displacement,
            "normal_dir": {
                "bc_xneg": self.orientation_dir((-1, 0, 0)),
                "bc_xpos": self.orientation_dir((+1, 0, 0)),
                "fault": self.orientation_dir((sindip, 0, cosdip)),
            },
            "horizontal_tangential_dir": {
                "bc_xneg": self.orientation_dir((0, -1, 0)),
                "bc_xpos": self.orientation_dir((0, +1, 0)),
            },
            "vertical_tangential_dir": {
                "bc_xneg": self.orientation_dir((0, 0, +1)),
                "bc_xpos": self.orientation_dir((0, 0, +1)),
            },
            "slip": self.slip,
            "traction_change": self.traction_change,
            "strike_dir": self.orientation_dir((0, +1, 0)),
            "up_dip_dir": self.orientation_dir((-cosdip, 0, sindip)),
        }
        return

    def getField(self, name, mesh_entity, pts):
        if isinstance(self.fields[name], dict):
            field = self.fields[name][mesh_entity](pts)
        else:
            field = self.fields[name](pts)
        return field

    def getMask(self, name, mesh_entity, pts):
        mask = None
        if name == "displacement":
            fdist = numpy.dot(pts-fault_center, fault_normal)
            index = numpy.argmin(numpy.abs(numpy.min(fdist)))
            mask = numpy.abs(fdist) < 1.0
        return mask
    
    def displacement(self, locs):
        """Compute displacement field at locations.
        """
        (npts, dim) = locs.shape

        fdist = numpy.dot(locs-fault_center, fault_normal)
        mask = numpy.zeros(npts)
        mask[fdist > 0.0] = +1.0
        mask[fdist < 0.0] = -1.0

        disp = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        disp[0,:, 0] = -0.5*numpy.cos(fault_dip) * mask
        disp[0,:, 1] = -1.0 * mask
        disp[0,:, 2] = 0.5*numpy.sin(fault_dip) * mask
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

    def slip(self, locs):
        """Compute slip field at locations.
        """
        (npts, dim) = locs.shape
        slip = numpy.zeros((1, npts, 3), dtype=numpy.float64)
        slip[:,:,1] = -2.0
        slip[:,:,2] = +1.0
        return slip

    def traction_change(self, locs):
        """Compute change in fault traction field at locations.
        """
        (npts, dim) = locs.shape
        traction_change = numpy.zeros((1, npts, 3), dtype=numpy.float64)
        return traction_change

    def orientation_dir(self, vector):
        def fn_dir(locs):
            (npts, dim) = locs.shape
            values = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
            for d in range(self.SPACE_DIM):
                values[:,:,d] = vector[d]
            return values
        return fn_dir


# End of file
