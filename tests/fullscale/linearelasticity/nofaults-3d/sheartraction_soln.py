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
# @file tests/fullscale/linearelasticity/nofaults-3d/sheartraction_soln.py
#
# @brief Analytical solution to shear displacement/traction problem.
#
# 3-D uniform shear test.
#
#             --->
#          ----------
#          |        |
#        | |        | ^
#        v |        | |
#          |        |
#          ----------
#             <--
#
# Dirichlet boundary conditions
# boundary_xneg: Ux(-6000,y,z) = a*y, Uy(-6000,y,z) = a*x, Uz=0
# boundary_yneg: Ux(x,-6000,z) = a*y, Uy(x,-6000,z) = a*y, Uz=0
# boundary_zneg: Uz=0
# Neumann boundary conditions
#   \tau_shear_horiz(x,0,z) = -2*mu*a
#   \tau_shear_horiz(+6000,y,z) = +2*mu*a

import numpy


# Physical properties
p_density = 2500.0
p_vs = 3000.0
p_vp = 5291.502622129181

p_mu = p_density * p_vs**2
p_lambda = p_density * p_vp**2 - 2 * p_mu

# Uniform stress field (plane strain)
sxx = 0.0
syy = 0.0
szz = 0.0
sxy = 5.0e+6
syz = 0.0
sxz = 0.0

# Uniform strain field
exx = 1.0 / (2 * p_mu) * (sxx - p_lambda / (3 * p_lambda + 2 * p_mu) * (sxx + syy + szz))
eyy = 1.0 / (2 * p_mu) * (syy - p_lambda / (3 * p_lambda + 2 * p_mu) * (sxx + syy + szz))
ezz = 1.0 / (2 * p_mu) * (szz - p_lambda / (3 * p_lambda + 2 * p_mu) * (sxx + syy + szz))

exy = 1.0 / (2 * p_mu) * (sxy)
eyz = 1.0 / (2 * p_mu) * (syz)
exz = 1.0 / (2 * p_mu) * (sxz)


# ----------------------------------------------------------------------
class AnalyticalSoln(object):
    """
    Analytical solution to shear problem.
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
            "initial_amplitude": {
                "bc_yneg": self.displacement,
                "bc_xneg": self.displacement,
                "bc_xpos": self.bc_xpos_traction,
                "bc_ypos": self.bc_ypos_traction,
                "bc_zneg": self.displacement,
            }
        }
        self.key = None
        return

    def getField(self, name, pts):
        if self.key is None:
            field = self.fields[name](pts)
        else:
            field = self.fields[name][self.key](pts)
        return field

    def displacement(self, locs):
        """
        Compute displacement field at locations.
        """
        (npts, dim) = locs.shape
        disp = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        disp[0,:, 0] = exx * locs[:, 0] + exy * locs[:, 1] + exz * locs[:, 2]
        disp[0,:, 1] = exy * locs[:, 0] + eyy * locs[:, 1] + eyz * locs[:, 2]
        disp[0,:, 2] = exz * locs[:, 0] + eyz * locs[:, 1] + ezz * locs[:, 2]
        return disp

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
        (npts, dim) = locs.shape
        strain = numpy.zeros((1, npts, self.TENSOR_SIZE), dtype=numpy.float64)
        strain[0,:, 0] = exx
        strain[0,:, 1] = eyy
        strain[0,:, 2] = ezz
        strain[0,:, 3] = exy
        strain[0,:, 4] = eyz
        strain[0,:, 5] = exz
        return strain

    def stress(self, locs):
        """
        Compute stress field at locations.
        """
        (npts, dim) = locs.shape
        stress = numpy.zeros((1, npts, self.TENSOR_SIZE), dtype=numpy.float64)
        stress[0,:, 0] = sxx
        stress[0,:, 1] = syy
        stress[0,:, 2] = szz
        stress[0,:, 3] = sxy
        stress[0,:, 4] = syz
        stress[0,:, 5] = sxz
        return stress

    def bc_xpos_traction(self, locs):
        """Compute initial traction at locations.
        """
        (npts, dim) = locs.shape
        traction = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        traction[0,:, 0] = sxy
        traction[0,:, 1] = 0.0
        traction[0,:, 2] = 0.0
        return traction

    def bc_ypos_traction(self, locs):
        """Compute initial traction at locations.
        """
        (npts, dim) = locs.shape
        traction = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        traction[0,:, 0] = -sxy
        traction[0,:, 1] = 0.0
        traction[0,:, 2] = 0.0
        return traction


# End of file
