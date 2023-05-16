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
# @file tests/fullscale/linearelasticity/faults-2d/sheablocks_soln.py
#
# @brief Analytical solution to shear displacement/traction problem.
#
# 2-D uniform shear test.
#
#              ---->
#          ------ -------
#          |    | |     | 
#        | |    | |     | ^
#        v |    | |     | |
#          |    | |     |
#          ------ -------
#              <----
#
# Dirichlet boundary conditions
#   Ux(-4000,y) = a*y, Uy(-4000,y) = a*x
#   Ux(x,-4000) = a*y, Uy(x,-4000) = a*y
# Neumann boundary conditions
#   \tau_shear(x,0) = -*mu*a
#   \tau_shear(+4000,y) = +*mu*a

import numpy


# Physical properties
p_density = 2500.0
p_vs = 3000.0
p_vp = 5291.502622129181

p_mu = p_density * p_vs**2
p_lambda = p_density * p_vp**2 - 2 * p_mu

# Uniform stress field (plane strain)
sxx = 0.0
sxy = 5.0e+6
syy = 0.0
szz = p_lambda / (2 * p_lambda + 2 * p_mu) * (sxx + syy)

# Uniform strain field
exx = 1.0 / (2 * p_mu) * (sxx - p_lambda / (3 * p_lambda + 2 * p_mu) * (sxx + syy + szz))
eyy = 1.0 / (2 * p_mu) * (syy - p_lambda / (3 * p_lambda + 2 * p_mu) * (sxx + syy + szz))
ezz = 1.0 / (2 * p_mu) * (szz - p_lambda / (3 * p_lambda + 2 * p_mu) * (sxx + syy + szz))

exy = 1.0 / (2 * p_mu) * (sxy)

# print(exx,eyy,exy,ezz)
# print(exx*p_lambda/(p_lambda+2*p_mu))


# ----------------------------------------------------------------------
class AnalyticalSoln(object):
    """Analytical solution to shear problem.
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
            "initial_amplitude": {
                "bc_xneg": self.displacement,
                "bc_xpos": self.displacement,
                "bc_yneg": self.bc_traction,
                "bc_ypos": self.bc_traction,
            },
            "slip": self.slip,
            "lagrange_multiplier_fault": self.lagrange_multiplier_fault,
        }
        return

    def getField(self, name, mesh_entity, pts):
        if name == "initial_amplitude":
            field = self.fields[name][mesh_entity](pts)
        else:
            field = self.fields[name](pts)
        return field

    def displacement(self, locs):
        """Compute displacement field at locations.
        """
        (npts, dim) = locs.shape
        disp = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        disp[0,:, 0] = exx * locs[:, 0] + exy * locs[:, 1]
        disp[0,:, 1] = eyy * locs[:, 1] + exy * locs[:, 0]
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

    def bc_traction(self, locs):
        """Compute initial traction at locations.
        """
        (npts, dim) = locs.shape
        traction = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        traction[0,:, 0] = -sxy
        traction[0,:, 1] = 0.0
        return traction

    def slip(self, locs):
        """Compute slip field.
        """
        (npts, dim) = locs.shape
        slip = numpy.zeros((1, npts, 2), dtype=numpy.float64)
        return slip

    def lagrange_multiplier_fault(self, locs):
        """Compute Lagrange multiplier on faults.
        """
        (npts, dim) = locs.shape
        traction = numpy.zeros((1, npts, 2), dtype=numpy.float64)
        traction[:,:,1] = sxy
        return traction


# End of file
