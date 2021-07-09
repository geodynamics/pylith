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
# @file tests/fullscale/linearelasticity/nofaults-2d/sheartraction_rate_soln.py
#
# @brief Analytical solution to time-dependent shear displacement/traction problem.
#
# 2-D uniform shear test.
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
#   Ux(-4000,y) = (a+v*H(t-t0))*y, Uy(-4000,y) = (a+v*H(t-t0))*x
#   Ux(x,-4000) = (a+v*H(t-t0))*y, Uy(x,-4000) = (a+v*H(t-t0))*y
# Neumann boundary conditions
#   \tau_shear(x,0) = -2*mu*(a+v*H(t-t0))
#   \tau_shear(+4000,y) = +2*mu*(a+v*H(t-t0))

import numpy

from pythia.pyre.units.time import year

# Physical properties
p_density = 2500.0
p_vs = 3000.0
p_vp = 5291.502622129181

p_mu = p_density * p_vs**2
p_lambda = p_density * p_vp**2 - 2 * p_mu

# Time steps
tsteps = numpy.arange(0.0, 5.01, 1.0)  # year

# Initial stress field (plane strain)
s0xx = 0.0
s0xy = 2.0e+6  # Pa
s0yy = 0.0
s0zz = p_lambda / (2 * p_lambda + 2 * p_mu) * (s0xx + s0yy)

# Rate of change of stress field
tR = 1.0  # year
sRxx = 0.0
sRxy = 0.5e+6  # Pa
sRyy = 0.0
sRzz = p_lambda / (2 * p_lambda + 2 * p_mu) * (sRxx + sRyy)


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
                "bc_yneg": self.bc_initial_displacement,
                "bc_xneg": self.bc_initial_displacement,
                "bc_xpos": self.bc_xpos_initial_traction,
                "bc_ypos": self.bc_ypos_initial_traction,
            },
            "rate_start_time": {
                "bc_yneg": self.bc_rate_time,
                "bc_xneg": self.bc_rate_time,
                "bc_xpos": self.bc_rate_time,
                "bc_ypos": self.bc_rate_time,
            },
            "rate_amplitude": {
                "bc_yneg": self.bc_velocity,
                "bc_xneg": self.bc_velocity,
                "bc_xpos": self.bc_xpos_rate_traction,
                "bc_ypos": self.bc_ypos_rate_traction,
            },
        }
        return

    def getField(self, name, mesh_entity, pts):
        if name in ["initial_amplitude", "rate_start_time", "rate_amplitude"]:
            field = self.fields[name][mesh_entity](pts)
        else:
            field = self.fields[name](pts)
        return field

    def displacement(self, locs):
        """Compute displacement field at locations.
        """
        strain = self.strain(locs)
        (ntpts, npts, tensorSize) = strain.shape
        disp = numpy.zeros((ntpts, npts, self.SPACE_DIM), dtype=numpy.float64)
        disp[:,:, 0] = strain[:,:, 0] * locs[:, 0] + strain[:,:, 3] * locs[:, 1]
        disp[:,:, 1] = strain[:,:, 3] * locs[:, 0] + strain[:,:, 1] * locs[:, 1]
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
        stress = self.stress(locs)
        sxx = stress[:,:, 0]
        syy = stress[:,:, 1]
        szz = stress[:,:, 2]
        sxy = stress[:,:, 3]
        strain = numpy.zeros(stress.shape, dtype=numpy.float64)
        strain[:,:, 0] = 1.0 / (2 * p_mu) * (sxx - p_lambda / (3 * p_lambda + 2 * p_mu) * (sxx + syy + szz))
        strain[:,:, 1] = 1.0 / (2 * p_mu) * (syy - p_lambda / (3 * p_lambda + 2 * p_mu) * (sxx + syy + szz))
        strain[:,:, 2] = 1.0 / (2 * p_mu) * (szz - p_lambda / (3 * p_lambda + 2 * p_mu) * (sxx + syy + szz))
        strain[:,:, 3] = 1.0 / (2 * p_mu) * (sxy)
        return strain

    def stress(self, locs):
        """Compute stress field at locations.
        """
        (npts, dim) = locs.shape
        ntpts = tsteps.shape[0]
        maskRate = tsteps >= tR
        sxx = s0xx + sRxx * maskRate * (tsteps - tR)
        syy = s0yy + sRyy * maskRate * (tsteps - tR)
        szz = s0zz + sRzz * maskRate * (tsteps - tR)
        sxy = s0xy + sRxy * maskRate * (tsteps - tR)

        ones = numpy.ones((1, npts), dtype=numpy.float64)
        stress = numpy.zeros((ntpts, npts, self.TENSOR_SIZE), dtype=numpy.float64)
        stress[:,:, 0] = numpy.dot(sxx.reshape((ntpts, 1)), ones)
        stress[:,:, 1] = numpy.dot(syy.reshape((ntpts, 1)), ones)
        stress[:,:, 2] = numpy.dot(szz.reshape((ntpts, 1)), ones)
        stress[:,:, 3] = numpy.dot(sxy.reshape((ntpts, 1)), ones)
        return stress

    def bc_initial_displacement(self, locs):
        """Compute initial displacement at locations.
        """
        exx = 1.0 / (2 * p_mu) * (s0xx - p_lambda / (3 * p_lambda + 2 * p_mu) * (s0xx + s0yy + s0zz))
        eyy = 1.0 / (2 * p_mu) * (s0yy - p_lambda / (3 * p_lambda + 2 * p_mu) * (s0xx + s0yy + s0zz))
        exy = 1.0 / (2 * p_mu) * (s0xy)
        (npts, dim) = locs.shape
        disp = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        disp[0,:, 0] = exx * locs[:, 0] + exy * locs[:, 1]
        disp[0,:, 1] = exy * locs[:, 0] + eyy * locs[:, 1]
        return disp

    def bc_velocity(self, locs):
        """Compute velocity at locations.
        """
        eRxx = 1.0 / (2 * p_mu) * (sRxx - p_lambda / (3 * p_lambda + 2 * p_mu) * (sRxx + sRyy + sRzz))
        eRyy = 1.0 / (2 * p_mu) * (sRyy - p_lambda / (3 * p_lambda + 2 * p_mu) * (sRxx + sRyy + sRzz))
        eRxy = 1.0 / (2 * p_mu) * (sRxy)
        (npts, dim) = locs.shape
        velocity = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        velocity[0,:, 0] = (eRxx * locs[:, 0] + eRxy * locs[:, 1]) / year.value
        velocity[0,:, 1] = (eRxy * locs[:, 0] + eRyy * locs[:, 1]) / year.value
        return velocity

    def bc_xpos_initial_traction(self, locs):
        """Compute initial traction at locations.
        """
        (npts, dim) = locs.shape
        traction = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        traction[0,:, 0] = s0xy
        traction[0,:, 1] = 0.0
        return traction

    def bc_ypos_initial_traction(self, locs):
        """Compute initial traction at locations.
        """
        (npts, dim) = locs.shape
        traction = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        traction[0,:, 0] = -s0xy
        traction[0,:, 1] = 0.0
        return traction

    def bc_rate_time(self, locs):
        """Compute rate start time at locations.
        """
        (npts, dim) = locs.shape
        return tR * numpy.ones((1, npts, 1), dtype=numpy.float64) * year.value

    def bc_xpos_rate_traction(self, locs):
        """Compute rate of change in traction on +x boundary at locations.
        """
        (npts, dim) = locs.shape
        traction = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        traction[0,:, 0] = sRxy / year.value
        traction[0,:, 1] = 0.0
        return traction

    def bc_ypos_rate_traction(self, locs):
        """Compute rate of change in traction on +y boundary at locations.
        """
        (npts, dim) = locs.shape
        traction = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        traction[0,:, 0] = -sRxy / year.value
        traction[0,:, 1] = 0.0
        return traction


# End of file
