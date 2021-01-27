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
# @file tests/fullscale/linearelasticity/nofaults-3d/sheartraction_rate_soln.py
#
# @brief Analytical solution to time-dependent shear displacement/traction problem.
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
#
# boundary_xneg:
#   Ux(-6000,y) = (a+v*H(t-t0))*y, Uy(-6000,y) = (a+v*H(t-t0))*x
#
# boundary_yneg:
#   Ux(x,-6000) = (a+v*H(t-t0))*y, Uy(x,-6000) = (a+v*H(t-t0))*y
#
# boundary_zneg:
#   Uz=0
#
# Neumann boundary conditions
#
# boundary_xpos:
#   \tau_shear(x,0) = -2*mu*(a+v*H(t-t0))
#
# boundary_ypos:
#   \tau_shear(+4000,y) = +2*mu*(a+v*H(t-t0))

import numpy

from pyre.units.time import year

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
s0yy = 0.0
s0zz = 0.0
s0xy = 2.0e+6  # Pa
s0yz = 0.0
s0xz = 0.0

# Rate of change of stress field
tR = 1.0  # year
sRxx = 0.0
sRyy = 0.0
sRzz = 0.0
sRxy = 0.5e+6  # Pa
sRyz = 0.0
sRxz = 0.0


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
                "bc_yneg": self.bc_initial_displacement,
                "bc_xneg": self.bc_initial_displacement,
                "bc_xpos": self.bc_xpos_initial_traction,
                "bc_ypos": self.bc_ypos_initial_traction,
                "bc_zneg": self.bc_initial_displacement,
            },
            "rate_start_time": {
                "bc_yneg": self.bc_rate_time,
                "bc_xneg": self.bc_rate_time,
                "bc_xpos": self.bc_rate_time,
                "bc_ypos": self.bc_rate_time,
                "bc_zneg": self.bc_rate_time,
            },
            "rate_amplitude": {
                "bc_yneg": self.bc_velocity,
                "bc_xneg": self.bc_velocity,
                "bc_xpos": self.bc_xpos_rate_traction,
                "bc_ypos": self.bc_ypos_rate_traction,
                "bc_zneg": self.bc_velocity,
            },
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
        strain = self.strain(locs)
        (ntpts, npts, tensorSize) = strain.shape
        disp = numpy.zeros((ntpts, npts, self.SPACE_DIM), dtype=numpy.float64)
        disp[:,:, 0] = strain[:,:, 0] * locs[:, 0] + strain[:,:, 3] * locs[:, 1]
        disp[:,:, 1] = strain[:,:, 3] * locs[:, 0] + strain[:,:, 1] * locs[:, 1]
        disp[:,:, 0] = strain[:,:, 0] * locs[:, 0] + strain[:,:, 3] * locs[:, 1] + strain[:,:, 5] * locs[:, 2]
        disp[:,:, 1] = strain[:,:, 3] * locs[:, 0] + strain[:,:, 1] * locs[:, 1] + strain[:,:, 4] * locs[:, 2]
        disp[:,:, 2] = strain[:,:, 5] * locs[:, 0] + strain[:,:, 4] * locs[:, 1] + strain[:,:, 2] * locs[:, 2]
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
        stress = self.stress(locs)
        sxx = stress[:,:, 0]
        syy = stress[:,:, 1]
        szz = stress[:,:, 2]
        sxy = stress[:,:, 3]
        syz = stress[:,:, 4]
        sxz = stress[:,:, 5]
        strain = numpy.zeros(stress.shape, dtype=numpy.float64)
        strain[:,:, 0] = 1.0 / (2 * p_mu) * (sxx - p_lambda / (3 * p_lambda + 2 * p_mu) * (sxx + syy + szz))
        strain[:,:, 1] = 1.0 / (2 * p_mu) * (syy - p_lambda / (3 * p_lambda + 2 * p_mu) * (sxx + syy + szz))
        strain[:,:, 2] = 1.0 / (2 * p_mu) * (szz - p_lambda / (3 * p_lambda + 2 * p_mu) * (sxx + syy + szz))
        strain[:,:, 3] = 1.0 / (2 * p_mu) * (sxy)
        strain[:,:, 4] = 1.0 / (2 * p_mu) * (syz)
        strain[:,:, 5] = 1.0 / (2 * p_mu) * (sxz)
        return strain

    def stress(self, locs):
        """
        Compute stress field at locations.
        """
        (npts, dim) = locs.shape
        ntpts = tsteps.shape[0]
        maskRate = tsteps >= tR
        sxx = s0xx + sRxx * maskRate * (tsteps - tR)
        syy = s0yy + sRyy * maskRate * (tsteps - tR)
        szz = s0zz + sRzz * maskRate * (tsteps - tR)
        sxy = s0xy + sRxy * maskRate * (tsteps - tR)
        syz = s0yz + sRyz * maskRate * (tsteps - tR)
        sxz = s0xz + sRxz * maskRate * (tsteps - tR)

        ones = numpy.ones((1, npts), dtype=numpy.float64)
        stress = numpy.zeros((ntpts, npts, self.TENSOR_SIZE), dtype=numpy.float64)
        stress[:,:, 0] = numpy.dot(sxx.reshape((ntpts, 1)), ones)
        stress[:,:, 1] = numpy.dot(syy.reshape((ntpts, 1)), ones)
        stress[:,:, 2] = numpy.dot(szz.reshape((ntpts, 1)), ones)
        stress[:,:, 3] = numpy.dot(sxy.reshape((ntpts, 1)), ones)
        stress[:,:, 4] = numpy.dot(syz.reshape((ntpts, 1)), ones)
        stress[:,:, 5] = numpy.dot(sxz.reshape((ntpts, 1)), ones)
        return stress

    def bc_initial_displacement(self, locs):
        """Compute initial displacement at locations.
        """
        exx = 1.0 / (2 * p_mu) * (s0xx - p_lambda / (3 * p_lambda + 2 * p_mu) * (s0xx + s0yy + s0zz))
        eyy = 1.0 / (2 * p_mu) * (s0yy - p_lambda / (3 * p_lambda + 2 * p_mu) * (s0xx + s0yy + s0zz))
        ezz = 1.0 / (2 * p_mu) * (s0zz - p_lambda / (3 * p_lambda + 2 * p_mu) * (s0xx + s0yy + s0zz))
        exy = 1.0 / (2 * p_mu) * (s0xy)
        eyz = 1.0 / (2 * p_mu) * (s0yz)
        exz = 1.0 / (2 * p_mu) * (s0xz)
        (npts, dim) = locs.shape
        disp = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        disp[0,:, 0] = exx * locs[:, 0] + exy * locs[:, 1] + exz * locs[:, 2]
        disp[0,:, 1] = exy * locs[:, 0] + eyy * locs[:, 1] + eyz * locs[:, 2]
        disp[0,:, 2] = exz * locs[:, 0] + eyz * locs[:, 1] + ezz * locs[:, 2]
        return disp

    def bc_velocity(self, locs):
        """Compute velocity at locations.
        """
        eRxx = 1.0 / (2 * p_mu) * (sRxx - p_lambda / (3 * p_lambda + 2 * p_mu) * (sRxx + sRyy + sRzz))
        eRyy = 1.0 / (2 * p_mu) * (sRyy - p_lambda / (3 * p_lambda + 2 * p_mu) * (sRxx + sRyy + sRzz))
        eRzz = 1.0 / (2 * p_mu) * (sRzz - p_lambda / (3 * p_lambda + 2 * p_mu) * (sRxx + sRyy + sRzz))
        eRxy = 1.0 / (2 * p_mu) * (sRxy)
        eRyz = 1.0 / (2 * p_mu) * (sRyz)
        eRxz = 1.0 / (2 * p_mu) * (sRxz)
        (npts, dim) = locs.shape
        velocity = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        velocity[0,:, 0] = (eRxx * locs[:, 0] + eRxy * locs[:, 1] + eRxz * locs[:, 2]) / year.value
        velocity[0,:, 1] = (eRxy * locs[:, 0] + eRyy * locs[:, 1] + eRyz * locs[:, 2]) / year.value
        velocity[0,:, 2] = (eRxz * locs[:, 0] + eRyz * locs[:, 1] + eRzz * locs[:, 2]) / year.value
        return velocity

    def bc_xpos_initial_traction(self, locs):
        """Compute initial traction at locations.
        """
        (npts, dim) = locs.shape
        traction = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        traction[0,:, 0] = s0xy
        traction[0,:, 1] = 0.0
        traction[0,:, 2] = 0.0
        return traction

    def bc_ypos_initial_traction(self, locs):
        """Compute initial traction at locations.
        """
        (npts, dim) = locs.shape
        traction = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        traction[0,:, 0] = -s0xy
        traction[0,:, 1] = 0.0
        traction[0,:, 2] = 0.0
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
        traction[0,:, 2] = 0.0
        return traction

    def bc_ypos_rate_traction(self, locs):
        """Compute rate of change in traction on +y boundary at locations.
        """
        (npts, dim) = locs.shape
        traction = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        traction[0,:, 0] = -sRxy / year.value
        traction[0,:, 1] = 0.0
        traction[0,:, 2] = 0.0
        return traction


# End of file
