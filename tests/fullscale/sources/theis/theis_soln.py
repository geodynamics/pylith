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
# @file tests/fullscale/poroelasticity/sources/theis/theis_soln.py
#
# @brief Analytical solution to Theis's problem.
#
#          Uy=0
#       ----------
#       |        |
# Ux=0  |        | Ux=0
#       |        |
#       |        |
#       ----------
#         +1 Pa
#
# Dirichlet boundary conditions
#   Ux(+-1,0) = 0
#   Uy(x,+1) = 0
# Neumann boundary conditions
#   \tau_normal(x,+1) = 1*Pa

import numpy

# Physical properties
rho_s = 2500  # kg / m**3
rho_f = 1000  # kg / m**3
mu_f = 0.001  # Pa*s
G = 3.0  # Pa
K_sg = 10.0  # Pa
K_fl = 2e9  # Pa
K_d = 4.0  # Pa
# K_u = 2.6941176470588233 # Pa
alpha = 1.0 # -
phi = 0.2  # -
M = 4.705882352941176# Pa
k = 1e-14  # m**2
grav = 9.80665 # m/s**2

ymax = 100.0  # m
ymin = 0.0  # m
xmax = 100.0  # m
xmin = 0.0  # m
P_0 = 20e6  # Pa

# Hydrology parameters
K = (k*grav*rho_f) / mu_f # m / s
b = 1.0 # m
S_s = (rho_f*grav) / K_fl # 1 / m
T = K*b # m**2 / s
Q = 0.02 # m**3 / s

# Height of column, m
L = ymax - ymin
H = xmax - xmin

M = 1.0 / (phi / K_fl + (alpha - phi) / K_sg)  # Pa
K_u = K_d + alpha * alpha * M  # Pa,      Cheng (B.5)
nu = (3.0 * K_d - 2.0 * G) / (2.0 * (3.0 * K_d + G))  # -,       Cheng (B.8)
nu_u = (3.0 * K_u - 2.0 * G) / (2.0 * (3.0 * K_u + G))  # -,       Cheng (B.9)
eta = (3.0 * alpha * G) / (3.0 * K_d + 4.0 * G)  # -,       Cheng (B.11)
S = (3.0 * K_u + 4.0 * G) / (M * (3.0 * K_d + 4.0 * G))  # Pa^{-1}, Cheng (B.14)
c = (k / mu_f) / S  # m^2 / s, Cheng (B.16)

# Time steps
ts = 0.0028666667  # sec
nts = 2
tsteps = numpy.arange(0.0, ts * nts, ts)  + ts # sec

# ----------------------------------------------------------------------


class AnalyticalSoln(object):
    """Analytical solution to Theis's problem
    """
    SPACE_DIM = 2
    TENSOR_SIZE = 4
    ITERATIONS = 16000

    def __init__(self):
        self.fields = {
            "displacement": self.displacement,
            "pressure": self.pressure,
            "porosity": self.porosity,
            "trace_strain": self.trace_strain,
            "solid_density": self.solid_density,
            "fluid_density": self.fluid_density,
            "fluid_viscosity": self.fluid_viscosity,
            "shear_modulus": self.shear_modulus,
            "drained_bulk_modulus": self.drained_bulk_modulus,
            "biot_coefficient": self.biot_coefficient,
            "biot_modulus": self.biot_modulus,
            "isotropic_permeability": self.isotropic_permeability,
            "initial_amplitude": {
                "x_neg": self.zero_vector,
                "x_pos": self.zero_vector,
                "y_pos_neu": self.y_pos_neu,
                "y_pos_dir": self.zero_scalar,
                "y_neg": self.zero_vector,
            }
        }
        return

    def getField(self, name, mesh_entity, pts):
        if name in "initial_amplitude":
            field = self.fields[name][mesh_entity](pts)
        else:
            field = self.fields[name](pts)
        return field

    def zero_scalar(self, locs):
        (npts, dim) = locs.shape
        return numpy.zeros((1, npts, 1), dtype=numpy.float64)

    def zero_vector(self, locs):
        (npts, dim) = locs.shape
        return numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)

    def solid_density(self, locs):
        """Compute solid_density field at locations.
        """
        (npts, dim) = locs.shape
        solid_density = rho_s * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return solid_density

    def fluid_density(self, locs):
        """Compute fluid density field at locations.
        """
        (npts, dim) = locs.shape
        fluid_density = rho_f * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return fluid_density

    def shear_modulus(self, locs):
        """Compute shear modulus field at locations.
        """
        (npts, dim) = locs.shape
        shear_modulus = G * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return shear_modulus

    def porosity(self, locs):
        """Compute porosity field at locations.
        """
        (npts, dim) = locs.shape
        porosity = phi * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return porosity

    def fluid_viscosity(self, locs):
        """Compute fluid_viscosity field at locations.
        """
        (npts, dim) = locs.shape
        fluid_viscosity = mu_f * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return fluid_viscosity

    def drained_bulk_modulus(self, locs):
        """Compute undrained bulk modulus field at locations.
        """
        (npts, dim) = locs.shape
        undrained_bulk_modulus = K_d * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return undrained_bulk_modulus

    def biot_coefficient(self, locs):
        """Compute biot coefficient field at locations.
        """
        (npts, dim) = locs.shape
        biot_coefficient = alpha * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return biot_coefficient

    def biot_modulus(self, locs):
        """Compute biot modulus field at locations.
        """
        (npts, dim) = locs.shape
        biot_modulus = M * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return biot_modulus

    def isotropic_permeability(self, locs):
        """Compute isotropic permeability field at locations.
        """
        (npts, dim) = locs.shape
        isotropic_permeability = k * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return isotropic_permeability

    def displacement(self, locs):
        """Compute displacement field at locations.
        """
        (npts, dim) = locs.shape
        ntpts = tsteps.shape[0]
        displacement = numpy.zeros((ntpts, npts, dim), dtype=numpy.float64)
        z = locs[:, 1]
        t_track = 0
        z_star = 1 - z / L

        # for t in tsteps:
        #     if t < 0.0:
        #         displacement[0, :, 1] = ((P_0 * L * (1.0 - 2.0 * nu_u)) / (2.0 * G * (1.0 - nu_u))) * (1.0 - z_star)
        #     else:
        #         t_star = (c * t) / ((2 * L)**2)
        #         displacement[t_track, :, 1] =  (((P_0 * L * (1.0 - 2.0 * nu_u)) / (2.0 * G * (1.0 - nu_u))) * (1.0 - z_star) + ((P_0 * L * (nu_u - nu)) / (2.0 * G * (1.0 - nu_u) * (1.0 - nu))) * self.F2(z_star, t_star))
        #     t_track += 1

        return displacement

    def pressure(self, locs):
        """Compute pressure field at locations.
        """
        (npts, dim) = locs.shape
        ntpts = tsteps.shape[0]
        pressure = numpy.zeros((ntpts, npts, 1), dtype=numpy.float64)
        x = locs[:, 0]
        y = locs[:, 1]

        t_track = 0

        for t in tsteps:
            u = (r**2 * (1 / M) / (4*t)
            s = -Q/(4*numpy.pi* (k/(mu_f*b)) ) * (numpy.euler_gamma + numpy.log(u))
            pressure[t_track, :, 0] = P_0 - s
            t_track += 1

        return pressure

    # Series functions


    def strain(self, locs):
        """Compute strain field at locations.
        """
        (npts, dim) = locs.shape
        ntpts = tsteps.shape[0]
        e_xx = 0.0
        e_yy = self.trace_strain(locs)
        e_zz = 0.0
        e_xy = 0.0

        strain = numpy.zeros((ntpts, npts, self.TENSOR_SIZE), dtype=numpy.float64)
        strain[:, :, 0] = exx
        strain[:, :, 1] = eyy
        strain[:, :, 2] = ezz
        strain[:, :, 3] = exy
        return strain

    def initial_pressure(self, locs):
        """Compute initial pressure at locations
        """
        (npts, dim) = locs.shape
        pressure = numpy.zeros((1, npts), dtype=numpy.float64)
        z = locs[:, 1]

        pressure[0, :] = P_0

        return pressure

    def initial_trace_strain(self, locs):
        """Compute initial trace strain field at locations.
        """
        (npts, dim) = locs.shape

        trace_strain = numpy.zeros((1, npts), dtype=numpy.float64)
        z = locs[:, 1]
        z_star = z / L

        trace_strain[0, :] = -(P_0 * (1.0 - 2.0 * nu_u)) / (2.0 * G * (1.0 - nu_u))

        return trace_strain


# End of file
