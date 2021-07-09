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
# @file tests/fullscale/poroelasticity/terzaghi/terzaghi_soln.py
#
# @brief Analytical solution to Terzaghi's problem.
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
mu_f = 1.0  # Pa*s
G = 3.0  # Pa
K_sg = 10.0  # Pa
K_fl = 8.0  # Pa
K_d = 4.0  # Pa
# K_u = 2.6941176470588233 # Pa
alpha = 0.6  # -
phi = 0.1  # -
# M = 4.705882352941176# Pa
k = 1.5  # m**2

ymax = 10.0  # m
ymin = 0.0  # m
xmax = 10.0  # m
xmin = 0.0  # m
P_0 = -1.0  # Pa

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
    """Analytical solution to Terzaghi's problem
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

        for t in tsteps:
            if t < 0.0:
                displacement[0, :, 1] = ((P_0 * L * (1.0 - 2.0 * nu_u)) / (2.0 * G * (1.0 - nu_u))) * (1.0 - z_star)
            else:
                t_star = (c * t) / ((2 * L)**2)
                displacement[t_track, :, 1] =  (((P_0 * L * (1.0 - 2.0 * nu_u)) / (2.0 * G * (1.0 - nu_u))) * (1.0 - z_star) + ((P_0 * L * (nu_u - nu)) / (2.0 * G * (1.0 - nu_u) * (1.0 - nu))) * self.F2(z_star, t_star))
            t_track += 1

        return displacement

    def pressure(self, locs):
        """Compute pressure field at locations.
        """
        (npts, dim) = locs.shape
        ntpts = tsteps.shape[0]
        pressure = numpy.zeros((ntpts, npts, 1), dtype=numpy.float64)
        z = locs[:, 1]
        t_track = 0

        for t in tsteps:
            z_star = 1 - z / L
            t_star = (c * t) / (4. * L**2)
            pressure[t_track, :, 0] = -((P_0 * eta) / (G * S)) * self.F1(z_star, t_star)
            t_track += 1

        return pressure

    def trace_strain(self, locs):
        """Compute trace strain field at locations.
        """
        (npts, dim) = locs.shape
        ntpts = tsteps.shape[0]
        trace_strain = numpy.zeros((ntpts, npts, 1), dtype=numpy.float64)
        z = locs[:, 1]
        t_track = 0

        for t in tsteps:
            z_star = z / L
            t_star = (c * t) / (4 * L**2)
            trace_strain[t_track, :, 0] = -((P_0 * L * (1.0 - 2.0 * nu_u)) / (2.0 * G * (1.0 - nu_u) * L)) \
                + ((P_0 * L * (nu_u - nu)) / (2.0 * G * (1.0 - nu_u) * (1.0 - nu))) * self.F3(z_star, t_star)
            t_track += 1

        return trace_strain

    # Series functions

    def F1(self, z_star, t_star):
        F1 = 0.
        for m in numpy.arange(1, 2 * self.ITERATIONS + 1, 2):
            F1 += 4. / (m * numpy.pi) * numpy.sin(0.5 * m * numpy.pi * z_star) * numpy.exp(-(m * numpy.pi)**2 * t_star)
        return F1

    def F2(self, z_star, t_star):
        F2 = 0.
        for m in numpy.arange(1, 2 * self.ITERATIONS + 1, 2):
            F2 += (8. / (m * numpy.pi)**2) * numpy.cos(0.5 * m * numpy.pi *
                                                       z_star) * (1. - numpy.exp(-(m * numpy.pi)**2 * t_star))
        return F2

    def F3(self, z_star, t_star):
        F3 = 0.
        for m in numpy.arange(1, 2 * self.ITERATIONS + 1, 2):
            F3 += (-4.0 / (m * numpy.pi * L)) * numpy.sin(0.5 * m * numpy.pi *
                                                          z_star) * (1.0 - numpy.exp(-(m * numpy.pi)**2 * t_star))
        return F3

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

    def stress(self, locs):
        """Compute stress field at locations.
        """
        (npts, dim) = locs.shape
        ntpts = tsteps.shape[0]
        poisson_ratio = (3 * K_d - 2 * G) / (2 * (3 * K_d + G))
        trace_strain = self.trace_strain(locs)
        pressure = self.pressure(locs)
        e_xx = 0.0
        e_yy = self.trace_strain(locs)
        e_xy = 0.0

        stress = numpy.zeros((ntpts, npts, self.TENSOR_SIZE), dtype=numpy.float64)
        stress[:, :, 0] = ((2 * G * poisson_ratio) / (1 - 2 * poisson_ratio)) * \
            trace_strain + 2 * G * e_xx - alpha * pressure
        stress[:, :, 1] = ((2 * G * poisson_ratio) / (1 - 2 * poisson_ratio)) * \
            trace_strain + 2 * G * e_yy - alpha * pressure
        stress[:, :, 2] = ((2 * G * poisson_ratio) / (1 - 2 * poisson_ratio)) * trace_strain - alpha * pressure
        stress[:, :, 3] = 2 * G * e_xy
        return stress

    def y_pos_neu(self, locs):
        """Compute initial traction at locations.
        """
        (npts, dim) = locs.shape
        traction = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        traction[:, :, 0] = 0.0
        traction[:, :, 1] = P_0
        return traction

    def initial_displacement(self, locs):
        """Compute initial displacement at locations
        """
        (npts, dim) = locs.shape
        displacement = numpy.zeros((1, npts, dim), dtype=numpy.float64)
        z = locs[:, 1]
        z_star = 1 - z / L

        displacement[0, :, 1] = ((P_0 * L * (1.0 - 2.0 * nu_u)) / (2.0 * G * (1.0 - nu_u))) * (1.0 - z_star)
        return displacement

    def initial_pressure(self, locs):
        """Compute initial pressure at locations
        """
        (npts, dim) = locs.shape
        pressure = numpy.zeros((1, npts), dtype=numpy.float64)
        z = locs[:, 1]

        pressure[0, :] = (-P_0 * eta) / (G * S)

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
