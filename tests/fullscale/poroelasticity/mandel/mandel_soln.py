# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================
# @file tests/fullscale/poroelasticity/mandel/mandel_soln.py
#
# @brief Analytical solution to Mandel's problem.
# Owing to the symmetry of the problem, we only need consider the quarter
# domain case.
#
# This is based on Cheng (Poroelasticity, Section 7.4) and its implementation in PETSc tutorial
# src/ts/tutorials/ex53.c.
#
# Notes:
# - The surface loading is impulsive, so we use a custom PETSc TS (pylith/utils/TSAdaptImpulse),
#   which uses a very small time step for the first step before using the user-specified time step.
# - The accuracy of the solution is poor in the first few time steps due to the impulsive loading
#   so we only check the last time step in the simulation, which is more accurate and select a
#   a tolerance appropriate for the discretization size.
# - We use a Dirichlet BC for displacement on the top surface, because the boundary value problem
#   assumes loading via a rigid plate (uniform displacement) with a total load that is not uniform.
#
#           -F
#        ----------
#        |        |
#  Ux=0  |        | P=0
#        |        |
#        |        |
#        ----------
#          Uy=0
#
# Dirichlet boundary conditions
#   Ux(0,y) = 0
#   Uy(x,0) = 0
# Neumann boundary conditions
#   \tau_normal(x,ymax) = -1*Pa

import numpy

# Physical properties
rho_s = 2500  # kg / m**3
rho_f = 1000  # kg / m**3
mu_f = 1.0e-3  # Pa*s
G = 30.0e+9  # Pa
K_sg = 100.0e+9  # Pa
K_fl = 80.0e+9  # Pa
K_d = 40.0e+9  # Pa
alpha = 0.6
phi = 0.1
k = 1.5e-13  # m**2

x_min = 0.0e+3  # m
x_max = 8.0e+3  # m
y_min = 0.0
y_max = 1.0e+3
l_x = x_max - x_min

P_0 = 10.0e+3  # Pa
F = l_x*P_0

dt0 = 1.0e-6 * 1.0e+8
dt = 5.0e+4
t_end = 100.0e+4
time = numpy.concatenate((numpy.array([dt0]), dt0+numpy.arange(dt, t_end+0.5*dt, dt)))
tsteps = numpy.array([time[-1]])


M = 1.0 / (phi / K_fl + (alpha - phi) / K_sg)  # Pa
K_u = K_d + alpha * alpha * M  # Pa,      Cheng (B.5)
nu = (3.0 * K_d - 2.0 * G) / (2.0 * (3.0 * K_d + G))  # -,       Cheng (B.8)
nu_u = (3.0 * K_u - 2.0 * G) / (2.0 * (3.0 * K_u + G))  # -,       Cheng (B.9)
eta = (3.0 * alpha * G) / (3.0 * K_d + 4.0 * G)  # -,       Cheng (B.11)
S = (3.0 * K_u + 4.0 * G) / (M * (3.0 * K_d + 4.0 * G))  # Pa^{-1}, Cheng (B.14)
c = (k / mu_f) / S  # m^2 / s, Cheng (B.16)
B = (3. * (nu_u - nu)) / (alpha * (1. - 2. * nu) * (1. + nu_u))


# ----------------------------------------------------------------------
class AnalyticalSoln(object):
    """Analytical solution to Mandel's problem
    """
    SPACE_DIM = 2
    TENSOR_SIZE = 4
    ITERATIONS = 300
    SMALL = 1e-25

    def __init__(self):
        self.fields = {
            "displacement": self.displacement,
            "pressure": self.pressure,
            "trace_strain": self.trace_strain,
            "porosity": self.porosity,
            "solid_density": self.solid_density,
            "fluid_density": self.fluid_density,
            "fluid_viscosity": self.fluid_viscosity,
            "shear_modulus": self.shear_modulus,
            "drained_bulk_modulus": self.drained_bulk_modulus,
            "biot_coefficient": self.biot_coefficient,
            "biot_modulus": self.biot_modulus,
            "isotropic_permeability": self.isotropic_permeability,
            "initial_amplitude": {
                "bc_xneg": self.zero_vector,
                "bc_xpos": self.zero_scalar,
                "bc_yneg": self.zero_vector,
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

    def porosity(self, locs):
        """Compute solid_density field at locations.
        """
        (npts, dim) = locs.shape
        porosity = phi * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return porosity

    def shear_modulus(self, locs):
        """Compute shear modulus field at locations.
        """
        (npts, dim) = locs.shape
        shear_modulus = G * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return shear_modulus

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
        drained_bulk_modulus = K_d * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return drained_bulk_modulus

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
        x = locs[:, 0]
        y = locs[:, 1]
        zeroArray = self.mandelZeros()

        n = numpy.arange(1, self.ITERATIONS+1, 1)
        a_n = zeroArray[n - 1]
        for i_t, t in enumerate(tsteps):
            A_s = numpy.sum((numpy.sin(a_n) * numpy.cos(a_n) / (a_n - numpy.sin(a_n) * numpy.cos(a_n))) * \
                numpy.exp(-1.0 * (a_n * a_n * c * t) / (l_x * l_x)))
            B_s = numpy.zeros(x.shape)
            for i_x, x_pt in enumerate(x):
                a_n = zeroArray[n - 1]
                B_s[i_x] = numpy.sum((numpy.cos(a_n) / (a_n - numpy.sin(a_n) * numpy.cos(a_n))) * \
                    numpy.sin((a_n * x_pt) / l_x) * numpy.exp(-1.0 * (a_n * a_n * c * t) / (l_x * l_x)))

            displacement[i_t, :, 0] = ((F * nu) / (2.0 * G * l_x) - (F * nu_u) / (G * l_x) * A_s) * x + F / G * B_s
            displacement[i_t, :, 1] = (-1 * (F * (1.0 - nu)) / (2.0 * G * l_x) + (F * (1 - nu_u)) / (G * l_x) * A_s) * y
        return displacement

    def pressure(self, locs):
        """Compute pressure field at locations.
        """
        (npts, dim) = locs.shape
        ntpts = tsteps.shape[0]
        pressure = numpy.zeros((ntpts, npts, 1), dtype=numpy.float64)
        x = locs[:, 0]
        zeroArray = self.mandelZeros()

        for i_t, t in enumerate(tsteps):
            n = numpy.arange(1, self.ITERATIONS+1, 1)
            a_n = zeroArray[n-1]
            p = numpy.zeros(x.shape)
            for i_x, x_pt in enumerate(x):
                p[i_x] = numpy.sum((numpy.sin(a_n) / (a_n - numpy.sin(a_n) * numpy.cos(a_n))) * \
                        (numpy.cos((a_n * x_pt) / l_x) - numpy.cos(a_n)) * numpy.exp(-1.0 * (a_n * a_n * c * t) / (l_x * l_x)))
            pressure[i_t, :, 0] = (2 * F * B * (1.0 + nu_u)) / (3.0 * l_x) * p
        return pressure

    def trace_strain(self, locs):
        """Compute trace strain field at locations.
        """
        (npts, dim) = locs.shape
        ntpts = tsteps.shape[0]
        trace_strain = numpy.zeros((ntpts, npts, 1), dtype=numpy.float64)
        x = locs[:, 0]
        zeroArray = self.mandelZeros()

        for i_t, t in enumerate(tsteps):
            n = numpy.arange(1, self.ITERATIONS+1, 1)
            a_n = zeroArray[n-1]
            eps_B = numpy.sum((numpy.exp((-1.0*a_n*a_n*c*t)/(l_x*l_x)) * numpy.sin(a_n)*numpy.cos(a_n)) / (a_n - numpy.sin(a_n)*numpy.cos(a_n)))
            eps_C = numpy.sum((numpy.exp((-1.0*a_n*a_n*c*t)/(a_n*a_n)) * numpy.sin(a_n)*numpy.cos(a_n)) / (a_n - numpy.sin(a_n)*numpy.cos(a_n)))
            eps_A = numpy.zeros(x.shape)
            for i_x, x_pt in enumerate(x):
                eps_A[i_x] = numpy.sum((a_n * numpy.exp((-1.0*a_n*a_n*c*t)/(l_x*l_x)) * numpy.cos(a_n)*numpy.cos( (a_n*x_pt)/l_x)) / (l_x * (a_n - numpy.sin(a_n)*numpy.cos(a_n))))

            trace_strain[i_t, :, 0] = (F/G)*eps_A + ( (F*nu)/(2.0*G*l_x)) - eps_B/(G*l_x) - (F*(1.0-nu))/(2.0*G*l_x) + eps_C/(G*l_x)
        return trace_strain

    # Series functions

    def mandelZeros(self):
        """Compute roots for analytical Mandel problem solutions
        """
        zeroArray = numpy.zeros(self.ITERATIONS)
        x0 = 0

        for i in numpy.arange(1, self.ITERATIONS + 1, 1):
            a1 = x0 + numpy.pi/4
            a2 = x0 + numpy.pi/2 - 10000*2.2204e-16
            am = a1
            for j in numpy.arange(0, self.ITERATIONS, 1):
                y1 = numpy.tan(a1) - ((1.0 - nu) / (nu_u - nu)) * a1
                y2 = numpy.tan(a2) - ((1.0 - nu) / (nu_u - nu)) * a2
                am = (a1 + a2) / 2.0
                ym = numpy.tan(am) - (1 - nu) / (nu_u - nu) * am
                if ((ym * y1) > 0):
                    a1 = am
                else:
                    a2 = am
                if (numpy.abs(y2) < self.SMALL):
                    am = a2
            zeroArray[i - 1] = am
            x0 += numpy.pi
        return zeroArray

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
        strain[:, :, 0] = e_xx
        strain[:, :, 1] = e_yy
        strain[:, :, 2] = e_zz
        strain[:, :, 3] = e_xy
        return strain

    def stress(self, locs):
        """Compute stress field at locations.
        """
        (npts, dim) = locs.shape
        ntpts = tsteps.shape[0]
        p_poisson_ratio = (3 * K_d - 2 * G) / (2 * (3 * K_d + G))
        trace_strain = self.trace_strain(locs)
        pressure = self.pressure(locs)
        e_xx = 0.0
        e_yy = self.trace_strain(locs)
        e_xy = 0.0

        stress = numpy.zeros((ntpts, npts, self.TENSOR_SIZE), dtype=numpy.float64)
        stress[:, :, 0] = ((2 * G * p_poisson_ratio) / (1 - 2 * p_poisson_ratio)) * \
            trace_strain + 2 * G * e_xx - alpha * pressure
        stress[:, :, 1] = ((2 * G * p_poisson_ratio) / (1 - 2 * p_poisson_ratio)) * \
            trace_strain + 2 * G * e_yy - alpha * pressure
        stress[:, :, 2] = ((2 * G * p_poisson_ratio) / (1 - 2 * p_poisson_ratio)) * trace_strain - alpha * pressure
        stress[:, :, 3] = 2 * G * e_xy
        return stress


# End of file
