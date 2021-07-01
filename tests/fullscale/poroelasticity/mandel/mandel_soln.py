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
# @file tests/fullscale/poroelasticity/mandel/mandel_soln.py
#
# @brief Analytical solution to Mandel's problem.
# Owing to the symmetry of the problem, we only need consider the quarter
# domain case.
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
mu_f = 1.0  # Pa*s
G = 3.0  # Pa
K_sg = 10.0  # Pa
K_fl = 8.0  # Pa
K_d = 4.0  # Pa
# K_u = 2.6941176470588233 # Pa
alpha = 0.6  # -
phi = 0.1
# M = 4.705882352941176 # Pa
k = 1.5  # m**2

xmin = 0.0  # m
xmax = 10.0  # m
ymin = 0.0  # m
ymax = 1.0  # m


# Height of column, m
a = (xmax - xmin)
b = (ymax - ymin)


vertical_stress = 1.0  # Pa
F = vertical_stress*a

M = 1.0 / (phi / K_fl + (alpha - phi) / K_sg)  # Pa
K_u = K_d + alpha * alpha * M  # Pa,      Cheng (B.5)
# K_d = K_u - alpha*alpha*M # Pa,      Cheng (B.5)
nu = (3.0 * K_d - 2.0 * G) / (2.0 * (3.0 * K_d + G))  # -,       Cheng (B.8)
nu_u = (3.0 * K_u - 2.0 * G) / (2.0 * (3.0 * K_u + G))  # -,       Cheng (B.9)
eta = (3.0 * alpha * G) / (3.0 * K_d + 4.0 * G)  # -,       Cheng (B.11)
S = (3.0 * K_u + 4.0 * G) / (M * (3.0 * K_d + 4.0 * G))  # Pa^{-1}, Cheng (B.14)
c = (k / mu_f) / S  # m^2 / s, Cheng (B.16)
B = (3. * (nu_u - nu)) / (alpha * (1. - 2. * nu) * (1. + nu_u))

# Time steps
ts = 0.0028666667  # sec
nts = 2
tsteps = numpy.arange(0.0, ts * nts, ts) + ts  # sec


# ----------------------------------------------------------------------
class AnalyticalSoln(object):
    """Analytical solution to Mandel's problem
    """
    SPACE_DIM = 2
    TENSOR_SIZE = 4
    ITERATIONS = 300
    EPS = 1e-25

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
            "undrained_bulk_modulus": self.undrained_bulk_modulus,
            "drained_bulk_modulus": self.drained_bulk_modulus,
            "biot_coefficient": self.biot_coefficient,
            "biot_modulus": self.biot_modulus,
            "isotropic_permeability": self.isotropic_permeability,
            "initial_amplitude": {
                "x_neg": self.zero_vector,
                "x_pos": self.zero_scalar,
                "y_neg": self.zero_vector,
                "y_pos": self.initial_displacement,
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

    def undrained_bulk_modulus(self, locs):
        """Compute undrained bulk modulus field at locations.
        """
        (npts, dim) = locs.shape
        undrained_bulk_modulus = K_u * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return undrained_bulk_modulus

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

    def y_pos(self, locs):
        """Compute initial traction at locations.

        :TODO: If this is the initial traction, then it should be a single time point (0).
        """
        (npts, dim) = locs.shape
        ntpts = tsteps.shape[0]
        t_track = 0

        displacement = numpy.zeros((ntpts, npts, self.SPACE_DIM), dtype=numpy.float64)
        displacement[:, :, 0] = 0.0
        for t in tsteps:
            displacement[t_track, :, 1] = F
            t_track += 1
        return traction

    def displacement(self, locs):
        """Compute displacement field at locations.
        """
        (npts, dim) = locs.shape
        ntpts = tsteps.shape[0]
        displacement = numpy.zeros((ntpts, npts, dim), dtype=numpy.float64)
        x = locs[:, 0]
        z = locs[:, 1]
        t_track = 0
        zeroArray = self.mandelZeros()

        for t in tsteps:
            A_x = 0.0
            B_x = 0.0

            for n in numpy.arange(1, self.ITERATIONS + 1, 1):
                a_n = zeroArray[n - 1]
                A_x += (numpy.sin(a_n) * numpy.cos(a_n) / (a_n - numpy.sin(a_n) * numpy.cos(a_n))) * \
                    numpy.exp(-1.0 * (a_n * a_n * c * t) / (a * a))
                B_x += (numpy.cos(a_n) / (a_n - numpy.sin(a_n) * numpy.cos(a_n))) * \
                    numpy.sin((a_n * x) / a) * numpy.exp(-1.0 * (a_n * a_n * c * t) / (a * a))

            displacement[t_track, :, 0] = ((F * nu) / (2.0 * G * a) - (F * nu_u) / (G * a) * A_x) * x + F / G * B_x
            displacement[t_track, :, 1] = (-1 * (F * (1.0 - nu)) / (2 * G * a) + (F * (1 - nu_u)) / (G * a) * A_x) * z
            t_track += 1

        return displacement

    def pressure(self, locs):
        """Compute pressure field at locations.
        """
        (npts, dim) = locs.shape
        ntpts = tsteps.shape[0]
        pressure = numpy.zeros((ntpts, npts, 1), dtype=numpy.float64)
        x = locs[:, 0]
        z = locs[:, 1]
        t_track = 0
        zeroArray = self.mandelZeros()

        for t in tsteps:

            if t == 0.0:
                pressure[t_track, :] = (1. / (3. * a)) * (B * (1. + nu_u)) * F
            else:
                p = 0.0
                for n in numpy.arange(1, self.ITERATIONS + 1, 1):
                    x_n = zeroArray[n - 1]
                    p += (numpy.sin(x_n) / (x_n - numpy.sin(x_n) * numpy.cos(x_n))) * \
                        (numpy.cos((x_n * x) / a) - numpy.cos(x_n)) * numpy.exp(-1.0 * (x_n * x_n * c * t) / (a * a))
                pressure[t_track, :, 0] = ((2.0 * (F * B * (1.0 + nu_u))) / (3.0 * a)) * p
            t_track += 1

        return pressure

    def trace_strain(self, locs):
        """Compute trace strain field at locations.
        """
        (npts, dim) = locs.shape
        ntpts = tsteps.shape[0]
        trace_strain = numpy.zeros((ntpts, npts, 1), dtype=numpy.float64)
        x = locs[:, 0]
        z = locs[:, 1]
        t_track = 0
        zeroArray = self.mandelZeros()

        for t in tsteps:

            eps_A = 0.0
            eps_B = 0.0
            eps_C = 0.0

            for i in numpy.arange(1, self.ITERATIONS+1,1):
                x_n = zeroArray[i-1]
                eps_A += (x_n * numpy.exp( (-1.0*x_n*x_n*c*t)/(a*a)) * numpy.cos(x_n)*numpy.cos( (x_n*x)/a)) / (a * (x_n - numpy.sin(x_n)*numpy.cos(x_n)))
                eps_B += ( numpy.exp( (-1.0*x_n*x_n*c*t)/(a*a)) * numpy.sin(x_n)*numpy.cos(x_n)) / (x_n - numpy.sin(x_n)*numpy.cos(x_n))
                eps_C += ( numpy.exp( (-1.0*x_n*x_n*c*t)/(x_n*x_n)) * numpy.sin(x_n)*numpy.cos(x_n)) / (x_n - numpy.sin(x_n)*numpy.cos(x_n))

            trace_strain[t_track,:,0] = (F/G)*eps_A + ( (F*nu)/(2.0*G*a)) - eps_B/(G*a) - (F*(1.0-nu))/(2/0*G*a) + eps_C/(G*a)
            t_track += 1

        return trace_strain

    # Series functions

    def mandelZeros(self):
        """Compute roots for analytical mandel problem solutions
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
                if (numpy.abs(y2) < self.EPS):
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
        p_poisson_ratio = (3 * p_K_d - 2 * p_G) / (2 * (3 * p_K_d + p_G))
        trace_strain = self.trace_strain(locs)
        pressure = self.pressure(locs)
        e_xx = 0.0
        e_yy = self.trace_strain(locs)
        e_xy = 0.0

        stress = numpy.zeros((ntpts, npts, self.TENSOR_SIZE), dtype=numpy.float64)
        stress[:, :, 0] = ((2 * p_G * p_poisson_ratio) / (1 - 2 * p_poisson_ratio)) * \
            trace_strain + 2 * p_G * e_xx - p_alpha * pressure
        stress[:, :, 1] = ((2 * p_G * p_poisson_ratio) / (1 - 2 * p_poisson_ratio)) * \
            trace_strain + 2 * p_G * e_yy - p_alpha * pressure
        stress[:, :, 2] = ((2 * p_G * p_poisson_ratio) / (1 - 2 * p_poisson_ratio)) * trace_strain - p_alpha * pressure
        stress[:, :, 3] = 2 * p_G * e_xy
        return stress

    def initial_traction(self, locs):
        """Compute traction at locations.

        :TODO: If this is the initial traction, then it should be a single time point (0).
        """
        (npts, dim) = locs.shape
        ntpts = tsteps.shape[0]
        traction = numpy.zeros((ntpts, npts, self.SPACE_DIM), dtype=numpy.float64)
        x = locs[:, 0]
        z = locs[:, 1]
        t_track = 0
        zeroArray = self.mandelZeros()

        for t in tsteps:

            sigma_zz_A = 0.0
            sigma_zz_B = 0.0

            for i in np.arange(1, self.ITERATIONS + 1, 1):
                x_n = zeroArray[i - 1]
                sigma_zz_A += (numpy.sin(x_n) / (x_n - numpy.sin(x_n) * numpy.cos(x_n))) * \
                    numpy.cos((x_n * x) / a) * numpy.exp(-1.0 * (x_n * x_n * c * t) / (a * a))
                sigma_zz_B += ((numpy.sin(x_n) * numpy.cos(x_n)) / (x_n - numpy.sin(x_n) *
                                                                    numpy.cos(x_n))) * numpy.exp(-1.0 * (x_n * x_n * c * t) / (a * a))

            traction[t_track, :, 0] = 0.0
            traction[t_track, :, 1] = -(F / a) - ((2.0 * F * (nu_u - nu)) / (a * (1.0 - nu))
                                                  ) * sigma_zz_A + ((2.0 * F) / a) * sigma_zz_B
            t_track += 1

        return traction

    def initial_displacement(self, locs):
        """Compute initial displacement at locations
        """
        (npts, dim) = locs.shape
        displacement = numpy.zeros((1, npts, dim), dtype=numpy.float64)
        x = locs[:, 0]
        z = locs[:, 1]

        displacement[0, :, 0] = 0.0  # (F*nu_u*x)/(2.*G*a)
        displacement[0, :, 1] = -1.0  # -1.*(F*(1.-nu_u)*z)/(2.*G*a)

        return displacement

    def initial_pressure(self, locs):
        """Compute initial pressure at locations
        """
        (npts, dim) = locs.shape
        pressure = numpy.zeros((1, npts), dtype=numpy.float64)
        x = locs[:, 0]
        z = locs[:, 1]
        t = 0.0

        pressure[0, :] = (1. / (3. * a)) * B * (1. + nu_u) * F

        return pressure

    def initial_trace_strain(self, locs):
        """Compute initial trace strain field at locations.
        """
        (npts, dim) = locs.shape
        zeroArray = self.mandelZeros()
        trace_strain = numpy.zeros((1, npts), dtype=numpy.float64)
        x = locs[:, 0]
        z = locs[:, 1]
        t = 0.0

        trace_strain[0, :] = 0.0

        return trace_strain

    def sigma_zz(self, locs):
        """Compute traction at locations.
        """
        (npts, dim) = locs.shape
        ntpts = tsteps.shape[0]
        traction = numpy.zeros((ntpts, npts, self.SPACE_DIM), dtype=numpy.float64)
        x = locs[:, 0]
        z = locs[:, 1]
        t_track = 0
        zeroArray = self.mandelZeros()

        for t in tsteps:

            sigma_zz_A = 0.0
            sigma_zz_B = 0.0

            for i in np.arange(1, self.ITERATIONS + 1, 1):
                x_n = zeroArray[i - 1]
                sigma_zz_A += (numpy.sin(x_n) / (x_n - numpy.sin(x_n) * numpy.cos(x_n))) * \
                    numpy.cos((x_n * x) / a) * numpy.exp(-1.0 * (x_n * x_n * c * t) / (a * a))
                sigma_zz_B += ((numpy.sin(x_n) * numpy.cos(x_n)) / (x_n - numpy.sin(x_n) *
                                                                    numpy.cos(x_n))) * numpy.exp(-1.0 * (x_n * x_n * c * t) / (a * a))

            traction[t_track, :, 0] = 0.0
            traction[t_track, :, 1] = -(F / a) - ((2.0 * F * (nu_u - nu)) / (a * (1.0 - nu))
                                                  ) * sigma_zz_A + ((2.0 * F) / a) * sigma_zz_B
            t_track += 1

        return traction


# End of file
