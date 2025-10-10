# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================
# @file tests/fullscale/poroelasticity/cryer/cryer_soln.py
#
# @brief Analytical solution to Cryer's problem.
# Owing to the symmetry of the problem, we only need consider the quarter
# domain case.
#
# This is based on Cheng (Poroelasticity, Section 7.10) and its implementation in PETSc tutorial
# src/ts/tutorials/ex53.c.
#
# Notes:
# - The traction loading is impulsive, so we use a custom PETSc TS (pylith/utils/TSAdaptImpulse),
#   which uses a very small time step for the first step before using the user-specified time step.
# - The accuracy of the solution is poor in the first few time steps due to the impulsive loading
#   so we only check the last time step in the simulation, which is more accurate and select a
#   a tolerance appropriate for the discretization size. 
#
# Dirichlet boundary conditions
#   Ux(0,y,z) = 0
#   Uy(x,0,z) = 0
#   Uy(x,y,0) = 0
# Neumann boundary conditions
#   \tau_normal(rmax) = -1*Pa

import numpy

# Physical properties
G = 30.0e+9
rho_s = 2500
rho_f = 1000
K_fl = 80.0e+9
K_sg = 100.0e+9
K_d = 40.0e+9
alpha = 0.6
phi = 0.1
k = 1.5e-13
mu_f = 0.001
P_0 = 10.0e+3
R_0 = 10.0e+3

M = 1.0 / ( phi / K_fl + (alpha - phi) /K_sg)
kappa = k/mu_f
K_u = K_d + alpha*alpha*M
S = (3*K_u + 4*G) / (M*(3*K_d + 4*G)) #(1/M) + ( (3*alpha*alpha) / (3*K_d + 4*G) )#
c = kappa / S
nu = (3*K_d - 2*G) / (2*(3*K_d + G))
nu_u = (3*K_u - 2*G) / (2*(3*K_u + G))
U_R_inf = -(P_0*R_0*(1-2*nu))/(2*G*(1+nu))
eta = (alpha*(1-2*nu))/(2*(1-nu))

dt0 = 1.0e-6 * 1.0e+8
dt = 5.0e+4
t_end = 100.0e+4
time = numpy.concatenate((numpy.array([dt0]), dt0+numpy.arange(dt, t_end+0.5*dt, dt)))
tsteps = numpy.array([time[-1]])

# ----------------------------------------------------------------------
class AnalyticalSoln(object):
    """
    Analytical solution to Cryer's problem
    """
    SPACE_DIM = 3
    TENSOR_SIZE = 4
    ITERATIONS = 100

    def __init__(self):
        self.fields = {
            "displacement": self.displacement,
            "pressure": self.pressure,
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
                "bc_yneg": self.zero_vector,
                "bc_zneg": self.zero_vector,
                "bc_shell_traction": self.surface_traction,
                "bc_shell_pressure": self.zero_scalar
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
        """
        Compute solid_density field at locations.
        """
        (npts, dim) = locs.shape
        solid_density = rho_s * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return solid_density

    def fluid_density(self, locs):
        """
        Compute fluid density field at locations.
        """
        (npts, dim) = locs.shape
        fluid_density = rho_f * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return fluid_density

    def porosity(self, locs):
        """
        Compute solid_density field at locations.
        """
        (npts, dim) = locs.shape
        porosity = phi * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return porosity

    def shear_modulus(self, locs):
        """
        Compute shear modulus field at locations.
        """
        (npts, dim) = locs.shape
        shear_modulus = G * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return shear_modulus

    def fluid_viscosity(self, locs):
        """
        Compute fluid_viscosity field at locations.
        """
        (npts, dim) = locs.shape
        fluid_viscosity = mu_f * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return fluid_viscosity

    def drained_bulk_modulus(self, locs):
        """
        Compute undrained bulk modulus field at locations.
        """
        (npts, dim) = locs.shape
        drained_bulk_modulus = K_d * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return drained_bulk_modulus

    def biot_coefficient(self, locs):
        """
        Compute biot coefficient field at locations.
        """
        (npts, dim) = locs.shape
        biot_coefficient = alpha * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return biot_coefficient

    def biot_modulus(self, locs):
        """
        Compute biot modulus field at locations.
        """
        (npts, dim) = locs.shape
        biot_modulus = M * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return biot_modulus

    def isotropic_permeability(self, locs):
        """
        Compute isotropic permeability field at locations.
        """
        (npts, dim) = locs.shape
        isotropic_permeability = k * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return isotropic_permeability

    def displacement(self, locs):
        """
        Compute displacement field at locations.
        """
        (npts, dim) = locs.shape
        ntpts = tsteps.shape[0]
        displacement = numpy.zeros((ntpts, npts, dim), dtype=numpy.float64)

        x_n = self.cryerZeros()
        R = numpy.sqrt(locs[:,0]*locs[:,0] + locs[:,1]*locs[:,1] + locs[:,2]*locs[:,2])
        theta = numpy.nan_to_num( numpy.arctan( numpy.nan_to_num( numpy.sqrt(locs[:,0]**2 + locs[:,1]**2) / locs[:,2] ) ) )
        R_star = R.reshape((R.size,1)) / R_0

        x_n.reshape((1,x_n.size))
        sqrt_xn = numpy.sqrt(x_n)

        E = (1-nu)**2*(1+nu_u)**2*x_n - 18*(1+nu)*(nu_u-nu)*(1-nu_u)

        for i_t, t in enumerate(tsteps):
            t_star = (c*t)/(R_0**2)
            r_exact_N =  R_star.ravel() - numpy.nan_to_num(numpy.sum(((12*(1 + nu)*(nu_u - nu)) / \
                                        ((1 - 2*nu)*E*R_star**2*x_n*numpy.sin(sqrt_xn)) ) * \
                                        (3*(nu_u - nu) * (numpy.sin(R_star*sqrt_xn) - R_star*sqrt_xn*numpy.cos(R_star*sqrt_xn)) + \
                                        (1 - nu)*(1 - 2*nu)*R_star**3*x_n*numpy.sin(sqrt_xn)) * \
                                        numpy.exp(-x_n*t_star),axis=1))

            displacement[i_t, :, 0] = (r_exact_N*U_R_inf)*numpy.cos(phi)*numpy.sin(theta)
            displacement[i_t, :, 1] = (r_exact_N*U_R_inf)*numpy.sin(phi)*numpy.sin(theta)
            displacement[i_t, :, 2] = (r_exact_N*U_R_inf)*numpy.cos(theta)

        return displacement

    def pressure(self, locs):
        """
        Compute pressure field at locations.
        """
        (npts, dim) = locs.shape
        ntpts = tsteps.shape[0]
        pressure = numpy.zeros((ntpts, npts, 1), dtype=numpy.float64)

        center = numpy.where(~locs.any(axis=1))[0]

        x_n = self.cryerZeros()
        sqrt_xn = numpy.sqrt(x_n)

        R = numpy.sqrt(locs[:,0]*locs[:,0] + locs[:,1]*locs[:,1] + locs[:,2]*locs[:,2])
        R_star = R.reshape([R.size,1]) / R_0
        x_n.reshape([1,x_n.size])

        E = (1-nu)**2 * (1+nu_u)**2 * x_n - 18*(1+nu)*(nu_u-nu)*(1-nu_u)

        for i_t, t in enumerate(tsteps):

            t_star = (c*t)/(R_0**2)
            pressure[i_t,:,0] = P_0*numpy.sum( ( (18*numpy.square(nu_u-nu) ) / (eta*E) ) * \
                                        ( (numpy.sin(R_star*sqrt_xn)) / (R_star*numpy.sin(sqrt_xn) ) - 1 ) * \
                                        numpy.exp(-x_n*t_star) , axis=1)

            pressure[i_t,center,0] = P_0*numpy.sum( ( (18*numpy.square(nu_u-nu) ) / (eta*E) ) * \
                                        ( sqrt_xn / (numpy.sin(sqrt_xn) ) - 1 ) * \
                                        numpy.exp(-x_n*t_star))
        return pressure


    # Series functions

    def cryerZeros(self):

        f = lambda x: numpy.tan(numpy.sqrt(x)) - (6*(nu_u - nu)*numpy.sqrt(x))/(6*(nu_u - nu) - (1 - nu)*(1 + nu_u)*x) # Compressible Constituents

        n_series = self.ITERATIONS
        a_n = numpy.zeros(n_series) # initializing roots array
        xtol = 1e-30
        rtol = 1e-15
        for i in range(1,n_series+1):
            a = numpy.square(i*numpy.pi) - (i+1)*numpy.pi
            b = numpy.square(i*numpy.pi) + (i+1)*numpy.pi
            # print('a: ',a)
            # print('b: ',b)
            f_c = 10
            f_c_old = 0
            rtol_flag = False
            c = 0
            it = 0
            while numpy.abs(f_c) > xtol and rtol_flag == False:
                c = (a + b) / 2
                f_c = f(c)

                # print('c: ',c)
                # print('f(c):',f_c)

                if f(a)*f_c < 0:
                    a = a
                    b = c
                elif f(b)*f_c < 0:
                    a = c
                    b = b
                else:
                    print('Bisection method failed')
                    # print('a: ',a)
                    # print('b: ',b)
                    break
                if numpy.abs(numpy.abs(f_c_old) - numpy.abs(f_c)) < rtol:
                    rtol_flag = True
                it += 1
                # print('f(c): ',f_c)
                # print('rtol: ',numpy.abs(numpy.abs(f_c_old) - numpy.abs(f_c)))
                # print('rtol flag: ',rtol_flag)
                f_c_old = f_c
            # print('n: ',i)
            # print('c: ',c)
            # print('f(c):',f_c)
            # print('iter: ',it)
            a_n[i-1] = c

        return(a_n)

    def surface_traction(self, locs):
        (npts, dim) = locs.shape
        traction = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        traction[:,:,-1] = -P_0

        return traction


# End of file
