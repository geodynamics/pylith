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
# @file tests/fullscale/viscoelasticity/nofaults-2d/axialstrainrate_genmaxwell_soln.py
#
# @brief Analytical solution to axial strain rate relaxation problem for a
# generalized Maxwell viscoelastic material.
#
# 2-D axial strain rate solution for linear generalized Maxwell viscoelastic material.
#
#             Uy=0
#          ----------
#          |        |
# Ux=0     |        |  Ux = U0 + V0*t
#          |        |
#          |        |
#          ----------
#            Uy=0
#
# Dirichlet boundary conditions
# Ux(-4000,y) = 0
# Ux(+4000,y) = U0 + V0*t
# Uy(x,-4000) = 0
# Uy(x,+4000) = 0
#

import numpy
import math
from pyre.units.time import year

# Physical properties.
p_density = 2500.0
p_vs = 3464.1016
p_vp = 6000.0
p_viscosity_1 = 9.46728e17
p_viscosity_2 = 4.73364e17
p_viscosity_3 = 1.893456e18
p_shear_ratio_1 = 0.25
p_shear_ratio_2 = 0.25
p_shear_ratio_3 = 0.25

# Stress and strain components, and number of Maxwell elements.
numMaxElements = 3
numComponents = 4

# Initial displacement (1 m) and velocity (2 m/year).
U0 = 1.0
V0 = 2.0/year.value

# Derived properties.
p_mu = p_density*p_vs*p_vs
p_lambda = p_density*p_vp*p_vp - 2.0*p_mu
p_youngs = p_mu*(3.0*p_lambda + 2.0*p_mu)/(p_lambda + p_mu)
p_poissons = 0.5*p_lambda/(p_lambda + p_mu)
p_shear_ratio_0 = 1.0 - p_shear_ratio_1 - p_shear_ratio_2 - p_shear_ratio_3
p_tau_1 = p_viscosity_1/(p_mu*p_shear_ratio_1)
p_tau_2 = p_viscosity_2/(p_mu*p_shear_ratio_2)
p_tau_3 = p_viscosity_3/(p_mu*p_shear_ratio_3)
p_shear_ratio = numpy.array([p_shear_ratio_1, p_shear_ratio_2, p_shear_ratio_3], dtype=numpy.float64)

# Time information.
dt = 0.025*year.value
startTime = dt
endTime = 0.5*year.value
numSteps = 20
t = numpy.linspace(startTime, endTime, num=numSteps, dtype=numpy.float64)

# Uniform strain rate field (plane strain).
e0 = U0/8000.0
edot0 = V0/8000.0
exx = e0 + edot0*t
eyy = numpy.zeros(numSteps, dtype=numpy.float64)
ezz = numpy.zeros(numSteps, dtype=numpy.float64)
exy = numpy.zeros(numSteps, dtype=numpy.float64)

# Deviatoric strains at t=0.
eMean0 = e0/3.0
eDev0 = numpy.zeros(numComponents, dtype=numpy.float64)
eDev0[0] = e0 - eMean0
eDev0[1] = -eMean0
eDev0[2] = -eMean0
eDev0[3] = 0.0

# Deviatoric strain rates.
eMeanDot0 = edot0/3.0
eDevDot0 = numpy.zeros(numComponents, dtype=numpy.float64)
eDevDot0[0] = edot0 - eMeanDot0
eDevDot0[1] = -eMeanDot0
eDevDot0[2] = -eMeanDot0
eDevDot0[3] = 0.0

# Deviatoric strains.
eMean = (exx + eyy + ezz)/3.0
eDev = numpy.zeros((numSteps, numComponents), dtype=numpy.float64)
eDev[:, 0] = exx - eMean
eDev[:, 1] = eyy - eMean
eDev[:, 2] = ezz - eMean
eDev[:, 3] = exy

# Loop over time steps.
eVis = numpy.zeros((numSteps, numComponents, numMaxElements), dtype=numpy.float64)
sDev = numpy.zeros((numSteps, numComponents), dtype=numpy.float64)

for timeStep in range(numSteps):
    timeFac1 = math.exp(-t[timeStep]/p_tau_1)
    timeFac2 = math.exp(-t[timeStep]/p_tau_2)
    timeFac3 = math.exp(-t[timeStep]/p_tau_3)
    # Viscous strains.
    eVis[timeStep,:, 0] = eDev0*timeFac1 + eDevDot0*p_tau_1*(1.0 - timeFac1)
    eVis[timeStep,:, 1] = eDev0*timeFac2 + eDevDot0*p_tau_2*(1.0 - timeFac2)
    eVis[timeStep,:, 2] = eDev0*timeFac3 + eDevDot0*p_tau_3*(1.0 - timeFac3)

    # Deviatoric stresses.
    sDev[timeStep,:] = p_shear_ratio_0*eDev[timeStep,:]
    for elementNum in range(numMaxElements):
        sDev[timeStep,:] += p_shear_ratio[elementNum]*eVis[timeStep,:, elementNum]

sDev *= 2.0*p_mu

# Total stresses.
sMean = eMean*(3.0*p_lambda + 2.0*p_mu)
sxx = sDev[:, 0] + sMean
syy = sDev[:, 1] + sMean
szz = sDev[:, 2] + sMean
sxy = sDev[:, 3]

# ----------------------------------------------------------------------
class AnalyticalSoln(object):
    """
    Analytical solution to axial extension problem.
    """
    SPACE_DIM = 2
    TENSOR_SIZE = 4

    def __init__(self):
        self.fields = {
            "displacement": self.displacement,
            "density": self.density,
            "shear_modulus": self.shear_modulus,
            "bulk_modulus": self.bulk_modulus,
            "shear_modulus_ratio": self.shear_modulus_ratio,
            "maxwell_time": self.maxwell_time,
            "cauchy_strain": self.strain,
            "cauchy_stress": self.stress,
            "viscous_strain": self.viscous_strain,
            "initial_amplitude": {
                "bc_xneg": self.initial_displacement,
                "bc_yneg": self.initial_displacement,
                "bc_xpos": self.initial_displacement,
                "bc_ypos": self.initial_displacement
                },
            "rate_start_time": {
                "bc_xneg": self.bc_rate_time,
                "bc_yneg": self.bc_rate_time,
                "bc_xpos": self.bc_rate_time,
                "bc_ypos": self.bc_rate_time
                },
            "rate_amplitude": {
                "bc_xneg": self.bc_velocity,
                "bc_yneg": self.bc_velocity,
                "bc_xpos": self.bc_velocity,
                "bc_ypos": self.bc_velocity
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
        disp = numpy.zeros((numSteps, npts, self.SPACE_DIM), dtype=numpy.float64)
        disp[:,:, 0] = numpy.dot(exx.reshape(numSteps, 1), (locs[:, 0] + 4000.0).reshape(1, npts))
        return disp

    def initial_displacement(self, locs):
        """
        Compute initial displacement field at locations.
        """
        (npts, dim) = locs.shape
        disp = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        disp[0,:, 0] = e0*(locs[:, 0] + 4000.0).reshape(1, npts)
        return disp

    def bc_velocity(self, locs):
        """
        Compute velocity field at locations.
        """
        (npts, dim) = locs.shape
        velocity = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        velocity[0,:, 0] = edot0*(locs[:, 0] + 4000.0).reshape(1, npts)
        return velocity

    def bc_rate_time(self, locs):
        """
        Compute rate start time at locations.
        """
        (npts, dim) = locs.shape
        rate_start_time = numpy.zeros((1, npts, 1), dtype=numpy.float64)
        return rate_start_time

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

    def maxwell_time(self, locs):
        """
        Compute Maxwell time field at locations.
        """
        (npts, dim) = locs.shape
        maxwell_time = numpy.zeros((1, npts, 3), dtype=numpy.float64)
        maxwell_time[0,:, 0] = p_tau_1
        maxwell_time[0,:, 1] = p_tau_2
        maxwell_time[0,:, 2] = p_tau_3
        return maxwell_time

    def shear_modulus_ratio(self, locs):
        """
        Compute shear modulus ratio field at locations.
        """
        (npts, dim) = locs.shape
        shear_modulus_ratio = numpy.zeros((1, npts, 3), dtype=numpy.float64)
        shear_modulus_ratio[0,:, 0] = p_shear_ratio_1
        shear_modulus_ratio[0,:, 1] = p_shear_ratio_2
        shear_modulus_ratio[0,:, 2] = p_shear_ratio_3
        return shear_modulus_ratio

    def strain(self, locs):
        """
        Compute strain field at locations.
        """
        (npts, dim) = locs.shape
        strain = numpy.zeros((numSteps, npts, self.TENSOR_SIZE), dtype=numpy.float64)
        strain[:,:, 0] = exx.reshape(numSteps, 1)
        strain[:,:, 1] = eyy.reshape(numSteps, 1)
        strain[:,:, 2] = ezz.reshape(numSteps, 1)
        strain[:,:, 3] = exy.reshape(numSteps, 1)
        return strain

    def stress(self, locs):
        """
        Compute stress field at locations.
        """
        (npts, dim) = locs.shape
        stress = numpy.zeros((numSteps, npts, self.TENSOR_SIZE), dtype=numpy.float64)
        stress[:,:, 0] = sxx.reshape(numSteps, 1)
        stress[:,:, 1] = syy.reshape(numSteps, 1)
        stress[:,:, 2] = szz.reshape(numSteps, 1)
        stress[:,:, 3] = sxy.reshape(numSteps, 1)
        return stress

    def viscous_strain(self, locs):
        """
        Compute viscous strain field at locations.
        """
        (npts, dim) = locs.shape
        viscous_strain = numpy.zeros((numSteps, npts, 3*self.TENSOR_SIZE), dtype=numpy.float64)
        viscous_strain[:,:, 0] = eVis[:, 0, 0].reshape(numSteps, 1)
        viscous_strain[:,:, 1] = eVis[:, 1, 0].reshape(numSteps, 1)
        viscous_strain[:,:, 2] = eVis[:, 2, 0].reshape(numSteps, 1)
        viscous_strain[:,:, 3] = eVis[:, 3, 0].reshape(numSteps, 1)
        viscous_strain[:,:, 4] = eVis[:, 0, 1].reshape(numSteps, 1)
        viscous_strain[:,:, 5] = eVis[:, 1, 1].reshape(numSteps, 1)
        viscous_strain[:,:, 6] = eVis[:, 2, 1].reshape(numSteps, 1)
        viscous_strain[:,:, 7] = eVis[:, 3, 1].reshape(numSteps, 1)
        viscous_strain[:,:, 8] = eVis[:, 0, 2].reshape(numSteps, 1)
        viscous_strain[:,:, 9] = eVis[:, 1, 2].reshape(numSteps, 1)
        viscous_strain[:,:, 10] = eVis[:, 2, 2].reshape(numSteps, 1)
        viscous_strain[:,:, 11] = eVis[:, 3, 2].reshape(numSteps, 1)
        return viscous_strain


# End of file
