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
# @file tests/fullscale/viscoelasticity/nofaults-3d/axialtraction_maxwell_soln.py
#
# @brief Analytical solution to axial traction problem for a Maxwell viscoelastic material.
#
# 3-D axial traction solution for linear Maxwell viscoelastic material.
#
#             Tz=T0
#          ----------
#          |        |
# Ux=0     |        |  Ux=0
#          |        |
#          |        |
#          ----------
#            Uz=0
#
# Dirichlet boundary conditions
# Ux(-4000,y,z) = 0
# Ux(+4000,y,z) = 0
# Uy(x,-4000,z) = 0
# Uy(x,+4000,z) = 0
# Uz(x,y,-8000) = 0
#
# Neumann boundary conditions
# Tz(x,y,0) = T0

import numpy

# Physical properties.
p_density = 2500.0
p_vs = 3464.1016
p_vp = 6000.0
p_viscosity = 9.46728e17

p_mu = p_density*p_vs*p_vs
p_lambda = p_density*p_vp*p_vp - 2.0*p_mu
p_youngs = p_mu*(3.0*p_lambda + 2.0*p_mu)/(p_lambda + p_mu)
p_poissons = 0.5*p_lambda/(p_lambda + p_mu)

# Time information.
year = 60.0*60.0*24.0*365.25
dt = 0.025*year
startTime = dt
endTime = 0.5*year
numSteps = 20
timeArray = numpy.linspace(startTime, endTime, num=numSteps, dtype=numpy.float64)

# Uniform stress field (plane strain).
T0 = -1.0e7
szz = T0*numpy.ones(numSteps, dtype=numpy.float64)
timeFac = numpy.exp(-p_youngs*timeArray/(6.0*p_viscosity*(1.0 - p_poissons)))
poisFac = (2.0*p_poissons - 1.0)/(1.0 - p_poissons)
sxx = T0*(1.0 + poisFac*timeFac)
syy = T0*(1.0 + poisFac*timeFac)
sxy = numpy.zeros(numSteps, dtype=numpy.float64)
syz = numpy.zeros(numSteps, dtype=numpy.float64)
sxz = numpy.zeros(numSteps, dtype=numpy.float64)

# Deviatoric stress.
meanStress = (sxx + syy + szz)/3.0
sDevxx = sxx - meanStress
sDevyy = syy - meanStress
sDevzz = szz - meanStress

# Uniform strain field.
exx = numpy.zeros(numSteps, dtype=numpy.float64)
eyy = numpy.zeros(numSteps, dtype=numpy.float64)
ezz = T0*(1.0 - 2.0*p_poissons)*(3.0 + 2.0*poisFac*timeFac)/p_youngs
exy = numpy.zeros(numSteps, dtype=numpy.float64)
eyz = numpy.zeros(numSteps, dtype=numpy.float64)
exz = numpy.zeros(numSteps, dtype=numpy.float64)

# outArray = numpy.column_stack((timeArray, syy, ezz))
# numpy.savetxt('axialtraction_maxwell_analytical.txt', outArray)

# Get viscous strains from deviatoric stress.
eVisxx = 0.5*sDevxx/p_mu
eVisyy = 0.5*sDevyy/p_mu
eViszz = 0.5*sDevzz/p_mu
eVisxy = 0.5*sxy/p_mu
eVisyz = 0.5*syz/p_mu
eVisxz = 0.5*sxz/p_mu

# ----------------------------------------------------------------------
class AnalyticalSoln(object):
    """Analytical solution to axial extension problem.
    """
    SPACE_DIM = 3
    TENSOR_SIZE = 6

    def __init__(self):
        self.fields = {
            "displacement": self.displacement,
            "density": self.density,
            "shear_modulus": self.shear_modulus,
            "bulk_modulus": self.bulk_modulus,
            "maxwell_time": self.maxwell_time,
            "cauchy_strain": self.strain,
            "cauchy_stress": self.stress,
            "viscous_strain": self.viscous_strain,
            "initial_amplitude": {
                "bc_xneg": self.initial_displacement,
                "bc_yneg": self.initial_displacement,
                "bc_zneg": self.initial_displacement,
                "bc_xpos": self.initial_displacement,
                "bc_ypos": self.initial_displacement,
                "bc_zpos": self.initial_traction
            }
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
        disp = numpy.zeros((numSteps, npts, self.SPACE_DIM), dtype=numpy.float64)
        disp[:,:, 2] = numpy.dot(ezz.reshape(numSteps, 1), (locs[:, 2] + 8000.0).reshape(1, npts))
        return disp

    def initial_displacement(self, locs):
        """Compute initial displacement field at locations.
        """
        (npts, dim) = locs.shape
        disp = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        return disp

    def initial_traction(self, locs):
        """Compute initial traction field at locations.
        """
        (npts, dim) = locs.shape
        traction = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        traction[:,:, 2] = T0
        return traction

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

    def maxwell_time(self, locs):
        """Compute Maxwell time field at locations.
        """
        (npts, dim) = locs.shape
        maxwell_time = p_viscosity * numpy.ones((1, npts, 1), dtype=numpy.float64)/p_mu
        return maxwell_time

    def strain(self, locs):
        """Compute strain field at locations.
        """
        (npts, dim) = locs.shape
        strain = numpy.zeros((numSteps, npts, self.TENSOR_SIZE), dtype=numpy.float64)
        strain[:,:, 0] = exx.reshape(numSteps, 1)
        strain[:,:, 1] = eyy.reshape(numSteps, 1)
        strain[:,:, 2] = ezz.reshape(numSteps, 1)
        strain[:,:, 3] = exy.reshape(numSteps, 1)
        strain[:,:, 4] = eyz.reshape(numSteps, 1)
        strain[:,:, 5] = exz.reshape(numSteps, 1)
        return strain

    def stress(self, locs):
        """Compute stress field at locations.
        """
        (npts, dim) = locs.shape
        stress = numpy.zeros((numSteps, npts, self.TENSOR_SIZE), dtype=numpy.float64)
        stress[:,:, 0] = sxx.reshape(numSteps, 1)
        stress[:,:, 1] = syy.reshape(numSteps, 1)
        stress[:,:, 2] = szz.reshape(numSteps, 1)
        stress[:,:, 3] = sxy.reshape(numSteps, 1)
        stress[:,:, 4] = syz.reshape(numSteps, 1)
        stress[:,:, 5] = sxz.reshape(numSteps, 1)
        return stress

    def viscous_strain(self, locs):
        """Compute viscous strain field at locations.
        """
        (npts, dim) = locs.shape
        viscous_strain = numpy.zeros((numSteps, npts, self.TENSOR_SIZE), dtype=numpy.float64)
        viscous_strain[:,:, 0] = eVisxx.reshape(numSteps, 1)
        viscous_strain[:,:, 1] = eVisyy.reshape(numSteps, 1)
        viscous_strain[:,:, 2] = eViszz.reshape(numSteps, 1)
        viscous_strain[:,:, 3] = eVisxy.reshape(numSteps, 1)
        viscous_strain[:,:, 4] = eVisyz.reshape(numSteps, 1)
        viscous_strain[:,:, 5] = eVisxz.reshape(numSteps, 1)
        return viscous_strain


# End of file
