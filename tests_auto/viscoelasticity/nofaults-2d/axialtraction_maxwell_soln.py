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
# @file tests_auto/viscoelasticity/nofaults-2d/axialtraction_maxwell_soln.py
#
# @brief Analytical solution to axial traction problem for a Maxwell viscoelastic material.
#
# 2-D axial traction solution for linear Maxwell viscoelastic material.
#
#             Uy=0
#          ----------
#          |        |
# Ux=0     |        |  Tx=T0
#          |        |
#          |        |
#          ----------
#            Uy=0
#
# Dirichlet boundary conditions
# Ux(-4000,y) = 0
# Uy(x,-4000) = 0
# Uy(x,+4000) = 0
#
# Neumann boundary conditions
# Tx(+4000,y) = T0

import numpy
import pdb

pdb.set_trace()

# Physical properties.
p_density = 2500.0
p_vs = 3464.1016
p_vp = 6000.0
p_viscosity = 9.46728e17

p_mu = p_density * p_vs**2
p_lambda = p_density * p_vp**2 - 2 * p_mu
p_youngs = p_mu*(3.0*p_lambda + 2.0*p_mu)/(p_lambda + p_mu)
p_poissons = 0.5*p_lambda/(p_lambda + p_mu)

# Time information.
year = 60.0*60.0*24.0*365.25
dt = 0.1*year
startTime = dt
endTime = 10.1*year
numSteps = 101
timeArray = numpy.linspace(startTime, endTime, num=numSteps, dtype=numpy.float64)

# Uniform stress field (plane strain).
T0 = -1.0e7
sxx = T0*numpy.ones(numSteps, dtype=numpy.float64)
timeFac = numpy.exp(-p_youngs*timeArray/(6.0*p_viscosity*(1.0 - p_poissons)))
poisFac = (2.0*p_poissons - 1.0)/(1.0 - p_poissons)
syy = T0*(1.0 + poisFac*timeFac)
szz = syy
sxy = numpy.zeros(numSteps, dtype=numpy.float64)

# Uniform strain field.
exx = T0*(1.0 - 2.0*p_poissons)*(3.0 + 2.0*poisFac*timeFac)/p_youngs
eyy = numpy.zeros(numSteps, dtype=numpy.float64)
ezz = numpy.zeros(numSteps, dtype=numpy.float64)
exy = numpy.zeros(numSteps, dtype=numpy.float64)

# Total deviatoric strain.
meanStrain = exx/3.0
eDevxx = exx - meanStrain
eDevyy = eyy - meanStrain
eDevzz = ezz - meanStrain
eDevxy = exy

# Elastic deviatoric strain.
eElasxx = T0*(1.0 - 2.0*p_poissons)*(3.0 + 2.0*poisFac)/p_youngs
meanElasStrain = eElasxx/3.0
eElasDevxx = eElasxx - meanElasStrain
eElasDevyy = -meanElasStrain
eElasDevzz = -meanElasStrain
eElasDevxy = 0.0

# Viscous strains.
eVisxx = eDevxx - eElasDevxx
eVisyy = eDevyy - eElasDevyy
eViszz = eDevzz - eElasDevzz
eVisxy = eDevxy - eElasDevxy

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
            "maxwell_time": self.maxwell_time,
            "cauchy_strain": self.strain,
            "cauchy_stress": self.stress,
            "viscous_strain": self.viscous_strain,
            "initial_amplitude": {
                "bc_xneg": self.initial_displacement,
                "bc_yneg": self.initial_displacement,
                "bc_ypos": self.initial_displacement,
                "bc_xpos": self.initial_traction
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
        disp[:, :, 0] = numpy.dot(exx.reshape(numSteps, 1), (locs[:, 0] + 4000.0).reshape(1, npts))
        return disp

    def initial_displacement(self, locs):
        """
        Compute initial displacement field at locations.
        """
        (npts, dim) = locs.shape
        disp = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        return disp

    def initial_traction(self, locs):
        """
        Compute initial traction field at locations.
        """
        (npts, dim) = locs.shape
        traction = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        traction[:,:,0] = T0
        return traction

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
        maxwell_time = p_viscosity * numpy.ones((1, npts, 1), dtype=numpy.float64)/p_mu
        return maxwell_time

    def strain(self, locs):
        """
        Compute strain field at locations.
        """
        (npts, dim) = locs.shape
        strain = numpy.zeros((numSteps, npts, self.TENSOR_SIZE), dtype=numpy.float64)
        strain[:, :, 0] = exx
        strain[:, :, 1] = eyy
        strain[:, :, 2] = ezz
        strain[:, :, 3] = exy
        return strain

    def stress(self, locs):
        """
        Compute stress field at locations.
        """
        (npts, dim) = locs.shape
        stress = numpy.zeros((numSteps, npts, self.TENSOR_SIZE), dtype=numpy.float64)
        stress[:, :, 0] = sxx
        stress[:, :, 1] = syy
        stress[:, :, 2] = szz
        stress[:, :, 3] = sxy
        return stress

    def viscous_strain(self, locs):
        """
        Compute viscous strain field at locations.
        """
        (npts, dim) = locs.shape
        viscous_strain = numpy.zeros((numSteps, npts, self.TENSOR_SIZE), dtype=numpy.float64)
        viscous_strain[:, :, 0] = eVisxx
        viscous_strain[:, :, 1] = eVisyy
        viscous_strain[:, :, 2] = eViszz
        viscous_strain[:, :, 3] = eVisxy
        return viscous_strain


# End of file
