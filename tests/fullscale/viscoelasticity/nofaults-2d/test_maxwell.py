#!/usr/bin/env nemesis

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
# @file tests/fullscale/viscoelasticity/nofaults-2d/axialtraction_maxwell_soln.py
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

import matplotlib.pyplot as plt
import numpy
import h5py
import pylab

# FE input file.
pylithInput = 'output/axialtraction_maxwell_quad_onecell-maxwell.h5'

# Physical properties.
p_density = 2500.0
p_vs = 3464.1016
p_vp = 6000.0
p_viscosity = 9.46728e17
p_mu = p_density * p_vs*p_vs
p_lambda = p_density * p_vp*p_vp - 2.0 * p_mu
p_youngs = p_mu*(3.0*p_lambda + 2.0*p_mu)/(p_lambda + p_mu)
p_poissons = 0.5*p_lambda/(p_lambda + p_mu)
year = 60.0*60.0*24.0*365.25

# Finite element results.
h5 = h5py.File(pylithInput, 'r')
time = h5['time'][:].flatten()
timeYears = time/year
numSteps = time.shape[0]
stress = h5['vertex_fields/cauchy_stress'][:]
syyFE = stress[:, 0, 1]
strain = h5['vertex_fields/cauchy_strain'][:]
exxFE = strain[:, 0, 0]
h5.close()

# Uniform stress field (plane strain).
T0 = -1.0e7
sxx = T0*numpy.ones(numSteps, dtype=numpy.float64)
timeFac = numpy.exp(-p_youngs*time/(6.0*p_viscosity*(1.0 - p_poissons)))
poisFac = (2.0*p_poissons - 1.0)/(1.0 - p_poissons)
syy = T0*(1.0 + poisFac*timeFac)
szz = syy
sxy = numpy.zeros(numSteps, dtype=numpy.float64)

# Deviatoric stress.
meanStress = (sxx + syy + szz)/3.0
sDevxx = sxx - meanStress
sDevyy = syy - meanStress
sDevzz = szz - meanStress

# Uniform strain field.
exx = T0*(1.0 - 2.0*p_poissons)*(3.0 + 2.0*poisFac*timeFac)/p_youngs
eyy = numpy.zeros(numSteps, dtype=numpy.float64)
ezz = numpy.zeros(numSteps, dtype=numpy.float64)
exy = numpy.zeros(numSteps, dtype=numpy.float64)

# Get viscous strains from deviatoric stress.
eVisxx = 0.5*sDevxx/p_mu
eVisyy = 0.5*sDevyy/p_mu
eViszz = 0.5*sDevzz/p_mu
eVisxy = 0.5*sxy/p_mu

# Total deviatoric strain.
meanStrain = exx/3.0
eDevxx = exx - meanStrain
eDevyy = eyy - meanStrain
eDevzz = ezz - meanStrain
eDevxy = exy

# Plot results.
# fig, ax1 = plt.subplots()
# color = 'tab:black'
# ax1.set_xlabel('time (years)')
# ax1.set_ylabel('stress_yy (Pa)', color=color)
# ax1.plot(timeYears, syy, 'k-')
# ax1.plot(timeFEYears, syyFE, 'r--')
# ax2 = ax1.twinx()

pylab.plot(timeYears, exx, 'k-', timeYears, exxFE, 'r--')
pylab.show()

eDiff = exx - exxFE
eDiffRel = eDiff/exx
pylab.plot(timeYears, eDiff, 'k-')
pylab.show()

pylab.plot(timeYears, eDiffRel, 'k-')
pylab.show()


# End of file
