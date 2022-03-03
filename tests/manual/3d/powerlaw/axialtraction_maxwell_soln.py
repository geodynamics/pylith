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
# @brief Analytical solution for single time step to axial traction problem for a Maxwell viscoelastic material.
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

#-------------------------------------------------------------------------------
def AnalyticalSoln(times, locs, params):
    """
    Compute stresses, strains, and state variables at time times.
    """
    p_density = params['p_density']
    p_vs = params['p_vs']
    p_vp = params['p_vp']
    p_viscosity = params['p_viscosity']
    p_power_law_exponent = params['p_power_law_exponent']
    p_power_law_reference_strain_rate = params['p_power_law_reference_strain_rate']
    p_power_law_reference_stress = params['p_power_law_reference_stress']
    p_mu = params['p_mu']
    p_lambda = params['p_lambda']
    p_youngs = params['p_youngs']
    p_poissons = params['p_poissons']
    T0 = params['T0']

    tensorSize = 6
    numSteps = times.shape[0]
    numPoints = locs.shape[0]
    points = locs.reshape(numPoints,3)
    stress = numpy.zeros((numSteps, tensorSize), dtype=numpy.float64)
    devStress = numpy.zeros((numSteps, tensorSize), dtype=numpy.float64)
    strain = numpy.zeros((numSteps, tensorSize), dtype=numpy.float64)
    devStrain = numpy.zeros((numSteps, tensorSize), dtype=numpy.float64)
    maxwellVisStrain = numpy.zeros((numSteps, tensorSize), dtype=numpy.float64)
    powerLawVisStrain = numpy.zeros((numSteps, tensorSize), dtype=numpy.float64)
    timeFac = numpy.exp(-p_youngs*times/(6.0*p_viscosity*(1.0 - p_poissons)))
    poisFac = (2.0*p_poissons - 1.0)/(1.0 - p_poissons)

    # Stress.
    stress[:,0] = T0*(1.0 + poisFac*timeFac)
    stress[:,1] = T0*(1.0 + poisFac*timeFac)
    stress[:,2] = T0

    # Deviatoric stress.
    meanStress = (stress[:,0] + stress[:,1] + stress[:,2])/3.0
    devStress[:,0] = stress[:,0] - meanStress
    devStress[:,1] = stress[:,1] - meanStress
    devStress[:,2] = stress[:,2] - meanStress

    # Uniform strain field.
    strain[:,2] = T0*(1.0 - 2.0*p_poissons)*(3.0 + 2.0*poisFac*timeFac)/p_youngs

    # Vertical displacement.
    dispz = numpy.dot(strain[:,2].reshape(numSteps, numPoints), (points[:,2] + 8000.0)).reshape(numSteps, numPoints)

    # Deviatoric strain.
    meanStrain = (strain[:,0] + strain[:,1] + strain[:,2])/3.0
    devStrain[:,0] = strain[:,0] - meanStrain
    devStrain[:,1] = strain[:,1] - meanStrain
    devStrain[:,2] = strain[:,2] - meanStrain

    # Maxwell viscous strain.
    maxwellVisStrain[:,0] = 0.5*devStress[:,0]/p_mu
    maxwellVisStrain[:,1] = 0.5*devStress[:,1]/p_mu
    maxwellVisStrain[:,2] = 0.5*devStress[:,2]/p_mu

    # Power-law viscous strain.
    powerLawVisStrain[:,0] = devStrain[:,0] - 0.5*devStress[:,0]/p_mu
    powerLawVisStrain[:,1] = devStrain[:,1] - 0.5*devStress[:,1]/p_mu
    powerLawVisStrain[:,2] = devStrain[:,2] - 0.5*devStress[:,2]/p_mu

    return (dispz, stress, devStress, strain, devStrain, maxwellVisStrain, powerLawVisStrain)


# End of file
