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
"""2-D axial traction solution for linear Maxwell viscoelastic material.

             Ty=T0
          ----------
          |        |
 Ux=0     |        |  Ux=0
          |        |
          |        |
          ----------
            Uy=0

 Dirichlet boundary conditions
 Ux(-4000,y,z) = 0
 Ux(+4000,y,z) = 0
 Uy(x,y,-8000) = 0

 Neumann boundary conditions
 Ty(x,y,0) = T0
"""

import numpy

#-------------------------------------------------------------------------------
p_density = 2500.0
p_vs = 3464.1016
p_vp = 6000.0
p_viscosity = 9.467279917257993e17
p_power_law_exponent = 1.0
p_power_law_reference_strain_rate = 1.0e-6
p_power_law_reference_stress = 2.0*p_viscosity*p_power_law_reference_strain_rate
p_mu = p_density*p_vs*p_vs
p_lambda = p_density*p_vp*p_vp - 2.0*p_mu
p_youngs = p_mu*(3.0*p_lambda + 2.0*p_mu)/(p_lambda + p_mu)
p_poissons = 0.5*p_lambda/(p_lambda + p_mu)

T0 = -1.0e+9


#-------------------------------------------------------------------------------
def AnalyticalSoln(times, locs):
    """
    Compute stresses, strains, and state variables at time times.
    """
    tensorSize = 4
    numSteps = times.shape[0]
    numPoints = locs.shape[0]
    points = locs.reshape(numPoints,2)
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
    stress[:,1] = T0
    stress[:,2] = T0*(1.0 + poisFac*timeFac)

    # Deviatoric stress.
    meanStress = (stress[:,0] + stress[:,1] + stress[:,2])/3.0
    devStress[:,0] = stress[:,0] - meanStress
    devStress[:,1] = stress[:,1] - meanStress
    devStress[:,2] = stress[:,1] - meanStress

    # Uniform strain field.
    strain[:,1] = T0*(1.0 - 2.0*p_poissons)*(3.0 + 2.0*poisFac*timeFac)/p_youngs

    # Vertical displacement.
    dispy = numpy.dot(strain[:,1].reshape(numSteps, numPoints), (points[:,1] + 8000.0)).reshape(numSteps, numPoints)

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
    powerLawVisStrain[:,2] = devStrain[:,2] - 0.5*devStress[:,1]/p_mu

    return (dispy, stress, devStress, strain, devStrain, maxwellVisStrain, powerLawVisStrain)


# End of file
