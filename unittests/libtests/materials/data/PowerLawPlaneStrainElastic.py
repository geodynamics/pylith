#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file unittests/libtests/materials/data/PowerLawPlaneStrainElastic.py

## @brief Python application for generating C++ data files for testing
## C++ PowerLawPlaneStrain object with elastic behavior.

from ElasticMaterialApp import ElasticMaterialApp

import numpy

# ----------------------------------------------------------------------
dimension = 2
numElasticConsts = 9
tensorSize = 3

# PowerLawPlaneStrainElastic class
class PowerLawPlaneStrainElastic(ElasticMaterialApp):
  """
  Python application for generating C++ data files for testing C++
  PowerLawPlaneStrain object with elastic behavior.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="powerlawplanestrainelastic"):
    """
    Constructor.
    """
    ElasticMaterialApp.__init__(self, name)

    # import pdb
    # pdb.set_trace()
    numLocs = 2

    self.dimension = dimension
    self.numLocs = numLocs

    self.dbPropertyValues = ["density", "vs", "vp",
                             "reference-strain-rate", "reference-stress",
                             "power-law-exponent"]
    self.numPropertyValues = numpy.array([1, 1, 1, 1, 1, 1], dtype=numpy.int32)

    self.dbStateVarValues = ["stress-zz-initial",
                             "viscous-strain-xx",
                             "viscous-strain-yy",
                             "viscous-strain-zz",
                             "viscous-strain-xy",
                             "stress4-xx",
                             "stress4-yy",
                             "stress4-zz",
                             "stress4-xy"
                             ]
    self.numStateVarValues = numpy.array([1, 4, 4], dtype=numpy.int32)

    densityA = 2500.0
    vsA = 3000.0
    vpA = vsA*3**0.5
    # Derive new values in based on previous value for power-law coefficient
    # and viscosity coefficient.
    powerLawCoeffA = 1.0/3.0e18
    refStrainRateA = 1.0e-6
    powerLawExponentA = 1.0
    strainA = [1.1e-4, 1.2e-4, 1.4e-4]
    initialStressA = [2.1e4, 2.2e4, 2.4e4]
    initialStrainA = [3.1e-4, 3.2e-4, 3.4e-4]
    muA = vsA*vsA*densityA
    lambdaA = vpA*vpA*densityA - 2.0*muA
    stressZZInitialA = numpy.array([1.0754e4], dtype=numpy.float64)

    viscosityCoeffA = (1.0/((3.0**0.5)**(powerLawExponentA + 1.0) \
                            * powerLawCoeffA))**(1.0/powerLawExponentA)
    refStressA = viscosityCoeffA * \
                 (2.0 * refStrainRateA) ** (1.0/powerLawExponentA)
    # refStressA = (refStrainRateA/powerLawCoeffA)**(1.0/powerLawExponentA)
    
    densityB = 2000.0
    vsB = 1200.0
    vpB = vsB*3**0.5
    powerLawCoeffB = 1.0/9.0e30
    refStrainRateB = 1.0e-6
    powerLawExponentB = 3.0
    strainB = [4.1e-4, 4.2e-4, 4.4e-4]
    initialStressB = [5.1e4, 5.2e4, 5.4e4]
    initialStrainB = [6.1e-4, 6.2e-4, 6.4e-4]
    muB = vsB*vsB*densityB
    lambdaB = vpB*vpB*densityB - 2.0*muB
    viscosityCoeffB = (1.0/((3.0**0.5)**(powerLawExponentB + 1.0) \
                            * powerLawCoeffB))**(1.0/powerLawExponentB)
    refStressB = viscosityCoeffB * \
                 (2.0 * refStrainRateB) ** (1.0/powerLawExponentB)
    # refStressB = (refStrainRateB/powerLawCoeffB)**(1.0/powerLawExponentB)
    stressZZInitialB = numpy.array([2.575e4], dtype=numpy.float64)

    self.lengthScale = 1.0e+3
    self.pressureScale = muA
    self.timeScale = 1.0
    self.densityScale = muA / (self.lengthScale / self.timeScale)**2
    self.strainRateScale = 1.0/self.timeScale

    self.dbProperties = numpy.array([ [densityA, vsA, vpA, \
                                       refStrainRateA, refStressA, \
                                       powerLawExponentA],
                                      [densityB, vsB, vpB, \
                                       refStrainRateB, refStressB, \
                                       powerLawExponentB] ], 
                                    dtype=numpy.float64)
    self.properties = numpy.array([ [densityA, muA, lambdaA, \
                                     refStrainRateA, refStressA, \
                                     powerLawExponentA],
                                    [densityB, muB, lambdaB, \
                                     refStrainRateB, refStressB, \
                                     powerLawExponentB] ],
                                     dtype=numpy.float64)

    # TEMPORARY, need to determine how to use initial state variables
    self.dbStateVars = numpy.zeros( (numLocs, 1 + 4 + 4),
                                    dtype=numpy.float64)
    self.dbStateVars[0,0] = stressZZInitialA
    self.dbStateVars[1,0] = stressZZInitialB

    self.stateVars = numpy.zeros( (numLocs, 1 + 4 + 4), dtype=numpy.float64)
    self.stateVars[0,0] = stressZZInitialA
    self.stateVars[1,0] = stressZZInitialB

    mu0 = self.pressureScale
    density0 = self.densityScale
    time0 = self.timeScale
    strainRate0 = self.strainRateScale
    self.propertiesNondim = \
        numpy.array([ [densityA/density0, muA/mu0, lambdaA/mu0, \
                       refStrainRateA/strainRate0, refStressA/mu0, \
                       powerLawExponentA],
                      [densityB/density0, muB/mu0, lambdaB/mu0, \
                       refStrainRateB/strainRate0, refStressB/mu0, \
                       powerLawExponentB] ],
                    dtype=numpy.float64)

    stressZZInitialANondim = stressZZInitialA/mu0
    stressZZInitialBNondim = stressZZInitialB/mu0
    
    self.stateVarsNondim = numpy.zeros( (numLocs, 1 + 4 + 4),
                                        dtype=numpy.float64)
    self.stateVarsNondim[0, 0] = stressZZInitialANondim
    self.stateVarsNondim[1, 0] = stressZZInitialBNondim

    self.initialStress = numpy.array([initialStressA,
                                      initialStressB],
                                    dtype=numpy.float64)
    self.initialStrain = numpy.array([initialStrainA,
                                      initialStrainB],
                                    dtype=numpy.float64)
    
    self.density = numpy.array([densityA,
                                densityB],
                               dtype=numpy.float64)

    self.strain = numpy.array([strainA,
                               strainB],
                               dtype=numpy.float64)
    
    self.stress = numpy.zeros( (numLocs, tensorSize), dtype=numpy.float64)
    self.elasticConsts = numpy.zeros( (self.numLocs, numElasticConsts), \
                                        dtype=numpy.float64)

    (self.elasticConsts[0,:], self.stress[0,:]) = \
        self._calcStress(strainA, muA, lambdaA, \
                           initialStressA, initialStrainA)
    (self.elasticConsts[1,:], self.stress[1,:]) = \
        self._calcStress(strainB, muB, lambdaB, \
                           initialStressB, initialStrainB)

    viscousStrainUpdated = numpy.zeros((numLocs, 4), dtype=numpy.float64)

    stressUpdated = numpy.zeros((numLocs, 4), dtype=numpy.float64)
    stressUpdated[0,0] = self.stress[0,0]
    stressUpdated[1,0] = self.stress[1,0]
    stressUpdated[0,1] = self.stress[0,1]
    stressUpdated[1,1] = self.stress[1,1]
    stressUpdated[0,3] = self.stress[0,2]
    stressUpdated[1,3] = self.stress[1,2]
    stressUpdated[0,2] = lambdaA * (strainA[0] + strainA[1]) + stressZZInitialA
    stressUpdated[1,2] = lambdaB * (strainB[0] + strainB[1]) + stressZZInitialB
    
    self.stateVarsUpdated = numpy.zeros( (numLocs, 1 + 4 + 4),
                                         dtype=numpy.float64)
    self.stateVarsUpdated[0,0] = stressZZInitialA
    self.stateVarsUpdated[1,0] = stressZZInitialB
    self.stateVarsUpdated[0,5:9] = stressUpdated[0,:]
    self.stateVarsUpdated[1,5:9] = stressUpdated[1,:]

    maxwellTimeA = self._getMaxwellTime(muA, refStrainRateA, refStressA,
                                        powerLawExponentA,
                                        stressUpdated[0,:])
    maxwellTimeB = self._getMaxwellTime(muB, refStrainRateB, refStressB,
                                        powerLawExponentB,
                                        stressUpdated[1,:])


    self.dtStableImplicit = 0.2 * min(maxwellTimeA, maxwellTimeB)
    self.dtStableExplicit = 1000.0 / vpA

    return


  def _scalarProduct(self, tensor1, tensor2):
    """
    Compute the scalar product of two tensors stored in vector form.
    """
    scalarProduct = tensor1[0] * tensor2[0] + \
                    tensor1[1] * tensor2[1] + \
                    tensor1[2] * tensor2[2] + \
                    2.0 * tensor1[3] * tensor2[3]
    return scalarProduct

    
  def _getMaxwellTime(self, mu, refStrainRate, refStress, powerLawExponent,
                      stress4):
    """
    Compute Maxwell time from stress4, viscosity coefficient, shear modulus, and
    power-law exponent.
    """
    meanStress = (stress4[0] + stress4[1] + stress4[2])/3.0
    devStress = numpy.array(stress4, dtype=numpy.float64)
    
    devStress[0] = stress4[0] - meanStress
    devStress[1] = stress4[1] - meanStress
    devStress[2] = stress4[2] - meanStress

    devStressProd = self._scalarProduct(devStress, devStress)
    effStress = (0.5 * devStressProd)**0.5
    maxwellTime = 1.0
    if (effStress != 0.0):
      maxwellTime = (refStress/effStress)**(powerLawExponent - 1.0) * \
                    (refStress/mu)/(refStrainRate * 6.0)

    return maxwellTime

  def _calcStress(self, strainV, muV, lambdaV, initialStressV, initialStrainV):
    """
    Compute stress and derivative of elasticity matrix.
    """
    C1111 = lambdaV + 2.0*muV
    C1122 = lambdaV
    C1112 = 0.0
    C2211 = lambdaV
    C2222 = lambdaV + 2.0*muV
    C2212 = 0.0
    C1211 = 0.0
    C1222 = 0.0
    C1212 = 2.0*muV
    elasticConsts = numpy.array([C1111, C1122, C1112,
                                 C2211, C2222, C2212,
                                 C1211, C1222, C1212], dtype=numpy.float64)

    strain = numpy.reshape(strainV, (tensorSize,1))
    initialStress = numpy.reshape(initialStressV, (tensorSize,1))
    initialStrain = numpy.reshape(initialStrainV, (tensorSize,1))
    elastic = numpy.array([ [C1111, C1122, C1112],
                            [C2211, C2222, C2212],
                            [C1211, C1222, C1212] ], dtype=numpy.float64)
    stress = numpy.dot(elastic, strain-initialStrain) + initialStress
    return (elasticConsts, numpy.ravel(stress))
  

# MAIN /////////////////////////////////////////////////////////////////
if __name__ == "__main__":

  app = PowerLawPlaneStrainElastic()
  app.run()


# End of file 
