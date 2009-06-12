#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# <LicenseText>
#
# ----------------------------------------------------------------------
#

## @file unittests/libtests/materials/data/PowerLaw3DElastic.py

## @brief Python application for generating C++ data files for testing
## C++ PowerLaw3D object with elastic behavior.

from ElasticMaterialApp import ElasticMaterialApp

import numpy

# ----------------------------------------------------------------------
dimension = 3
numElasticConsts = 21
tensorSize = 6

# PowerLaw3DElastic class
class PowerLaw3DElastic(ElasticMaterialApp):
  """
  Python application for generating C++ data files for testing C++
  PowerLaw3D object with elastic behavior.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="powerlaw3delastic"):
    """
    Constructor.
    """
    ElasticMaterialApp.__init__(self, name)

    numLocs = 2

    self.dimension = dimension
    self.numLocs = numLocs

    self.dbPropertyValues = ["density", "vs", "vp",
                             "viscosity_coeff", "power_law_exponent"]
    self.numPropertyValues = numpy.array([1, 1, 1, 1, 1], dtype=numpy.int32)

    self.dbStateVarValues = ["viscous-strain-xx",
                             "viscous-strain-yy",
                             "viscous-strain-zz",
                             "viscous-strain-xy",
                             "viscous-strain-yz",
                             "viscous-strain-xz",
                             "stress-xx",
                             "stress-yy",
                             "stress-zz",
                             "stress-xy",
                             "stress-yz",
                             "stress-xz",
                             ]
    self.numStateVarValues = numpy.array([6, 6], dtype=numpy.int32)

    densityA = 2500.0
    vsA = 3000.0
    vpA = vsA*3**0.5
    viscosityCoeffA = 1.0e18
    powerLawExponentA = 1.0
    strainA = [1.1e-4, 1.2e-4, 1.3e-4, 1.4e-4, 1.5e-4, 1.6e-4]
    initialStressA = [2.1e4, 2.2e4, 2.3e4, 2.4e4, 2.5e4, 2.6e4]
    initialStrainA = [3.1e-4, 3.2e-4, 3.3e-4, 3.4e-4, 3.5e-4, 3.6e-4]
    muA = vsA*vsA*densityA
    lambdaA = vpA*vpA*densityA - 2.0*muA
    
    densityB = 2000.0
    vsB = 1200.0
    vpB = vsB*3**0.5
    viscosityCoeffB = 1.0e10
    powerLawExponentB = 3.0
    strainB = [4.1e-4, 4.2e-4, 4.3e-4, 4.4e-4, 4.5e-4, 4.6e-4]
    initialStressB = [5.1e4, 5.2e4, 5.3e4, 5.4e4, 5.5e4, 5.6e4]
    initialStrainB = [6.1e-4, 6.2e-4, 6.3e-4, 6.4e-4, 6.5e-4, 6.6e-4]
    muB = vsB*vsB*densityB
    lambdaB = vpB*vpB*densityB - 2.0*muB

    self.lengthScale = 1.0e+3
    self.pressureScale = muA
    self.timeScale = 1.0
    self.densityScale = 1.0e+3
    self.viscosityCoeffScaleA = \
                  (self.pressureScale**(1.0/powerLawExponentA))/self.timeScale
    self.viscosityCoeffScaleB = \
                  (self.pressureScale**(1.0/powerLawExponentB))/self.timeScale

    self.dbProperties = numpy.array([ [densityA, vsA, vpA, \
                                       viscosityCoeffA, powerLawExponentA],
                                      [densityB, vsB, vpB, \
                                       viscosityCoeffB, powerLawExponentB] ], 
                                    dtype=numpy.float64)
    self.properties = numpy.array([ [densityA, muA, lambdaA, \
                                     viscosityCoeffA, powerLawExponentA],
                                    [densityB, muB, lambdaB, \
                                     viscosityCoeffB, powerLawExponentB] ],
                                     dtype=numpy.float64)

    # TEMPORARY, need to determine how to use initial state variables
    self.dbStateVars = numpy.zeros( (numLocs, tensorSize+tensorSize),
                                    dtype=numpy.float64)
    self.stateVars = numpy.zeros( (numLocs, tensorSize+tensorSize),
                                  dtype=numpy.float64)

    mu0 = self.pressureScale
    density0 = self.densityScale
    time0 = self.timeScale
    viscosityCoeff0 = self.viscosityCoeffScaleA
    viscosityCoeff1 = self.viscosityCoeffScaleB
    self.propertiesNondim = \
        numpy.array([ [densityA/density0, muA/mu0, lambdaA/mu0, \
                       viscosityCoeffA/viscosityCoeff0, powerLawExponentA],
                      [densityB/density0, muB/mu0, lambdaB/mu0, \
                       viscosityCoeffB/viscosityCoeff1, powerLawExponentB] ],
                    dtype=numpy.float64)

    self.stateVarsNondim = self.stateVars # no scaling

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

    maxwellTimeA = self._getMaxwellTime(muA, viscosityCoeffA, \
                                        powerLawExponentA, self.stress[0,:])
    maxwellTimeB = self._getMaxwellTime(muB, viscosityCoeffB, \
                                        powerLawExponentB, self.stress[1,:])

    viscousStrainUpdated = numpy.zeros((numLocs, tensorSize),
                                       dtype=numpy.float64)
    stressUpdated = self.stress
    
    self.stateVarsUpdated = numpy.array( [viscousStrainUpdated[0,:],
                                          stressUpdated[0,:],
                                         viscousStrainUpdated[1,:],
                                          stressUpdated[1,:]],
                                         dtype=numpy.float64)

    self.dtStableImplicit = 0.1*min(maxwellTimeA, maxwellTimeB)
    return


  def _scalarProduct(self, tensor1, tensor2):
    """
    Compute the scalar product of two tensors stored in vector form.
    """
    scalarProduct = tensor1[0] * tensor2[0] + \
                    tensor1[1] * tensor2[1] + \
                    tensor1[2] * tensor2[2] + \
                    2.0 * (tensor1[3] * tensor2[3] + \
                           tensor1[4] * tensor2[4] + \
                           tensor1[5] * tensor2[5])
    return scalarProduct

    
  def _getMaxwellTime(self, mu, viscosityCoeff, powerLawExponent, stress):
    """
    Compute Maxwell time from stress, viscosity coefficient, shear modulus, and
    power-law exponent.
    """
    meanStress = (stress[0] + stress[1] + stress[2])/3.0
    devStress = numpy.array(stress, dtype=numpy.float64)
    
    devStress[0] = stress[0] - meanStress
    devStress[1] = stress[1] - meanStress
    devStress[2] = stress[2] - meanStress

    devStressProd = self._scalarProduct(devStress, devStress)
    effStress = (0.5 * devStressProd)**0.5
    maxwellTime = 1.0
    if (effStress != 0.0):
      maxwellTime = (viscosityCoeff/effStress)**(powerLawExponent - 1.0) * \
                    (viscosityCoeff/mu)

    return maxwellTime

  def _calcStress(self, strainV, muV, lambdaV, initialStressV, initialStrainV):
    """
    Compute stress and derivative of elasticity matrix.
    """
    C1111 = lambdaV + 2.0*muV
    C1122 = lambdaV
    C1133 = lambdaV
    C1112 = 0.0
    C1123 = 0.0
    C1113 = 0.0
    C2222 = lambdaV + 2.0*muV
    C2233 = lambdaV
    C2212 = 0.0
    C2223 = 0.0
    C2213 = 0.0
    C3333 = lambdaV + 2.0*muV
    C3312 = 0.0
    C3323 = 0.0
    C3313 = 0.0
    C1212 = 2.0*muV
    C1223 = 0.0
    C1213 = 0.0
    C2323 = 2.0*muV
    C2313 = 0.0
    C1313 = 2.0*muV
    elasticConsts = numpy.array([C1111, C1122, C1133, C1112, C1123, C1113,
                                 C2222, C2233, C2212, C2223, C2213,
                                 C3333, C3312, C3323, C3313,
                                 C1212, C1223, C1213,
                                 C2323, C2313,
                                 C1313], dtype=numpy.float64)

    strain = numpy.reshape(strainV, (6,1))
    initialStress = numpy.reshape(initialStressV, (tensorSize,1))
    initialStrain = numpy.reshape(initialStrainV, (tensorSize,1))
    elastic = numpy.array([ [C1111, C1122, C1133, C1112, C1123, C1113],
                            [C1122, C2222, C2233, C2212, C2223, C2213],
                            [C1133, C2233, C3333, C3312, C3323, C3313],
                            [C1112, C2212, C3312, C1212, C1223, C1213],
                            [C1123, C2223, C3323, C1223, C2323, C2313],
                            [C1113, C2213, C3313, C1213, C2313, C1313] ],
                          dtype=numpy.float64)
    stress = numpy.dot(elastic, strain-initialStrain) + initialStress
    return (elasticConsts, numpy.ravel(stress))
  

# MAIN /////////////////////////////////////////////////////////////////
if __name__ == "__main__":

  app = PowerLaw3DElastic()
  app.run()


# End of file 
