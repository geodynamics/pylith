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

## @file unittests/libtests/materials/data/PowerLawPlaneStrainTimeDep.py

## @brief Python application for generating C++ data files for testing
## C++ PowerLawPlaneStrain object with viscoelastic behavior.

from ElasticMaterialApp import ElasticMaterialApp

import numpy

# ----------------------------------------------------------------------
dimension = 2
numElasticConsts = 9
tensorSize = 3

# PowerLawPlaneStrainTimeDep class
class PowerLawPlaneStrainTimeDep(ElasticMaterialApp):
  """
  Python application for generating C++ data files for testing C++
  PowerLawPlaneStrain object using viscoelastic behavior.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="powerlawplanestraintimedep"):
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
    self.propertyValues = ["density", "mu", "lambda",
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

    self.alpha = 0.5
    self.dt = 2.0e5

    densityA = 2500.0
    vsA = 3000.0
    vpA = vsA*3**0.5
    powerLawCoeffA = 1.0/3.0e18
    refStrainRateA = 1.0e-6
    powerLawExponentA = 1.0
    strainA = [1.1e-4, 1.2e-4, 1.4e-4]
    initialStressA = [2.1e4, 2.2e4, 2.4e4]
    initialStrainA = [3.6e-5, 3.5e-5, 3.3e-5]
    muA = vsA*vsA*densityA
    lambdaA = vpA*vpA*densityA - 2.0*muA
    viscosityCoeffA = (1.0/((3.0**0.5)**(powerLawExponentA + 1.0) \
                            * powerLawCoeffA))**(1.0/powerLawExponentA)
    refStressA = viscosityCoeffA * \
                 (2.0 * refStrainRateA) ** (1.0/powerLawExponentA)
    stressZZInitialA = 2.3e4

    densityB = 2000.0
    vsB = 1200.0
    vpB = vsB*3**0.5
    powerLawCoeffB = 1.0/9.0e36
    refStrainRateB = 1.0e-6
    powerLawExponentB = 3.0
    strainB = [4.1e-4, 4.2e-4, 4.4e-4]
    initialStressB = [5.1e4, 5.2e4, 5.4e4]
    initialStrainB = [6.1e-5, 6.2e-5, 6.6e-5]
    muB = vsB*vsB*densityB
    lambdaB = vpB*vpB*densityB - 2.0*muB
    viscosityCoeffB = (1.0/((3.0**0.5)**(powerLawExponentB + 1.0) \
                            * powerLawCoeffB))**(1.0/powerLawExponentB)
    refStressB = viscosityCoeffB * \
                 (2.0 * refStrainRateB) ** (1.0/powerLawExponentB)
    stressZZInitialB = 5.3e4

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
    self.dbStateVars = numpy.zeros( (numLocs, tensorSize),
                                    dtype=numpy.float64)
    self.dbStateVars[0, 0] = stressZZInitialA
    self.dbStateVars[1, 0] = stressZZInitialB

    mu0 = self.pressureScale
    density0 = self.densityScale
    time0 = self.timeScale
    strainRate0 = self.strainRateScale
    self.propertiesNondim = \
                          numpy.array([ [densityA/density0, muA/mu0, \
                                         lambdaA/mu0, \
                                         refStrainRateA/strainRate0, \
                                         refStressA/mu0, \
                                         powerLawExponentA], \
                                        [densityB/density0, muB/mu0, \
                                         lambdaB/mu0, \
                                         refStrainRateB/strainRate0, \
                                         refStressB/mu0, \
                                         powerLawExponentB] ], \
                                      dtype=numpy.float64)

    self.initialStress = numpy.array([initialStressA,
                                      initialStressB],
                                    dtype=numpy.float64)
    self.initialStrain = numpy.array([initialStrainA,
                                      initialStrainB],
                                    dtype=numpy.float64)
    
    self.density = numpy.array([densityA,
                                densityB],
                               dtype=numpy.float64)

    # Define state variables
    visStrainA = numpy.array([4.1e-5, 4.2e-5, 4.3e-5, 4.4e-5],
                             dtype=numpy.float64)
    visStrainB = numpy.array([1.1e-5, 1.2e-5, 1.3e-5, 1.4e-5],
                             dtype=numpy.float64)
    stress4A = numpy.array([3.1e4, 3.2e4, 3.3e4, 3.4e4], dtype=numpy.float64)
    stress4B = numpy.array([5.1e4, 5.2e4, 5.3e4, 5.4e4], dtype=numpy.float64)
    stressNondimA = stress4A/mu0
    stressNondimB = stress4B/mu0
    stressZZInitialANondim = stressZZInitialA/mu0
    stressZZInitialBNondim = stressZZInitialB/mu0

    stateVarsA = numpy.concatenate(([stressZZInitialA], visStrainA, stress4A))
    stateVarsB = numpy.concatenate(([stressZZInitialB], visStrainB, stress4B))
    self.stateVars = numpy.array((stateVarsA, stateVarsB), dtype=numpy.float64)

    stateVarsNondimA = numpy.concatenate(([stressZZInitialANondim], visStrainA,
                                          stressNondimA))
    stateVarsNondimB = numpy.concatenate(([stressZZInitialBNondim], visStrainB,
                                          stressNondimB))
    self.stateVarsNondim = numpy.array((stateVarsNondimA, stateVarsNondimB),
                                       dtype=numpy.float64)
    
    self.strain = numpy.array([strainA, strainB],
                               dtype=numpy.float64)
    self.stress = numpy.zeros( (numLocs, tensorSize), dtype=numpy.float64)
    self.stateVarsUpdated = numpy.zeros( (numLocs, 1 + 4 + 4),
                                         dtype=numpy.float64)
    self.elasticConsts = numpy.zeros( (self.numLocs, numElasticConsts),
                                      dtype=numpy.float64)

    (self.elasticConsts[0,:], self.stress[0,:], self.stateVarsUpdated[0,:]) = \
                              self._calcStress(strainA, 
                                               muA, lambdaA, refStrainRateA,
                                               refStressA,
                                               powerLawExponentA,
                                               visStrainA, stress4A,
                                               initialStressA, initialStrainA,
                                               stateVarsA)
    (self.elasticConsts[1,:], self.stress[1,:], self.stateVarsUpdated[1,:]) = \
                              self._calcStress(strainB, 
                                               muB, lambdaB, refStrainRateB, 
                                               refStressB,
                                               powerLawExponentB,
                                               visStrainB, stress4B,
                                               initialStressB, initialStrainB,
                                               stateVarsB)

    # Use state variables to compute Maxwell times (and stable time step size).
    maxwellTimeA = self._getMaxwellTime(muA, refStrainRateA, refStressA, \
                                        powerLawExponentA,
                                        self.stateVarsUpdated[0,5:])

    maxwellTimeB = self._getMaxwellTime(muB, refStrainRateB, refStressB, \
                                        powerLawExponentB,
                                        self.stateVarsUpdated[1,5:])

    self.dtStableImplicit = 0.2 * min(maxwellTimeA, maxwellTimeB)
    self.dtStableExplicit = 1000.0 / vpA

    return


  def _bracket(self, effStressInitialGuess, ae, b, c, d, alpha, dt, effStressT,
               powerLawExponentV, refStrainRateV, refStressV):
    """
    Function to bracket the effective stress.
    """
    maxIterations = 50
    bracketFactor = 1.6

    x1 = 0.0
    x2 = 0.0

    if effStressInitialGuess > 0.0:
      x1 = effStressInitialGuess - 0.5 * effStressInitialGuess
      x2 = effStressInitialGuess + 0.5 * effStressInitialGuess
    else:
      x1 = 500.0
      x2 = 1500.0

    funcValue1 = self._effStressFunc(x1, ae, b, c, d, alpha, dt,
                                     effStressT, powerLawExponentV,
                                     refStrainRateV, refStressV)
    funcValue2 = self._effStressFunc(x2, ae, b, c, d, alpha, dt,
                                     effStressT, powerLawExponentV,
                                     refStrainRateV, refStressV)

    iteration = 0
    bracketed = False

    while iteration < maxIterations:
      if (funcValue1 * funcValue2) < 0.0:
        bracketed = True
        break
      if abs(funcValue1) < abs(funcValue2):
        x1 += bracketFactor * (x1 - x2)
        x1 = max(x1, 0.0)

        funcValue1 = self._effStressFunc(x1, ae, b, c, d, alpha, dt,
                                         effStressT, powerLawExponentV,
                                         refStrainRateV, refStressV)
      else:
        x2 += bracketFactor * (x1 - x2)
        x2 = max(x2, 0.0)

        funcValue2 = self._effStressFunc(x2, ae, b, c, d, alpha, dt,
                                         effStressT, powerLawExponentV,
                                         refStrainRateV, refStressV)
      iteration += 1

    if bracketed == False:
      raise RuntimeError("Unable to bracket root.")

    return x1, x2
    
    
  def _getMaxwellTime(self, mu, refStrainRate, refStress, powerLawExponent,
                      stress4):
    """
    Compute Maxwell time from stress, reference stress and strain rate, shear
    modulus, and power-law exponent.
    """
    meanStress = (stress4[0] + stress4[1] + stress4[2])/3.0
    devStress = numpy.array(stress4, dtype=numpy.float64)
    
    devStress[0] = stress4[0] - meanStress
    devStress[1] = stress4[1] - meanStress
    devStress[2] = stress4[2] - meanStress
    devStress[3] = stress4[3]

    devStressProd = self._scalarProduct(devStress, devStress)
    effStress = (0.5 * devStressProd)**0.5
    maxwellTime = 1.0
    if (effStress != 0.0):
      maxwellTime = (refStress/effStress)**(powerLawExponent - 1.0) * \
                    (refStress/mu)/(refStrainRate * 6.0)

    return maxwellTime


  def _scalarProduct(self, tensor1, tensor2):
    """
    Compute the scalar product of two tensors stored in vector form.
    """
    scalarProduct = tensor1[0] * tensor2[0] + \
                    tensor1[1] * tensor2[1] + \
                    tensor1[2] * tensor2[2] + \
                    2.0 * tensor1[3] * tensor2[3]
    return scalarProduct

    
  def _calcStressComponent(self, strainVal, strainComp, stressComp, strainTpdt,
                           muV, lambdaV, refStrainRateV, refStressV,
                           powerLawExponentV, visStrainT, stressT,
                           initialStress, initialStrain, stateVars):
    """
    Function to compute a particular stress component as a function of a
    strain component.
    """
    strainTest = numpy.array(strainTpdt, dtype=numpy.float64)
    strainTest[strainComp] = strainVal
    stressTpdt, visStrainTpdt = self._computeStress(strainTest, muV, lambdaV,
                                                    refStrainRateV,
                                                    refStressV,
                                                    powerLawExponentV,
                                                    visStrainT,
                                                    stressT,
                                                    initialStress,
                                                    initialStrain,
                                                    stateVars)
    return stressTpdt[stressComp]


  def _effStressFunc(self, effStressTpdt, ae, b, c, d, alpha, dt, effStressT,
                     powerLawExponentV, refStrainRateV, refStressV):
    """
    Function to compute effective stress function for a given effective stress.
    """

    factor1 = 1.0 - alpha
    effStressTau = factor1 * effStressT + alpha * effStressTpdt
    gammaTau = refStrainRateV * (effStressTau/refStressV)** \
               (powerLawExponentV - 1.0) / refStressV
    a = ae + alpha * dt * gammaTau
    effStressFunc = a * a * effStressTpdt * effStressTpdt - b + \
                    c * gammaTau - d * d * gammaTau * gammaTau

    return effStressFunc

    
  def _computeStress(self, strainTpdt, muV, lambdaV, refStrainRateV, refStressV,
                     powerLawExponentV, visStrainT, stressT,
                     initialStress, initialStrain, stateVars):
    """
    Function to compute stresses and viscous strains using the
    effective stress function algorithm.
    """
    import scipy.optimize
    
    # Constants
    mu2 = 2.0 * muV
    lamPlusMu = lambdaV + muV
    bulkModulus = lambdaV + 2.0 * muV/3.0
    ae = 1.0/mu2
    timeFac = self.dt * (1.0 - self.alpha)
    diag = numpy.array([1.0, 1.0, 1.0, 0.0], dtype=numpy.float64)

    # Initial stress values
    initialStress4 = numpy.array([initialStress[0], initialStress[1],
                                  stateVars[0], initialStress[2]],
                                 dtype=numpy.float64)
    meanStressInitial = (initialStress4[0] + initialStress4[1] +
                         initialStress4[3])/3.0
    devStressInitial = initialStress4 - meanStressInitial * diag
    stressInvar2Initial = 0.5 * self._scalarProduct(devStressInitial,
                                                    devStressInitial)
    
    # Initial strain values
    initialStrain4 = numpy.array([initialStrain[0], initialStrain[1],
                                 0.0, initialStrain[2]], dtype=numpy.float64)
    meanStrainInitial = (initialStrain[0] + initialStrain[1])/3.0

    # Values for current time step
    strainTpdt4 = numpy.array([strainTpdt[0], strainTpdt[1],
                                 0.0, strainTpdt[2]], dtype=numpy.float64)
    meanStrainTpdt = (strainTpdt[0] + strainTpdt[1])/3.0 - meanStrainInitial
    meanStressTpdt = 3.0 * bulkModulus * meanStrainTpdt
    strainPPTpdt = strainTpdt4 - meanStrainTpdt * diag - \
                   visStrainT - initialStrain4
    strainPPInvar2Tpdt = 0.5 * self._scalarProduct(strainPPTpdt, strainPPTpdt)

    # Values for previous time step
    meanStressT = (stressT[0] + stressT[1] + stressT[2])/3.0
    devStressT = stressT - diag * meanStressT
    stressInvar2T = 0.5 * self._scalarProduct(devStressT, devStressT)
    effStressT = stressInvar2T**0.5

    # Finish defining parameters needed for root-finding algorithm.
    b = strainPPInvar2Tpdt + \
        ae * self._scalarProduct(strainPPTpdt, devStressInitial) + \
        ae * ae * stressInvar2Initial
    c = (self._scalarProduct(strainPPTpdt, devStressT) + \
         ae * self._scalarProduct(devStressT, devStressInitial))  * timeFac
    d = timeFac * effStressT

    # Bracket the root
    effStressInitialGuess = effStressT

    x1, x2 = self._bracket(effStressInitialGuess, ae, b, c, d, self.alpha,
                           self.dt, effStressT, powerLawExponentV,
                           refStrainRateV, refStressV)

    # Find the root using Brent's method (from scipy)
    rootTolerance = 1.0e-14
    effStressTpdt = scipy.optimize.brentq(self._effStressFunc, x1, x2,
                                          args=(ae, b, c, d, self.alpha,
                                                self.dt, effStressT,
                                                powerLawExponentV,
                                                refStrainRateV, refStressV),
                                          xtol=rootTolerance)
    
    # Compute stresses from the effective stress.
    effStressTau = (1.0 - self.alpha) * effStressT + self.alpha * effStressTpdt
    gammaTau = refStrainRateV * ((effStressTau/refStressV)** \
                      (powerLawExponentV - 1.0)) / refStressV
    factor1 = 1.0/(ae + self.alpha * self.dt * gammaTau)
    factor2 = timeFac * gammaTau
    devStressTpdt = 0.0
    stressTpdt = numpy.zeros( (4), dtype=numpy.float64)
    visStrainTpdt = numpy.zeros( (4), dtype=numpy.float64)

    for iComp in range(4):
      devStressTpdt = factor1 * (strainPPTpdt[iComp] - \
                                 factor2 * devStressT[iComp] + \
                                 ae * devStressInitial[iComp])
      stressTpdt[iComp] = devStressTpdt + diag[iComp] * \
                          (meanStressTpdt + meanStressInitial)
      devStressTau = (1.0 - self.alpha) * devStressT[iComp] + \
                     self.alpha * devStressTpdt
      deltaVisStrain = self.dt * gammaTau * devStressTau
      visStrainTpdt[iComp] = visStrainT[iComp] + deltaVisStrain

    return stressTpdt, visStrainTpdt

  
  def _calcStress(self, strainV, muV, lambdaV, refStrainRateV, refStressV,
                  powerLawExponentV,visStrainV, stressV,
                  initialStressV, initialStrainV, stateVars):
    """
    Compute stress, updated state variables and derivative of elasticity matrix.
    This assumes behavior is always viscoelastic.
    """
    import scipy.misc

    # Define some numpy arrays
    strainTpdt = numpy.array(strainV, dtype=numpy.float64)
    visStrainT = numpy.array(visStrainV, dtype=numpy.float64)
    stressT = numpy.array(stressV, dtype=numpy.float64)
    initialStress = numpy.array(initialStressV, dtype=numpy.float64)
    initialStrain = numpy.array(initialStrainV, dtype=numpy.float64)

    stressTpdt, visStrainTpdt = self._computeStress(strainTpdt, muV, lambdaV,
                                                    refStrainRateV, refStressV,
                                                    powerLawExponentV,
                                                    visStrainT, stressT,
                                                    initialStress,
                                                    initialStrain,
                                                    stateVars)

    stateVarsUpdated = numpy.concatenate(([stateVars[0]], visStrainTpdt,
                                          stressTpdt))

    # Compute components of tangent constitutive matrix using numerical
    # derivatives.
    derivDx = 1.0e-12
    derivOrder = 3
    elasticConstsList = []
    components = [0, 1, 3]
    stress3 = stressTpdt[components]

    for stressComp in components:
      for strainComp in range(tensorSize):
        dStressDStrain = scipy.misc.derivative(self._calcStressComponent,
                                               strainTpdt[strainComp],
                                               dx=derivDx,
                                               args=(strainComp,
                                                     stressComp,
                                                     strainTpdt, muV, lambdaV,
                                                     refStrainRateV, refStressV,
                                                     powerLawExponentV,
                                                     visStrainT,
                                                     stressT, initialStress,
                                                     initialStrain,
                                                     stateVars),
                                               order=derivOrder)
        elasticConstsList.append(dStressDStrain)

    elasticConsts = numpy.array(elasticConstsList, dtype=numpy.float64)

    return (elasticConsts, numpy.ravel(stress3),
            numpy.ravel(stateVarsUpdated))
  

# MAIN /////////////////////////////////////////////////////////////////
if __name__ == "__main__":

  app = PowerLawPlaneStrainTimeDep()
  app.run()


# End of file 
