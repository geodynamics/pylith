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

## @file unittests/libtests/materials/data/PowerLaw3DTimeDep.py

## @brief Python application for generating C++ data files for testing
## C++ PowerLaw3D object with viscoelastic behavior.

from ElasticMaterialApp import ElasticMaterialApp

import numpy

# ----------------------------------------------------------------------
dimension = 3
numElasticConsts = 36
tensorSize = 6

# PowerLaw3DTimeDep class
class PowerLaw3DTimeDep(ElasticMaterialApp):
  """
  Python application for generating C++ data files for testing C++
  PowerLaw3D object using viscoelastic behavior.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="powerlaw3dtimedep"):
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
    self.stateVarValues = ["viscous-strain", "stress"]
    self.numStateVarValues = numpy.array([6, 6], dtype=numpy.int32)

    self.alpha = 0.5
    self.dt = 2.0e5

    densityA = 2500.0
    vsA = 3000.0
    vpA = vsA*3**0.5
    powerLawCoeffA = 1.0/3.0e18
    refStrainRateA = 1.0e-6
    powerLawExponentA = 1.0
    strainA = [1.1e-4, 1.2e-4, 1.3e-4, 1.4e-4, 1.5e-4, 1.6e-4]
    initialStressA = [2.1e4, 2.2e4, 2.3e4, 2.4e4, 2.5e4, 2.6e4]
    initialStrainA = [3.6e-5, 3.5e-5, 3.4e-5, 3.3e-5, 3.2e-5, 3.1e-5]
    muA = vsA*vsA*densityA
    lambdaA = vpA*vpA*densityA - 2.0*muA
    viscosityCoeffA = (1.0/((3.0**0.5)**(powerLawExponentA + 1.0) \
                            * powerLawCoeffA))**(1.0/powerLawExponentA)
    refStressA = viscosityCoeffA * \
                 (2.0 * refStrainRateA) ** (1.0/powerLawExponentA)

    densityB = 2000.0
    vsB = 1200.0
    vpB = vsB*3**0.5
    powerLawCoeffB = 1.0/9.0e36
    refStrainRateB = 1.0e-6
    powerLawExponentB = 3.0
    strainB = [4.1e-4, 4.2e-4, 4.3e-4, 4.4e-4, 4.5e-4, 4.6e-4]
    initialStressB = [5.1e4, 5.2e4, 5.3e4, 5.4e4, 5.5e4, 5.6e4]
    initialStrainB = [6.1e-5, 6.2e-5, 6.3e-5, 6.6e-5, 6.5e-5, 6.4e-5]
    muB = vsB*vsB*densityB
    lambdaB = vpB*vpB*densityB - 2.0*muB
    viscosityCoeffB = (1.0/((3.0**0.5)**(powerLawExponentB + 1.0) \
                            * powerLawCoeffB))**(1.0/powerLawExponentB)
    refStressB = viscosityCoeffB * \
                 (2.0 * refStrainRateB) ** (1.0/powerLawExponentB)

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
    visStrainA = [4.1e-5, 4.2e-5, 4.3e-5, 4.4e-5, 4.5e-5, 4.6e-5]
    visStrainB = [1.1e-5, 1.2e-5, 1.3e-5, 1.4e-5, 1.5e-5, 1.6e-5]
    stressA = [3.1e4, 3.2e4, 3.3e4, 3.4e4, 3.5e4, 3.6e4]
    stressB = [5.1e4, 5.2e4, 5.3e4, 5.4e4, 5.5e4, 5.6e4]
    stressNondimA = numpy.array(stressA)/mu0
    stressNondimB = numpy.array(stressB)/mu0

    self.stateVars = numpy.array([ [visStrainA, stressA],
                                   [visStrainB, stressB] ],
                                 dtype=numpy.float64)
    self.stateVarsNondim = numpy.array([ [visStrainA, stressNondimA],
                                         [visStrainB, stressNondimB] ],
                                       dtype=numpy.float64)
    
    self.strain = numpy.array([strainA, strainB],
                               dtype=numpy.float64)
    self.stress = numpy.zeros( (numLocs, tensorSize), dtype=numpy.float64)
    self.stateVarsUpdated = numpy.zeros( (numLocs, tensorSize + tensorSize),
                                         dtype=numpy.float64)
    self.elasticConsts = numpy.zeros( (self.numLocs, numElasticConsts),
                                      dtype=numpy.float64)

    (self.elasticConsts[0,:], self.stress[0,:], self.stateVarsUpdated[0,:]) = \
                              self._calcStress(strainA, 
                                               muA, lambdaA, refStrainRateA,
                                               refStressA,
                                               powerLawExponentA,
                                               visStrainA, stressA,
                                               initialStressA, initialStrainA)
    (self.elasticConsts[1,:], self.stress[1,:], self.stateVarsUpdated[1,:]) = \
                              self._calcStress(strainB, 
                                               muB, lambdaB, refStrainRateB, 
                                               refStressB,
                                               powerLawExponentB,
                                               visStrainB, stressB,
                                               initialStressB, initialStrainB)

    # Use state variables to compute Maxwell times (and stable time step size).
    maxwellTimeA = self._getMaxwellTime(muA, refStrainRateA, refStressA, \
                                        powerLawExponentA, stressA)

    maxwellTimeB = self._getMaxwellTime(muB, refStrainRateB, refStressB, \
                                        powerLawExponentB, stressB)

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
                      stress):
    """
    Compute Maxwell time from stress, reference stress and strain rate, shear
    modulus, and power-law exponent.
    """
    meanStress = (stress[0] + stress[1] + stress[2])/3.0
    devStress = numpy.array(stress, dtype=numpy.float64)
    
    devStress[0] = stress[0] - meanStress
    devStress[1] = stress[1] - meanStress
    devStress[2] = stress[2] - meanStress
    devStress[3] = stress[3]
    devStress[4] = stress[4]
    devStress[5] = stress[5]

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
                    2.0 * (tensor1[3] * tensor2[3] + \
                           tensor1[4] * tensor2[4] + \
                           tensor1[5] * tensor2[5])
    return scalarProduct

    
  def _calcStressComponent(self, strainVal, strainComp, stressComp, strainTpdt,
                           muV, lambdaV, refStrainRateV, refStressV,
                           powerLawExponentV, visStrainT, stressT,
                           initialStress, initialStrain):
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
                                                    initialStrain)
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
                     initialStress, initialStrain):
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
    diag = numpy.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0], dtype=numpy.float64)

    # Initial stress values
    meanStressInitial = (initialStress[0] + initialStress[1] +
                         initialStress[2])/3.0
    devStressInitial = initialStress - meanStressInitial * diag
    stressInvar2Initial = 0.5 * self._scalarProduct(devStressInitial,
                                                    devStressInitial)
    
    # Initial strain values
    meanStrainInitial = (initialStrain[0] + initialStrain[1] +
                         initialStrain[2])/3.0

    # Values for current time step
    meanStrainTpdt = (strainTpdt[0] + strainTpdt[1] + strainTpdt[2])/3.0 - \
                     meanStrainInitial
    meanStressTpdt = 3.0 * bulkModulus * meanStrainTpdt
    strainPPTpdt = strainTpdt - meanStrainTpdt * diag - \
                   visStrainT - initialStrain
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
    stressTpdt = numpy.zeros( (tensorSize), dtype=numpy.float64)
    visStrainTpdt = numpy.zeros( (tensorSize), dtype=numpy.float64)

    for iComp in range(tensorSize):
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
                  initialStressV, initialStrainV):
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
                                                    initialStrain)

    stateVarsUpdated = numpy.array( [visStrainTpdt, stressTpdt],
                                    dtype=numpy.float64)

    # Compute components of tangent constitutive matrix using numerical
    # derivatives.
    derivDx = 1.0e-12
    derivOrder = 3
    elasticConstsList = []

    for stressComp in range(tensorSize):
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
                                                     initialStrain),
                                               order=derivOrder)
        elasticConstsList.append(dStressDStrain)

    elasticConsts = numpy.array(elasticConstsList, dtype=numpy.float64)

    return (elasticConsts, numpy.ravel(stressTpdt),
            numpy.ravel(stateVarsUpdated))
  

# MAIN /////////////////////////////////////////////////////////////////
if __name__ == "__main__":

  app = PowerLaw3DTimeDep()
  app.run()


# End of file 
