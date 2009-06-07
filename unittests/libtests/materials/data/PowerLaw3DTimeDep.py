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

## @file unittests/libtests/materials/data/PowerLaw3DTimeDep.py

## @brief Python application for generating C++ data files for testing
## C++ PowerLaw3D object with viscoelastic behavior.

from ElasticMaterialApp import ElasticMaterialApp

import numpy

# ----------------------------------------------------------------------
dimension = 3
numElasticConsts = 21
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

    numLocs = 2

    self.dimension = dimension
    self.numLocs = numLocs

    self.dbPropertyValues = ["density", "vs", "vp",
                             "viscosity_coeff", "power_law_exponent"]
    self.propertyValues = ["density", "mu", "lambda",
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
    self.stateVarValues = ["viscous-strain", "stress"]
    self.numStateVarValues = numpy.array([6, 6], dtype=numpy.int32)

    self.alpha = 0.5
    self.dt = 2.0e5

    densityA = 2500.0
    vsA = 3000.0
    vpA = vsA*3**0.5
    viscosityCoeffA = 1.0e18
    powerLawExponentA = 1.0
    strainA = [1.1e-4, 1.2e-4, 1.3e-4, 1.4e-4, 1.5e-4, 1.6e-4]
    initialStressA = [2.1e4, 2.2e4, 2.3e4, 2.4e4, 2.5e4, 2.6e4]
    initialStrainA = [3.6e-4, 3.5e-4, 3.4e-4, 3.3e-4, 3.2e-4, 3.1e-4]
    muA = vsA*vsA*densityA
    lambdaA = vpA*vpA*densityA - 2.0*muA

    densityB = 2000.0
    vsB = 1200.0
    vpB = vsB*3**0.5
    viscosityCoeffB = 1.2e16
    powerLawExponentB = 3.0
    strainB = [4.1e-4, 4.2e-4, 4.3e-4, 4.4e-4, 4.5e-4, 4.6e-4]
    initialStressB = [5.1e4, 5.2e4, 5.3e4, 5.4e4, 5.5e4, 5.6e4]
    initialStrainB = [6.1e-4, 6.2e-4, 6.3e-4, 6.6e-4, 6.5e-4, 6.4e-4]
    muB = vsB*vsB*densityB
    lambdaB = vpB*vpB*densityB - 2.0*muB

    # TEMPORARY, need to determine how to use initial strain
    initialStrainA = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    initialStrainB = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    
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
    self.dbStateVars = numpy.zeros( (numLocs, tensorSize),
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
    visStrainA = [4.1e-5, 4.2e-5 4.3e-5, 4.4e-5, 4.5e-5, 4.6e-5]
    visStrainB = [1.1e-5, 1.2e-5 1.3e-5, 1.4e-5, 1.5e-5, 1.6e-5]
    stressA = [3.1e4, 3.2e4, 3.3e4, 3.4e4, 3.5e4, 3.6e4]
    stressB = [5.1e4, 5.2e4, 5.3e4, 5.4e4, 5.5e4, 5.6e4]
    stressNondimA = stressA/mu0
    stressNondimB = stressB/mu0

    self.stateVars = numpy.array([ [visStrainA, stressA],
                                   [visStrainB, stressB] ],
                                 dtype=numpy.float64)
    self.stateVarsNondim = numpy.array([ [visStrainA, stressNondimA],
                                         [visStrainB, stressNondimB] ],
                                       dtype=numpy.float64)
    
    self.strain = numpy.array([strainA, strainB],
                               dtype=numpy.float64)
    self.stress = numpy.zeros( (numLocs, tensorSize), dtype=numpy.float64)
    self.elasticConsts = numpy.zeros( (self.numLocs, numElasticConsts),
                                      dtype=numpy.float64)

    (self.elasticConsts[0,:], self.stress[0,:]) = \
                              self._calcStress(strainA, 
                                               muA, lambdaA, viscosityCoeffA,
                                               powerLawExponentA,
                                               visStrainA, stressA,
                                               initialStressA, initialStrainA)
    (self.elasticConsts[1,:], self.stress[1,:]) = \
                              self._calcStress(strainB, 
                                               muB, lambdaB, viscosityCoeffB, 
                                               visStrainB, stressB,
                                               initialStressB, initialStrainB)

    maxwellTimeA = self._getMaxwellTime(muA, viscosityCoeffA, \
                                        powerLawExponentA, self.stress[0,:])
    maxwellTimeB = self._getMaxwellTime(muB, viscosityCoeffB, \
                                        powerLawExponentB, self.stress[1,:])
    self.dtStableImplicit = 0.1*min(maxwellTimeA, maxwellTimeB)
    return


  def _bracket(effStressInitialGuess, effStressParamsTuple):
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

    funcValue1 = self._effStressFunc(x1, effStressParamsTuple)
    funcValue2 = self._effStressFunc(x2, effStressParamsTuple)

    iteration = 0
    bracketed = False

    while iteration < maxIterations:
      if (funcValue1 * funcValue2) < 0.0:
        bracketed = True
        break
      if abs(funcValue1) < abs(funcValue2):
        x1 += bracketFactor * (x1 - x2)
        x1 = max(x1, 0.0)
        funcValue1 = self._effStressFunc(x1, effStressParamsTuple)
      else:
        x2 += bracketFactor * (x1 - x2)
        x2 = max(x2, 0.0)
        funcValue2 = self._effStressFunc(x2, effStressParamsTuple)
      iteration += 1

    if bracketed == False:
      raise RuntimeError("Unable to bracket root.")

    return x1, x2
    
    
  def _getMaxwellTime(self, mu, viscosityCoeff, powerLawExponent, stress):
    """
    Compute Maxwell time from stress, viscosity coefficient, shear modulus, and
    power-law exponent.
    """
    meanStress = (stress[0] + stress[1] + stress[2])/3.0
    devStress = stress
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
      maxwellTime = (viscosityCoeff/effStress)**(powerLawExponent - 1.0) * \
                    (viscosityCoeff/mu)

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
                           mu, lambda, viscosityCoeff,
                           powerLawExponent, visStrainT, stressT,
                           initialStress, initialStrain):
    """
    Function to compute a particular stress component as a function of a
    strain component.
    """
    strainTpdt[strainComp] = strainVal
    stressTpdt = self._computeStress(strainTpdt, mu, lambda, viscosityCoeff,
                                    powerLawExponent, visStrainT, stressT,
                                    initialStress, initialStrain)
    return stressTpdt[stressComp]

    
  def _computeStress(strainTpdt, muV, lambdaV, viscosityCoeffV,
                     powerLawExponentV, visStrainT, stressT,
                     initialStress, initialStrain):
    """
    Function to compute stresses using the effective stress function algorithm.
    """
    import scipy
    
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
    
    # Initial strain values
    meanStrainInitial = (initialStrain[0] + initialStrain[1] +
                         initialStrain[2])/3.0

    # Values for current time step
    meanStrainTpdt = (strainTpdt[0] + strainTpdt[1] + strainTpdt[2])/3.0
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

    # Define parameters for effective stress function
    effStressParamsTuple = (ae, b, c, d, self.alpha, self.dt, effStressT,
                            powerLawExpV, viscosityCoeffV)

    # Bracket the root
    effStressInitialGuess = effStressT
    x1, x2 = self._bracket(effStressInitialGuess, effStressParamsTuple)

    # Find the root using Brent's method (from scipy)
    effStressTpdt = scipy.optimize.brentq(self._effStressFunc, x1, x2,
                                          args=effStressParamsTuple)
    
    # Compute stresses from the effective stress.
    effStressTau = (1.0 - self.alpha) * effStressT + self.alpha * effStressTpdt
    gammaTau = 0.5 * ((effStressTau/viscosityCoeffV)**(powerLawExpV - 1.0)) / \
               viscosityCoeff
    factor1 = 1.0/(ae + self.alpha * self.dt * gammaTau)
    factor2 = timeFac * gammTau
    devStressTpdt = 0.0
    stressTpdt = numpy.zeros( (tensorSize), dtype=numpy.float64)

    for iComp in range(tensorSize):
      devStressTpdt = factor1 * (strainPPTpdt[iComp] - \
                                 factor2 * devStressT[iComp] + \
                                 ae * devStressInitial[iComp])
      stressTpdt[iComp] = devStressTpdt + diag[iComp] * \
                          (meanStressTpdt + meanStressInitial)

    return stressTpdt, effStressParamsTuple

  
  def _calcStress(self, strainV, muV, lambdaV, viscosityCoeffV,
                  powerLawExponentV,visStrainV, stressV,
                  initialStressV, initialStrainV):
    """
    Compute stress and derivative of elasticity matrix.
    This assumes behavior is always viscoelastic.
    """
    import scipy

    # Define some numpy arrays
    strainTpdt = numpy.array(strainV, dtype=numpy.float64)
    visStrainT = numpy.array(visStrainV, dtype=numpy.float64)
    stressT = numpy.array(stressV, dtype=numpy.float64)
    initialStress = numpy.array(initialStressV, dtype=numpy.float64)
    initialStrain = numpy.array(initialStrainV, dtype=numpy.float64)

    stressTpdt = self._computeStress(strainTpdt, muV, lambdaV,
                                     viscosityCoeffV, powerLawExponentV,
                                     visStrainT, stressT,
                                     initialStress, initialStrain)

    # Compute components of tangent constitutive matrix using numerical
    # derivatives.
    calcStressCompParamsList = [0, 0, strainTpdt, muV, lambdaV, viscosityCoeffV,
                                powerLawExponentV, visStrainT, stressT,
                                initialStress, initialStrain]
    derivDx = 1.0
    derivOrder = 5
    elasticConstsList = []

    for stressComponent in range(tensorSize):
      for strainComponent in range(stressComponent, tensorSize):
        calcStressCompParamsList[0] = strainComponent
        calcStressCompParamsList[1] = stressComponent
        calcStressCompParamsTuple = tuple(calcStressCompParamsList)
        dStressDStrain = scipy.misc.derivative(self._calcStressComponent,
                                               strainTpdt[strainComp],
                                               dx=derivDx,
                                               args=calcStressCompParamsTuple,
                                               order=derivOrder)
        elasticConstsList.append(dStressDStrain)

                                               
    
      
    

    C1111 = bulkModulus + 4.0*visFac
    C1122 = bulkModulus - 2.0*visFac
    C1133 = C1122
    C1112 = 0.0
    C1123 = 0.0
    C1113 = 0.0
    C2222 = C1111
    C2233 = C1122
    C2212 = 0.0
    C2223 = 0.0
    C2213 = 0.0
    C3333 = C1111
    C3312 = 0.0
    C3323 = 0.0
    C3313 = 0.0
    C1212 = 6.0*visFac
    C1223 = 0.0
    C1213 = 0.0
    C2323 = C1212
    C2313 = 0.0
    C1313 = C1212
    elasticConsts = numpy.array([C1111, C1122, C1133, C1112, C1123, C1113,
                                 C2222, C2233, C2212, C2223, C2213,
                                 C3333, C3312, C3323, C3313,
                                 C1212, C1223, C1213,
                                 C2323, C2313,
                                 C1313], dtype=numpy.float64)

    stressV = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    expFac = math.exp(-self.dt/maxwellTimeV)
    print "expFac:",expFac
    print "visStrain",visStrainV
    elasFac = 2.0*muV
    devStrainTpdt = 0.0
    devStrainT = 0.0
    devStressTpdt = 0.0
    visStrain = 0.0
    for iComp in range(tensorSize):
      devStrainTpdt = strainV[iComp] - diag[iComp]*meanStrainTpdt
      devStrainT = totalStrainR[iComp] - diag[iComp]*meanStrainT
      visStrain = expFac*visStrainV[iComp] + dq*(devStrainTpdt - devStrainT)
      devStressTpdt = elasFac*visStrain
      stressV[iComp] = diag[iComp]*meanStressTpdt + devStressTpdt + \
                       initialStressV[iComp]
      
    stress = numpy.reshape(stressV, (tensorSize,1))
    return (elasticConsts, numpy.ravel(stress))
  

# MAIN /////////////////////////////////////////////////////////////////
if __name__ == "__main__":

  app = PowerLaw3DTimeDep()
  app.run()


# End of file 
