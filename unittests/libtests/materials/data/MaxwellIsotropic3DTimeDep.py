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

## @file unittests/libtests/materials/data/MaxwellIsotropic3DTimeDep.py

## @brief Python application for generating C++ data files for testing
## C++ MaxwellIsotropic3D object with viscoelastic behavior.

from ElasticMaterialApp import ElasticMaterialApp

import numpy

# ----------------------------------------------------------------------
dimension = 3
numElasticConsts = 21
tensorSize = 6

# MaxwellIsotropic3DTimeDep class
class MaxwellIsotropic3DTimeDep(ElasticMaterialApp):
  """
  Python application for generating C++ data files for testing C++
  MaxwellIsotropic3D object using viscoelastic behavior.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="maxwellisotropic3dtimedep"):
    """
    Constructor.
    """
    ElasticMaterialApp.__init__(self, name)

    numLocs = 2

    self.dimension = dimension
    self.numLocs = numLocs

    self.dbPropertyValues = ["density", "vs", "vp", "viscosity"]
    self.propertyValues = ["density", "mu", "lambda", "maxwellTime"]
    self.numPropertyValues = numpy.array([1, 1, 1, 1], dtype=numpy.int32)

    self.dbStateVarValues = ["total-strain-xx",
                             "total-strain-yy",
                             "total-strain-zz",
                             "total-strain-xy",
                             "total-strain-yz",
                             "total-strain-xz",
                             "viscous-strain-xx",
                             "viscous-strain-yy",
                             "viscous-strain-zz",
                             "viscous-strain-xy",
                             "viscous-strain-yz",
                             "viscous-strain-xz",
                             ]
    self.stateVarValues = ["total-strain", "viscous-strain"]
    self.numStateVarValues = numpy.array([6, 6], dtype=numpy.int32)

    self.dt = 2.0e5

    densityA = 2500.0
    vsA = 3000.0
    vpA = vsA*3**0.5
    viscosityA = 1.0e18
    strainA = [1.1e-4, 1.2e-4, 1.3e-4, 1.4e-4, 1.5e-4, 1.6e-4]
    initialStressA = [2.1e4, 2.2e4, 2.3e4, 2.4e4, 2.5e4, 2.6e4]
    #initialStrainA = [3.6e-4, 3.5e-4, 3.4e-4, 3.3e-4, 3.2e-4, 3.1e-4]
    muA = vsA*vsA*densityA
    lambdaA = vpA*vpA*densityA - 2.0*muA
    maxwellTimeA = viscosityA / muA
    meanStrainA = (strainA[1] + strainA[2] + strainA[3])/3.0

    densityB = 2000.0
    vsB = 1200.0
    vpB = vsB*3**0.5
    viscosityB = 1.0e19
    strainB = [1.2e-4, 2.3e-4, 3.4e-4, 4.5e-4, 5.6e-4, 6.7e-4]
    strainB = [4.1e-4, 4.2e-4, 4.3e-4, 4.4e-4, 4.5e-4, 4.6e-4]
    initialStressB = [5.1e4, 5.2e4, 5.3e4, 5.4e4, 5.5e4, 5.6e4]
    #initialStrainB = [6.1e-4, 6.2e-4, 6.3e-4, 6.6e-4, 6.5e-4, 6.4e-4]
    muB = vsB*vsB*densityB
    lambdaB = vpB*vpB*densityB - 2.0*muB
    maxwellTimeB = viscosityB / muB
    meanStrainB = (strainB[1] + strainB[2] + strainB[3])/3.0

    # TEMPORARY, need to determine how to use initial strain
    initialStrainA = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    initialStrainB = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    
    diag = numpy.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0],
                       dtype=numpy.float64)

    self.lengthScale = 1.0e+3
    self.pressureScale = muA
    self.timeScale = 1.0
    self.densityScale = 1.0e+3

    self.dbProperties = numpy.array([ [densityA, vsA, vpA, viscosityA],
                                      [densityB, vsB, vpB, viscosityB] ], 
                                    dtype=numpy.float64)
    self.properties = numpy.array([ [densityA, muA, lambdaA, maxwellTimeA],
                                    [densityB, muB, lambdaB, maxwellTimeB] ],
                                     dtype=numpy.float64)

    # TEMPORARY, need to determine how to use initial state variables
    self.dbStateVars = numpy.zeros( (numLocs, tensorSize),
                                    dtype=numpy.float64)

    mu0 = self.pressureScale
    density0 = self.densityScale
    time0 = self.timeScale
    self.propertiesNondim = \
        numpy.array([ [densityA/density0, muA/mu0, lambdaA/mu0, maxwellTimeA/time0],
                      [densityB/density0, muB/mu0, lambdaB/mu0, maxwellTimeB/time0] ],
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

    # Simplest approach for now is to assume this is the first step
    # after the elastic solution.  In that case, both the total strain
    # from the last step (total_strain) and the total viscous strain
    # (viscous_strain) are defined by the assigned elastic strain.
    totalStrainA = strainA
    totalStrainB = strainB
    viscousStrainA = numpy.array(strainA) - diag*meanStrainA
    viscousStrainB = numpy.array(strainB) - diag*meanStrainB
    self.stateVars = numpy.array([ [viscousStrainA, totalStrainA],
                                   [viscousStrainB, totalStrainB] ],
                                 dtype=numpy.float64)
    self.stateVarsNondim = self.stateVars # no scaling
    
    self.strain = numpy.array([strainA, strainB],
                               dtype=numpy.float64)
    self.stress = numpy.zeros( (numLocs, tensorSize), dtype=numpy.float64)
    self.elasticConsts = numpy.zeros( (self.numLocs, numElasticConsts),
                                      dtype=numpy.float64)

    (self.elasticConsts[0,:], self.stress[0,:]) = \
                              self._calcStress(strainA, muA, lambdaA,
                                               maxwellTimeA, 
                                               totalStrainA, viscousStrainA,
                                               initialStressA, initialStrainA)
    (self.elasticConsts[1,:], self.stress[1,:]) = \
                              self._calcStress(strainB, muB, lambdaB,
                                               maxwellTimeB, 
                                               totalStrainB, viscousStrainB,
                                               initialStressB, initialStrainB)

    self.dtStableImplicit = 0.1*min(maxwellTimeA, maxwellTimeB)
    return


  def _calcStress(self, strainV, muV, lambdaV, maxwellTimeV, totalStrainV,
                  viscousStrainV, initialStressV, initialStrainV):
    """
    Compute stress and derivative of elasticity matrix.
    This assumes behavior is always viscoelastic.
    """
    import math
    
    bulkModulus = lambdaV + 2.0 * muV/3.0
    diag = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0]

    totalStrainR = numpy.array(totalStrainV) - numpy.array(initialStrainV)
    print totalStrainV
    print initialStrainV
    print totalStrainR

    traceStrainT = totalStrainR[0] + totalStrainR[1] + totalStrainR[2]
    traceStrainTpdt = strainV[0] + strainV[1] + strainV[2]
    meanStrainT = traceStrainT / 3.0
    meanStrainTpdt = traceStrainTpdt / 3.0
    meanStressTpdt = bulkModulus * traceStrainTpdt
    timeFrac = 1.0e-5
    numTerms = 5
    dq = 0.0
    if maxwellTimeV < timeFrac*self.dt:
      fSign = 1.0
      factorial = 1.0
      fraction = 1.0
      dq = 1.0
      for iTerm in range(2, numTerms + 1):
        factorial *= iTerm
        fSign *= -1.0
        fraction *= self.dt/maxwellTimeV
        dq += fSign*fraction/factorial
    else:
      dq = maxwellTimeV*(1.0-math.exp(-self.dt/maxwellTimeV))/self.dt

    visFac = muV*dq/3.0

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
    print "viscousStrain",viscousStrainV
    elasFac = 2.0*muV
    devStrainTpdt = 0.0
    devStrainT = 0.0
    devStressTpdt = 0.0
    viscousStrain = 0.0
    for iComp in range(tensorSize):
      devStrainTpdt = strainV[iComp] - diag[iComp]*meanStrainTpdt
      devStrainT = totalStrainR[iComp] - diag[iComp]*meanStrainT
      viscousStrain = expFac*viscousStrainV[iComp] + dq*(devStrainTpdt - devStrainT)
      devStressTpdt = elasFac*viscousStrain
      stressV[iComp] = diag[iComp]*meanStressTpdt + devStressTpdt + \
                       initialStressV[iComp]
      
    stress = numpy.reshape(stressV, (tensorSize,1))
    return (elasticConsts, numpy.ravel(stress))
  

# MAIN /////////////////////////////////////////////////////////////////
if __name__ == "__main__":

  app = MaxwellIsotropic3DTimeDep()
  app.run()


# End of file 
