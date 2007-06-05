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

    self.dimension = 3

    self.numDBValues = 4
    self.dbValues = ["density", "vs", "vp" , "viscosity"]
    self.numParameters = 6
    self.numParamValues = [1, 1, 1, 1, 6, 6]
    self.parameterNames = ["density", "mu", "lambda", "maxwellTime", "strainT", "visStrain"]

    self.dt = 2.0e5

    densityA = 2500.0
    vsA = 3000.0
    vpA = vsA*3**0.5
    viscosityA = 1.0e18
    strainA = [1.1e-4, 2.2e-4, 3.3e-4, 4.4e-4, 5.5e-4, 6.6e-4]
    meanStrainA = (strainA[1] + strainA[2] + strainA[3])/3.0
    
    densityB = 2000.0
    vsB = 1200.0
    vpB = vsB*3**0.5
    viscosityB = 1.0e19
    strainB = [1.2e-4, 2.3e-4, 3.4e-4, 4.5e-4, 5.6e-4, 6.7e-4]
    meanStrainB = (strainB[1] + strainB[2] + strainB[3])/3.0

    diag = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0]

    self.dbData = numpy.array([ [densityA, vsA, vpA, viscosityA],
                                [densityB, vsB, vpB, viscosityB] ],
                              dtype=numpy.float64)
    muA = vsA*vsA*densityA
    lambdaA = vpA*vpA*densityA - 2.0*muA
    maxwellTimeA = viscosityA/muA
    visStrainA = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    muB = vsB*vsB*densityB
    lambdaB = vpB*vpB*densityB - 2.0*muB
    maxwellTimeB = viscosityB/muB
    visStrainB = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    # Simplest approach for now is to assume this is the first step after the elastic solution.
    # In that case, both the total strain from the last step (strainT) and the total viscous
    # strain (visStrain) are defined by the assigned elastic strain.
    strainTA = strainA[:]
    strainTB = strainB[:]
    for i in range(6):
      visStrainA[i] = strainA[i] - diag[i] * meanStrainA
      visStrainB[i] = strainB[i] - diag[i] * meanStrainB

    vecParamsA = numpy.hstack((strainTA, visStrainA))
    vecParamsB = numpy.hstack((strainTB, visStrainB))
    vecParams = numpy.vstack((vecParamsA, vecParamsB))
    scalarParams = numpy.array([ [densityA, muA, lambdaA, maxwellTimeA],
                                       [densityB, muB, lambdaB, maxwellTimeB] ],
                                     dtype=numpy.float64)
    self.parameterData = numpy.hstack((scalarParams, vecParams))
    
    self.numLocs = 2
    numElasticConsts = 21
    self.density = numpy.array([densityA, densityB],
                               dtype=numpy.float64)

    self.strain = numpy.array([strainA, strainB],
                               dtype=numpy.float64)
    self.stress = numpy.zeros( (self.numLocs, 6), dtype=numpy.float64)
    self.elasticConsts = numpy.zeros( (self.numLocs, numElasticConsts),
                                      dtype=numpy.float64)

    (self.elasticConsts[0,:], self.stress[0,:]) = \
                              self._calcStress(strainA, muA, lambdaA, maxwellTimeA, strainTA, visStrainA)
    (self.elasticConsts[1,:], self.stress[1,:]) = \
                              self._calcStress(strainB, muB, lambdaB, maxwellTimeB, strainTB, visStrainB)
    return


  def _calcStress(self, strainV, muV, lambdaV, maxwellTimeV, strainTV, visStrainV):
    """
    Compute stress and derivative of elasticity matrix.
    This assumes behavior is always viscoelastic.
    """
    import math
    
    bulkModulus = lambdaV + 2.0 * muV/3.0

    diag = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0]
    traceStrainT = strainTV[0] + strainTV[1] + strainTV[2]
    traceStrainTpdt = strainV[0] + strainV[1] + strainV[2]
    meanStrainT = traceStrainT/3.0
    meanStrainTpdt = traceStrainTpdt/3.0
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
    C1212 = 3.0*visFac
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
    elasFac = 2.0*muV
    devStrainTpdt = 0.0
    devStrainT = 0.0
    devStressTpdt = 0.0
    visStrain = 0.0
    for iComp in range(6):
      devStrainTpdt = strainV[iComp] - diag[iComp]*meanStrainTpdt
      devStrainT = strainTV[iComp] - diag[iComp]*meanStrainT
      visStrain = expFac*visStrainV[iComp] + dq*(devStrainTpdt - devStrainT)
      devStressTpdt = elasFac*visStrain
      stressV[iComp] = diag[iComp]*meanStressTpdt + devStressTpdt
      
    stress = numpy.reshape(stressV, (6,1))
    return (elasticConsts, numpy.ravel(stress))
  

# MAIN /////////////////////////////////////////////////////////////////
if __name__ == "__main__":

  app = MaxwellIsotropic3DTimeDep()
  app.run()


# End of file 
