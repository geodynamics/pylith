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

## @file unittests/libtests/materials/data/GenMaxwellIsotropic3DTimeDep.py

## @brief Python application for generating C++ data files for testing
## C++ GenMaxwellIsotropic3D object with viscoelastic behavior.

from ElasticMaterialApp import ElasticMaterialApp

import numpy

# ----------------------------------------------------------------------
# GenMaxwellIsotropic3DTimeDep class
class GenMaxwellIsotropic3DTimeDep(ElasticMaterialApp):
  """
  Python application for generating C++ data files for testing C++
  GenMaxwellIsotropic3D object using viscoelastic behavior.
  """

  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="genmaxwellisotropic3dtimedep"):
    """
    Constructor.
    """
    ElasticMaterialApp.__init__(self, name)

    # import pdb
    # pdb.set_trace()

    self.dimension = 3
    self.tensorSize = 6

    self.numMaxwellModels = 3
    self.numDBValues = 3+2*self.numMaxwellModels
    self.dbValues = ["density", "vs", "vp" ,
                     "shear_ratio_1", "shear_ratio_2", "shear_ratio_3",
                     "viscosity_1", "viscosity_2", "viscosity_3"]
    self.numParameters = 7
    self.numParamValues = [1, 1, 1,
                           self.numMaxwellModels, self.numMaxwellModels,
                           self.tensorSize,
                           self.tensorSize*self.numMaxwellModels]
    self.parameterNames = ["density", "mu", "lambda",
                           "shearRatio", "maxwellTime",
                           "strainT", "visStrain"]

    self.dt = 2.0e5

    densityA = 2500.0
    vsA = 3000.0
    vpA = vsA*3**0.5
    shearRatioA = [0.5, 0.1, 0.2]
    viscosityA = [1.0e18, 1.0e17, 1.0e19]
    visStrainA = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    strainA = [1.1e-4, 2.2e-4, 3.3e-4, 4.4e-4, 5.5e-4, 6.6e-4]
    meanStrainA = (strainA[0] + strainA[1] + strainA[2])/3.0
    elasDataA = [densityA, vsA, vpA]
    matDbA = elasDataA + shearRatioA + viscosityA
    muA = vsA*vsA*densityA
    lambdaA = vpA*vpA*densityA - 2.0*muA
    
    densityB = 2000.0
    vsB = 1200.0
    vpB = vsB*3**0.5
    shearRatioB = [0.2, 0.2, 0.2]
    viscosityB = [1.0e18, 1.0e19, 1.0e20]
    visStrainB = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    strainB = [1.2e-4, 2.3e-4, 3.4e-4, 4.5e-4, 5.6e-4, 6.7e-4]
    meanStrainB = (strainB[0] + strainB[1] + strainB[2])/3.0
    elasDataB = [densityB, vsB, vpB]
    matDbB = elasDataB + shearRatioB + viscosityB
    muB = vsB*vsB*densityB
    lambdaB = vpB*vpB*densityB - 2.0*muB

    matDb = matDbA + matDbB
    self.dbData = numpy.array(matDb, dtype=numpy.float64)

    diag = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0]

    # Simplest approach for now is to assume this is the first step after the
    # elastic solution. In that case, both the total strain from the last step
    # (strainT) and the total viscous strain (visStrain) are defined by the
    # assigned elastic strain.
    strainTA = strainA[:]
    strainTB = strainB[:]

    maxwellTimeA = [0.0, 0.0, 0.0]
    maxwellTimeB = [0.0, 0.0, 0.0]
    for model in range(self.numMaxwellModels):
      if shearRatioA[model] != 0.0:
        maxwellTimeA[model] = viscosityA[model]/(muA*shearRatioA[model])
        for i in range(self.tensorSize):
          visStrainA[i+self.tensorSize*model] = strainA[i] - diag[i] * meanStrainA
      if shearRatioB[model] != 0.0:
        maxwellTimeB[model] = viscosityB[model]/(muB*shearRatioB[model])
        for i in range(self.tensorSize):
          visStrainB[i+self.tensorSize*model] = strainB[i] - diag[i] * meanStrainB
        
    vecParamsA = numpy.hstack((shearRatioA, maxwellTimeA, strainTA, visStrainA))
    vecParamsB = numpy.hstack((shearRatioB, maxwellTimeB, strainTB, visStrainB))
    vecParams = numpy.vstack((vecParamsA, vecParamsB))
    scalarParams = numpy.array([ [densityA, muA, lambdaA],
                                 [densityB, muB, lambdaB] ],
                               dtype=numpy.float64)
    self.parameterData = numpy.hstack((scalarParams, vecParams))
    
    self.numLocs = 2
    numElasticConsts = 21
    self.density = numpy.array([densityA, densityB], dtype=numpy.float64)

    self.strain = numpy.array([strainA, strainB], dtype=numpy.float64)
    self.stress = numpy.zeros( (self.numLocs, 6), dtype=numpy.float64)
    self.elasticConsts = numpy.zeros( (self.numLocs, numElasticConsts),
                                      dtype=numpy.float64)

    (self.elasticConsts[0,:], self.stress[0,:]) = \
                              self._calcStress(strainA, muA, lambdaA,
                                               shearRatioA, maxwellTimeA,
                                               strainTA, visStrainA)
    (self.elasticConsts[1,:], self.stress[1,:]) = \
                              self._calcStress(strainB, muB, lambdaB,
                                               shearRatioB, maxwellTimeB,
                                               strainTB, visStrainB)
    self.dtStableImplicit = 0.1*min(min(maxwellTimeA),
                                    min(maxwellTimeB))
    return


  def _calcStress(self, strainV, muV, lambdaV,
                  shearRatioV, maxwellTimeV, strainTV, visStrainV):
    """
    Compute stress and derivative of elasticity matrix.
    This assumes behavior is always viscoelastic.
    """
    import math
    # import pdb
    
    bulkModulus = lambdaV + 2.0 * muV/3.0

    diag = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0]
    traceStrainT = strainTV[0] + strainTV[1] + strainTV[2]
    traceStrainTpdt = strainV[0] + strainV[1] + strainV[2]
    meanStrainT = traceStrainT/3.0
    meanStrainTpdt = traceStrainTpdt/3.0
    meanStressTpdt = bulkModulus * traceStrainTpdt

    timeFrac = 1.0e-10
    numTerms = 5
    visFrac = 0.0
    visFac = 0.0
    dq = [0.0, 0.0, 0.0]
    for model in range(self.numMaxwellModels):
      visFrac += shearRatioV[model]
      if shearRatioV[model] != 0.0:
        if self.dt < timeFrac*maxwellTimeV[model]:
          fSign = 1.0
          factorial = 1.0
          fraction = 1.0
          dq[model] = 1.0
          for iTerm in range(2, numTerms + 1):
            factorial *= iTerm
            fSign *= -1.0
            fraction *= self.dt/maxwellTimeV[model]
            dq[model] += fSign*fraction/factorial
        elif maxwellTimeV[model] < timeFrac*self.dt:
          dq[model] = maxwellTimeV[model]/dt
        else:
          dq[model] = maxwellTimeV[model] * \
          (1.0-math.exp(-self.dt/maxwellTimeV[model]))/self.dt
        visFac += shearRatioV[model] * dq[model]

    elasFrac = 1.0 - visFrac
    shearFac = muV*(elasFrac + visFac)/3.0

    C1111 = bulkModulus + 4.0*shearFac
    C1122 = bulkModulus - 2.0*shearFac
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
    C1212 = 6.0*shearFac
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

    elasFac = 2.0*muV
    devStrainTpdt = 0.0
    devStrainT = 0.0
    deltaStrain = 0.0
    devStressTpdt = 0.0
    visStrain = 0.0
    for iComp in range(self.tensorSize):
      devStrainTpdt = strainV[iComp] - diag[iComp]*meanStrainTpdt
      devStrainT = strainTV[iComp] - diag[iComp]*meanStrainT
      deltaStrain = devStrainTpdt - devStrainT
      devStressTpdt = elasFrac*devStrainTpdt
      for model in range(self.numMaxwellModels):
        if shearRatioV[model] != 0.0:
          visStrain = math.exp(-self.dt/maxwellTimeV[model]) * \
                      visStrainV[iComp + self.tensorSize * model] + \
                      dq[model] * deltaStrain
          devStressTpdt += shearRatioV[model] * visStrain
      devStressTpdt = elasFac * devStressTpdt
      stressV[iComp] = diag[iComp]*meanStressTpdt + devStressTpdt
      
    stress = numpy.reshape(stressV, (6,1))
    return (elasticConsts, numpy.ravel(stress))
  

# MAIN /////////////////////////////////////////////////////////////////
if __name__ == "__main__":

  app = GenMaxwellIsotropic3DTimeDep()
  app.run()


# End of file 
