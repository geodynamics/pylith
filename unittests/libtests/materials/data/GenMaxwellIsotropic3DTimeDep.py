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
dimension = 3
numElasticConsts = 21
tensorSize = 6
numMaxwellModels = 3

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

    numLocs = 2
    self.dt = 2.0e5
    self.dimension = dimension
    self.numLocs = numLocs

    self.dbPropertyValues = ["density", "vs", "vp",
                             "shear_ratio_1", "shear_ratio_2", "shear_ratio_3",
                             "viscosity_1", "viscosity_2", "viscosity_3"]
    self.numPropertyValues = numpy.array([1, 1, 1, 1, 1, 1, 1, 1, 1], dtype=numpy.int32)

    self.dbStateVarValues = ["total-strain-xx",
                             "total-strain-yy",
                             "total-strain-zz",
                             "total-strain-xy",
                             "total-strain-yz",
                             "total-strain-xz",
                             "viscous-strain-1-xx",
                             "viscous-strain-1-yy",
                             "viscous-strain-1-zz",
                             "viscous-strain-1-xy",
                             "viscous-strain-1-yz",
                             "viscous-strain-1-xz",
                             "viscous-strain-2-xx",
                             "viscous-strain-2-yy",
                             "viscous-strain-2-zz",
                             "viscous-strain-2-xy",
                             "viscous-strain-2-yz",
                             "viscous-strain-2-xz",
                             "viscous-strain-3-xx",
                             "viscous-strain-3-yy",
                             "viscous-strain-3-zz",
                             "viscous-strain-3-xy",
                             "viscous-strain-3-yz",
                             "viscous-strain-3-xz",
                             ]
    self.numStateVarValues = numpy.array([tensorSize]*(1+numMaxwellModels),
                                         dtype=numpy.int32)

    densityA = 2500.0
    vsA = 3000.0
    vpA = vsA*3**0.5
    shearRatioA = [0.5, 0.1, 0.2]
    viscosityA = [1.0e18, 1.0e17, 1.0e19]
    strainA = [1.1e-4, 2.2e-4, 3.3e-4, 4.4e-4, 5.5e-4, 6.6e-4]
    initialStressA = [2.1e4, 2.2e4, 2.3e4, 2.4e4, 2.5e4, 2.6e4]
    #initialStrainA = [3.1e-4, 3.2e-4, 3.3e-4, 3.4e-4, 3.5e-4, 3.6e-4]
    muA = vsA*vsA*densityA
    lambdaA = vpA*vpA*densityA - 2.0*muA

    densityB = 2000.0
    vsB = 1200.0
    vpB = vsB*3**0.5
    shearRatioB = [0.2, 0.2, 0.2]
    viscosityB = [1.0e18, 1.0e19, 1.0e20]
    strainB = [1.2e-4, 2.3e-4, 3.4e-4, 4.5e-4, 5.6e-4, 6.7e-4]
    initialStressB = [5.1e4, 5.2e4, 5.3e4, 5.4e4, 5.5e4, 5.6e4]
    #initialStrainB = [6.1e-4, 6.2e-4, 6.3e-4, 6.4e-4, 6.5e-4, 6.6e-4]
    muB = vsB*vsB*densityB
    lambdaB = vpB*vpB*densityB - 2.0*muB

    # TEMPORARY, need to determine how to use initial strain
    initialStrainA = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    initialStrainB = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    
    # Simplest approach for now is to assume this is the first step after the
    # elastic solution. In that case, both the total strain from the last step
    # (strainT) and the total viscous strain (visStrain) are defined by the
    # assigned elastic strain.
    diag = numpy.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0], dtype=numpy.float64)
    strainTA = numpy.array(strainA)
    strainTB = numpy.array(strainB)
    meanStrainA = (strainA[0] + strainA[1] + strainA[2])/3.0
    meanStrainB = (strainB[0] + strainB[1] + strainB[2])/3.0

    maxwellTimeA = [0.0, 0.0, 0.0]
    maxwellTimeB = [0.0, 0.0, 0.0]
    visStrainA = numpy.zeros( (numMaxwellModels, tensorSize), dtype=numpy.float64)
    visStrainB = numpy.zeros( (numMaxwellModels, tensorSize), dtype=numpy.float64)
    for imodel in xrange(numMaxwellModels):
      if shearRatioA[imodel] != 0.0:
        maxwellTimeA[imodel] = viscosityA[imodel]/(muA*shearRatioA[imodel])
        visStrainA[imodel,:] = strainA[:] - diag[:] * meanStrainA
      if shearRatioB[imodel] != 0.0:
        maxwellTimeB[imodel] = viscosityB[imodel]/(muB*shearRatioB[imodel])
        visStrainB[imodel,:] = strainB[:] - diag[:] * meanStrainB

    self.lengthScale = 1.0e+3
    self.pressureScale = muA
    self.timeScale = 1.0
    self.densityScale = 1.0e+3

    propA = [densityA, vsA, vpA] + shearRatioA + viscosityA
    propB = [densityB, vsB, vpB] + shearRatioB + viscosityB
    self.dbProperties = numpy.array([propA, propB], dtype=numpy.float64)
    propA = [densityA, muA, lambdaA] + shearRatioA + maxwellTimeA
    propB = [densityB, muB, lambdaB] + shearRatioB + maxwellTimeB
    self.properties = numpy.array([propA, propB], dtype=numpy.float64)

    # TEMPORARY, need to determine how to use initial state variables
    self.dbStateVars = numpy.zeros( (numLocs, tensorSize+numMaxwellModels*tensorSize),
                                    dtype=numpy.float64)
    self.stateVars = numpy.zeros( (numLocs, tensorSize+numMaxwellModels*tensorSize),
                                  dtype=numpy.float64)
    self.stateVars[0,0:tensorSize] = strainTA
    self.stateVars[0,tensorSize:(1+numMaxwellModels)*tensorSize] = visStrainA.ravel()
    self.stateVars[1,0:tensorSize] = strainTB
    self.stateVars[1,tensorSize:(1+numMaxwellModels)*tensorSize] = visStrainB.ravel()

    mu0 = self.pressureScale
    density0 = self.densityScale
    time0 = self.timeScale
    self.propertiesNondim = \
        numpy.array([ [densityA/density0, muA/mu0, lambdaA/mu0,
                       shearRatioA[0], shearRatioA[1], shearRatioA[2],
                       maxwellTimeA[0]/time0, maxwellTimeA[1]/time0, maxwellTimeA[2]/time0],
                      [densityB/density0, muB/mu0, lambdaB/mu0,
                       shearRatioB[0], shearRatioB[1], shearRatioB[2],
                       maxwellTimeB[0]/time0, maxwellTimeB[1]/time0, maxwellTimeB[2]/time0] ],
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
    self.elasticConsts = numpy.zeros( (numLocs, numElasticConsts), \
                                        dtype=numpy.float64)

    (self.elasticConsts[0,:], self.stress[0,:]) = \
        self._calcStress(strainA, muA, lambdaA, shearRatioA, maxwellTimeA,
                         strainTA, visStrainA,
                         initialStressA, initialStrainA)
    (self.elasticConsts[1,:], self.stress[1,:]) = \
        self._calcStress(strainB, muB, lambdaB, shearRatioB, maxwellTimeB,
                         strainTB, visStrainB,
                         initialStressB, initialStrainB)
    self.dtStableImplicit = 0.1*min(min(maxwellTimeA), min(maxwellTimeB))

    return


  def _calcStress(self, strainV, muV, lambdaV, shearRatioV, maxwellTimeV,
                  strainTV, visStrainV, initialStressV, initialStrainV):
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
    for imodel in range(numMaxwellModels):
      visFrac += shearRatioV[imodel]
      if shearRatioV[imodel] != 0.0:
        if self.dt < timeFrac*maxwellTimeV[imodel]:
          fSign = 1.0
          factorial = 1.0
          fraction = 1.0
          dq[imodel] = 1.0
          for iTerm in range(2, numTerms + 1):
            factorial *= iTerm
            fSign *= -1.0
            fraction *= self.dt/maxwellTimeV[imodel]
            dq[imodel] += fSign*fraction/factorial
        elif maxwellTimeV[imodel] < timeFrac*self.dt:
          dq[imodel] = maxwellTimeV[imodel]/dt
        else:
          dq[imodel] = maxwellTimeV[imodel] * \
          (1.0-math.exp(-self.dt/maxwellTimeV[imodel]))/self.dt
        visFac += shearRatioV[imodel] * dq[imodel]

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
    for iComp in xrange(tensorSize):
      devStrainTpdt = strainV[iComp] - diag[iComp]*meanStrainTpdt
      devStrainT = strainTV[iComp] - diag[iComp]*meanStrainT
      deltaStrain = devStrainTpdt - devStrainT
      devStressTpdt = elasFrac*devStrainTpdt
      for imodel in range(numMaxwellModels):
        if shearRatioV[imodel] != 0.0:
          visStrain = math.exp(-self.dt/maxwellTimeV[imodel]) * \
                      visStrainV[imodel,iComp] + \
                      dq[imodel] * deltaStrain
          devStressTpdt += shearRatioV[imodel] * visStrain
      devStressTpdt = elasFac * devStressTpdt
      stressV[iComp] = diag[iComp]*meanStressTpdt + devStressTpdt + \
                       initialStressV[iComp]
      
    stress = numpy.reshape(stressV, (6,1))
    return (elasticConsts, numpy.ravel(stress))
  

# MAIN /////////////////////////////////////////////////////////////////
if __name__ == "__main__":

  app = GenMaxwellIsotropic3DTimeDep()
  app.run()


# End of file 
