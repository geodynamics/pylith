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

## @file unittests/libtests/materials/data/GenMaxwellIsotropic3DElastic.py

## @brief Python application for generating C++ data files for testing
## C++ GenMaxwellIsotropic3D object with elastic behavior.

from ElasticMaterialApp import ElasticMaterialApp

import numpy

# ----------------------------------------------------------------------
dimension = 3
numElasticConsts = 21
tensorSize = 6
numMaxwellModels = 3

# GenMaxwellIsotropic3DElastic class
class GenMaxwellIsotropic3DElastic(ElasticMaterialApp):
  """
  Python application for generating C++ data files for testing C++
  GenMaxwellIsotropic3D object with elastic behavior.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="genmaxwellisotropic3delastic"):
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
    initialStrainA = [3.1e-4, 3.2e-4, 3.3e-4, 3.4e-4, 3.5e-4, 3.6e-4]
    muA = vsA*vsA*densityA
    lambdaA = vpA*vpA*densityA - 2.0*muA

    densityB = 2000.0
    vsB = 1200.0
    vpB = vsB*3**0.5
    shearRatioB = [0.2, 0.2, 0.2]
    viscosityB = [1.0e18, 1.0e19, 1.0e20]
    strainB = [1.2e-4, 2.3e-4, 3.4e-4, 4.5e-4, 5.6e-4, 6.7e-4]
    initialStressB = [5.1e4, 5.2e4, 5.3e4, 5.4e4, 5.5e4, 5.6e4]
    initialStrainB = [6.1e-4, 6.2e-4, 6.3e-4, 6.4e-4, 6.5e-4, 6.6e-4]
    muB = vsB*vsB*densityB
    lambdaB = vpB*vpB*densityB - 2.0*muB

    maxwellTimeA = [0.0, 0.0, 0.0]
    maxwellTimeB = [0.0, 0.0, 0.0]
    for i in xrange(numMaxwellModels):
      if shearRatioA[i] != 0.0:
        maxwellTimeA[i] = viscosityA[i]/(muA*shearRatioA[i])
      if shearRatioB[i] != 0.0:
        maxwellTimeB[i] = viscosityB[i]/(muB*shearRatioB[i])

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
    self.elasticConsts = numpy.zeros( (self.numLocs, numElasticConsts), \
                                        dtype=numpy.float64)

    (self.elasticConsts[0,:], self.stress[0,:]) = \
        self._calcStress(strainA, muA, lambdaA, \
                           initialStressA, initialStrainA)
    (self.elasticConsts[1,:], self.stress[1,:]) = \
        self._calcStress(strainB, muB, lambdaB, \
                           initialStressB, initialStrainB)
    self.dtStableImplicit = 0.1*min(min(maxwellTimeA), min(maxwellTimeB))

    return


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

  app = GenMaxwellIsotropic3DElastic()
  app.run()


# End of file 
