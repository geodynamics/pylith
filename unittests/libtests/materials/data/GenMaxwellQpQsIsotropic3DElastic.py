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

## @file unittests/libtests/materials/data/GenMaxwellQpQsIsotropic3DElastic.py

## @brief Python application for generating C++ data files for testing
## C++ GenMaxwellQpQsIsotropic3D object with elastic behavior.

from ElasticMaterialApp import ElasticMaterialApp

import numpy

# ----------------------------------------------------------------------
dimension = 3
numElasticConsts = 36
tensorSize = 6
numMaxwellModels = 3

# GenMaxwellQpQsIsotropic3DElastic class
class GenMaxwellQpQsIsotropic3DElastic(ElasticMaterialApp):
  """
  Python application for generating C++ data files for testing C++
  GenMaxwellQpQsIsotropic3D object with elastic behavior.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="genmaxwellbulkisotropic3delastic"):
    """
    Constructor.
    """
    ElasticMaterialApp.__init__(self, name)

    numLocs = 2
    self.dt = 2.0e5
    self.dimension = dimension
    self.numLocs = numLocs

    self.dbPropertyValues = ["density", "vs", "vp",
                             "shear-ratio-1", 
                             "shear-ratio-2", 
                             "shear-ratio-3",
                             "shear-viscosity-1", 
                             "shear-viscosity-2", 
                             "shear-viscosity-3",
                             "bulk-ratio-1", 
                             "bulk-ratio-2", 
                             "bulk-ratio-3",
                             "bulk-viscosity-1", 
                             "bulk-viscosity-2", 
                             "bulk-viscosity-3",
                             ]
    self.numPropertyValues = numpy.array([1, 1, 1,
                                          1, 1, 1,
                                          1, 1, 1,
                                          1, 1, 1,
                                          1, 1, 1], dtype=numpy.int32)

    self.dbStateVarValues = ["total-strain-xx",
                             "total-strain-yy",
                             "total-strain-zz",
                             "total-strain-xy",
                             "total-strain-yz",
                             "total-strain-xz",
                             "viscous-deviatoric-strain-1-xx",
                             "viscous-deviatoric-strain-1-yy",
                             "viscous-deviatoric-strain-1-zz",
                             "viscous-deviatoric-strain-1-xy",
                             "viscous-deviatoric-strain-1-yz",
                             "viscous-deviatoric-strain-1-xz",
                             "viscous-deviatoric-strain-2-xx",
                             "viscous-deviatoric-strain-2-yy",
                             "viscous-deviatoric-strain-2-zz",
                             "viscous-deviatoric-strain-2-xy",
                             "viscous-deviatoric-strain-2-yz",
                             "viscous-deviatoric-strain-2-xz",
                             "viscous-deviatoric-strain-3-xx",
                             "viscous-deviatoric-strain-3-yy",
                             "viscous-deviatoric-strain-3-zz",
                             "viscous-deviatoric-strain-3-xy",
                             "viscous-deviatoric-strain-3-yz",
                             "viscous-deviatoric-strain-3-xz",
                             "viscous-mean-strain-1",
                             "viscous-mean-strain-2",
                             "viscous-mean-strain-3",
                             ]

    self.numStateVarValues = numpy.array([1, 1, 1, 1, 1, 1,
                                          1, 1, 1, 1, 1, 1,
                                          1, 1, 1, 1, 1, 1,
                                          1, 1, 1, 1, 1, 1,
                                          1, 1, 1], dtype=numpy.int32)

    densityA = 2500.0
    vsA = 3000.0
    vpA = vsA*3**0.5
    shearRatioA = [0.5, 0.1, 0.2]
    shearViscosityA = [1.0e+18, 1.0e+17, 1.0e+19]
    bulkRatioA = [0.4, 0.3, 0.1] 
    bulkViscosityA = [2.0e+18, 2.0e+17, 2.0e+19]
    strainA = [1.1e-4, 2.2e-4, 3.3e-4, 4.4e-4, 5.5e-4, 6.6e-4]
    initialStressA = [2.1e4, 2.2e4, 2.3e4, 2.4e4, 2.5e4, 2.6e4]
    initialStrainA = [3.1e-4, 3.2e-4, 3.3e-4, 3.4e-4, 3.5e-4, 3.6e-4]
    muA = vsA*vsA*densityA
    lambdaA = vpA*vpA*densityA - 2.0*muA
    kA = lambdaA + 2/3.0*muA

    densityB = 2000.0
    vsB = 1200.0
    vpB = vsB*3**0.5
    shearRatioB = [0.2, 0.2, 0.2]
    shearViscosityB = [1.0e18, 1.0e19, 1.0e20]
    bulkRatioB = [0.2, 0.2, 0.2] 
    bulkViscosityB = [2.0e18, 2.0e19, 2.0e20]
    strainB = [1.2e-4, 2.3e-4, 3.4e-4, 4.5e-4, 5.6e-4, 6.7e-4]
    initialStressB = [5.1e4, 5.2e4, 5.3e4, 5.4e4, 5.5e4, 5.6e4]
    initialStrainB = [6.1e-4, 6.2e-4, 6.3e-4, 6.4e-4, 6.5e-4, 6.6e-4]
    muB = vsB*vsB*densityB
    lambdaB = vpB*vpB*densityB - 2.0*muB
    kB = lambdaB + 2/3.0*muB

    maxwellTimeA = [0.0, 0.0, 0.0]
    maxwellTimeB = [0.0, 0.0, 0.0]
    for i in xrange(numMaxwellModels):
      if shearRatioA[i] != 0.0:
        maxwellTimeA[i] = shearViscosityA[i]/muA
      if shearRatioB[i] != 0.0:
        maxwellTimeB[i] = shearViscosityB[i]/muB

    maxwellTimeBulkA = [0.0, 0.0, 0.0]
    maxwellTimeBulkB = [0.0, 0.0, 0.0]
    for i in xrange(numMaxwellModels):
      if bulkRatioA[i] != 0.0:
        maxwellTimeBulkA[i] = bulkViscosityA[i]/kA
      if bulkRatioB[i] != 0.0:
        maxwellTimeBulkB[i] = bulkViscosityB[i]/kB

    self.lengthScale = 1.0e+3
    self.pressureScale = muA
    self.timeScale = 1.0
    self.densityScale = muA / (self.lengthScale / self.timeScale)**2

    propA = [densityA, vsA, vpA] + shearRatioA + shearViscosityA + bulkRatioA + bulkViscosityA
    propB = [densityB, vsB, vpB] + shearRatioB + shearViscosityB + bulkRatioB + bulkViscosityB
    self.dbProperties = numpy.array([propA, propB], dtype=numpy.float64)
    propA = [densityA, muA, kA] + shearRatioA + maxwellTimeA + bulkRatioA + maxwellTimeBulkA
    propB = [densityB, muB, kB] + shearRatioB + maxwellTimeB + bulkRatioB + maxwellTimeBulkB
    self.properties = numpy.array([propA, propB], dtype=numpy.float64)

    # TEMPORARY, need to determine how to use initial state variables
    self.dbStateVars = numpy.zeros( (numLocs, tensorSize+numMaxwellModels*tensorSize + numMaxwellModels),
                                    dtype=numpy.float64)
    self.stateVars = numpy.zeros( (numLocs, tensorSize+numMaxwellModels*tensorSize + numMaxwellModels),
                                  dtype=numpy.float64)

    mu0 = self.pressureScale
    density0 = self.densityScale
    time0 = self.timeScale
    self.propertiesNondim = \
        numpy.array([ [densityA/density0, muA/mu0, kA/mu0,
                       shearRatioA[0], shearRatioA[1], shearRatioA[2],
                       maxwellTimeA[0]/time0, maxwellTimeA[1]/time0, maxwellTimeA[2]/time0,
                       bulkRatioA[0], bulkRatioA[1], bulkRatioA[2],
                       maxwellTimeBulkA[0]/time0, maxwellTimeBulkA[1]/time0, maxwellTimeBulkA[2]/time0],
                      [densityB/density0, muB/mu0, kB/mu0,
                       shearRatioB[0], shearRatioB[1], shearRatioB[2],
                       maxwellTimeB[0]/time0, maxwellTimeB[1]/time0, maxwellTimeB[2]/time0, 
                       bulkRatioB[0], bulkRatioB[1], bulkRatioB[2],
                       maxwellTimeBulkB[0]/time0, maxwellTimeBulkB[1]/time0, maxwellTimeBulkB[2]/time0] ],
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
    self.dtStableImplicit = 0.2*min(min(maxwellTimeA), min(maxwellTimeB))
    self.dtStableExplicit = 1000.0 / vpA

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
    C2211 = lambdaV
    C2222 = lambdaV + 2.0*muV
    C2233 = lambdaV
    C2212 = 0.0
    C2223 = 0.0
    C2213 = 0.0
    C3311 = lambdaV
    C3322 = lambdaV
    C3333 = lambdaV + 2.0*muV
    C3312 = 0.0
    C3323 = 0.0
    C3313 = 0.0
    C1211 = 0.0
    C1222 = 0.0
    C1233 = 0.0
    C1212 = 2.0*muV
    C1223 = 0.0
    C1213 = 0.0
    C2311 = 0.0
    C2322 = 0.0
    C2333 = 0.0
    C2312 = 0.0
    C2323 = 2.0*muV
    C2313 = 0.0
    C1311 = 0.0
    C1322 = 0.0
    C1333 = 0.0
    C1312 = 0.0
    C1323 = 0.0
    C1313 = 2.0*muV
    elasticConsts = numpy.array([C1111, C1122, C1133, C1112, C1123, C1113,
                                 C2211, C2222, C2233, C2212, C2223, C2213,
                                 C3311, C3322, C3333, C3312, C3323, C3313,
                                 C1211, C1222, C1233, C1212, C1223, C1213,
                                 C2311, C2322, C2333, C2312, C2323, C2313,
                                 C1311, C1322, C1333, C1312, C1323, C1313],
				 dtype=numpy.float64)

    strain = numpy.reshape(strainV, (6,1))
    initialStress = numpy.reshape(initialStressV, (tensorSize,1))
    initialStrain = numpy.reshape(initialStrainV, (tensorSize,1))
    elastic = numpy.array([ [C1111, C1122, C1133, C1112, C1123, C1113],
                            [C2211, C2222, C2233, C2212, C2223, C2213],
                            [C3311, C3322, C3333, C3312, C3323, C3313],
                            [C1211, C1222, C1233, C1212, C1223, C1213],
                            [C2311, C2322, C2333, C2312, C2323, C2313],
                            [C1311, C1322, C1333, C1312, C1323, C1313] ],
                          dtype=numpy.float64)
    stress = numpy.dot(elastic, strain-initialStrain) + initialStress
    return (elasticConsts, numpy.ravel(stress))
  

# MAIN /////////////////////////////////////////////////////////////////
if __name__ == "__main__":

  app = GenMaxwellQpQsIsotropic3DElastic()
  app.run()


# End of file 
