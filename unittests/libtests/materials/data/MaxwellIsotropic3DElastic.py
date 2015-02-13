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

## @file unittests/libtests/materials/data/MaxwellIsotropic3DElastic.py

## @brief Python application for generating C++ data files for testing
## C++ MaxwellIsotropic3D object with elastic behavior.

from ElasticMaterialApp import ElasticMaterialApp

import numpy

# ----------------------------------------------------------------------
dimension = 3
numElasticConsts = 36
tensorSize = 6

# MaxwellIsotropic3DElastic class
class MaxwellIsotropic3DElastic(ElasticMaterialApp):
  """
  Python application for generating C++ data files for testing C++
  MaxwellIsotropic3D object with elastic behavior.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="maxwellisotropic3delastic"):
    """
    Constructor.
    """
    ElasticMaterialApp.__init__(self, name)

    numLocs = 2

    self.dimension = dimension
    self.numLocs = numLocs

    self.dbPropertyValues = ["density", "vs", "vp", "viscosity"]
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
    self.numStateVarValues = numpy.array([6, 6], dtype=numpy.int32)

    densityA = 2500.0
    vsA = 3000.0
    vpA = vsA*3**0.5
    viscosityA = 1.0e18
    strainA = [1.1e-4, 1.2e-4, 1.3e-4, 1.4e-4, 1.5e-4, 1.6e-4]
    initialStressA = [2.1e4, 2.2e4, 2.3e4, 2.4e4, 2.5e4, 2.6e4]
    initialStrainA = [3.1e-5, 3.2e-5, 3.3e-5, 3.4e-5, 3.5e-5, 3.6e-5]
    muA = vsA*vsA*densityA
    lambdaA = vpA*vpA*densityA - 2.0*muA
    maxwellTimeA = viscosityA / muA
    
    densityB = 2000.0
    vsB = 1200.0
    vpB = vsB*3**0.5
    viscosityB = 1.0e18
    strainB = [4.1e-4, 4.2e-4, 4.3e-4, 4.4e-4, 4.5e-4, 4.6e-4]
    initialStressB = [5.1e4, 5.2e4, 5.3e4, 5.4e4, 5.5e4, 5.6e4]
    initialStrainB = [6.1e-5, 6.2e-5, 6.3e-5, 6.4e-5, 6.5e-5, 6.6e-5]
    muB = vsB*vsB*densityB
    lambdaB = vpB*vpB*densityB - 2.0*muB
    maxwellTimeB = viscosityB / muB

    self.lengthScale = 1.0e+3
    self.pressureScale = muA
    self.timeScale = 1.0
    self.densityScale = muA / (self.lengthScale / self.timeScale)**2

    self.dbProperties = numpy.array([ [densityA, vsA, vpA, viscosityA],
                                      [densityB, vsB, vpB, viscosityB] ], 
                                    dtype=numpy.float64)
    self.properties = numpy.array([ [densityA, muA, lambdaA, maxwellTimeA],
                                    [densityB, muB, lambdaB, maxwellTimeB] ],
                                     dtype=numpy.float64)

    # TEMPORARY, need to determine how to use initial state variables
    self.dbStateVars = numpy.zeros( (numLocs, tensorSize+tensorSize),
                                    dtype=numpy.float64)
    self.stateVars = numpy.zeros( (numLocs, tensorSize+tensorSize),
                                  dtype=numpy.float64)

    mu0 = self.pressureScale
    density0 = self.densityScale
    time0 = self.timeScale
    self.propertiesNondim = \
        numpy.array([ [densityA/density0, muA/mu0, lambdaA/mu0, maxwellTimeA/time0],
                      [densityB/density0, muB/mu0, lambdaB/mu0, maxwellTimeB/time0] ],
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
    self.stateVarsUpdated = numpy.zeros( (numLocs, tensorSize + tensorSize), \
                                         dtype=numpy.float64)

    (self.elasticConsts[0,:], self.stress[0,:], self.stateVarsUpdated[0,:]) = \
        self._calcStress(strainA, muA, lambdaA, \
                           initialStressA, initialStrainA)
    (self.elasticConsts[1,:], self.stress[1,:], self.stateVarsUpdated[1,:]) = \
        self._calcStress(strainB, muB, lambdaB, \
                           initialStressB, initialStrainB)
    self.dtStableImplicit = 0.2*min(maxwellTimeA, maxwellTimeB)
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

    initialStress = numpy.reshape(initialStressV, (tensorSize,1))
    initialStrain = numpy.reshape(initialStrainV, (tensorSize,1))
    strain = numpy.reshape(strainV, (6,1)) - initialStrain
    elastic = numpy.array([ [C1111, C1122, C1133, C1112, C1123, C1113],
                            [C2211, C2222, C2233, C2212, C2223, C2213],
                            [C3311, C3322, C3333, C3312, C3323, C3313],
                            [C1211, C1222, C1233, C1212, C1223, C1213],
                            [C2311, C2322, C2333, C2312, C2323, C2313],
                            [C1311, C1322, C1333, C1312, C1323, C1313] ],
                          dtype=numpy.float64)
    stress = numpy.dot(elastic, strain) + initialStress
    meanStrain = (strain[0] + strain[1] + strain[2])/3.0
    strainVec = numpy.array(strainV, dtype=numpy.float64)
    viscousStrain = [strain[0] - meanStrain,
                     strain[1] - meanStrain,
                     strain[2] - meanStrain,
                     strain[3],
                     strain[4],
                     strain[5]]
    stateVarsUpdated = numpy.array( [strainVec, viscousStrain],
                                    dtype=numpy.float64)
    return (elasticConsts, numpy.ravel(stress), numpy.ravel(stateVarsUpdated))
  

# MAIN /////////////////////////////////////////////////////////////////
if __name__ == "__main__":

  app = MaxwellIsotropic3DElastic()
  app.run()


# End of file 
