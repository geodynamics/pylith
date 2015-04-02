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

## @file unittests/libtests/materials/data/ElasticPlaneStrain.py

## @brief Python application for generating C++ data files for testing
## C++ ElasticPlaneStrain object.

from ElasticMaterialApp import ElasticMaterialApp

import numpy

# ----------------------------------------------------------------------
dimension = 2
numElasticConsts = 9
tensorSize = 3

# ElasticPlaneStrain class
class ElasticPlaneStrain(ElasticMaterialApp):
  """
  Python application for generating C++ data files for testing C++
  ElasticPlaneStrain object.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="elasticplanestrain"):
    """
    Constructor.
    """
    ElasticMaterialApp.__init__(self, name)

    numLocs = 2

    self.dimension = dimension
    self.numLocs = numLocs
    
    self.dbPropertyValues = ["density", "vs", "vp"]    
    self.propertyValues = ["density", "mu", "lambda"]
    self.numPropertyValues = numpy.array([1, 1, 1], dtype=numpy.int32)

    densityA = 2500.0
    vsA = 3000.0
    vpA = vsA*3**0.5
    strainA = [1.1e-4, 1.2e-4, 1.3e-4]
    initialStressA = [2.1e4, 2.2e4, 2.3e4]
    initialStrainA = [3.1e-4, 3.2e-4, 3.3e-4]
    muA = vsA*vsA*densityA
    lambdaA = vpA*vpA*densityA - 2.0*muA
    
    densityB = 2000.0
    vsB = 1200.0
    vpB = vsB*3**0.5
    strainB = [4.1e-4, 4.2e-4, 4.3e-4]
    initialStressB = [5.1e4, 5.2e4, 5.3e4]
    initialStrainB = [6.1e-4, 6.2e-4, 6.3e-4]
    muB = vsB*vsB*densityB
    lambdaB = vpB*vpB*densityB - 2.0*muB

    self.lengthScale = 1.0e+3
    self.pressureScale = muA
    self.timeScale = 1.0
    self.densityScale = muA / (self.lengthScale / self.timeScale)**2
    
    self.dbProperties = numpy.array([ [densityA, vsA, vpA],
                                      [densityB, vsB, vpB] ], 
                                    dtype=numpy.float64)
    self.properties = numpy.array([ [densityA, muA, lambdaA],                                                       [densityB, muB, lambdaB] ],
                                     dtype=numpy.float64)

    mu0 = self.pressureScale
    density0 = self.densityScale
    self.propertiesNondim = \
        numpy.array([ [densityA/density0, muA/mu0, lambdaA/mu0],
                      [densityB/density0, muB/mu0, lambdaB/mu0] ],
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

    self.strain = numpy.array([strainA,
                               strainB],
                               dtype=numpy.float64)
    
    stress = numpy.zeros( (numLocs, tensorSize), dtype=numpy.float64)
    elasticConsts = numpy.zeros( (numLocs, numElasticConsts),
                                 dtype=numpy.float64)

    (elasticConsts[0,:], stress[0,:]) = \
        self._calcStress(strainA, densityA, muA, lambdaA,
                         initialStressA, initialStrainA)
    (elasticConsts[1,:], stress[1,:]) = \
        self._calcStress(strainB, densityB, muB, lambdaB,
                         initialStressB, initialStrainB)

    self.stress = stress
    self.elasticConsts = elasticConsts

    self.dtStableExplicit = 1000.0 / vpA

    return


  def _calcStress(self, strainV, densityV, muV, lambdaV,
                  initialStressV, initialStrainV):
    """
    Compute stress and derivative of elasticity matrix.
    """
    C1111 = 4*muV * (lambdaV + muV) / (lambdaV + 2*muV)
    C1122 = 2*muV*lambdaV / (lambdaV + 2*muV)
    C1112 = 0.0
    C2211 = 2*muV*lambdaV / (lambdaV + 2*muV)
    C2222 = 4*muV * (lambdaV + muV) / (lambdaV + 2*muV)
    C2212 = 0.0
    C1211 = 0.0
    C1222 = 0.0
    C1212 = 2.0*muV
    elasticConsts = numpy.array([C1111, C1122, C1112,
                                 C2211, C2222, C2212,
                                 C1211, C1222, C1212], dtype=numpy.float64)

    strain = numpy.reshape(strainV, (tensorSize,1))
    initialStress = numpy.reshape(initialStressV, (tensorSize,1))
    initialStrain = numpy.reshape(initialStrainV, (tensorSize,1))
    elastic = numpy.array([ [C1111, C1122, C1112],
                            [C2211, C2222, C2212],
                            [C1211, C1222, C1212] ],
                          dtype=numpy.float64)
    stress = numpy.dot(elastic,strain-initialStrain) + initialStress
    return (elasticConsts, numpy.ravel(stress))
  

# MAIN /////////////////////////////////////////////////////////////////
if __name__ == "__main__":

  app = ElasticPlaneStrain()
  app.run()


# End of file 
