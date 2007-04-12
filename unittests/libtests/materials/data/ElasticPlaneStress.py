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

## @file unittests/libtests/materials/data/ElasticPlaneStrain.py

## @brief Python application for generating C++ data files for testing
## C++ ElasticPlaneStrain object.

from ElasticMaterialApp import ElasticMaterialApp

import numpy

# ----------------------------------------------------------------------
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

    self.dimension = 2

    self.numDBValues = 3
    self.dbValues = ["density", "vs", "vp"]
    self.numParameters = 3
    self.parameterNames = ["density", "mu", "lambda"]

    densityA = 2500.0
    vsA = 3000.0
    vpA = vsA*3**0.5
    strainA = [1.1e-4, 2.2e-4, 3.3e-4]
    
    densityB = 2000.0
    vsB = 1200.0
    vpB = vsB*3**0.5
    strainB = [1.2e-4, 2.3e-4, 3.4e-4]

    self.dbData = numpy.array([ [densityA, vsA, vpA],
                                [densityB, vsB, vpB] ],
                              dtype=numpy.float64)
    muA = vsA*vsA*densityA
    lambdaA = vpA*vpA*densityA - 2.0*muA
    muB = vsB*vsB*densityB
    lambdaB = vpB*vpB*densityB - 2.0*muB
    self.parameterData = numpy.array([ [densityA, muA, lambdaA],
                                       [densityB, muB, lambdaB] ],
                                     dtype=numpy.float64)
    
    self.numLocs = 2
    numElasticConsts = 6
    self.density = numpy.array([densityA, densityB],
                               dtype=numpy.float64)

    self.strain = numpy.array([strainA, strainB],
                               dtype=numpy.float64)
    self.stress = numpy.zeros( (self.numLocs, 3), dtype=numpy.float64)
    self.elasticConsts = numpy.zeros( (self.numLocs, numElasticConsts),
                                      dtype=numpy.float64)

    (self.elasticConsts[0,:], self.stress[0,:]) = \
                              self._calcStress(strainA, densityA, muA, lambdaA)
    (self.elasticConsts[1,:], self.stress[1,:]) = \
                              self._calcStress(strainB, densityB, muB, lambdaB)
    return


  def _calcStress(self, strainV, densityV, muV, lambdaV):
    """
    Compute stress and derivative of elasticity matrix.
    """
    C1111 = 4*muV * (lambdaV + muV) / (lambdaV + 2*muV)
    C1122 = 2*muV*lambdaV / (lambdaV + 2*muV)
    C1112 = 0.0
    C2222 = 4*muV * (lambdaV + muV) / (lambdaV + 2*muV)
    C2212 = 0.0
    C1212 = 2.0*muV
    elasticConsts = numpy.array([C1111, C1122, C1112,
                                 C2222, C2212,
                                 C1212], dtype=numpy.float64)

    strain = numpy.reshape(strainV, (3,1))
    elastic = numpy.array([ [C1111, C1122, C1112],
                            [C1122, C2222, C2212],
                            [C1112, C2212, C1212] ],
                          dtype=numpy.float64)
    stress = numpy.dot(elastic, strain)
    return (elasticConsts, numpy.ravel(stress))
  

# MAIN /////////////////////////////////////////////////////////////////
if __name__ == "__main__":

  app = ElasticPlaneStrain()
  app.run()


# End of file 
