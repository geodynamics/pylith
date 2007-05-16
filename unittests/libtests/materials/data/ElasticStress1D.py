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

## @file unittests/libtests/materials/data/ElasticStress1D.py

## @brief Python application for generating C++ data files for testing
## C++ ElasticStress1D object.

from ElasticMaterialApp import ElasticMaterialApp

import numpy

# ----------------------------------------------------------------------
# ElasticStress1D class
class ElasticStress1D(ElasticMaterialApp):
  """
  Python application for generating C++ data files for testing C++
  ElasticStress1D object.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="elasticstress1d"):
    """
    Constructor.
    """
    ElasticMaterialApp.__init__(self, name)

    self.dimension = 1

    self.numDBValues = 3
    self.dbValues = ["density", "vs", "vp"]
    self.numParameters = 3
    self.numParamValues = [1, 1, 1]
    self.parameterNames = ["density", "mu", "lambda"]

    densityA = 2500.0
    vsA = 3000.0
    vpA = vsA*3**0.5
    strainA = [1.1e-4]
    
    densityB = 2000.0
    vsB = 1200.0
    vpB = vsB*3**0.5
    strainB = [1.2e-4]

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
    numElasticConsts = 1
    self.density = numpy.array([densityA, densityB],
                               dtype=numpy.float64)

    self.strain = numpy.array([strainA, strainB],
                               dtype=numpy.float64)
    self.stress = numpy.zeros( (self.numLocs, 1), dtype=numpy.float64)
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
    C1111 = muV * (3*lambdaV + 2*muV) / (lambdaV + muV)
    elasticConsts = numpy.array([C1111], dtype=numpy.float64)

    strain = numpy.reshape(strainV, (1,1))
    elastic = numpy.array([ [C1111] ],
                          dtype=numpy.float64)
    stress = numpy.dot(elastic, strain)
    return (elasticConsts, numpy.ravel(stress))
  

# MAIN /////////////////////////////////////////////////////////////////
if __name__ == "__main__":

  app = ElasticStress1D()
  app.run()


# End of file 
