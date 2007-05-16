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

## @file unittests/libtests/materials/data/ElasticStrain1D.py

## @brief Python application for generating C++ data files for testing
## C++ ElasticStrain1D object.

from ElasticMaterialApp import ElasticMaterialApp

import numpy

# ----------------------------------------------------------------------
# ElasticStrain1D class
class ElasticStrain1D(ElasticMaterialApp):
  """
  Python application for generating C++ data files for testing C++
  ElasticStrain1D object.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="elasticstrain1d"):
    """
    Constructor.
    """
    ElasticMaterialApp.__init__(self, name)

    self.dimension = 1

    self.numDBValues = 2
    self.dbValues = ["density", "vp"]
    self.numParameters = 2
    self.numParamValues = [1, 1]
    self.parameterNames = ["density", "lambda2mu"]

    densityA = 2500.0
    vpA = 5000.0
    strainA = [1.1e-4]
    
    densityB = 2000.0
    vpB = 3000.0
    strainB = [1.2e-4]

    self.dbData = numpy.array([ [densityA, vpA],
                                [densityB, vpB] ],
                              dtype=numpy.float64)
    lambda2muA = vpA*vpA*densityA
    lambda2muB = vpB*vpB*densityB
    self.parameterData = numpy.array([ [densityA, lambda2muA],
                                       [densityB, lambda2muB] ],
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
                              self._calcStress(strainA, densityA, lambda2muA)
    (self.elasticConsts[1,:], self.stress[1,:]) = \
                              self._calcStress(strainB, densityB, lambda2muB)
    return


  def _calcStress(self, strainV, densityV, lambda2muV):
    """
    Compute stress and derivative of elasticity matrix.
    """
    C1111 = lambda2muV
    elasticConsts = numpy.array([C1111], dtype=numpy.float64)

    strain = numpy.reshape(strainV, (1,1))
    elastic = numpy.array([ [C1111] ],
                          dtype=numpy.float64)
    stress = numpy.dot(elastic, strain)
    return (elasticConsts, numpy.ravel(stress))
  

# MAIN /////////////////////////////////////////////////////////////////
if __name__ == "__main__":

  app = ElasticStrain1D()
  app.run()


# End of file 
