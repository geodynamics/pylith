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

## @file unittests/libtests/materials/data/ElasticIsotropic3D.py

## @brief Python application for generating C++ data files for testing
## C++ ElasticIsotropic3D object.

from ElasticMaterialApp import ElasticMaterialApp

import numpy

# ----------------------------------------------------------------------
# ElasticIsotropic3D class
class ElasticIsotropic3D(ElasticMaterialApp):
  """
  Python application for generating C++ data files for testing C++
  ElasticIsotropic3D object.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="elasticisotropic3d"):
    """
    Constructor.
    """
    ElasticMaterialApp.__init__(self, name)

    self.numDBValues = 3
    self.dbValues = ["density", "vs", "vp"]
    self.numParameters = 3
    self.parameterNames = ["density", "mu", "lambda"]

    densityA = 2500.0
    vsA = 3000.0
    vpA = vsA*3**0.5
    densityB = 2000.0
    vsB = 1200.0
    vpB = vsB*3**0.5
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
    numElasticConsts = 21
    self.density = numpy.array([densityA, densityB],
                               dtype=numpy.float64)
    self.elasticConsts = numpy.zeros( (self.numLocs, numElasticConsts),
                                      dtype=numpy.float64)
    self.elasticConsts[0,:] = self._calcElasticConsts(densityA, muA, lambdaA)
    self.elasticConsts[1,:] = self._calcElasticConsts(densityB, muB, lambdaB)
    return


  def _calcElasticConsts(self, densityV, muV, lambdaV):
    """
    Compute elastic constants.
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
    C1212 = muV
    C1223 = 0.0
    C1213 = 0.0
    C2323 = muV
    C2313 = 0.0
    C1313 = muV
    elasticConsts = numpy.array([C1111, C1122, C1133, C1112, C1123, C1113,
                                 C2222, C2233, C2212, C2223, C2213,
                                 C3333, C3312, C3323, C3313,
                                 C1212, C1223, C1213,
                                 C2323, C2313,
                                 C1313], dtype=numpy.float64)
    return elasticConsts
  

# MAIN /////////////////////////////////////////////////////////////////
if __name__ == "__main__":

  app = ElasticIsotropic3D()
  app.run()


# End of file 
