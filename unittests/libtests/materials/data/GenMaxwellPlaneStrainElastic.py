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

## @file unittests/libtests/materials/data/GenMaxwellPlaneStrainElastic.py

## @brief Python application for generating C++ data files for testing
## C++ GenMaxwellPlaneStrain object with elastic behavior.

from ElasticMaterialApp import ElasticMaterialApp

import numpy

# ----------------------------------------------------------------------
dimension = 2
numElasticConsts = 9
tensorSize = 3
numMaxwellModels = 3

# GenMaxwellPlaneStrainElastic class
class GenMaxwellPlaneStrainElastic(ElasticMaterialApp):
  """
  Python application for generating C++ data files for testing C++
  GenMaxwellPlaneStrain object with elastic behavior.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="genmaxwellplanestrainelastic"):
    """
    Constructor.
    """
    ElasticMaterialApp.__init__(self, name)

    # import pdb
    # pdb.set_trace()
    
    numLocs = 2

    self.dimension = dimension
    self.numLocs = numLocs

    self.dbPropertyValues = ["density", "vs", "vp",
                             "shear-ratio-1", "shear-ratio-2", "shear-ratio-3",
                             "viscosity-1", "viscosity-2", "viscosity-3"]
    self.numPropertyValues = numpy.array([1, 1, 1, 1, 1, 1, 1, 1, 1],
                                         dtype=numpy.int32)

    self.dbStateVarValues = ["stress-zz-initial",
                             "total-strain-xx",
                             "total-strain-yy",
                             "total-strain-xy",
                             "viscous-strain-1-xx",
                             "viscous-strain-1-yy",
                             "viscous-strain-1-zz",
                             "viscous-strain-1-xy",
                             "viscous-strain-2-xx",
                             "viscous-strain-2-yy",
                             "viscous-strain-2-zz",
                             "viscous-strain-2-xy",
                             "viscous-strain-3-xx",
                             "viscous-strain-3-yy",
                             "viscous-strain-3-zz",
                             "viscous-strain-3-xy"
                             ]
    self.numStateVarValues = numpy.array([1, tensorSize, 4, 4, 4],
                                         dtype=numpy.int32)

    densityA = 2500.0
    vsA = 3000.0
    vpA = vsA*3**0.5
    shearRatioA = [0.5, 0.1, 0.2]
    viscosityA = [1.0e18, 1.0e17, 1.0e19]
    strainA = [1.1e-4, 1.2e-4, 1.4e-4]
    initialStressA = [2.1e4, 2.2e4, 2.4e4]
    initialStrainA = [3.1e-5, 3.2e-5, 3.4e-5]
    muA = vsA*vsA*densityA
    lambdaA = vpA*vpA*densityA - 2.0*muA
    stressInitialZZA = numpy.array([1.5e4], dtype=numpy.float64)


    densityB = 2000.0
    vsB = 1200.0
    vpB = vsB*3**0.5
    shearRatioB = [0.2, 0.2, 0.2]
    viscosityB = [1.0e18, 1.0e19, 1.0e20]
    strainB = [4.1e-4, 4.2e-4, 4.4e-4]
    initialStressB = [5.1e4, 5.2e4, 5.4e4]
    initialStrainB = [6.1e-5, 6.2e-5, 6.4e-5]
    muB = vsB*vsB*densityB
    lambdaB = vpB*vpB*densityB - 2.0*muB
    stressInitialZZB = numpy.array([4.5e4], dtype=numpy.float64)

    maxwellTimeA = [1.0e30, 1.0e30, 1.0e30]
    maxwellTimeB = [1.0e30, 1.0e30, 1.0e30]
    for i in xrange(numMaxwellModels):
      if shearRatioA[i] != 0.0:
        maxwellTimeA[i] = viscosityA[i]/(muA*shearRatioA[i])
      if shearRatioB[i] != 0.0:
        maxwellTimeB[i] = viscosityB[i]/(muB*shearRatioB[i])

    self.lengthScale = 1.0e+3
    self.pressureScale = muA
    self.timeScale = 1.0
    self.densityScale = muA / (self.lengthScale / self.timeScale)**2

    propA = [densityA, vsA, vpA] + shearRatioA + viscosityA
    propB = [densityB, vsB, vpB] + shearRatioB + viscosityB
    self.dbProperties = numpy.array([propA, propB], dtype=numpy.float64)
    propA = [densityA, muA, lambdaA] + shearRatioA + maxwellTimeA
    propB = [densityB, muB, lambdaB] + shearRatioB + maxwellTimeB
    self.properties = numpy.array([propA, propB], dtype=numpy.float64)

    # TEMPORARY, need to determine how to use initial state variables
    # At present, only the first (stressInitialZZ) is being used.
    self.dbStateVars = numpy.zeros( (numLocs,
                                     1 + tensorSize + 4 * numMaxwellModels),
                                    dtype=numpy.float64)
    self.dbStateVars[0, 0] = stressInitialZZA
    self.dbStateVars[1, 0] = stressInitialZZB
    
    self.stateVars = numpy.zeros( (numLocs,
                                   1 + tensorSize + 4 * numMaxwellModels),
                                  dtype=numpy.float64)
    self.stateVars[0, 0] = stressInitialZZA
    self.stateVars[1, 0] = stressInitialZZB

    mu0 = self.pressureScale
    density0 = self.densityScale
    time0 = self.timeScale
    self.propertiesNondim = \
        numpy.array([ [densityA/density0, muA/mu0, lambdaA/mu0,
                       shearRatioA[0], shearRatioA[1], shearRatioA[2],
                       maxwellTimeA[0]/time0, maxwellTimeA[1]/time0,
                       maxwellTimeA[2]/time0],
                      [densityB/density0, muB/mu0, lambdaB/mu0,
                       shearRatioB[0], shearRatioB[1], shearRatioB[2],
                       maxwellTimeB[0]/time0, maxwellTimeB[1]/time0,
                       maxwellTimeB[2]/time0] ],
                    dtype=numpy.float64)

    stressInitialZZANondim = stressInitialZZA/mu0
    stressInitialZZBNondim = stressInitialZZB/mu0
    
    self.stateVarsNondim = numpy.zeros( (numLocs,
                                         1 + tensorSize + 4 * numMaxwellModels),
                                        dtype=numpy.float64)

    self.stateVarsNondim[0, 0] = stressInitialZZANondim
    self.stateVarsNondim[1, 0] = stressInitialZZBNondim

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
    self.stateVarsUpdated = numpy.zeros((self.numLocs,
                                         1 + tensorSize + 4 * numMaxwellModels),
                                        dtype=numpy.float64)

    (self.elasticConsts[0,:], self.stress[0,:], self.stateVarsUpdated[0,:]) = \
                              self._calcStress(strainA, muA, lambdaA,
                                               initialStressA, initialStrainA,
                                               self.stateVars[0,:])
    (self.elasticConsts[1,:], self.stress[1,:], self.stateVarsUpdated[1,:]) = \
                              self._calcStress(strainB, muB, lambdaB,
                                               initialStressB, initialStrainB,
                                               self.stateVars[1,:])
    self.dtStableImplicit = 0.2*min(min(maxwellTimeA), min(maxwellTimeB))
    self.dtStableExplicit = 1000.0 / vpA

    return


  def _calcStress(self, strainV, muV, lambdaV, initialStressV, initialStrainV,
                  stateVars):
    """
    Compute stress and derivative of elasticity matrix.
    """
    C1111 = lambdaV + 2.0 * muV
    C1122 = lambdaV
    C1112 = 0.0
    C2211 = lambdaV
    C2222 = lambdaV + 2.0 * muV
    C2212 = 0.0
    C1211 = 0.0
    C1222 = 0.0
    C1212 = 2.0 * muV
    elasticConsts = numpy.array([C1111, C1122, C1112,
                                 C2211, C2222, C2212,
                                 C1211, C1222, C1212],
                                dtype=numpy.float64)

    initialStress = numpy.reshape(initialStressV, (tensorSize,1))
    initialStrain = numpy.reshape(initialStrainV, (tensorSize,1))
    stressZZInitial = numpy.array([stateVars[0]], dtype=numpy.float64)
    strain = numpy.reshape(strainV, (tensorSize,1)) - initialStrain
    
    elastic = numpy.array([ [C1111, C1122, C1112],
                            [C2211, C2222, C2212],
                            [C1211, C1222, C1212] ],
                          dtype=numpy.float64)
    stress = numpy.dot(elastic, strain) + initialStress
    meanStrain = (strain[0,0] + strain[1,0])/3.0
    viscousStrain = numpy.array([strain[0,0] - meanStrain,
                                 strain[1,0] - meanStrain,
                                 -meanStrain,
                                 strain[2,0]],
                                dtype=numpy.float64)
    viscousStrainVec = numpy.ravel(viscousStrain)
    strainVec = numpy.array(strainV, dtype=numpy.float64)

    stateVarsUpdated = numpy.concatenate((stressZZInitial, strainVec,
                                          viscousStrainVec, viscousStrainVec,
                                          viscousStrainVec))
    return (elasticConsts, numpy.ravel(stress), numpy.ravel(stateVarsUpdated))
  

# MAIN /////////////////////////////////////////////////////////////////
if __name__ == "__main__":

  app = GenMaxwellPlaneStrainElastic()
  app.run()


# End of file 
