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

## @file unittests/libtests/materials/data/DruckerPragerPlaneStrainElastic.py

## @brief Python application for generating C++ data files for testing
## C++ DruckerPragerPlaneStrain object with elastic behavior.

from ElasticMaterialApp import ElasticMaterialApp

import numpy
import math

# ----------------------------------------------------------------------
dimension = 2
numElasticConsts = 9
tensorSize = 3

# DruckerPragerPlaneStrainElastic class
class DruckerPragerPlaneStrainElastic(ElasticMaterialApp):
  """
  Python application for generating C++ data files for testing C++
  DruckerPragerPlaneStrain object with elastic behavior.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="druckerpragerplanestrainelastic"):
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
                             "friction-angle", "cohesion",
                             "dilatation-angle"]
    self.numPropertyValues = numpy.array([1, 1, 1, 1, 1, 1], dtype=numpy.int32)

    self.dbStateVarValues = ["stress-zz-initial",
                             "plastic-strain-xx",
                             "plastic-strain-yy",
                             "plastic-strain-zz",
                             "plastic-strain-xy"
                             ]
    self.numStateVarValues = numpy.array([1, 4], dtype=numpy.int32)

    densityA = 2500.0
    vsA = 3000.0
    vpA = vsA*3**0.5
    # First case has different values for friction angle and dilatation angle.
    frictionAngleA = math.radians(30.0)
    dilatationAngleA = math.radians(20.0)
    cohesionA = 3.0e5
    strainA = [1.1e-4, 1.2e-4, 1.4e-4]
    initialStressA = [2.1e4, 2.2e4, 2.4e4]
    initialStrainA = [3.1e-4, 3.2e-4, 3.4e-4]
    muA = vsA*vsA*densityA
    lambdaA = vpA*vpA*densityA - 2.0*muA
    stressZZInitialA = 1.075e+4

    denomFrictionA = math.sqrt(3.0) * (3.0 - math.sin(frictionAngleA))
    denomDilatationA = math.sqrt(3.0) * (3.0 - math.sin(dilatationAngleA))
    alphaYieldA = 2.0 * math.sin(frictionAngleA)/denomFrictionA
    betaA = 6.0 * cohesionA * math.cos(frictionAngleA)/denomFrictionA
    alphaFlowA = 2.0 * math.sin(dilatationAngleA)/denomDilatationA
    
    densityB = 2000.0
    vsB = 1200.0
    vpB = vsB*3**0.5
    # Second case has same values for friction angle and dilatation angle.
    frictionAngleB = math.radians(25.0)
    dilatationAngleB = math.radians(25.0)
    cohesionB = 1.0e5
    strainB = [4.1e-4, 4.2e-4, 4.4e-4]
    initialStressB = [5.1e4, 5.2e4, 5.4e4]
    initialStrainB = [6.1e-4, 6.2e-4, 6.4e-4]
    muB = vsB*vsB*densityB
    lambdaB = vpB*vpB*densityB - 2.0*muB
    denomFrictionB = math.sqrt(3.0) * (3.0 - math.sin(frictionAngleB))
    denomDilatationB = math.sqrt(3.0) * (3.0 - math.sin(dilatationAngleB))
    alphaYieldB = 2.0 * math.sin(frictionAngleB)/denomFrictionB
    betaB = 6.0 * cohesionB * math.cos(frictionAngleB)/denomFrictionB
    alphaFlowB = 2.0 * math.sin(dilatationAngleB)/denomDilatationB
    stressZZInitialB = 2.575e+4

    self.lengthScale = 1.0e+3
    self.pressureScale = muA
    self.timeScale = 1.0
    self.densityScale = muA / (self.lengthScale / self.timeScale)**2

    self.dbProperties = numpy.array([ [densityA, vsA, vpA, \
                                       frictionAngleA, cohesionA, \
                                       dilatationAngleA],
                                      [densityB, vsB, vpB, \
                                       frictionAngleB, cohesionB, \
                                       dilatationAngleB] ], 
                                    dtype=numpy.float64)
    self.properties = numpy.array([ [densityA, muA, lambdaA, \
                                     alphaYieldA, betaA, \
                                     alphaFlowA],
                                    [densityB, muB, lambdaB, \
                                     alphaYieldB, betaB, \
                                     alphaFlowB] ],
                                     dtype=numpy.float64)

    # TEMPORARY, need to determine how to use initial state variables
    self.dbStateVars = numpy.zeros( (numLocs, 1 + 4), dtype=numpy.float64)
    self.dbStateVars[0, 0] = stressZZInitialA
    self.dbStateVars[1, 0] = stressZZInitialB

    self.stateVars = numpy.zeros( (numLocs, 1 + 4), dtype=numpy.float64)
    self.stateVars[0, 0] = stressZZInitialA
    self.stateVars[1, 0] = stressZZInitialB

    mu0 = self.pressureScale
    density0 = self.densityScale
    self.propertiesNondim = \
        numpy.array([ [densityA/density0, muA/mu0, lambdaA/mu0, \
                       alphaYieldA, betaA/mu0, \
                       alphaFlowA],
                      [densityB/density0, muB/mu0, lambdaB/mu0, \
                       alphaYieldB, betaB/mu0, \
                       alphaFlowB] ],
                    dtype=numpy.float64)

    stressZZInitialANondim = stressZZInitialA/mu0
    stressZZInitialBNondim = stressZZInitialB/mu0

    self.stateVarsNondim = numpy.zeros( (numLocs, 1 + 4), dtype=numpy.float64)
    self.stateVarsNondim[0, 0] = stressZZInitialANondim
    self.stateVarsNondim[1, 0] = stressZZInitialBNondim
    
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

    self.dtStableImplicit = 1.0e10
    self.dtStableExplicit = 1000.0 / vpA

    stateVarsUpdatedA = numpy.array([stressZZInitialA, 0.0, 0.0, 0.0, 0.0],
                                    dtype=numpy.float64)
    stateVarsUpdatedB = numpy.array([stressZZInitialB, 0.0, 0.0, 0.0, 0.0],
                                    dtype=numpy.float64)

    self.stateVarsUpdated = numpy.array( [stateVarsUpdatedA, stateVarsUpdatedB],
                                         dtype=numpy.float64)

    return


  def _calcStress(self, strainV, muV, lambdaV, initialStressV, initialStrainV):
    """
    Compute stress and derivative of elasticity matrix.
    """
    C1111 = lambdaV + 2.0*muV
    C1122 = lambdaV
    C1112 = 0.0
    C2211 = lambdaV
    C2222 = lambdaV + 2.0*muV
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
                            [C1211, C1222, C1212] ], dtype=numpy.float64)
    stress = numpy.dot(elastic, strain-initialStrain) + initialStress
    return (elasticConsts, numpy.ravel(stress))
  

# MAIN /////////////////////////////////////////////////////////////////
if __name__ == "__main__":

  app = DruckerPragerPlaneStrainElastic()
  app.run()


# End of file 
