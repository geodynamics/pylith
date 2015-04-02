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

## @file unittests/libtests/materials/data/DruckerPragerPlaneStrainTimeDep.py

## @brief Python application for generating C++ data files for testing
## C++ DruckerPragerPlaneStrain object with time dependent behavior.

from ElasticMaterialApp import ElasticMaterialApp

import numpy
import math

# ----------------------------------------------------------------------
dimension = 2
numElasticConsts = 9
tensorSize = 3

# DruckerPragerPlaneStrainTimeDep class
class DruckerPragerPlaneStrainTimeDep(ElasticMaterialApp):
  """
  Python application for generating C++ data files for testing C++
  DruckerPragerPlaneStrain object with time dependent behavior.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="druckerpragerplanestraintimedep"):
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
    self.propertyValues = ["density", "mu", "lambda",
                           "alpha_yield", "beta", "alpha_flow"]
    self.numPropertyValues = numpy.array([1, 1, 1, 1, 1, 1], dtype=numpy.int32)

    self.dbStateVarValues = ["stress-zz-initial",
                             "plastic-strain-xx",
                             "plastic-strain-yy",
                             "plastic-strain-zz",
                             "plastic-strain-xy"
                             ]
    self.numStateVarValues = numpy.array([1, 4], dtype=numpy.int32)

    self.dt = 2.0e5

    densityA = 2500.0
    vsA = 3000.0
    vpA = vsA*3**0.5
    # First case has same values for friction angle and dilatation angle.
    frictionAngleA = math.radians(30.0)
    dilatationAngleA = math.radians(20.0)
    cohesionA = 3.0e5
    strainA = [-2.1e-4, 1.2e-4, 1.1e-5]
    initialStressA = [2.1e4, 2.2e4, 2.4e4]
    initialStrainA = [3.6e-5, 3.5e-5, 3.3e-5]
    muA = vsA*vsA*densityA
    lambdaA = vpA*vpA*densityA - 2.0*muA
    stressZZInitialA = 2.3e4

    denomFrictionA = math.sqrt(3.0) * (3.0 - math.sin(frictionAngleA))
    denomDilatationA = math.sqrt(3.0) * (3.0 - math.sin(dilatationAngleA))
    alphaYieldA = 2.0 * math.sin(frictionAngleA)/denomFrictionA
    betaA = 6.0 * cohesionA * math.cos(frictionAngleA)/denomFrictionA
    alphaFlowA = 2.0 * math.sin(dilatationAngleA)/denomDilatationA
    
    densityB = 2000.0
    vsB = 1200.0
    vpB = vsB*3**0.5
    # Second case has different values for friction angle and dilatation angle.
    frictionAngleB = math.radians(25.0)
    dilatationAngleB = math.radians(25.0)
    cohesionB = 1.0e4
    strainB = [4.1e-4, 4.2e-4, 1.4e-4]
    initialStressB = [5.6e4, 5.5e4, 5.3e4]
    initialStrainB = [6.6e-5, 6.5e-5, 6.2e-5]
    muB = vsB*vsB*densityB
    lambdaB = vpB*vpB*densityB - 2.0*muB
    stressZZInitialB = 5.4e4

    denomFrictionB = math.sqrt(3.0) * (3.0 - math.sin(frictionAngleB))
    denomDilatationB = math.sqrt(3.0) * (3.0 - math.sin(dilatationAngleB))
    alphaYieldB = 2.0 * math.sin(frictionAngleB)/denomFrictionB
    betaB = 6.0 * cohesionB * math.cos(frictionAngleB)/denomFrictionB
    alphaFlowB = 2.0 * math.sin(dilatationAngleB)/denomDilatationB

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

    stateVarsA = [stressZZInitialA, 4.1e-5, 4.2e-5, 4.4e-5, 4.5e-5]
    stateVarsB = [stressZZInitialB, 1.1e-5, 1.2e-5, 1.4e-5, 1.5e-5]
    plasStrainA = stateVarsA[1:]
    plasStrainB = stateVarsB[1:]
    self.stateVars = numpy.array([[stateVarsA], [stateVarsB] ],
                                 dtype=numpy.float64)

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
    self.stateVarsNondim = numpy.array([[stateVarsA], [stateVarsB] ],
                                       dtype=numpy.float64)
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
    self.stateVarsUpdated = numpy.zeros( (numLocs, 5),
                                         dtype=numpy.float64)
    self.elasticConsts = numpy.zeros( (self.numLocs, numElasticConsts), \
                                        dtype=numpy.float64)

    (self.elasticConsts[0,:], self.stress[0,:], self.stateVarsUpdated[0,:]) = \
                              self._calcStress(strainA, muA, lambdaA,
                                               alphaYieldA, betaA, alphaFlowA,
                                               plasStrainA,
                                               initialStressA, initialStrainA,
                                               stateVarsA)
    (self.elasticConsts[1,:], self.stress[1,:], self.stateVarsUpdated[1,:]) = \
                              self._calcStress(strainB, muB, lambdaB,
                                               alphaYieldB, betaB, alphaFlowB,
                                               plasStrainB,
                                               initialStressB, initialStrainB,
                                               stateVarsB)

    self.dtStableImplicit = 1.0e+99
    self.dtStableExplicit = 1000.0 / vpA

    return


  def _scalarProduct(self, tensor1, tensor2):
    """
    Compute the scalar product of two tensors stored in vector form (length 4).
    """
    scalarProduct = tensor1[0] * tensor2[0] + \
                    tensor1[1] * tensor2[1] + \
                    tensor1[2] * tensor2[2] + \
                    2.0 * tensor1[3] * tensor2[3]
    return scalarProduct


  def _calcStressComponent(self, strainVal, strainComp, stressComp, strainTpdt,
                           muV, lambdaV, alphaYieldV, betaV, alphaFlowV,
                           plasStrainT, initialStress, initialStrain,
                           stateVars):
    """
    Function to compute a particular stress component as a function of a
    strain component.
    """
    strainTest = numpy.array(strainTpdt, dtype=numpy.float64)
    strainTest[strainComp] = strainVal
    stressTpdt, visStrainTpdt = self._computeStress(strainTest, muV, lambdaV,
                                                    alphaYieldV, betaV,
                                                    alphaFlowV,
                                                    plasStrainT,
                                                    initialStress,
                                                    initialStrain,
                                                    stateVars)
    return stressTpdt[stressComp]


  def _computeStress(self, strainTpdt, muV, lambdaV, alphaYieldV, betaV,
                     alphaFlowV, plasStrainT, initialStress, initialStrain,
                     stateVars):
    """
    Function to compute stresses and plastic strains.
    """
    
    # Constants
    mu2 = 2.0 * muV
    lamPlusMu = lambdaV + muV
    bulkModulus = lambdaV + mu2/3.0
    ae = 1.0/mu2
    am = 1.0/(3.0 * bulkModulus)
    diag = numpy.array([1.0, 1.0, 1.0, 0.0], dtype=numpy.float64)

    # Values from previous time step
    meanPlasStrainT = (plasStrainT[0] + plasStrainT[1] + plasStrainT[2])/3.0
    devPlasStrainT = plasStrainT - meanPlasStrainT * diag
    
    # Initial stress values
    initialStress4 = numpy.array([initialStress[0], initialStress[1],
                                 stateVars[0], initialStress[2]],
                                dtype=numpy.float64)
    meanStressInitial = (initialStress[0] + initialStress[1] + stateVars[0])/3.0
    devStressInitial = initialStress4 - meanStressInitial * diag
    
    # Initial strain values
    initialStrain4 = numpy.array([initialStrain[0], initialStrain[1],
                                 0.0, initialStrain[2]], dtype=numpy.float64)
    meanStrainInitial = (initialStrain[0] + initialStrain[1])/3.0
    devStrainInitial = initialStrain4 - meanStrainInitial * diag

    # Values for current time step
    strainTpdt4 = numpy.array([strainTpdt[0], strainTpdt[1],
                                 0.0, strainTpdt[2]], dtype=numpy.float64)
    meanStrainTpdt = (strainTpdt[0] + strainTpdt[1])/3.0
    meanStrainPPTpdt = meanStrainTpdt - meanPlasStrainT - meanStrainInitial
    strainPPTpdt = strainTpdt4 - meanStrainTpdt * diag - \
                   devPlasStrainT - devStrainInitial

    # Compute trial elastic stresses and yield function to see if yield should
    # occur.
    trialDevStress = strainPPTpdt/ae + devStressInitial
    trialMeanStress = meanStrainPPTpdt/am + meanStressInitial
    yieldFunction = 3.0 * alphaYieldV * trialMeanStress + \
                    math.sqrt(0.5 * self._scalarProduct(trialDevStress,
                                                        trialDevStress)) - \
                                                        betaV

    # If yield function is greater than zero, compute elastoplastic stress.
    if (yieldFunction >= 0.0):
      devStressInitialProd = self._scalarProduct(devStressInitial,
                                                 devStressInitial)
      strainPPTpdtProd = self._scalarProduct(strainPPTpdt, strainPPTpdt)
      d = math.sqrt(ae * ae * devStressInitialProd + \
                    2.0 * ae * self._scalarProduct(devStressInitial, \
                                                   strainPPTpdt) + \
                    strainPPTpdtProd)
      dFac = math.sqrt(2.0) * d
      testMult = 2.0 * ae * am * \
                 (3.0 * alphaYieldV * (meanStrainPPTpdt/am + meanStressInitial) + \
                  d/(math.sqrt(2.0) * ae) - betaV)/ \
                  (6.0 * alphaYieldV * alphaFlowV * ae + am)
      plasticMult = min(testMult, dFac)
      deltaMeanPlasticStrain = plasticMult * alphaFlowV
      meanStressTpdt = (meanStrainPPTpdt - deltaMeanPlasticStrain)/am + \
                       meanStressInitial

      deltaDevPlasticStrain = 0.0
      stressTpdt4 = numpy.zeros( (4), dtype=numpy.float64)
      stressTpdt = numpy.zeros( (tensorSize), dtype=numpy.float64)
      plasStrainTpdt = numpy.zeros( (4), dtype=numpy.float64)

      for iComp in range(4):
        deltaDevPlasticStrain = plasticMult * (strainPPTpdt[iComp] + \
                                               ae * devStressInitial[iComp])/ \
                                               dFac
        devStressTpdt = (strainPPTpdt[iComp] - deltaDevPlasticStrain)/ae + \
                        devStressInitial[iComp]
        stressTpdt4[iComp] = devStressTpdt + diag[iComp] * meanStressTpdt
        plasStrainTpdt[iComp] = plasStrainT[iComp] + deltaDevPlasticStrain + \
                                diag[iComp] * deltaMeanPlasticStrain
      
      stressTpdt[0] = stressTpdt4[0]
      stressTpdt[1] = stressTpdt4[1]
      stressTpdt[2] = stressTpdt4[3]

    return stressTpdt, plasStrainTpdt


  def _calcStress(self, strainV, muV, lambdaV, alphaYieldV, betaV, alphaFlowV,
                  plasStrainV, initialStressV, initialStrainV, stateVars):
    """
    Compute stress, updated state variables and derivative of elasticity matrix.
    """
    import scipy.misc

    # Define some numpy arrays
    strainTpdt = numpy.array(strainV, dtype=numpy.float64)
    plasStrainT = numpy.array(plasStrainV, dtype=numpy.float64)
    initialStress = numpy.array(initialStressV, dtype=numpy.float64)
    initialStrain = numpy.array(initialStrainV, dtype=numpy.float64)
    stressTpdt, plasStrainTpdt = self._computeStress(strainTpdt, muV, lambdaV,
                                                     alphaYieldV, betaV,
                                                     alphaFlowV, plasStrainT,
                                                     initialStress,
                                                     initialStrain,
                                                     stateVars)
    stateVarsUpdated = numpy.array( (stateVars[0], plasStrainTpdt[0],
                                     plasStrainTpdt[1], plasStrainTpdt[2],
                                     plasStrainTpdt[3]))
    # Compute components of tangent constitutive matrix using numerical
    # derivatives.
    derivDx = 1.0e-12
    derivOrder = 3
    elasticConstsList = []

    for stressComp in range(tensorSize):
      for strainComp in range(tensorSize):
        dStressDStrain = scipy.misc.derivative(self._calcStressComponent,
                                               strainTpdt[strainComp],
                                               dx=derivDx,
                                               args=(strainComp,
                                                     stressComp,
                                                     strainTpdt, muV, lambdaV,
                                                     alphaYieldV, betaV,
                                                     alphaFlowV,
                                                     plasStrainT,
                                                     initialStress,
                                                     initialStrain,
                                                     stateVars),
                                               order=derivOrder)
        elasticConstsList.append(dStressDStrain)

    elasticConsts = numpy.array(elasticConstsList, dtype=numpy.float64)

    return (elasticConsts, numpy.ravel(stressTpdt),
            numpy.ravel(stateVarsUpdated))
 

# MAIN /////////////////////////////////////////////////////////////////
if __name__ == "__main__":

  app = DruckerPragerPlaneStrainTimeDep()
  app.run()


# End of file 
