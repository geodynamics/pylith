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
# Copyright (c) 2010-2014 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file unittests/libtests/materials/data/DruckerPrager3DTimeDep.py

## @brief Python application for generating C++ data files for testing
## C++ DruckerPrager3D object with time dependent behavior.

from ElasticMaterialApp import ElasticMaterialApp

import numpy
import math

# ----------------------------------------------------------------------
dimension = 3
numElasticConsts = 36
tensorSize = 6

# DruckerPrager3DTimeDep class
class DruckerPrager3DTimeDep(ElasticMaterialApp):
  """
  Python application for generating C++ data files for testing C++
  DruckerPrager3D object with time dependent behavior.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="druckerprager3dtimedep"):
    """
    Constructor.
    """
    ElasticMaterialApp.__init__(self, name)

    numLocs = 2

    self.dimension = dimension
    self.numLocs = numLocs

    self.dbPropertyValues = ["density", "vs", "vp",
                             "friction-angle", "cohesion",
                             "dilatation-angle"]
    self.propertyValues = ["density", "mu", "lambda",
                           "alpha_yield", "beta", "alpha_flow"]
    self.numPropertyValues = numpy.array([1, 1, 1, 1, 1, 1], dtype=numpy.int32)

    self.dbStateVarValues = ["plastic-strain-xx",
                             "plastic-strain-yy",
                             "plastic-strain-zz",
                             "plastic-strain-xy",
                             "plastic-strain-yz",
                             "plastic-strain-xz"
                             ]
    self.stateVarValues = ["plastic-strain"]
    self.numStateVarValues = numpy.array([6], dtype=numpy.int32)

    self.dt = 2.0e5

    densityA = 2500.0
    vsA = 3000.0
    vpA = vsA*3**0.5
    # First case has same values for friction angle and dilatation angle.
    frictionAngleA = math.radians(30.0)
    dilatationAngleA = math.radians(20.0)
    cohesionA = 3.0e5
    strainA = [-2.1e-4, 1.2e-4, 1.3e-4, 1.1e-5, 1.1e-5, 1.1e-5]
    initialStressA = [2.1e4, 2.2e4, 2.3e4, 2.4e4, 2.5e4, 2.6e4]
    initialStrainA = [3.6e-5, 3.5e-5, 3.4e-5, 3.3e-5, 3.2e-5, 3.1e-5]
    # initialStressA = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    # initialStrainA = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    muA = vsA*vsA*densityA
    lambdaA = vpA*vpA*densityA - 2.0*muA

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
    # frictionAngleB = 0.0
    # dilatationAngleB = 0.0
    cohesionB = 1.0e4
    strainB = [4.1e-4, 4.2e-4, 4.3e-4, 1.4e-4, 1.5e-4, 1.6e-4]
    initialStressB = [5.6e4, 5.5e4, 5.4e4, 5.3e4, 5.2e4, 5.1e4]
    initialStrainB = [6.6e-5, 6.5e-5, 6.4e-5, 6.3e-5, 6.2e-5, 6.1e-5]
    # initialStressB = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    # initialStrainB = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    muB = vsB*vsB*densityB
    lambdaB = vpB*vpB*densityB - 2.0*muB
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
    self.dbStateVars = numpy.zeros( (numLocs, tensorSize), dtype=numpy.float64)

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

    self.initialStress = numpy.array([initialStressA,
                                      initialStressB],
                                    dtype=numpy.float64)
    self.initialStrain = numpy.array([initialStrainA,
                                      initialStrainB],
                                    dtype=numpy.float64)
    
    self.density = numpy.array([densityA,
                                densityB],
                               dtype=numpy.float64)

    # Define state variables
    plasStrainA = [4.1e-5, 4.2e-5, 4.3e-5, 4.4e-5, 4.5e-5, 4.6e-5]
    plasStrainB = [1.1e-5, 1.2e-5, 1.3e-5, 1.4e-5, 1.5e-5, 1.6e-5]
    # plasStrainA = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    # plasStrainB = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    self.stateVars = numpy.array([[plasStrainA], [plasStrainB] ],
                                 dtype=numpy.float64)
    self.stateVarsNondim = self.stateVars # No scaling.
    self.strain = numpy.array([strainA,
                               strainB],
                               dtype=numpy.float64)
    
    self.stress = numpy.zeros( (numLocs, tensorSize), dtype=numpy.float64)
    self.stateVarsUpdated = numpy.zeros( (numLocs, tensorSize),
                                         dtype=numpy.float64)
    self.elasticConsts = numpy.zeros( (self.numLocs, numElasticConsts), \
                                        dtype=numpy.float64)

    (self.elasticConsts[0,:], self.stress[0,:], self.stateVarsUpdated[0,:]) = \
                              self._calcStress(strainA, muA, lambdaA,
                                               alphaYieldA, betaA, alphaFlowA,
                                               plasStrainA,
                                               initialStressA, initialStrainA)
    (self.elasticConsts[1,:], self.stress[1,:], self.stateVarsUpdated[1,:]) = \
                              self._calcStress(strainB, muB, lambdaB,
                                               alphaYieldB, betaB, alphaFlowB,
                                               plasStrainB,
                                               initialStressB, initialStrainB)

    self.dtStableImplicit = 1.0e+99
    self.dtStableExplicit = 1000.0 / vpA

    return


  def _scalarProduct(self, tensor1, tensor2):
    """
    Compute the scalar product of two tensors stored in vector form.
    """
    scalarProduct = tensor1[0] * tensor2[0] + \
                    tensor1[1] * tensor2[1] + \
                    tensor1[2] * tensor2[2] + \
                    2.0 * (tensor1[3] * tensor2[3] + \
                           tensor1[4] * tensor2[4] + \
                           tensor1[5] * tensor2[5])
    return scalarProduct


  def _calcStressComponent(self, strainVal, strainComp, stressComp, strainTpdt,
                           muV, lambdaV, alphaYieldV, betaV, alphaFlowV,
                           plasStrainT, initialStress, initialStrain):
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
                                                    initialStrain)
    return stressTpdt[stressComp]


  def _computeStress(self, strainTpdt, muV, lambdaV, alphaYieldV, betaV,
                     alphaFlowV, plasStrainT, initialStress, initialStrain):
    """
    Function to compute stresses and plastic strains.
    """
    
    # Constants
    mu2 = 2.0 * muV
    lamPlusMu = lambdaV + muV
    bulkModulus = lambdaV + mu2/3.0
    ae = 1.0/mu2
    am = 1.0/(3.0 * bulkModulus)
    diag = numpy.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0], dtype=numpy.float64)

    # Values from previous time step
    meanPlasStrainT = (plasStrainT[0] + plasStrainT[1] + plasStrainT[2])/3.0
    devPlasStrainT = plasStrainT - meanPlasStrainT * diag
    
    # Initial stress values
    meanStressInitial = (initialStress[0] + initialStress[1] +
                         initialStress[2])/3.0
    devStressInitial = initialStress - meanStressInitial * diag
    
    # Initial strain values
    meanStrainInitial = (initialStrain[0] + initialStrain[1] +
                         initialStrain[2])/3.0
    devStrainInitial = initialStrain - meanStrainInitial * diag

    # Values for current time step
    meanStrainTpdt = (strainTpdt[0] + strainTpdt[1] + strainTpdt[2])/3.0
    meanStrainPPTpdt = meanStrainTpdt - meanPlasStrainT - meanStrainInitial
    strainPPTpdt = strainTpdt - meanStrainTpdt * diag - \
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
      stressTpdt = numpy.zeros( (tensorSize), dtype=numpy.float64)
      plasStrainTpdt = numpy.zeros( (tensorSize), dtype=numpy.float64)

      for iComp in range(tensorSize):
        deltaDevPlasticStrain = plasticMult * (strainPPTpdt[iComp] + \
                                               ae * devStressInitial[iComp])/ \
                                               dFac
        devStressTpdt = (strainPPTpdt[iComp] - deltaDevPlasticStrain)/ae + \
                        devStressInitial[iComp]
        stressTpdt[iComp] = devStressTpdt + diag[iComp] * meanStressTpdt
        plasStrainTpdt[iComp] = plasStrainT[iComp] + deltaDevPlasticStrain + \
                                diag[iComp] * deltaMeanPlasticStrain
      
    return stressTpdt, plasStrainTpdt


  def _calcStress(self, strainV, muV, lambdaV, alphaYieldV, betaV, alphaFlowV,
                  plasStrainV, initialStressV, initialStrainV):
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
                                                     initialStrain)
    stateVarsUpdated = numpy.array( [plasStrainTpdt], dtype=numpy.float64)
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
                                                     initialStrain),
                                               order=derivOrder)
        elasticConstsList.append(dStressDStrain)

    elasticConsts = numpy.array(elasticConstsList, dtype=numpy.float64)

    return (elasticConsts, numpy.ravel(stressTpdt),
            numpy.ravel(stateVarsUpdated))
 

# MAIN /////////////////////////////////////////////////////////////////
if __name__ == "__main__":

  app = DruckerPrager3DTimeDep()
  app.run()


# End of file 
