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

## @file unittests/libtests/materials/data/GenMaxwellPlaneStrainTimeDep.py

## @brief Python application for generating C++ data files for testing
## C++ GenMaxwellPlaneStrain object with viscoelastic behavior.

from ElasticMaterialApp import ElasticMaterialApp

import numpy
import math

# ----------------------------------------------------------------------
dimension = 2
numElasticConsts = 9
tensorSize = 3
numMaxwellModels = 3

# GenMaxwellPlaneStrainTimeDep class
class GenMaxwellPlaneStrainTimeDep(ElasticMaterialApp):
  """
  Python application for generating C++ data files for testing C++
  GenMaxwellPlaneStrain object using viscoelastic behavior.
  """

  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="genmaxwellplanestraintimedep"):
    """
    Constructor.
    """
    ElasticMaterialApp.__init__(self, name)

    # import pdb
    # pdb.set_trace()

    numLocs = 2
    self.dt = 2.0e5
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
    initialStrainA = [3.6e-5, 3.5e-5, 3.3e-5]
    muA = vsA*vsA*densityA
    lambdaA = vpA*vpA*densityA - 2.0*muA
    meanStrainA = (strainA[0] + strainA[1])/3.0
    stressInitialZZA = numpy.array([2.0e4], dtype=numpy.float64)

    densityB = 2000.0
    vsB = 1200.0
    vpB = vsB*3**0.5
    shearRatioB = [0.2, 0.2, 0.2]
    viscosityB = [1.0e18, 1.0e19, 1.0e20]
    strainB = [4.1e-4, 4.2e-4, 4.4e-4]
    initialStressB = [5.1e4, 5.2e4, 5.4e4]
    initialStrainB = [6.1e-5, 6.2e-5, 6.6e-5]
    muB = vsB*vsB*densityB
    lambdaB = vpB*vpB*densityB - 2.0*muB
    meanStrainB = (strainB[0] + strainB[1])/3.0
    stressInitialZZB = numpy.array([5.0e4], dtype=numpy.float64)

    # Compute Maxwell time and viscous strain.
    diag = numpy.array([1.0, 1.0, 1.0, 0.0], dtype=numpy.float64)
    strainTA = [strainA[0], strainA[1], 0.0, strainA[2]]
    strainTB = [strainB[0], strainB[1], 0.0, strainB[2]]
    
    maxwellTimeA = [1.0e30, 1.0e30, 1.0e30]
    maxwellTimeB = [1.0e30, 1.0e30, 1.0e30]
    visStrainA = numpy.zeros( (numMaxwellModels, 4), dtype=numpy.float64)
    visStrainB = numpy.zeros( (numMaxwellModels, 4), dtype=numpy.float64)
    for imodel in xrange(numMaxwellModels):
      if shearRatioA[imodel] != 0.0:
        maxwellTimeA[imodel] = viscosityA[imodel]/(muA*shearRatioA[imodel])
        visStrainA[imodel,:] = numpy.array(strainTA, dtype=numpy.float64) - \
                               diag * meanStrainA
      if shearRatioB[imodel] != 0.0:
        maxwellTimeB[imodel] = viscosityB[imodel]/(muB*shearRatioB[imodel])
        visStrainB[imodel,:] = numpy.array(strainTB, dtype=numpy.float64) - \
                               diag * meanStrainB

    dbPropA = [densityA, vsA, vpA] + shearRatioA + viscosityA
    dbPropB = [densityB, vsB, vpB] + shearRatioB + viscosityB
    self.dbProperties = numpy.array([dbPropA, dbPropB], dtype=numpy.float64)
    propA = [densityA, muA, lambdaA] + shearRatioA + maxwellTimeA
    propB = [densityB, muB, lambdaB] + shearRatioB + maxwellTimeB
    self.properties = numpy.array([propA, propB], dtype=numpy.float64)

    # TEMPORARY, need to determine how to use initial state variables
    # At present, only the first (stressInitialZZ) is being used.
    self.dbStateVars = numpy.zeros((numLocs,
                                    1 + tensorSize + 4 * numMaxwellModels),
                                   dtype=numpy.float64)
    self.dbStateVars[0, 0] = stressInitialZZA
    self.dbStateVars[1, 0] = stressInitialZZB


    self.lengthScale = 1.0e+3
    self.pressureScale = muA
    self.timeScale = 1.0
    self.densityScale = muA / (self.lengthScale / self.timeScale)**2

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

    self.initialStress = numpy.array([initialStressA,
                                      initialStressB],
                                    dtype=numpy.float64)
    self.initialStrain = numpy.array([initialStrainA,
                                      initialStrainB],
                                    dtype=numpy.float64)
    
    self.density = numpy.array([densityA,
                                densityB],
                               dtype=numpy.float64)

    # Simplest approach for now is to assume this is the first step
    # after the elastic solution.  In that case, both the total strain
    # from the last step (total_strain) and the total viscous strain
    # (viscous_strain) are defined by the assigned elastic strain.
    # Revised approach.  For a better test, I am setting the total strain
    # for the current time step to be equal to the strain from the previous
    # time step plus a constant amount.
    totalStrainA = [strainA[0] + 1.0e-5,
                    strainA[1] + 1.0e-5,
                    strainA[2] + 1.0e-5]
    totalStrainB = [strainB[0] + 1.0e-5,
                    strainB[1] + 1.0e-5,
                    strainB[2] + 1.0e-5]
    
    visStrainVecA = numpy.ravel(visStrainA)
    visStrainVecB = numpy.ravel(visStrainB)
    strainVecA = numpy.array(strainA, dtype=numpy.float64)
    strainVecB = numpy.array(strainB, dtype=numpy.float64)
    stateVarsA = numpy.concatenate((stressInitialZZA, strainVecA,
                                    visStrainVecA))
    stateVarsB = numpy.concatenate((stressInitialZZB, strainVecB,
                                    visStrainVecB))
                                   
    self.stateVars = numpy.array([stateVarsA, stateVarsB], dtype=numpy.float64)

    stressInitialZZANondim = stressInitialZZA/mu0
    stressInitialZZBNondim = stressInitialZZB/mu0
    self.stateVarsNondim = numpy.zeros( (numLocs,
                                         1 + tensorSize + 4 * numMaxwellModels),
                                        dtype=numpy.float64)
    self.stateVarsNondim[:] = self.stateVars[:] # no scaling
    # Scale stressInitialZZ
    self.stateVarsNondim[0, 0] = stressInitialZZANondim
    self.stateVarsNondim[1, 0] = stressInitialZZBNondim

    self.strain = numpy.array([totalStrainA, totalStrainB], dtype=numpy.float64)
    self.stress = numpy.zeros( (numLocs, tensorSize), dtype=numpy.float64)
    self.stateVarsUpdated = numpy.zeros((numLocs,
                                         1 + tensorSize + 4 * numMaxwellModels),
                                        dtype=numpy.float64)
    self.elasticConsts = numpy.zeros((numLocs, numElasticConsts), \
                                     dtype=numpy.float64)

    (self.elasticConsts[0,:], self.stress[0,:], self.stateVarsUpdated[0,:]) = \
                              self._calcStress(strainA, muA, lambdaA,
                                               shearRatioA, maxwellTimeA,
                                               totalStrainA, visStrainA,
                                               initialStressA, initialStrainA,
                                               stateVarsA)
    (self.elasticConsts[1,:], self.stress[1,:], self.stateVarsUpdated[1,:]) = \
                              self._calcStress(strainB, muB, lambdaB,
                                               shearRatioB, maxwellTimeB,
                                               totalStrainB, visStrainB,
                                               initialStressB, initialStrainB,
                                               stateVarsB)
    self.dtStableImplicit = 0.2*min(min(maxwellTimeA), min(maxwellTimeB))
    self.dtStableExplicit = 1000.0 / vpA

    return


  def _calcStressComponent(self, strainVal, strainComp, stressComp, strainTpdt,
                           muV, lambdaV, shearRatioV, maxwellTimeV, dqV,
                           totalStrainT, viscousStrainT,
                           initialStress, initialStrain,
                           stateVars):
    """
    Function to compute a particular stress component as a function of a
    strain component.
    """
    strainTest = numpy.array(strainTpdt, dtype=numpy.float64)
    strainTest[strainComp] = strainVal
    stressTpdt, viscousStrainTpdt = self._computeStress(strainTest,
                                                        muV,
                                                        lambdaV,
                                                        shearRatioV,
                                                        maxwellTimeV,
                                                        dqV,
                                                        totalStrainT,
                                                        viscousStrainT,
                                                        initialStress,
                                                        initialStrain,
                                                        stateVars)
    return stressTpdt[stressComp]


  def _computeStress(self, strainTpdt, muV, lambdaV, shearRatioV, maxwellTimeV,
                     dqV, strainT, viscousStrainT,
                     initialStress, initialStrain, stateVars):
    """
    Function to compute stresses and viscous strains for a given strain.
    """
    
    bulkModulus = lambdaV + 2.0 * muV / 3.0
    diag = numpy.array([1.0, 1.0, 1.0, 0.0], dtype=numpy.float64)

    # Initial stresses and strains
    meanStrainInitial = (initialStrain[0] + initialStrain[1]) / 3.0
    meanStressInitial = \
    (initialStress[0] + initialStress[1] + stateVars[0]) / 3.0

    initialStressV = numpy.array([initialStress[0],
                                  initialStress[1],
                                  stateVars[0],
                                  initialStress[2]], dtype=numpy.float64)
    initialStrainV = numpy.array([initialStrain[0],
                                  initialStrain[1],
                                  0.0,
                                  initialStrain[2]], dtype=numpy.float64)

    devStrainInitial = initialStrainV - diag * meanStrainInitial
    devStressInitial = initialStressV - diag * meanStressInitial

    # Strains from previous time step
    strainTV = numpy.array([strainT[0], strainT[1], 0.0, strainT[2]],
                           dtype=numpy.float64)
    meanStrainT = (strainT[0] + strainT[1]) / 3.0
    devStrainT = strainTV - diag * meanStrainT - devStrainInitial

    # Strains and mean stress for current time step
    strainTpdtV = numpy.array([strainTpdt[0], strainTpdt[1], 0.0,
                               strainTpdt[2]], dtype=numpy.float64)
    meanStrainTpdt = (strainTpdt[0] + strainTpdt[1]) / 3.0
    devStrainTpdt = strainTpdtV - diag * meanStrainTpdt - devStrainInitial
    meanStressTpdt = 3.0 * bulkModulus * \
                     (meanStrainTpdt - meanStrainInitial) + meanStressInitial

    # Compute viscous factors
    visFrac = 0.0
    visFac = 0.0
    expFac = numpy.zeros((numMaxwellModels), dtype=numpy.float64)
    for model in range(numMaxwellModels):
      visFrac += shearRatioV[model]
      visFac += shearRatioV[model] * dqV[model]
      expFac[model] = math.exp(-self.dt/maxwellTimeV[model])
                      

    elasFrac = 1.0 - visFrac
    stressTpdt = numpy.zeros( (4), dtype=numpy.float64)
    viscousStrainTpdt = numpy.zeros( (3, 4), dtype=numpy.float64)
    elasFac = 2.0 * muV

    # Compute viscous strain and stress
    for iComp in range(4):
      deltaStrain = devStrainTpdt[iComp] - devStrainT[iComp]
      devStressTpdt = elasFrac * devStrainTpdt[iComp]
      for model in range(numMaxwellModels):
        if shearRatioV[model] != 0.0:
          viscousStrainTpdt[model, iComp] = expFac[model] * \
                                            viscousStrainT[model, iComp] + \
                                            dqV[model] * deltaStrain
          devStressTpdt += shearRatioV[model] * viscousStrainTpdt[model, iComp]
      devStressTpdt = elasFac * devStressTpdt + devStressInitial[iComp]
      stressTpdt[iComp] = diag[iComp] * meanStressTpdt + devStressTpdt
    
    stressTpdtV = numpy.array([stressTpdt[0], stressTpdt[1], stressTpdt[3]],
                              dtype=numpy.float64)
      
    return stressTpdtV, viscousStrainTpdt

  
  def _computeViscousFactor(self, maxwellTime):
    """
    Compute viscous strain factor for a given Maxwell time.
    """
    
    timeFrac = 1.0e-5
    numTerms = 5
    dq = 0.0
    if maxwellTime < timeFrac*self.dt:
      fSign = 1.0
      factorial = 1.0
      fraction = 1.0
      dq = 1.0
      for iTerm in range(2, numTerms + 1):
        factorial *= iTerm
        fSign *= -1.0
        fraction *= self.dt/maxwellTime
        dq += fSign*fraction/factorial
    else:
      dq = maxwellTime*(1.0-math.exp(-self.dt/maxwellTime))/self.dt

    return dq

  
  def _calcStress(self, strainV, muV, lambdaV, shearRatioV, maxwellTimeV,
                  totalStrainV, viscousStrainV, initialStressV, initialStrainV,
                  stateVars):
    """
    Compute stress and derivative of elasticity matrix.
    This assumes behavior is always viscoelastic.
    """
    import scipy.misc
    
    # Define some numpy arrays
    strainTpdt = numpy.array(totalStrainV, dtype=numpy.float64)
    strainT = numpy.array(strainV, dtype=numpy.float64)
    viscousStrainT = numpy.array(viscousStrainV, dtype=numpy.float64)
    initialStress = numpy.array(initialStressV, dtype=numpy.float64)
    initialStrain = numpy.array(initialStrainV, dtype=numpy.float64)
    stressZZInitial = numpy.array([stateVars[0]], dtype=numpy.float64)

    # Compute viscous factors for each Maxwell model
    dqV = numpy.zeros(numMaxwellModels, dtype=numpy.float64)
    for model in range(numMaxwellModels):
      dqV[model] = self._computeViscousFactor(maxwellTimeV[model])

    # Compute current stress and viscous strain
    stressTpdt, viscousStrainTpdt = self._computeStress(strainTpdt,
                                                        muV,
                                                        lambdaV,
                                                        shearRatioV,
                                                        maxwellTimeV,
                                                        dqV,
                                                        strainT,
                                                        viscousStrainT,
                                                        initialStress,
                                                        initialStrain,
                                                        stateVars)
    # Form updated state variables
    strainTpdtVec = numpy.ravel(strainTpdt)
    viscousStrainTpdtVec = numpy.ravel(viscousStrainTpdt)
    stateVarsUpdated = numpy.concatenate((stressZZInitial, strainTpdtVec,
                                          viscousStrainTpdtVec))

    # Compute components of tangent constitutive matrix using numerical
    # derivatives.
    derivDx = 1.0e-12
    derivOrder = 3
    elasticConstsList = []

    for stressComp in range(tensorSize):
      for strainComp in range(tensorSize):
        dStressDStrain = stressComp + strainComp
        dStressDStrain = scipy.misc.derivative(self._calcStressComponent,
                                               strainTpdt[strainComp],
                                               dx=derivDx,
                                               args=(strainComp,
                                                     stressComp,
                                                     strainTpdt,
                                                     muV,
                                                     lambdaV,
                                                     shearRatioV,
                                                     maxwellTimeV,
                                                     dqV,
                                                     strainT,
                                                     viscousStrainT,
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

  app = GenMaxwellPlaneStrainTimeDep()
  app.run()


# End of file 
