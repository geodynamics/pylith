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
# Copyright (c) 2010 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file unittests/libtests/materials/data/MaxwellPlaneStrainTimeDep.py

## @brief Python application for generating C++ data files for testing
## C++ MaxwellPlaneStrain object with viscoelastic behavior.

from ElasticMaterialApp import ElasticMaterialApp

import numpy

# ----------------------------------------------------------------------
dimension = 2
numElasticConsts = 9
tensorSize = 3

# MaxwellPlaneStrainTimeDep class
class MaxwellPlaneStrainTimeDep(ElasticMaterialApp):
  """
  Python application for generating C++ data files for testing C++
  MaxwellPlaneStrain object using viscoelastic behavior.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="maxwellplanestraintimedep"):
    """
    Constructor.
    """
    ElasticMaterialApp.__init__(self, name)

    import pdb
    pdb.set_trace()

    numLocs = 2

    self.dimension = dimension
    self.numLocs = numLocs

    self.dbPropertyValues = ["density", "vs", "vp", "viscosity"]
    self.propertyValues = ["density", "mu", "lambda", "maxwellTime"]
    self.numPropertyValues = numpy.array([1, 1, 1, 1], dtype=numpy.int32)

    self.dbStateVarValues = ["stress-zz-initial",
                             "total-strain-xx",
                             "total-strain-yy",
                             "total-strain-xy",
                             "viscous-strain-xx",
                             "viscous-strain-yy",
                             "viscous-strain-zz",
                             "viscous-strain-xy"
                             ]
    self.stateVarValues = ["stress-zz-initial","total-strain", "viscous-strain"]
    self.numStateVarValues = numpy.array([1, 3, 4], dtype=numpy.int32)

    self.dt = 2.0e5

    densityA = 2500.0
    vsA = 3000.0
    vpA = vsA*3**0.5
    viscosityA = 1.0e18
    strainA = [1.1e-4, 1.2e-4, 1.4e-4]
    initialStressA = [2.1e4, 2.2e4, 2.4e4]
    initialStrainA = [3.6e-5, 3.5e-5, 3.3e-5]
    muA = vsA*vsA*densityA
    lambdaA = vpA*vpA*densityA - 2.0*muA
    maxwellTimeA = viscosityA / muA
    meanStrainA = (strainA[0] + strainA[1])/3.0
    stressInitialZZA = numpy.array([2.0e4], dtype=numpy.float64)

    densityB = 2000.0
    vsB = 1200.0
    vpB = vsB*3**0.5
    viscosityB = 1.0e19
    strainB = [4.1e-4, 4.2e-4, 4.4e-4]
    initialStressB = [5.1e4, 5.2e4, 5.4e4]
    initialStrainB = [6.1e-5, 6.2e-5, 6.6e-5]
    muB = vsB*vsB*densityB
    lambdaB = vpB*vpB*densityB - 2.0*muB
    maxwellTimeB = viscosityB / muB
    meanStrainB = (strainB[0] + strainB[1])/3.0
    stressInitialZZB = numpy.array([5.0e4], dtype=numpy.float64)

    diag = numpy.array([1.0, 1.0, 1.0, 0.0],
                       dtype=numpy.float64)

    strainTA = [strainA[0], strainA[1], 0.0, strainA[2]]
    strainTB = [strainB[0], strainB[1], 0.0, strainB[2]]

    self.lengthScale = 1.0e+3
    self.pressureScale = muA
    self.timeScale = 1.0
    self.densityScale = 1.0e+3

    self.dbProperties = numpy.array([ [densityA, vsA, vpA, viscosityA],
                                      [densityB, vsB, vpB, viscosityB] ], 
                                    dtype=numpy.float64)
    self.properties = numpy.array([ [densityA, muA, lambdaA, maxwellTimeA],
                                    [densityB, muB, lambdaB, maxwellTimeB] ],
                                     dtype=numpy.float64)

    # TEMPORARY, need to determine how to use initial state variables
    # At present, only the first (stressInitialZZ) is being used.
    self.dbStateVars = numpy.zeros( (numLocs, tensorSize+5),
                                    dtype=numpy.float64)
    self.dbStateVars[0, 0] = stressInitialZZA
    self.dbStateVars[1, 0] = stressInitialZZB

    mu0 = self.pressureScale
    density0 = self.densityScale
    time0 = self.timeScale
    self.propertiesNondim = \
        numpy.array([ [densityA/density0, muA/mu0, lambdaA/mu0, \
                       maxwellTimeA/time0],
                      [densityB/density0, muB/mu0, lambdaB/mu0, \
                       maxwellTimeB/time0] ],
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
    viscousStrainA = numpy.array(strainTA) - diag * meanStrainA
    viscousStrainB = numpy.array(strainTB) - diag * meanStrainB
    viscousStrainVecA = numpy.ravel(viscousStrainA)
    viscousStrainVecB = numpy.ravel(viscousStrainB)
    strainVecA = numpy.array(strainA, dtype=numpy.float64)
    strainVecB = numpy.array(strainB, dtype=numpy.float64)
    stateVarsA = numpy.concatenate((stressInitialZZA, strainVecA,
                                    viscousStrainVecA))
    stateVarsB = numpy.concatenate((stressInitialZZB, strainVecB,
                                    viscousStrainVecB))
    self.stateVars = numpy.array([ stateVarsA,
                                   stateVarsB ],
                                 dtype=numpy.float64)
    stressInitialZZANondim = stressInitialZZA/mu0
    stressInitialZZBNondim = stressInitialZZB/mu0
    self.stateVarsNondim = numpy.zeros( (numLocs, tensorSize + 5),
                                        dtype=numpy.float64)
    self.stateVarsNondim[:] = self.stateVars[:] # no scaling
    # Scale stressInitialZZ
    self.stateVarsNondim[0, 0] = stressInitialZZANondim
    self.stateVarsNondim[1, 0] = stressInitialZZBNondim
    
    self.strain = numpy.array([totalStrainA, totalStrainB],
                               dtype=numpy.float64)
    self.stress = numpy.zeros( (numLocs, tensorSize), dtype=numpy.float64)
    self.stateVarsUpdated = numpy.zeros( (numLocs, tensorSize + 5),
                                         dtype=numpy.float64)
    self.elasticConsts = numpy.zeros( (numLocs, numElasticConsts),
                                      dtype=numpy.float64)

    (self.elasticConsts[0,:], self.stress[0,:], self.stateVarsUpdated[0,:]) = \
                     self._calcStress(strainA, 
                                      muA, lambdaA, maxwellTimeA,
                                      totalStrainA, viscousStrainA,
                                      initialStressA, initialStrainA,
                                      stateVarsA)
    (self.elasticConsts[1,:], self.stress[1,:], self.stateVarsUpdated[1,:]) = \
                              self._calcStress(strainB, 
                                               muB, lambdaB, maxwellTimeB, 
                                               totalStrainB, viscousStrainB,
                                               initialStressB, initialStrainB,
                                               stateVarsB)

    self.stateVarsUpdated[0, 0] = stressInitialZZA
    self.stateVarsUpdated[1, 0] = stressInitialZZB

    self.dtStableImplicit = 0.2*min(maxwellTimeA, maxwellTimeB)
    return


  def _calcStressComponent(self, strainVal, strainComp, stressComp, strainTpdt,
                           muV, lambdaV, maxwellTimeV, dqV, totalStrainT,
                           viscousStrainT, initialStress, initialStrain,
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
                                                        maxwellTimeV,
                                                        dqV,
                                                        totalStrainT,
                                                        viscousStrainT,
                                                        initialStress,
                                                        initialStrain,
                                                        stateVars)
    return stressTpdt[stressComp]


  def _computeStress(self, strainTpdt, muV, lambdaV, maxwellTimeV, dqV,
                     strainT, viscousStrainT,
                     initialStress, initialStrain, stateVars):
    """
    Function to compute stresses and viscous strains for a given strain.
    """
    import math
    
    bulkModulus = lambdaV + 2.0 * muV / 3.0
    diag = [1.0, 1.0, 0.0]

    # Initial stresses and strains
    meanStrainInitial = \
    (initialStrain[0] + initialStrain[1]) / 3.0
    meanStressInitial = \
    (initialStress[0] + initialStress[1] + stateVars[0]) / 3.0

    devStrainInitial = initialStrain - numpy.array(diag) * meanStrainInitial
    devStressInitial = initialStress - numpy.array(diag) * meanStressInitial

    meanStrainT = (strainT[0] + strainT[1]) / 3.0
    meanStrainTpdt = (strainTpdt[0] + strainTpdt[1]) / 3.0
    meanStressTpdt = 3.0 * bulkModulus * \
                     (meanStrainTpdt - meanStrainInitial) + meanStressInitial

    stressTpdt = numpy.zeros( (tensorSize), dtype=numpy.float64)
    viscousStrainTpdt = numpy.zeros( (4), dtype=numpy.float64)

    expFac = math.exp(-self.dt/maxwellTimeV)
    elasFac = 2.0 * muV

    devStrainTpdt11 = strainTpdt[0] - meanStrainTpdt
    devStrainTpdt22 = strainTpdt[1] - meanStrainTpdt
    devStrainTpdt33 = -meanStrainTpdt
    devStrainTpdt12 = strainTpdt[2]

    devStrainT11 = strainT[0] - meanStrainT
    devStrainT22 = strainT[1] - meanStrainT
    devStrainT33 = -meanStrainT
    devStrainT12 = strainT[2]

    viscousStrainTpdt = [expFac * viscousStrainT[0] + \
                        dqV * (devStrainTpdt11 - devStrainT11),
                        expFac * viscousStrainT[1] + \
                        dqV * (devStrainTpdt22 - devStrainT22),
                        expFac * viscousStrainT[2] + \
                        dqV * (devStrainTpdt33 - devStrainT33),
                        expFac * viscousStrainT[3] + \
                        dqV * (devStrainTpdt12 - devStrainT12)]

    devStressTpdt11 = elasFac * \
                      (viscousStrainTpdt[0] - devStrainInitial[0]) + \
                      devStressInitial[0]
    devStressTpdt22 = elasFac * \
                      (viscousStrainTpdt[1] - devStrainInitial[1]) + \
                      devStressInitial[1]
    devStressTpdt33 = elasFac * viscousStrainTpdt[2] + stateVars[0]
    devStressTpdt12 = elasFac * \
                      (viscousStrainTpdt[3] - devStrainInitial[2]) + \
                      devStressInitial[2]

    stressTpdt = [meanStressTpdt + devStressTpdt11,
                  meanStressTpdt + devStressTpdt22,
                  devStressTpdt12]
      
    return stressTpdt, viscousStrainTpdt

                                                        
  def _computeViscousFactor(self, maxwellTime):
    """
    Compute viscous strain factor for a given Maxwell time.
    """
    import math
    
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

  
  def _calcStress(self, strainV, muV, lambdaV, maxwellTimeV, totalStrainV,
                  viscousStrainV, initialStressV, initialStrainV, stateVars):
    """
    Compute stress, derivative of elasticity matrix, and updated state
    variables. This assumes behavior is always viscoelastic.
    """
    import scipy.misc

    # Define some numpy arrays
    strainTpdt = numpy.array(totalStrainV, dtype=numpy.float64)
    strainT = numpy.array(strainV, dtype=numpy.float64)
    viscousStrainT = numpy.array(viscousStrainV, dtype=numpy.float64)
    initialStress = numpy.array(initialStressV, dtype=numpy.float64)
    initialStrain = numpy.array(initialStrainV, dtype=numpy.float64)
    stressZZInitial = numpy.array([stateVars[0]], dtype=numpy.float64)
    
    # Compute current stress and viscous strain.
    dqV = self._computeViscousFactor(maxwellTimeV)
    stressTpdt, viscousStrainTpdt = self._computeStress(strainTpdt,
                                                        muV,
                                                        lambdaV,
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

  app = MaxwellPlaneStrainTimeDep()
  app.run()


# End of file 
