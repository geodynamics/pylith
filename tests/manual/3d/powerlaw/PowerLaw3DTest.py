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
# Copyright (c) 2010-2017 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file tests/manual/3d/powerlaw/PowerLaw3DTest.py

## @brief Python script to test power-law implementation.

import numpy
# import pdb
# pdb.set_trace()

# ----------------------------------------------------------------------
dimension = 3
numElasticConsts = 36
tensorSize = 6
numLocs = 4
year = 60.0*60.0*24.0*365.25
debugOutput = True
debugging = True

# Examples:
# axialstrain_powerlaw_n1_norefstate
# axialstrain_powerlaw_norefstate
# axialtraction_powerlaw_norefstate
# shearstrain_powerlaw_norefstate

# Parameter values.
alpha = 0.5
dt = numpy.array([0.025, 0.05, 0.1, 0.1], dtype=numpy.float64)
density = numpy.array([2500.0, 2500.0, 2500.0, 2500.0], dtype=numpy.float64)
vs = numpy.array([3464.1016, 3464.1016, 3464.1016, 3464.1016], dtype=numpy.float64)
vp = numpy.array([6000.0, 6000.0, 6000.0, 6000.0], dtype=numpy.float64)
powerLawReferenceStrainRate = numpy.array([1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6], dtype=numpy.float64)
powerLawReferenceStress = numpy.array([5680367950354.795, 1.798919e+10, 1.798919e+10, 1.798919e+10], dtype=numpy.float64)
powerLawExponent = numpy.array([1.00, 3.50, 3.50, 3.50], dtype=numpy.float64)
totalStrain = numpy.array([[0.00, 0.00, 0.00025, 0.00, 0.00, 0.00],
                           [0.00, 0.00, 0.00025, 0.00, 0.00, 0.00],
                           [0.00, 0.00, -0.0111111, 0.00, 0.00, 0.00],
                           [0.00, 0.00, 0.00, 0.00, 0.000318649, 0.000238987]], dtype=numpy.float64)
initialStress = numpy.zeros((numLocs, tensorSize), dtype=numpy.float64)
initialStrain = numpy.zeros((numLocs, tensorSize), dtype=numpy.float64)
visStrainT = numpy.zeros((numLocs, tensorSize), dtype=numpy.float64)
stressT = numpy.zeros((numLocs, tensorSize), dtype=numpy.float64)

# Normalization.
lengthScale = numpy.array([1000.0, 1000.0, 1000.0, 1000.0], dtype=numpy.float64)
timeScale = year*numpy.array([1.0, 600000.0, 50.0, 50000.0], dtype=numpy.float64)
pressureScale = numpy.array([3.0e10, 3.0e10, 3.0e10, 3.0e10], dtype=numpy.float64)
densityScale = numpy.array([2.98765e19, 1.07555e+31, 7.46912e+22, 7.46912e+28], dtype=numpy.float64)

# Derived parameters.
strainRateScale = 1.0/timeScale
mu = vs*vs*density
lam = vp*vp*density - 2.0*mu
                           
# Normalized properties.
densityND = density/densityScale
muND = mu/pressureScale
lamND = lam/pressureScale
powerLawReferenceStrainRateND = powerLawReferenceStrainRate/strainRateScale
powerLawReferenceStressND = powerLawReferenceStress/pressureScale
initialStressND = (initialStress.transpose()/pressureScale).transpose()
stressTND = (stressT.transpose()/pressureScale).transpose()

# Arrays to hold solutions.
stress = numpy.zeros( (numLocs, tensorSize), dtype=numpy.float64)
stressUpdated = numpy.zeros( (numLocs, tensorSize), dtype=numpy.float64)
visStrainUpdated = numpy.zeros( (numLocs, tensorSize), dtype=numpy.float64)
elasticConsts = numpy.zeros( (numLocs, tensorSize, tensorSize), dtype=numpy.float64)
maxwellTime = numpy.zeros(numLocs, dtype=numpy.float64)
dtStableImplicit = numpy.zeros(numLocs, dtype=numpy.float64)
diagVec = numpy.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0], dtype=numpy.float64)


def bracket(effStressInitialGuess, ae, b, c, d, alpha, dt, effStressT, powerLawExponentV, refStrainRateV, refStressV):
    """
    Function to bracket the effective stress.
    """
    maxIterations = 50
    bracketFactor = 1.6

    x1 = 0.0
    x2 = 0.0

    if effStressInitialGuess > 0.0:
        x1 = effStressInitialGuess - 0.5 * effStressInitialGuess
        x2 = effStressInitialGuess + 0.5 * effStressInitialGuess
    else:
        x1 = 500.0
        x2 = 1500.0

    funcValue1 = effStressFunc(x1, ae, b, c, d, alpha, dt, effStressT, powerLawExponentV, refStrainRateV, refStressV)
    funcValue2 = effStressFunc(x2, ae, b, c, d, alpha, dt, effStressT, powerLawExponentV, refStrainRateV, refStressV)

    iteration = 0
    bracketed = False

    while iteration < maxIterations:
        if (funcValue1 * funcValue2) < 0.0:
            bracketed = True
            break
        if abs(funcValue1) < abs(funcValue2):
            x1 += bracketFactor * (x1 - x2)
            x1 = max(x1, 0.0)
            funcValue1 = effStressFunc(x1, ae, b, c, d, alpha, dt, effStressT, powerLawExponentV, refStrainRateV, refStressV)
        else:
            x2 += bracketFactor * (x1 - x2)
            x2 = max(x2, 0.0)
            funcValue2 = effStressFunc(x2, ae, b, c, d, alpha, dt, effStressT, powerLawExponentV, refStrainRateV, refStressV)

        iteration += 1

    if bracketed == False:
        raise RuntimeError("Unable to bracket root.")

    return x1, x2
    
    
def getMaxwellTime(mu, refStrainRate, refStress, powerLawExponent, stress):
    """
    Compute Maxwell time from stress, reference stress and strain rate, shear
    modulus, and power-law exponent.
    """
    meanStress = (stress[0] + stress[1] + stress[2])/3.0
    devStress = numpy.array(stress, dtype=numpy.float64)
    
    devStress[0] = stress[0] - meanStress
    devStress[1] = stress[1] - meanStress
    devStress[2] = stress[2] - meanStress
    devStress[3] = stress[3]
    devStress[4] = stress[4]
    devStress[5] = stress[5]

    devStressProd = scalarProduct(devStress, devStress)
    effStress = (0.5 * devStressProd)**0.5
    maxwellTime = 1.0
    if (effStress != 0.0):
        maxwellTime = (refStress/effStress)**(powerLawExponent - 1.0) * (refStress/mu)/(refStrainRate * 6.0)

    return maxwellTime


def scalarProduct(tensor1, tensor2):
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

    
def calcStressComponent(strainVal, strainComp, stressComp, strainTpdt, muV, lambdaV, refStrainRateV, refStressV,
                        powerLawExponentV, visStrainT, stressT, initialStress, initialStrain, alpha, dt, debugOutput):
    """
    Function to compute a particular stress component as a function of a
    strain component.
    """
    strainTest = numpy.array(strainTpdt, dtype=numpy.float64)
    strainTest[strainComp] = strainVal
    stressTpdt, visStrainTpdt = computeStress(strainTest, muV, lambdaV, refStrainRateV, refStressV, powerLawExponentV,
                                              visStrainT, stressT, initialStress, initialStrain, alpha, dt, debugOutput)
    return stressTpdt[stressComp]


def effStressFunc(effStressTpdt, ae, b, c, d, alpha, dt, effStressT, powerLawExponentV, refStrainRateV, refStressV):
    """
    Function to compute effective stress function for a given effective stress.
    """

    factor1 = 1.0 - alpha
    effStressTau = factor1 * effStressT + alpha * effStressTpdt
    gammaTau = refStrainRateV * (effStressTau/refStressV)** (powerLawExponentV - 1.0) / refStressV
    a = ae + alpha * dt * gammaTau
    effStressFunc = a * a * effStressTpdt * effStressTpdt - b + c * gammaTau - d * d * gammaTau * gammaTau

    return effStressFunc

    
def computeStress(strainTpdt, muV, lambdaV, refStrainRateV, refStressV, powerLawExponentV, visStrainT, stressT,
                  initialStress, initialStrain, alpha, dt, debugOutput):
    """
    Function to compute stresses and viscous strains using the effective stress function algorithm.
    """
    import scipy.optimize
    
    # Constants
    mu2 = 2.0 * muV
    lamPlusMu = lambdaV + muV
    bulkModulus = lambdaV + 2.0 * muV/3.0
    ae = 1.0/mu2
    timeFac = dt * (1.0 - alpha)
    diag = numpy.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0], dtype=numpy.float64)

    # Initial stress values
    meanStressInitial = (initialStress[0] + initialStress[1] + initialStress[2])/3.0
    devStressInitial = initialStress - meanStressInitial * diag
    stressInvar2Initial = 0.5 * scalarProduct(devStressInitial, devStressInitial)
    
    # Initial strain values
    meanStrainInitial = (initialStrain[0] + initialStrain[1] + initialStrain[2])/3.0

    # Values for current time step
    meanStrainTpdt = (strainTpdt[0] + strainTpdt[1] + strainTpdt[2])/3.0 - meanStrainInitial
    meanStressTpdt = 3.0 * bulkModulus * meanStrainTpdt
    strainPPTpdt = strainTpdt - meanStrainTpdt * diag - visStrainT - initialStrain
    strainPPInvar2Tpdt = 0.5 * scalarProduct(strainPPTpdt, strainPPTpdt)

    # Values for previous time step
    meanStressT = (stressT[0] + stressT[1] + stressT[2])/3.0
    devStressT = stressT - diag * meanStressT
    stressInvar2T = 0.5 * scalarProduct(devStressT, devStressT)
    effStressT = stressInvar2T**0.5

    # Finish defining parameters needed for root-finding algorithm.
    b = strainPPInvar2Tpdt + ae*scalarProduct(strainPPTpdt, devStressInitial) + ae*ae*stressInvar2Initial
    c = (scalarProduct(strainPPTpdt, devStressT) + ae*scalarProduct(devStressT, devStressInitial))*timeFac
    d = timeFac*effStressT

    # Bracket the root
    effStressInitialGuess = effStressT

    x1, x2 = bracket(effStressInitialGuess, ae, b, c, d, alpha, dt, effStressT, powerLawExponentV, refStrainRateV, refStressV)

    # Find the root using Brent's method (from scipy)
    rootTolerance = 1.0e-14
    effStressTpdt = scipy.optimize.brentq(effStressFunc, x1, x2, args=(ae, b, c, d, alpha, dt, effStressT, powerLawExponentV,
                                                                       refStrainRateV, refStressV), xtol=rootTolerance)
    
    # Compute stresses from the effective stress.
    effStressTau = (1.0 - alpha) * effStressT + alpha * effStressTpdt
    gammaTau = refStrainRateV * ((effStressTau/refStressV)**(powerLawExponentV - 1.0)) / refStressV
    factor1 = 1.0/(ae + alpha * dt * gammaTau)
    factor2 = timeFac * gammaTau
    devStressTpdt = 0.0
    stressTpdt = numpy.zeros( (tensorSize), dtype=numpy.float64)
    visStrainTpdt = numpy.zeros( (tensorSize), dtype=numpy.float64)

    printOutput = debugging and debugOutput
    if printOutput:
        print('  b:  %#16.12g' % b)
        print('  c:  %#16.12g' % c)
        print('  d:  %#16.12g' % d)
        print('  meanStressTpdt: %#16.12g' % meanStressTpdt)
        print('  effStressTpdt:  %#16.12g' % effStressTpdt)
        print('  effStressTau:   %#16.12g' % effStressTau)
        print('  effStressT:     %#16.12g' % effStressT)
        print('  gammaTau:       %#16.12g' % gammaTau)
        print('  factor1:        %#16.12g' % factor1)
        print('  factor2:        %#16.12g' % factor2)
      
    for iComp in range(tensorSize):
        devStressTpdt = factor1 * (strainPPTpdt[iComp] - factor2 * devStressT[iComp] + ae * devStressInitial[iComp])
        stressTpdt[iComp] = devStressTpdt + diag[iComp] * (meanStressTpdt + meanStressInitial)
        devStressTau = (1.0 - alpha) * devStressT[iComp] + alpha * devStressTpdt
        deltaVisStrain = dt * gammaTau * devStressTau
        visStrainTpdt[iComp] = visStrainT[iComp] + deltaVisStrain

    return stressTpdt, visStrainTpdt

  
def calcStress(strainV, muV, lambdaV, refStrainRateV, refStressV, powerLawExponentV,visStrainV, stressV,
               initialStressV, initialStrainV, alpha, dt):
    """
    Compute stress, updated state variables and derivative of elasticity matrix.
    This assumes behavior is always viscoelastic.
    """
    import scipy.misc

    # Define some numpy arrays
    strainTpdt = numpy.array(strainV, dtype=numpy.float64)
    visStrainT = numpy.array(visStrainV, dtype=numpy.float64)
    stressT = numpy.array(stressV, dtype=numpy.float64)
    initialStress = numpy.array(initialStressV, dtype=numpy.float64)
    initialStrain = numpy.array(initialStrainV, dtype=numpy.float64)
    debugOutput = True

    stressTpdt, visStrainTpdt = computeStress(strainTpdt, muV, lambdaV, refStrainRateV, refStressV, powerLawExponentV, visStrainT, stressT,
                                              initialStress, initialStrain, alpha, dt, debugOutput)

    # Compute components of tangent constitutive matrix using numerical
    # derivatives.
    # derivDx = 1.0e-16
    cutOff = 1.0e-10
    derivOrder = 3
    testStrain = numpy.abs(strainTpdt)
    aboveCutoff = numpy.where(testStrain > cutOff)
    minStrain = numpy.amin(testStrain[aboveCutoff])
    derivDx = max(1.0e-2*minStrain, 1.0e-12)
    debugOutput = False

    if debugging:
        print('  derivDx:  %#16.12g' % derivDx)

    elasticConsts = numpy.zeros((tensorSize, tensorSize), dtype=numpy.float64)
    for stressComp in range(tensorSize):
        for strainComp in range(tensorSize):
            elasticConsts[stressComp, strainComp] = \
                                      scipy.misc.derivative(calcStressComponent, strainTpdt[strainComp], dx=derivDx,
                                                            args=(strainComp, stressComp, strainTpdt, muV, lambdaV, refStrainRateV,
                                                                  refStressV, powerLawExponentV, visStrainT, stressT, initialStress,
                                                                  initialStrain, alpha, dt, debugOutput), order=derivOrder)


    return (elasticConsts, numpy.ravel(stressTpdt), numpy.ravel(stressTpdt), numpy.ravel(visStrainTpdt))
  

# MAIN /////////////////////////////////////////////////////////////////

for locNum in range(numLocs):
    print("")
    print("Simulation number %d:" % locNum)
    print("")
    (elasticConsts[locNum,:,:], stress[locNum,:], stressUpdated[locNum,:], visStrainUpdated[locNum,:]) = \
                                calcStress(totalStrain[locNum,:], muND[locNum], lamND[locNum], powerLawReferenceStrainRateND[locNum],
                                           powerLawReferenceStressND[locNum], powerLawExponent[locNum], visStrainT[locNum,:],
                                           stressTND[locNum,:], initialStressND[locNum,:], initialStrain[locNum,:], alpha, dt[locNum])
    meanStress = numpy.sum(stress[locNum, 0:3])/3.0
    devStress = stress[locNum,:] - meanStress*diagVec
    print("    Stress:")
    print(stress[locNum, :])
    print("    Stress updated:")
    print(stressUpdated[locNum, :])
    print("    Deviatoric stress:")
    print(devStress)
    print("    Viscous strain updated:")
    print(visStrainUpdated[locNum, :])
    print("    Total strain:")
    print(totalStrain[locNum, :])
    print("    Jacobian:")
    print(elasticConsts[locNum, :,:])

    # Use state variables to compute Maxwell times (and stable time step size).
    maxwellTime[locNum] = getMaxwellTime(muND[locNum], powerLawReferenceStrainRateND[locNum], powerLawReferenceStressND[locNum],
                                         powerLawExponent[locNum], stress[locNum,:])
    dtStableImplicit[locNum] = 0.2 * maxwellTime[locNum]
    print("    Maxwell time:")
    print(maxwellTime[locNum])


# End of file 
