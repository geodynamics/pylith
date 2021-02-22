#!/usr/bin/env python

import numpy
import h5py
import matplotlib.pyplot as pyplot
from pythia.pyre.units.time import year
from pythia.pyre.units.unitparser import parser

stepSizes = ["dt01", "dt02", "dt05", "dt10"]

# Time and stress/strain normalization information.
uparser = parser()
timeUnits = uparser.parse('1.0*year').value
strainUnits = 1.0e-6
stressUnits = uparser.parse('1.0*MPa').value

# Line styles.
stylePlastic = ['k-', 'ro', 'bs', 'gD']
styleElastic = ['k--', 'r--', 'b--', 'g--']

# -----------------------------------------------------------------------------


def computeStressInvar(fileName):
    """
    Function to compute second deviatoric stress invariant.
    """

    h5 = h5py.File(fileName, "r", driver="sec2")
    times = h5['time'][:, 0, 0]
    timeYears = times / timeUnits
    numSteps = times.shape[0]
    stressInvarMean = numpy.zeros(numSteps, dtype=numpy.float64)

    stress = h5['cell_fields/stress'][:,:,:] / stressUnits
    h5.close()
    pres = (stress[:,:, 0] + stress[:,:, 1] + stress[:,:, 2]) / 3.0
    s11 = stress[:,:, 0] - pres
    s22 = stress[:,:, 1] - pres
    s33 = stress[:,:, 2] - pres
    s12 = stress[:,:, 3]
    s23 = stress[:,:, 4]
    s13 = stress[:,:, 5]

    stressInvar = numpy.sqrt(0.5 * (s11 * s11 + s22 * s22 + s33 * s33 +
                                    2.0 * (s12 * s12 + s13 * s13 + s23 * s23)))

    print("File %s:" % fileName)
    for stepNum in range(numSteps):
        meanVal = numpy.mean(stressInvar[stepNum,:])
        minVal = numpy.amin(stressInvar[stepNum,:])
        maxVal = numpy.amax(stressInvar[stepNum,:])
        stdDev = numpy.std(stressInvar[stepNum,:])
        print("Step number:  %d" % stepNum)
        print("Stress values in MPa:")
        print("  Mean of stress invariant:                %g" % meanVal)
        print("  Minimum of stress invariant:             %g" % minVal)
        print("  Maximum of stress invariant:             %g" % maxVal)
        print("  Standard deviation of stress invariant   %g" % stdDev)
        stressInvarMean[stepNum] = meanVal

    return (stressInvarMean, timeYears)


def computeStrainInvar(fileName):
    """
    Function to compute second plastic strain invariant.
    """

    h5 = h5py.File(fileName, "r", driver="sec2")
    times = h5['time'][:, 0, 0]
    timeYears = times / timeUnits
    numSteps = times.shape[0]
    strainInvarMean = numpy.zeros(numSteps, dtype=numpy.float64)

    plasStrain = h5['cell_fields/plastic_strain'][:,:,:] / strainUnits
    h5.close()
    e11 = plasStrain[:,:, 0]
    e22 = plasStrain[:,:, 1]
    e33 = plasStrain[:,:, 2]
    e12 = plasStrain[:,:, 3]
    e23 = plasStrain[:,:, 4]
    e13 = plasStrain[:,:, 5]

    strainInvar = numpy.sqrt(0.5 * (e11 * e11 + e22 * e22 + e33 * e33 +
                                    2.0 * (e12 * e12 + e13 * e13 + e23 * e23)))

    print("File %s:" % fileName)
    for stepNum in range(numSteps):
        meanVal = numpy.mean(strainInvar[stepNum,:])
        minVal = numpy.amin(strainInvar[stepNum,:])
        maxVal = numpy.amax(strainInvar[stepNum,:])
        stdDev = numpy.std(strainInvar[stepNum,:])
        print("Step number:  %d" % stepNum)
        print("Strain values in microStrain:")
        print("  Mean of plastic strain invariant:                %g" % meanVal)
        print("  Minimum of plastic strain invariant:             %g" % minVal)
        print("  Maximum of plastic strain invariant:             %g" % maxVal)
        print("  Standard deviation of plastic strain invariant:  %g" % stdDev)
        strainInvarMean[stepNum] = meanVal

    return (strainInvarMean, timeYears)


# -----------------------------------------------------------------------------
numStepSizes = len(stepSizes)

pyplot.figure(1)
ax1 = pyplot.subplot(121)
ax2 = pyplot.subplot(122)

for istep, stepSize in enumerate(stepSizes):

    # Elastic
    filename = "output/elastic_%s-statevars.h5" % stepSize
    (stressInvarElastic, timeElastic) = computeStressInvar(filename)

    # Plastic
    filename = "output/plastic_%s-statevars.h5" % stepSize
    (stressInvarPlastic, timePlastic) = computeStressInvar(filename)
    (strainInvarPlastic, timePlastic) = computeStrainInvar(filename)

    ax1.plot(timePlastic, strainInvarPlastic, stylePlastic[istep])

    ax2.plot(timeElastic, stressInvarElastic, styleElastic[istep],
             timePlastic, stressInvarPlastic, stylePlastic[istep])


ax1.set_xlabel('Time (years)')
ax1.set_ylabel('Plastic strain invariant (microStrain)')

ax2.set_xlabel('Time (years)')
ax2.set_ylabel('Stress invariant (MPa)')

pyplot.show()
