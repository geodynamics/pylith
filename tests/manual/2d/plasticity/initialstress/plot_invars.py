#!/usr/bin/env python

import numpy
import h5py
import matplotlib.pyplot as pyplot

from pythia.pyre.units.time import year
from pythia.pyre.units.unitparser import parser

# Input files.
stepSizes = ["dt01", "dt02", "dt05", "dt10"]


# Time and stress/strain normalization information.
uparser = parser()
timeUnits = uparser.parse('1.0*year').value
strainUnits = 1.0e-6
stressUnits = uparser.parse("1.0*MPa").value

# Line styles.
stylePlastic = ['k-', 'ro', 'bs', 'gD']
styleElastic = ['k--', 'r--', 'b--', 'g--']

# -----------------------------------------------------------------------------


def computeStressInvar(fileName):
    """
    Function to compute second deviatoric stress invariant.
    Exclude this for now because we need stress_zz to do it correctly.
    """

    h5 = h5py.File(fileName, "r", driver="sec2")
    time = h5['time'][:]
    stress = h5['cell_fields/stress'][:] / stressUnits
    h5.close()

    timeYears = time / timeUnits
    numSteps = time.shape[0]
    stressInvarMean = numpy.zeros(numSteps, dtype=numpy.float64)

    pres = (stress[:,:, 0] + stress[:,:, 1] + stress[:,:, 2]) / 3.0
    s11 = stress[:,:, 0] - pres
    s22 = stress[:,:, 1] - pres
    s12 = stress[:,:, 2]

    stressInvar = numpy.sqrt(0.5 * (s11 * s11 + s22 * s22 + 2.0 * (s12 * s12)))

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

    strainInvar = numpy.sqrt(0.5 * (e11 * e11 + e22 * e22 + e33 * e33 +
                                    2.0 * (e12 * e12)))

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

# :TODO: Add stress computations, including effects of stress_zz.

pyplot.figure(1)
ax = pyplot.subplot(111)

for istep, stepSize in enumerate(stepSizes):
    filename = "output/plastic_%s-statevars.h5" % stepSize
    (strainInvarMean, time) = computeStrainInvar(filename)
    ax.plot(time, strainInvarMean, stylePlastic[istep])

ax.set_xlabel('Time (years)')
ax.set_ylabel('Plastic strain invariant (microStrain)')

pyplot.show()
