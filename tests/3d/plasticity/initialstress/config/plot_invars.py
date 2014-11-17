#!/usr/bin/env python

import numpy
import h5py
import matplotlib.pyplot as plt
from pyre.units.time import year
from pyre.units.unitparser import parser
# import pdb

# pdb.set_trace()

# Input files.
elasPref = "../output/elastic_sshear_"
plasPref = "../output/plastic_sshear_"
elasSuff = "-elastic.h5"
plasSuff = "-plastic.h5"
stepSizes = ["dt01", "dt02", "dt05", "dt10"]

# Time and stress/strain normalization information.
totalTime = 20.0
uparser = parser()
timeUnits = uparser.parse('1.0*year').value
strainUnits = 1.0e-6
stressUnits = 1.0e6

# Line styles.
plasLines = ['k-', 'ro', 'bs', 'gD']
elasLines = ['k--', 'r--', 'b--', 'g--']

#-----------------------------------------------------------------------------
def _computeStressInvar(fileName):
  """
  Function to compute second deviatoric stress invariant.
  """

  h5 = h5py.File(fileName, "r", driver="sec2")
  times = h5['time'][:,0,0]
  timeYears = times/timeUnits
  numSteps = times.shape[0]
  stressInvarMean = numpy.zeros(numSteps, dtype=numpy.float64)

  stress = h5['cell_fields/stress'][:,:,:]/stressUnits
  h5.close()
  pres = (stress[:,:,0] + stress[:,:,1] + stress[:,:,2])/3.0
  s11 = stress[:,:,0] - pres
  s22 = stress[:,:,1] - pres
  s33 = stress[:,:,2] - pres
  s12 = stress[:,:,3]
  s23 = stress[:,:,4]
  s13 = stress[:,:,5]

  stressInvar = numpy.sqrt(0.5 * (s11 * s11 + s22 * s22 + s33 * s33 + \
                                  2.0 * (s12 * s12 + s13 * s13 + s23 * s23)))

  print "File %s:" % fileName
  for stepNum in range(numSteps):
    meanVal = numpy.mean(stressInvar[stepNum,:])
    minVal = numpy.amin(stressInvar[stepNum,:])
    maxVal = numpy.amax(stressInvar[stepNum,:])
    stdDev = numpy.std(stressInvar[stepNum,:])
    print "Step number:  %d" % stepNum
    print "Stress values in MPa:"
    print "  Mean of stress invariant:                %g" % meanVal
    print "  Minimum of stress invariant:             %g" % minVal
    print "  Maximum of stress invariant:             %g" % maxVal
    print "  Standard deviation of stress invariant   %g" % stdDev
    stressInvarMean[stepNum] = meanVal

  return (stressInvarMean, timeYears)


def _computeStrainInvar(fileName):
  """
  Function to compute second plastic strain invariant.
  """

  h5 = h5py.File(fileName, "r", driver="sec2")
  times = h5['time'][:,0,0]
  timeYears = times/timeUnits
  numSteps = times.shape[0]
  strainInvarMean = numpy.zeros(numSteps, dtype=numpy.float64)

  plasStrain = h5['cell_fields/plastic_strain'][:,:,:]/strainUnits
  h5.close()
  e11 = plasStrain[:,:,0]
  e22 = plasStrain[:,:,1]
  e33 = plasStrain[:,:,2]
  e12 = plasStrain[:,:,3]
  e23 = plasStrain[:,:,4]
  e13 = plasStrain[:,:,5]

  strainInvar = numpy.sqrt(0.5 * (e11 * e11 + e22 * e22 + e33 * e33 + \
                                  2.0 * (e12 * e12 + e13 * e13 + e23 * e23)))

  print "File %s:" % fileName
  for stepNum in range(numSteps):
    meanVal = numpy.mean(strainInvar[stepNum,:])
    minVal = numpy.amin(strainInvar[stepNum,:])
    maxVal = numpy.amax(strainInvar[stepNum,:])
    stdDev = numpy.std(strainInvar[stepNum,:])
    print "Step number:  %d" % stepNum
    print "Strain values in microStrain:"
    print "  Mean of plastic strain invariant:                %g" % meanVal
    print "  Minimum of plastic strain invariant:             %g" % minVal
    print "  Maximum of plastic strain invariant:             %g" % maxVal
    print "  Standard deviation of plastic strain invariant:  %g" % stdDev
    strainInvarMean[stepNum] = meanVal

  return (strainInvarMean, timeYears)
  

#-----------------------------------------------------------------------------
numStepSizes = len(stepSizes)

pStrainTimes = []
pStrains = []
pStressTimes = []
eStressTimes = []
pStresses = []
eStresses = []

for stepSize in range(numStepSizes):
  eStressFile = elasPref + stepSizes[stepSize] + elasSuff
  (eStress, eStressTime) = _computeStressInvar(eStressFile)
  pStressFile = plasPref + stepSizes[stepSize] + plasSuff
  (pStress, pStressTime) = _computeStressInvar(pStressFile)
  pStrainFile = plasPref + stepSizes[stepSize] + plasSuff
  (pStrain, pStrainTime) = _computeStrainInvar(pStrainFile)
  pStrainTimes.append(pStrainTime)
  pStrains.append(pStrain)
  pStressTimes.append(pStressTime)
  pStresses.append(pStress)
  eStressTimes.append(eStressTime)
  eStresses.append(eStress)

plt.figure(1)
plt.subplot(121)
plt.plot(pStrainTimes[0], pStrains[0], plasLines[0],
         pStrainTimes[1], pStrains[1], plasLines[1],
         pStrainTimes[2], pStrains[2], plasLines[2],
         pStrainTimes[3], pStrains[3], plasLines[3])
plt.xlabel('Time (years)')
plt.ylabel('Plastic strain invariant (microStrain)')

plt.subplot(122)
plt.plot(pStressTimes[0], pStresses[0], plasLines[0],
         pStressTimes[1], pStresses[1], plasLines[1],
         pStressTimes[2], pStresses[2], plasLines[2],
         pStressTimes[3], pStresses[3], plasLines[3],
         eStressTimes[0], eStresses[0], elasLines[1],
         eStressTimes[1], eStresses[1], elasLines[1],
         eStressTimes[2], eStresses[2], elasLines[2],
         eStressTimes[3], eStresses[3], elasLines[3])
plt.xlabel('Time (years)')
plt.ylabel('Stress invariant (MPa)')

plt.show()
