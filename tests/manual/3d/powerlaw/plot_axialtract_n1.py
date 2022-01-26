#!/usr/bin/env python

# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2018 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#
# @file tests/manual/3d/powerlaw/plot_axialtract_n1.py
#
# @brief Plot power-law solutions using different time step sizes and compare to analytical Maxwell.
#

import numpy
import h5py
import matplotlib.pyplot as plt
from axialtraction_maxwell_soln import AnalyticalSoln

# Year in seconds.
year = 60.0*60.0*24.0*365.25

# Input files.
stepSizes = numpy.array([0.01, 0.02, 0.05, 0.1, 0.2], dtype=numpy.float64)
stepSizesStr = ['0.01', '0.02', '0.05', '0.1', '0.2']
lineDefs = ['r+', 'g+', 'b+', 'c+', 'm+', 'y+']
lineDefAnl = 'k-'
legendAnl = 'Analytical'
numSims = len(stepSizes)
basePrefix = 'axialtraction_powerlaw_n1_norefstate_dt'
baseSuffix = '_onecell_hex'

    
def getVars(fileName):
    """
    Get time, sxx, ezz, and dispz from given file.
    """
    h5 = h5py.File(fileName, 'r')
    time = h5['time'][:].flatten()
    timeYears = time/year
    locs = h5['geometry/vertices'][:]
    stress = h5['vertex_fields/cauchy_stress'][:]
    strain = h5['vertex_fields/cauchy_strain'][:]
    disp = h5['vertex_fields/displacement'][:]
    sxx = stress[:,0,0]
    ezz = strain[:,0,2]
    dispz = disp[:,0,2]

    h5.close()

    return (timeYears, sxx, ezz, dispz, locs)


def scanLogfile(fileName):
    """
    Get Jacobian info from log file.
    """
    matchStr = '||J - Jfd||_F/||J||_F ='
    endStr = ','
    f = open(fileName, 'r')
    lines = reversed(f.readlines())
    lastLine = None
    for line in lines:
        if (matchStr in line):
            lastLine = line
            break
    val = float((lastLine.split(matchStr))[1].split(endStr)[0])

    return val
    
# Create subplots and loop over simulations.
fig, a = plt.subplots(2,2)
jacobianDiff = numpy.zeros(numSims, dtype=numpy.float64)

for simNum in range(numSims):
    baseName = basePrefix + stepSizesStr[simNum] + baseSuffix
    h5File = 'output/' + baseName + '-viscomat.h5'
    logFile = baseName + '.log'
    (timeYears, sxx, ezz, dispz, locs) = getVars(h5File)
    timeSecs = year*timeYears
    jacobianDiff[simNum] = scanLogfile(logFile)
    if (simNum == 0):
        (dispzAnl, stressAnl, devStressAnl, strainAnl, devStrainAnl,
         maxwellVisStrain, powerLawVisStrain) = AnalyticalSoln(timeSecs, locs[0,:].reshape(1,3))
        a[0][0].plot(timeYears, stressAnl[:,0], lineDefAnl, label=legendAnl)
        a[1][0].plot(timeYears, strainAnl[:,2], lineDefAnl, label=legendAnl)
        a[0][1].plot(timeYears, dispzAnl, lineDefAnl, label=legendAnl)
    a[0][0].plot(timeYears, sxx, lineDefs[simNum], label="dt={0}".format(stepSizesStr[simNum]))
    a[1][0].plot(timeYears, ezz, lineDefs[simNum], label="dt={0}".format(stepSizesStr[simNum]))
    a[0][1].plot(timeYears, dispz, lineDefs[simNum], label="dt={0}".format(stepSizesStr[simNum]))

a[1][1].loglog(stepSizes, jacobianDiff, 'k+-')
a[0][0].set_xlabel('Time (years)')
a[0][0].set_ylabel('Stress_xx (Pa)')
a[0][0].legend(loc="upper right")
a[1][0].set_xlabel('Time (years)')
a[1][0].set_ylabel('Strain_zz')
a[1][0].legend(loc="upper right")
a[0][1].set_xlabel('Time (years)')
a[0][1].set_ylabel('Displacement_z (m)')
a[0][1].legend(loc="upper right")
a[1][1].set_xlabel('Time step size (years)')
a[1][1].set_ylabel('Jacobian difference')
plt.show()

# End of file
