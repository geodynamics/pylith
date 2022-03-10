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
# @file tests/manual/2d/powerlaw/plot_sheartract.py
#
# @brief Plot power-law solutions using different time step sizes.
#

import numpy
import h5py
import matplotlib.pyplot as plt

# Year in seconds.
year = 60.0*60.0*24.0*365.25

# Input files.
lineDefs = ['k-', 'r+', 'g+', 'b+', 'c+', 'm+']
basePrefix = 'sheartraction_powerlaw_dt'

def getVars(fileName):
    """
    Get time, second deviatoric stress invariant, exy, dispy, and vertex coords from given file.
    """
    h5 = h5py.File(fileName, 'r')
    time = h5['time'][:].flatten()
    timeYears = time/year
    locs = h5['geometry/vertices'][:]
    devStress = h5['vertex_fields/deviatoric_stress'][:,0,:]
    strain = h5['vertex_fields/cauchy_strain'][:,0,:]
    disp = h5['vertex_fields/displacement'][:,0,:]
    stressInvar = numpy.sqrt(devStress[:,0]*devStress[:,0] + devStress[:,1]*devStress[:,1] + devStress[:,2]*devStress[:,2] + \
                             2.0*devStress[:,3]*devStress[:,3])
    exy = strain[:,3]
    dispy = disp[:,1]

    h5.close()

    return (timeYears, stressInvar, exy, dispy, locs)


def scanLogfile(fileName):
    """
    Get Jacobian info from log file.
    """
    matchStr = '||J - Jfd||_F/||J||_F ='
    endStr = ','
    lastLine = None
    val = None
    try:
        f = open(fileName, 'r')
        lines = reversed(f.readlines())
        for line in lines:
            if (matchStr in line):
                lastLine = line
                break
        if (lastLine):
            val = float((lastLine.split(matchStr))[1].split(endStr)[0])
    except FileNotFoundError:
        pass

    return val

    
def run(axs, colNum, stepSizes, useJacobian):
    """
    Create subplots and loop over simulations.
    """
    numSims = len(stepSizes)
    jacobianDiff = numpy.zeros(numSims, dtype=numpy.float64)
    jacobianInfo = False

    for simNum in range(numSims):
        dt = stepSizes[simNum]
        dtStr = repr(dt)
        baseName = basePrefix + dtStr
        h5File = 'output/' + baseName + '-viscomat.h5'
        logFile = baseName + '.log'
        (timeYears, stressInvar, exy, dispy, locs) = getVars(h5File)
        jacobianDiff[simNum] = scanLogfile(logFile)
        if (jacobianDiff[simNum]):
            jacobianInfo = True
        axs[0][colNum].plot(timeYears, stressInvar, lineDefs[simNum], label="dt={0}".format(dtStr))
        axs[1][colNum].plot(timeYears, exy, lineDefs[simNum], label="dt={0}".format(dtStr))
        axs[2][colNum].plot(timeYears, dispy, lineDefs[simNum], label="dt={0}".format(dtStr))

    axs[0][colNum].set_xlabel('Time (years)')
    axs[0][colNum].set_ylabel('2nd Dev Stress Invar (Pa)')
    axs[0][colNum].legend(loc="upper right")
    axs[0][colNum].set_title('Shear traction')
    axs[1][colNum].set_xlabel('Time (years)')
    axs[1][colNum].set_ylabel('Strain_xy')
    axs[1][colNum].legend(loc="upper right")
    axs[2][colNum].set_xlabel('Time (years)')
    axs[2][colNum].set_ylabel('Displacement_y (m)')
    axs[2][colNum].legend(loc="upper right")
    jacobianPlot = useJacobian and jacobianInfo
    if (jacobianPlot):
        axs[3][colNum].loglog(stepSizes, jacobianDiff, 'k+-')
        axs[3][colNum].set_xlabel('Time step size (years)')
        axs[3][colNum].set_ylabel('Jacobian difference')

    return

# End of file
