#!/usr/bin/env nemesis

# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2022 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------
#
# @file tests/fullscale/viscoelasticity/nofaults-2d/axialtraction_maxwell_soln.py
#
# @brief Simple test of gravitational stresses with updated coordinates and a tilted top surface.
#
# 2-D axial traction solution for linear Maxwell viscoelastic material.
#
#             Uy=0
#          ----------
#          |        |
# Ux=0     |        |  Tx=T0
#          |        |
#          |        |
#          ----------
#            Uy=0
#
# Dirichlet boundary conditions
# Ux(-4000,y) = 0
# Uy(x,-4000) = 0
# Uy(x,+4000) = 0
#
# Neumann boundary conditions
# Tx(+4000,y) = T0

import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import numpy
import h5py
import pylab
import pdb
pdb.set_trace()

# FE input file.
pylithInput = 'output/gravity_tilt_surf_quad-boundary.h5'

# Properties.
year = 60.0*60.0*24.0*365.25

# Finite element results.
h5 = h5py.File(pylithInput, 'r')
coords = h5['geometry/vertices'][:]
cellH5 = h5['topology/cells']
cellDim = cellH5.attrs['cell_dim']
cells = numpy.array(cellH5[:], dtype=numpy.int64)
time = h5['time'][:].flatten()
timeYears = time/year
numSteps = time.shape[0]
disp = h5['vertex_fields/displacement'][:]
h5.close()

# Velocity field.
dispDiff = numpy.diff(disp, axis=0)
timeDiff = numpy.diff(time)
velocity = dispDiff/timeDiff.reshape(numSteps - 1,1,1)

# Undeformed coordinates and sorted arrays.
coordsOrig = coords - disp[0,:,:]
sortInds = numpy.argsort(coordsOrig[:,0])
coordsOrigSort = coordsOrig[sortInds,:]
dispSort = disp[:,sortInds,:]
velocitySort = velocity[:,sortInds,:]

# Level surface plot.
xLevel = coordsOrigSort[:,0]
yLevel = numpy.zeros_like(xLevel)
yMean = numpy.mean(coordsOrigSort[:,1]) + numpy.mean(dispSort[0,:,1])

# Displacement plot.
color = cm.rainbow(numpy.linspace(0, 1, numSteps))
for stepNum in range(numSteps):
    pylab.plot(coordsOrigSort[:,0], dispSort[stepNum,:,1], c=color[stepNum,:])
pylab.plot(xLevel, yMean*numpy.ones_like(xLevel), 'k-')
pylab.show()

# Velocity plot.
for stepNum in range(numSteps - 1):
    pylab.plot(coordsOrigSort[:,0], velocitySort[stepNum,:,1], c=color[stepNum,:])
pylab.show()

# End of file
