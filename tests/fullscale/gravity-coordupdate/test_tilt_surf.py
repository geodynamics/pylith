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

# FE input file.
pylithInput = 'output/grav_tilt_surf_norefstate_quad-boundary.h5'

# Properties.
year = 60.0*60.0*24.0*365.25

# Finite element results.
h5 = h5py.File(pylithInput, 'r')
coords = h5['geometry/vertices'][:]
cellH5 = h5['viz/topology/cells']
cellDim = cellH5.attrs['cell_dim']
cells = numpy.array(cellH5[:], dtype=numpy.int64)
time = h5['time'][:].flatten()
timeYears = time/year
numSteps = time.shape[0]
dispDiff = h5['vertex_fields/displacement'][:]
h5.close()

# Velocity field.
# dispDiff = numpy.diff(disp, axis=0)
dt0 = time[1] - time[0]
timeDiff = numpy.diff(time, prepend=-dt0)
velocity = dispDiff/timeDiff.reshape(numSteps,1,1)
velocityScaled = velocity*year

# Undeformed coordinates and sorted arrays.
# coordsOrig = coords - disp[0,:,:]
sortInds = numpy.argsort(coords[:,0])
coordsSort = coords[sortInds,:]
dispSort = dispDiff[:,sortInds,:]
velocitySort = velocity[:,sortInds,:]
velocityScaledSort = velocityScaled[:,sortInds,:]

# Plots.
fig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(12,4))
xLevel = coordsSort[:,0]
yLevel = numpy.zeros_like(xLevel)
yMean = numpy.mean(coordsSort[:,1])

# Deformed coordinates plot.
color = cm.rainbow(numpy.linspace(0, 1, numSteps))
coordsDeformed = coordsSort.copy()
for stepNum in range(numSteps):
    coordsDeformed += dispSort[stepNum,:,:]
    ax1.plot(coordsDeformed[:,0], coordsDeformed[:,1], c=color[stepNum,:])

ax1.plot(coordsSort[:,0], yMean*numpy.ones_like(xLevel), 'k-')
ax1.set_title('Upper surface position')
ax1.set(xlabel='X position (m)', ylabel='Y position (m)')

# Displacement plot.
for stepNum in range(numSteps):
    ax2.plot(coordsSort[:,0], dispSort[stepNum,:,1], c=color[stepNum,:])
ax2.set_title('Upper surface displacement increment')
ax2.set(xlabel='X position (m)', ylabel='Y displacement (m)')

# Velocity plot.
for stepNum in range(numSteps - 1):
    ax3.plot(coordsSort[:,0], velocitySort[stepNum,:,1], c=color[stepNum,:])
ax3.set_title('Upper surface velocity')
ax3.set(xlabel='X position (m)', ylabel='Y velocity (m/year)')

plt.show()

# End of file
