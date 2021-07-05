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

## @file tests/manual/3d/powerlaw/CylinderTest.py

## @brief Python script to test power-law implementation for steady-state solution
##        of a pressurized cylinder.

import math
import numpy
import h5py
import netCDF4
from pylith.meshio.Xdmf import Xdmf
import pdb
pdb.set_trace()

from pylith.testing.FullTestApp import FullTestCase

from cylinderpres_soln import AnalyticalSoln
import cylinderpres_gendb

# ----------------------------------------------------------------------
# Filenames.
h5Prefix = 'output/cylinder_pres_powerlaw_'
meshPrefix = 'meshes/mesh_cylinder_'
spatialdbPrefix = 'mat_powerlaw_cylinder_pres_refstate_'
spatialdbSuffix = '.spatialdb'
cellTypes = ['tet', 'hex']
meshSuffix = '.exo'
dispSuffix = '-domain.h5'
stressSuffix = '-viscomat.h5'
configFiles = [['cylinder_pres_powerlaw_refstate.cfg', 'cylinder_pres_powerlaw_tet_refstate.cfg'],
               ['cylinder_pres_powerlaw_refstate.cfg', 'cylinder_pres_powerlaw_hex_refstate.cfg']]
numTests = len(cellTypes)

# Tolerances.
velAbsTol = 1.0e-5
stressAbsTol = 1.0e-5
velRelTol = 1.0e-5
stressRelTol = 1.0e-5

# ----------------------------------------------------------------------
# Loop over cell types.
for testNum in range(numTests):
    appName = cellTypes[testNum]
    meshFile = meshPrefix + cellTypes[testNum] + meshSuffix
    spatialdbFile = spatialdbPrefix + cellTypes[testNum] + spatialdbSuffix
    stressFile = h5Prefix + cellTypes[testNum] + stressSuffix
    dispFile = h5Prefix + cellTypes[testNum] + dispSuffix

    exodus = netCDF4.Dataset(meshFile, 'r')
    try:
        x = exodus.variables['coordx'][:]
        y = exodus.variables['coordy'][:]
        z = exodus.variables['coordz'][:]
        vertices = numpy.column_stack((x,y,z))
    except:
        vertices = exodus.variables['coord'][:].transpose()

    cylinderpres_gendb.generate_refstate_db(vertices, spatialdbFile)
    FullTestCase.run_pylith(appName, configFiles[testNum])

    # Read stress info from HDF5 file.
    h5Stress = h5py.File(stressFile, 'r')
    coords = h5Stress['geometry/vertices'][:]
    connect = numpy.array(h5Stress['topology/cells'][:], dtype=numpy.int64)
    cellCoords = coords[connect, :]
    cellCenters = numpy.mean(cellCoords, axis=1)
    stressNum = h5Stress['cell_fields/cauchy_stress'][-1,:,:]
    h5Stress.close()

    # Read displacement info from HDF5 file.
    h5Disp = h5py.File(dispFile, 'r')
    time = h5Disp['time'][:,0,0]
    dt = time[-1] - time[-2]
    dispTnMinus1 = h5Disp['vertex_fields/displacement'][-2,:,:]
    dispTn = h5Disp['vertex_fields/displacement'][-1,:,:]
    dispIncrNum = dispTn - dispTnMinus1
    h5Disp.close()

    # Compute analytical solution.
    soln = AnalyticalSoln()
    dispIncrAnl = soln.displacement_incr(coords, dt)
    stressAnl = soln.stress(coords)[0,:,:]

    # Compute difference.
    dispDiff = dispIncrAnl - dispIncrNum
    stressDiff = stressAnl - stressNum

# End of file 
