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

## @file tests/manual/3d/powerlaw/radial_bc.py

## @brief Python script to generate displacement BC for cylinder test.

import math
import numpy
import h5py
from spatialdata.geocoords.CSCart import CSCart
from spatialdata.spatialdb.SimpleIOAscii import createWriter

# ----------------------------------------------------------------------
# Filenames.
inFile = 'output/cylinder_pres_powerlaw_hex-viscomat_info.h5'
outFile = 'cylinder_disp_bc.spatialdb'

# Radial displacement value.
ur = -10.0

# Get coordinates.
h5 = h5py.File(inFile, 'r')
coords = h5['geometry/vertices'][:]
h5.close()

# Compute Cartesian displacements.
angs = numpy.arctan2(coords[:,1], coords[:,0])
ca = numpy.cos(angs)
sa = numpy.sin(angs)

ux = ur*ca
uy = ur*sa
uz = numpy.zeros_like(ux)

# Write spatial database.
cs = CSCart()
cs.inventory.spaceDim = 3
cs._configure()
data = {
    'points': coords,
    'coordsys': cs,
    'data_dim': 3,
    'values': [
    {'name': 'initial_amplitude_x',
     'units': 'm',
     'data': ux},
    {'name': 'initial_amplitude_y',
     'units': 'm',
     'data': uy},
    {'name': 'initial_amplitude_z',
     'units': 'm',
     'data': uz},
    ]}

io = createWriter(outFile)
io.write(data)

# End of file 
