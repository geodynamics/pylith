#!/usr/bin/env nemesis
# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================
# PREREQUISITES: numpy, netCDF4, spatialdata (PyLith)

import netCDF4
import numpy
import sys
filenameExodus = "mesh_cellsize.exo"
filenameDB = "matprops.spatialdb"
minPeriod = 10.0

# ======================================================================

# ----------------------------------------------------------------------
# Cell size based on minimum wavelength with Vs from spatial database


def getCellSizeDB(points):

    from spatialdata.geocoords.CSCart import CSCart
    from spatialdata.spatialdb.SimpleDB import SimpleDB
    from spatialdata.spatialdb.SimpleIOAscii import SimpleIOAscii

    # Coordinate system for mesh (must match coordsys in ExodusII file)
    cs = CSCart()
    cs._configure()

    # Spatial database with physical properties (Vs)
    dbIO = SimpleIOAscii()
    dbIO.inventory.filename = filenameDB
    dbIO._configure()
    db = SimpleDB()
    db.inventory.iohandler = dbIO
    db.inventory.description = "Physical properties"
    db.inventory.queryType = "linear"
    db._configure()

    (npoints, spacedim) = points.shape

    # Query database
    db.open()
    db.setQueryValues(["vs"])
    data = numpy.zeros((npoints, 1), dtype=numpy.float64)
    err = numpy.zeros((npoints,), dtype=numpy.int32)
    db.multiquery(data, err, points, cs)
    db.close()

    vs = data[:, 0]
    cellSize = minPeriod * vs / 10.0
    return cellSize


# ----------------------------------------------------------------------
# Cell size based on analytical function of vertex coordinates.
def getCellSizeFn(points):
    """Cell size is based on distance from a target and grows at a
    geometric rate.
    """
    # Coordinates of target
    target = (5.0e+3, -10.0e+3, -10.0e+3)

    # Compute distance from target
    dist = ((points[:, 0] - target[0])**2 +
            (points[:, 1] - target[1])**2 +
            (points[:, 2] - target[2])**2)**0.5
    bias_factor = 1.05  # Geometric rate
    dxStart = 1.0e+3  # Discretization size at target
    npts = numpy.ceil(numpy.log(1 - dist / dxStart * (1 - bias_factor)) / numpy.log(bias_factor))
    cellSize = dxStart * bias_factor**npts
    return cellSize


# ----------------------------------------------------------------------
# Get coordinates of points from ExodusII file.
exodus = netCDF4.Dataset(filenameExodus, 'a')
try:
    x = exodus.variables['coordx'][:]
    y = exodus.variables['coordy'][:]
    z = exodus.variables['coordz'][:]
    points = numpy.column_stack((x,y,z))
except:
    points = exodus.variables['coord'][:].transpose()

cellSizeDB = getCellSizeDB(points)
cellSizeFn = getCellSizeFn(points)

# Add cell size info to ExodusII file
if not 'num_nod_var' in exodus.dimensions.keys():
    exodus.createDimension('num_nod_var', 2)

    name_nod_var = exodus.createVariable('name_nod_var', 'S1',
                                         ('num_nod_var', 'len_string',))
    name_nod_var[0,:] = netCDF4.stringtoarr("cell_size_db", 33)
    name_nod_var[1,:] = netCDF4.stringtoarr("cell_size_fn", 33)

    vals_nod_var1 = exodus.createVariable('vals_nod_var1', numpy.float64,
                                          ('time_step', 'num_nodes',))
    vals_nod_var2 = exodus.createVariable('vals_nod_var2', numpy.float64,
                                          ('time_step', 'num_nodes',))


time_whole = exodus.variables['time_whole']
time_whole[0] = 0.0
vals_nod_var1 = exodus.variables['vals_nod_var1']
vals_nod_var2 = exodus.variables['vals_nod_var2']
vals_nod_var1[0, :] = cellSizeDB.transpose()
vals_nod_var2[0, :] = cellSizeFn.transpose()

exodus.close()


# End of file
