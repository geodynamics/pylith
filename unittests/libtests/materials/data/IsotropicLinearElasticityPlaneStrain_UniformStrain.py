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

## @file unittests/libtests/materials/data/IsotropicLinearElasticityPlaneStrain_UniformStrain.py

## @brief Python application for generating spatial database files for
## testing IsotropicLinearElasticityPlaneStrain via Method of
## Manufactured Solutions for uniform strain.

# ----------------------------------------------------------------------
# Domain
XLIM = (-4.0e+3, +4.0e+3)
YLIM = XLIM
DX = XLIM[1]-XLIM[0]

# Material properties
DENSITY = 2500.0
VS = 3000.0
VP = 2**0.5 * VS

# ----------------------------------------------------------------------
import numpy

from spatialdata.spatialdb.SimpleGridAscii import SimpleGridAscii
from spatialdata.geocoords.CSCart import CSCart

x = numpy.arange(XLIM[0], XLIM[1]+0.1*DX, DX, dtype=numpy.float64)
y = numpy.arange(YLIM[0], YLIM[1]+0.1*DX, DX, dtype=numpy.float64)
xgrid, ygrid = numpy.meshgrid(x, y)
points = numpy.vstack((xgrid.ravel(), ygrid.ravel())).transpose()
npts = points.shape[0]

density = DENSITY*numpy.ones((npts,))
vs = VS*numpy.ones((npts,))
vp = VP*numpy.ones((npts,))

# Create coordinate system for spatial database
cs = CSCart()
cs._configure()
cs.setSpaceDim(2)

# ----------------------------------------------------------------------
def generateAuxFields():
    writer = SimpleGridAscii()
    writer.inventory.filename = "IsotropicLinearElasticityPlaneStrain_UniformStrain_aux.spatialdb"
    writer._configure()
    writer.write({'points': points,
                  'x': x,
                  'y': y,
                  'coordsys': cs,
                  'data_dim': 2,
                  'values': [{'name': "vs", 'units': "m/s", 'data': vs},
                             {'name': "vp", 'units': "m/s", 'data': vp},
                             {'name': "density", 'units': "kg/m**3", 'data': density},
                           ]})
    return


# ----------------------------------------------------------------------
def generateSolution():
    EXX = 0.1
    EYY = 0.25
    EXY = 0.3
    
    disp = numpy.zeros((npts, 2))
    disp[:,0] = EXX*points[:,0] + EXY*points[:,1]
    disp[:,1] = EXY*points[:,0] + EYY*points[:,1]

    disp_dot = 0*disp

    # Create writer for spatial database file
    writer = SimpleGridAscii()
    writer.inventory.filename = "IsotropicLinearElasticityPlaneStrain_UniformStrain_soln.spatialdb"
    writer._configure()
    writer.write({'points': points,
                  'x': x,
                  'y': y,
                  'coordsys': cs,
                  'data_dim': 2,
                  'values': [{'name': "displacement_x", 'units': "m", 'data': disp[:,0]},
                             {'name': "displacement_y", 'units': "m", 'data': disp[:,1]},
                             {'name': "displacement_dot_x", 'units': "m/s", 'data': disp_dot[:,0]},
                             {'name': "displacement_dot_y", 'units': "m/s", 'data': disp_dot[:,1]},
                           ]})
    
    return


# ----------------------------------------------------------------------
def generatePerturbation():
    PERT_DX = 500.0
    PERT_AMPLITUDE = 1.0e-2
    
    x = numpy.arange(XLIM[0], XLIM[1]+0.1*PERT_DX, PERT_DX, dtype=numpy.float64)
    y = numpy.arange(YLIM[0], YLIM[1]+0.1*PERT_DX, PERT_DX, dtype=numpy.float64)
    xgrid, ygrid = numpy.meshgrid(x, y)
    points = numpy.vstack((xgrid.ravel(), ygrid.ravel())).transpose()
    npts = points.shape[0]

    disp = PERT_AMPLITUDE*(numpy.random.rand(npts,2)-0.5)
    disp_dot = 0*disp

    # Create writer for spatial database file
    writer = SimpleGridAscii()
    writer.inventory.filename = "IsotropicLinearElasticityPlaneStrain_UniformStrain_pert.spatialdb"
    writer._configure()
    writer.write({'points': points,
                  'x': x,
                  'y': y,
                  'coordsys': cs,
                  'data_dim': 2,
                  'values': [{'name': "displacement_x", 'units': "m", 'data': disp[:,0]},
                             {'name': "displacement_y", 'units': "m", 'data': disp[:,1]},
                             {'name': "displacement_dot_x", 'units': "m/s", 'data': disp_dot[:,0]},
                             {'name': "displacement_dot_y", 'units': "m/s", 'data': disp_dot[:,1]},
                           ]})
    
    return


# ======================================================================
def generate():
    generateAuxFields()
    generateSolution()
    generatePerturbation()


# MAIN /////////////////////////////////////////////////////////////////
if __name__ == "__main__":
    generate()
  
# End of file 
