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
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file unittests/libtests/materials/data/IsotropicLinearElasticityPlaneStrain_Tri3.py

## @brief Python application for generating spatial database files for
## testing IsotropicLinearElasticityPlaneStrain via Method of
## Manufactured Solutions.

import numpy

from spatialdata.spatialdb.SimpleIOAscii import SimpleIOAscii
from spatialdata.geocoords.CSCart import CSCart

points = numpy.array([[-4.0e+3, -4.0e+3],
                      [-4.0e+3, +4.0e+3],
                      [+4.0e+3, -4.0e+3],
                      [+4.0e+3, +4.0e+3]], dtype=numpy.float64)

npts = points.shape[0]

density = 2500.0*numpy.ones((npts,))
vs = 3000.0*numpy.ones((npts,))
vp = 3**0.5*3000.0*numpy.ones((npts,))

# Create coordinate system for spatial database
cs = CSCart()
cs._configure()
cs.setSpaceDim(2)

# ----------------------------------------------------------------------
def generateAuxFields():
    writer = SimpleIOAscii()
    writer.inventory.filename = "IsotropicLinearElasticityPlaneStrain_UniStrain_aux.spatialdb"
    writer._configure()
    writer.write({'points': points,
                  'coordsys': cs,
                  'data_dim': 2,
                  'values': [{'name': "vs", 'units': "m/s", 'data': vs},
                             {'name': "vp", 'units': "m/s", 'data': vp},
                             {'name': "density", 'units': "kg/m**3", 'data': density},
                           ]})
    return


# ----------------------------------------------------------------------
def generateSolution():
    exx0 = 0.1
    eyy0 = 0.25
    exy0 = 0.3
    
    exxR = 0.4
    eyyR = 0.7
    exyR = -0.2
    
    t = 0.1
    exx = exx0 + exxR*t
    eyy = eyy0 + eyyR*t
    exy = exy0 + exyR*t
    disp = numpy.zeros((npts, 2))
    disp[:,0] = exx*points[:,0] + exy*points[:,1]
    disp[:,1] = exy*points[:,0] + eyy*points[:,1]

    vel = numpy.zeros(disp.shape)
    vel[:,0] = exxR*points[:,0] + exyR*points[:,1]
    vel[:,1] = exyR*points[:,0] + eyyR*points[:,1]
    vel *= 0.0 # :TEMPORARY: :KLUDGE:

    disp_dot = vel

    vel_dot = numpy.zeros(vel.shape)

    # Create writer for spatial database file
    writer = SimpleIOAscii()
    writer.inventory.filename = "IsotropicLinearElasticityPlaneStrain_UniStrain_soln.spatialdb"
    writer._configure()
    writer.write({'points': points,
                  'coordsys': cs,
                  'data_dim': 2,
                  'values': [{'name': "displacement_x", 'units': "m", 'data': disp[:,0]},
                             {'name': "displacement_y", 'units': "m", 'data': disp[:,1]},
                             {'name': "velocity_x", 'units': "m/s", 'data': vel[:,0]},
                             {'name': "velocity_y", 'units': "m/s", 'data': vel[:,1]},
                             {'name': "displacement_dot_x", 'units': "m/s", 'data': disp_dot[:,0]},
                             {'name': "displacement_dot_y", 'units': "m/s", 'data': disp_dot[:,1]},
                             {'name': "velocity_dot_x", 'units': "m/s**2", 'data': vel_dot[:,0]},
                             {'name': "velocity_dot_y", 'units': "m/s**2", 'data': vel_dot[:,1]},
                           ]})
    
    return


# ----------------------------------------------------------------------
def generatePerturbation():
    dt = 0.05

    x = numpy.arange(-4.0e+3, 4.01e+3, 0.5e+3, dtype=numpy.float64)
    y = numpy.arange(-4.0e+3, 4.01e+3, 0.5e+3, dtype=numpy.float64)
    numX = x.shape[0]
    numY = y.shape[0]
    points = numpy.zeros((numX*numY,2), dtype=numpy.float64)
    npts = numX*numY
    for iY in xrange(numY):
        points[iY*numX:(iY+1)*numX,0] = x
        points[iY*numX:(iY+1)*numX,1] = y[iY]

    disp = 0.1*(numpy.random.rand(npts,2)-0.5)

    vel = 1.0/dt * disp
    vel *= 0.0 # :TEMPORARY: :KLUDGE:
  
    disp_dot = vel

    vel_dot = numpy.zeros(vel.shape)

    # Create writer for spatial database file
    writer = SimpleIOAscii()
    writer.inventory.filename = "IsotropicLinearElasticityPlaneStrain_UniStrain_pert.spatialdb"
    writer._configure()
    writer.write({'points': points,
                  'coordsys': cs,
                  'data_dim': 2,
                  'values': [{'name': "displacement_x", 'units': "m", 'data': disp[:,0]},
                             {'name': "displacement_y", 'units': "m", 'data': disp[:,1]},
                             {'name': "velocity_x", 'units': "m/s", 'data': vel[:,0]},
                             {'name': "velocity_y", 'units': "m/s", 'data': vel[:,1]},
                             {'name': "displacement_dot_x", 'units': "m/s", 'data': disp_dot[:,0]},
                             {'name': "displacement_dot_y", 'units': "m/s", 'data': disp_dot[:,1]},
                             {'name': "velocity_dot_x", 'units': "m/s**2", 'data': vel_dot[:,0]},
                             {'name': "velocity_dot_y", 'units': "m/s**2", 'data': vel_dot[:,1]},
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
