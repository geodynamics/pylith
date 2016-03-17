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

points = numpy.array([[-4.0, -4.0],
                      [-4.0, +4.0],
                      [+4.0, -4.0],
                      [+4.0, +4.0]], dtype=numpy.float64)

npts = points.shape[0]

density = 2500.0*numpy.ones((npts,))
vs = 3000.0*numpy.ones((npts,))
vp = 3**0.5*vs*numpy.ones((npts,))
modulus_mu = density*vs**2
modulus_lambda = density*vp**2 - 2.0*modulus_mu

# Create coordinate system for spatial database
cs = CSCart()
cs._configure()
cs.setSpaceDim(2)

# Auxiliary Fields
class AuxFields(object):
  """
  Python class for generation spatial database with auxiliary fields.
  """

  @staticmethod
  def generate():
    writer = SimpleIOAscii()
    writer.inventory.filename = "IsotropicLinearElasticityPlaneStrain_UniStrain_aux.spatialdb"
    writer._configure()
    writer.write({'points': points,
                  'coordsys': cs,
                  'data_dim': 2,
                  'values': [{'name': "vs", 'units': "m/s", 'data': vs},
                             {'name': "vp", 'units': "m/s", 'data': vp},
                             {'name': "density", 'units': "kg/m**3", 'data': density},
                           ]}
               )


# Solution @ t1
class Solution1(object):
  """
  Python class for generation spatial database with solution at t1.
  """
  exx = 0.1
  eyy = 0.2
  exy = 0.3
  t = 1.0

  @staticmethod
  def generate():

    disp = numpy.zeros((npts, 2))
    vel = numpy.zeros(disp.shape)
    disp[:,0] = Solution1.exx*points[:,0] + Solution1.exy*points[:,1]
    disp[:,1] = Solution1.exy*points[:,0] + Solution1.eyy*points[:,1]
  
    disp_dot = numpy.zeros(disp.shape)
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
                           ]}
               )
    
    return


# Solution @ t2
class Solution2(object):
  """
  Python class for generation spatial database with test solution
  field at t2 (not a solution to the problem).
  """

  @staticmethod
  def generate():
    import numpy
    x = numpy.arange(-4.0, 4.01, 0.1, dtype=numpy.float64)
    y = numpy.arange(-4.0, 4.01, 0.1, dtype=numpy.float64)
    numX = x.shape[0]
    numY = y.shape[0]
    points = numpy.zeros((numX*numY,2), dtype=numpy.float64)
    npts = numX*numY
    for iY in xrange(numY):
      points[iY*numX:(iY+1)*numX,0] = x
      points[iY*numX:(iY+1)*numX,1] = y[iY]

    import numpy.random
    disp = 0.1*(numpy.random.rand(npts,2)-0.5)
    vel = numpy.zeros(disp.shape)
    print points
  
    disp_dot = numpy.zeros(disp.shape)
    vel_dot = numpy.zeros(vel.shape)

    # Create writer for spatial database file
    writer = SimpleIOAscii()
    writer.inventory.filename = "IsotropicLinearElasticityPlaneStrain_Random_soln.spatialdb"
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
                           ]}
               )
    
    return


# ======================================================================
def generate():
  AuxFields.generate()
  Solution1.generate()
  Solution2.generate()


# MAIN /////////////////////////////////////////////////////////////////
if __name__ == "__main__":
  generate()
  
# End of file 
