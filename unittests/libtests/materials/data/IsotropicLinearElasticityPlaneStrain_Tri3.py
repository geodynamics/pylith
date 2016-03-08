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


#points = numpy.array([[-4.0, -4.0],
#                      [-4.0, +4.0],
#                      [+4.0, -4.0],
#                      [+4.0, +4.0]], dtype=numpy.float64)
points = numpy.array([[-0.5, 0.0],
                      [+0.5, 0.0]], dtype=numpy.float64)

exx = 0.1
eyy = 0.2
exy = 0.3

npts = points.shape[0]


density = 2500.0*numpy.ones((npts,))
vs = 3000.0*numpy.ones((npts,))
vp = 3**0.5*vs*numpy.ones((npts,))
modulus_mu = density*vs**2
modulus_lambda = density*vp**2 - 2.0*modulus_mu

# Create coordinate system for spatial database
from spatialdata.geocoords.CSCart import CSCart
cs = CSCart()
cs._configure()
cs.setSpaceDim(2)



# GenerateApp class
class GenerateApp(object):
  """
  Python application for generating spatial database for residual calculation.
  """


  def generateResidualDB(self):
    disp = numpy.zeros((npts, 2))
    vel = numpy.zeros(disp.shape)
    disp[:,0] = exx*points[:,0] + exy*points[:,1]
    disp[:,1] = exy*points[:,0] + eyy*points[:,1]
  
    disp_dot = numpy.zeros(disp.shape)
    vel_dot = numpy.zeros(vel.shape)

    # Create writer for spatial database file
    from spatialdata.spatialdb.SimpleIOAscii import SimpleIOAscii
    writer = SimpleIOAscii()
    writer.inventory.filename = "matinitialize.spatialdb"
    writer._configure()
    writer.write({'points': points,
                  'coordsys': cs,
                  'data_dim': 1,
                  'values': [{'name': "vs", 'units': "m/s", 'data': vs},
                             {'name': "vp", 'units': "m/s", 'data': vp},
                             {'name': "density", 'units': "kg/m**3", 'data': density},
                             {'name': "displacement_x", 'units': "m", 'data': disp[:,0]},
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
  app = GenerateApp()
  app.generateResidualDB()


# MAIN /////////////////////////////////////////////////////////////////
if __name__ == "__main__":
  generate()
  
# End of file 
