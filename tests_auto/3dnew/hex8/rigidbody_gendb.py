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
# Copyright (c) 2010-2013 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file tests/3dnew/hex8/rigidbody_gendb.py
##
## @brief Python script to generate spatial database with displacement
## boundary conditions for the rigid body motion test.

import numpy

class GenerateDB(object):
  """
  Python object to generate spatial database with displacement
  boundary conditions for the rigid body motion test.
  """
  
  def __init__(self):
    """
    Constructor.
    """
    return
  
  
  def run(self):
    """
    Generate the database.
    """
    # Domain
    x = numpy.arange(-4000.0, 4000.1, 1000.0)
    y = numpy.arange(-4000.0, 4000.1, 1000.0)
    z = numpy.arange(-6000.0, 0000.1, 1000.0)
    nptsX = x.shape[0]
    nptsY = y.shape[0]
    nptsZ = z.shape[0]

    xx = x * numpy.ones( (nptsY*nptsZ, 1), dtype=numpy.float64)
    yy0 = y * numpy.ones( (nptsZ, 1), dtype=numpy.float64)
    yy = yy0.ravel() * numpy.ones( (nptsX, 1), dtype=numpy.float64)
    zz = z * numpy.ones( (nptsX*nptsY, 1), dtype=numpy.float64)
    xyz = numpy.zeros( (nptsX*nptsY*nptsZ, 3), dtype=numpy.float64)
    xyz[:,0] = numpy.ravel(xx)
    xyz[:,1] = numpy.ravel(numpy.transpose(yy))
    xyz[:,2] = numpy.ravel(numpy.transpose(zz))

    from rigidbody_soln import AnalyticalSoln
    soln = AnalyticalSoln()
    disp = soln.displacement(xyz)

    from spatialdata.geocoords.CSCart import CSCart
    cs = CSCart()
    cs.inventory.spaceDim = 3
    cs._configure()
    data = {'locs': xyz,
            'coordsys': cs,
            'data_dim': 2,
            'values': [{'name': "displacement-x",
                        'units': "m",
                        'data': numpy.ravel(disp[:,0])},
                       {'name': "displacement-y",
                        'units': "m",
                        'data': numpy.ravel(disp[:,1])},
                       {'name': "displacement-z",
                        'units': "m",
                        'data': numpy.ravel(disp[:,2])}]}

    from spatialdata.spatialdb.SimpleIOAscii import SimpleIOAscii
    io = SimpleIOAscii()
    io.inventory.filename = "rigidbody_disp.spatialdb"
    io._configure()
    io.write(data)
    return

# ======================================================================
if __name__ == "__main__":
  app = GenerateDB()
  app.run()


# End of file 
