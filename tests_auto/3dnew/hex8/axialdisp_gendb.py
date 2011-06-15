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
# Copyright (c) 2010-2011 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file tests/3d/hex8/axialdisp_gendb.py
##
## @brief Python script to generate spatial database with displacement
## boundary conditions for the axial displacement test.

import numpy

class GenerateDB(object):
  """
  Python object to generate spatial database with displacement
  boundary conditions for the axial displacement test.
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
    z = numpy.arange(-6000.0,    0.1, 1000.0)
    nxpts = x.shape[0]
    nypts = y.shape[0]
    nzpts = z.shape[0]

    xx = numpy.tile(x, (1,nypts*nzpts))
    yy = numpy.tile(y, (nxpts,nzpts)).transpose()
    zz = numpy.tile(z, (nxpts*nypts,1)).transpose()

    xyz = numpy.zeros( (nxpts*nypts*nzpts, 3), dtype=numpy.float64)
    xyz[:,0] = numpy.ravel(xx)
    xyz[:,1] = numpy.ravel(yy)
    xyz[:,2] = numpy.ravel(zz)

    from axialdisp_soln import AnalyticalSoln
    soln = AnalyticalSoln()
    disp = soln.displacement(xyz)

    from spatialdata.geocoords.CSCart import CSCart
    cs = CSCart()
    cs.inventory.spaceDim = 3
    cs._configure()
    data = {'points': xyz,
            'coordsys': cs,
            'data_dim': 3,
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
    io.inventory.filename = "axial_disp.spatialdb"
    io._configure()
    io.write(data)
    return


# ======================================================================
if __name__ == "__main__":
  app = GenerateDB()
  app.run()


# End of file 
