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

## @file tests/2d/tri3/sheardisp_gendb.py
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
    npts = x.shape[0]

    xx = x * numpy.ones( (npts, 1), dtype=numpy.float64)
    yy = y * numpy.ones( (npts, 1), dtype=numpy.float64)
    xy = numpy.zeros( (npts**2, 2), dtype=numpy.float64)
    xy[:,0] = numpy.ravel(xx)
    xy[:,1] = numpy.ravel(numpy.transpose(yy))

    from sheardisp_soln import AnalyticalSoln
    soln = AnalyticalSoln()
    disp = soln.displacement(xy)

    from spatialdata.geocoords.CSCart import CSCart
    cs = CSCart()
    cs.inventory.spaceDim = 2
    cs._configure()
    data = {'points': xy,
            'coordsys': cs,
            'data_dim': 2,
            'values': [{'name': "displacement-x",
                        'units': "m",
                        'data': numpy.ravel(disp[0,:,0])},
                       {'name': "displacement-y",
                        'units': "m",
                        'data': numpy.ravel(disp[0,:,1])}]}

    from spatialdata.spatialdb.SimpleIOAscii import SimpleIOAscii
    io = SimpleIOAscii()
    io.inventory.filename = "shear_disp.spatialdb"
    io._configure()
    io.write(data)
    return

# ======================================================================
if __name__ == "__main__":
  app = GenerateDB()
  app.run()


# End of file 
