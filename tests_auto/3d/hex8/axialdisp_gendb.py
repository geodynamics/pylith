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

## @file tests/2d/tri3/axialdisp_gendb.py
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
    from axialdisp_soln import AnalyticalSoln
    soln = AnalyticalSoln()

    from spatialdata.geocoords.CSCart import CSCart
    cs = CSCart()
    cs.inventory.spaceDim = 3
    cs._configure()

    from spatialdata.spatialdb.SimpleIOAscii import SimpleIOAscii
    io = SimpleIOAscii()

    for component in ["x","y","z"]:
      if component == "x":
        xyz = numpy.array([[-40.0e+3, 0.0, 0.0],
                           [+40.0e+3, 0.0, 0.0]], dtype=numpy.float64)
        ic = 0
      elif component == "y":
        xyz = numpy.array([[0.0, -40.0e+3, 0.0],
                           [0.0, +40.0e+3, 0.0]], dtype=numpy.float64)
        ic = 1
      elif component == "z":
        xyz = numpy.array([[0.0, 0.0, -40.0e+3],
                           [0.0, 0.0,   0.0]], dtype=numpy.float64)
        ic = 2
      disp = soln.displacement(xyz)

      io.inventory.filename = "axial_disp%s.spatialdb" % component
      io._configure()
      data = {'points': xyz,
              'coordsys': cs,
              'data_dim': 1,
              'values': [{'name': "displacement-%s" % component,
                          'units': "m",
                          'data': numpy.ravel(disp[0,:,ic])},
                         ]}
      io.write(data)
    
    return

# ======================================================================
if __name__ == "__main__":
  app = GenerateDB()
  app.run()


# End of file 
