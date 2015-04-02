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

## @file tests/2d/slipdir/genspatialdb.py
##
## @brief Python script to generate spatial database with displacement
## boundary conditions.

import numpy

# ======================================================================
class GenerateDB(object):
  """
  Python object to generate spatial database with displacement
  boundary conditions.
  """
  
  def __init__(self):
    """
    Constructor.
    """
    self.soln = None
    self.filename = None
    return
  
  
  def run(self):
    """
    Generate the database.
    """
    # Domain
    x = numpy.arange(-4000.0, 4000.1, 500.0)
    y = numpy.arange(-4000.0, 4000.1, 500.0)
    npts = x.shape[0]

    xx = x * numpy.ones( (npts, 1), dtype=numpy.float64)
    yy = y * numpy.ones( (npts, 1), dtype=numpy.float64)
    xy = numpy.zeros( (npts**2, 2), dtype=numpy.float64)
    xy[:,0] = numpy.ravel(xx)
    xy[:,1] = numpy.ravel(numpy.transpose(yy))

    disp = self.soln.displacement(xy)

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
    io.inventory.filename = self.filename
    io._configure()
    io.write(data)
    return

# ======================================================================
class GenDBFaultX(GenerateDB):
  """
  Python object to generate spatial database with displacement
  boundary conditions for the faultx test.
  """
  
  def __init__(self):
    """
    Constructor.
    """
    from solution import SolnFaultX
    self.soln = SolnFaultX()
    self.filename = "faultx_disp.spatialdb"
    return
  
  
# ======================================================================
class GenDBFaultY(GenerateDB):
  """
  Python object to generate spatial database with displacement
  boundary conditions for the faulty test.
  """
  
  def __init__(self):
    """
    Constructor.
    """
    from solution import SolnFaultY
    self.soln = SolnFaultY()
    self.filename = "faulty_disp.spatialdb"
    return
  
  
# ======================================================================
class GenDBFaultXYP(GenerateDB):
  """
  Python object to generate spatial database with displacement
  boundary conditions for the faultxyp test.
  """
  
  def __init__(self):
    """
    Constructor.
    """
    from solution import SolnFaultXYP
    self.soln = SolnFaultXYP()
    self.filename = "faultxyp_disp.spatialdb"
    return
  
  
# ======================================================================
class GenDBFaultXYN(GenerateDB):
  """
  Python object to generate spatial database with displacement
  boundary conditions for the faultxyn test.
  """
  
  def __init__(self):
    """
    Constructor.
    """
    from solution import SolnFaultXYN
    self.soln = SolnFaultXYN()
    self.filename = "faultxyn_disp.spatialdb"
    return
  
  
# End of file 
