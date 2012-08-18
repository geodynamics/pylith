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
# Copyright (c) 2010-2012 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @brief Generate HDF5 data files for eqinfo tests.

import h5py
import numpy

# ----------------------------------------------------------------------
class TestData(object):
  """
  Abstract base class for test data.
  """
  
  def __init__(self):
    self.vertices = None
    self.cells = None
    self.slip = None
    self.time = None
    self.filename = None
    return

  def write(self):
    h5 = h5py.File(self.filename, "w", driver='sec2')
    h5.create_dataset('geometry/vertices', data=self.vertices)
    h5.create_dataset('topology/cells', data=self.cells)
    h5.create_dataset('time', data=self.time)
    slip = h5.create_dataset('vertex_fields/slip', data=self.slip)
    slip.attrs['vector_field_type'] = 'vector'
    h5.close()
    return


# ----------------------------------------------------------------------
class TestDataLineA(TestData):
  """
  Test data for fault as line.
  """
  
  def __init__(self):
    TestData.__init__(self)
    self.vertices = numpy.array([[-1.0, 0.0],
                                 [ 0.0, 0.0],
                                 [ 0.0, +1.5]], 
                                dtype=numpy.float64)
    self.cells = numpy.array([[0, 1],
                              [1, 2]], 
                             dtype=numpy.int32)
    self.slip = numpy.array([
        # t=0
        [[0.4, 0.0],
         [1.0, 0.0],
         [0.8, 0.0]],
        # t=1
        [[0.0, 0.0],
         [0.0, 0.0],
         [0.8, 0.0]]],
                            dtype=numpy.float64)
    self.time = numpy.array([0.0, 1.0])
    self.filename = "line_one.h5"


# ======================================================================
data = TestDataLineA()
data.write()

#data = TestDataTri()
#data.write()

#data = TestDataQuad()
#data.write()

# End of file 
