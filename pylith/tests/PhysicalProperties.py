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
# Copyright (c) 2010-2016 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file pylith/tests/PhysicalProperties.py
##
## @brief Check physical properties output from PyLith.

import numpy
import h5py

def check_properties(testcase, filename, mesh, properties):
  """
  Check properties.
  """
  h5 = h5py.File(filename, "r", driver="sec2")
  
  # Check cells
  cells = h5['topology/cells'][:]
  (ncells, ncorners) = cells.shape
  testcase.assertEqual(mesh['ncells'], ncells)
  testcase.assertEqual(mesh['ncorners'], ncorners)

  # Check vertices
  vertices = h5['geometry/vertices'][:]
  (nvertices, spaceDim) = vertices.shape
  testcase.assertEqual(mesh['nvertices'], nvertices)
  testcase.assertEqual(mesh['spaceDim'], spaceDim)

  # Check physical properties
  tolerance = 1.0e-6

  for name in properties.keys():
    istep = 0
    icomp = 0

    propertyE = properties[name][istep,:,icomp]
    property = h5['cell_fields/%s' % name][istep,:,icomp]

    okay = numpy.zeros((ncells,), dtype=numpy.bool)

    maskR = numpy.abs(propertyE > tolerance)
    ratio = numpy.abs(1.0 - property[maskR]/propertyE[maskR])
    if len(ratio) > 0:
      okay[maskR] = ratio < tolerance

    maskD = ~maskR
    diff = numpy.abs(property[maskD] - propertyE[maskD])
    if len(diff) > 0:
      okay[maskD] = diff < tolerance

    if numpy.sum(okay) != ncells:
      print "Error in values for physical property '%s'." % name
      print "Expected values:",propertyE
      print "Output values:",property
      print "Expected values (not okay): ",propertyE[~okay]
      print "Computed values (not okay): ",property[~okay]
      print "Coordinates (not okay): ",vertices[~okay,:]
      h5.close()
    testcase.assertEqual(ncells, numpy.sum(okay))

  h5.close()
  return


# End of file
