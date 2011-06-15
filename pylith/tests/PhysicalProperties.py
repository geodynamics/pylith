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

## @file pylith/tests/PhysicalProperties.py
##
## @brief Check physical properties output from PyLith.

import numpy

def check_properties(testcase, filename, mesh, properties):
  """
  Check properties.
  """
  data = testcase.reader.read(filename)
  
  # Check cells
  (ncells, ncorners) = data['cells'].shape
  testcase.assertEqual(mesh['ncells'], ncells)
  testcase.assertEqual(mesh['ncorners'], ncorners)

  # Check vertices
  (nvertices, spaceDim) = data['vertices'].shape
  testcase.assertEqual(mesh['nvertices'], nvertices)
  testcase.assertEqual(mesh['spaceDim'], spaceDim)

  # Check physical properties
  tolerance = 1.0e-6

  for name in properties.keys():
    propertyE = properties[name]
    property = data['cell_fields'][name]
    ratio = numpy.abs(1.0 - property[:]/propertyE[:,0])
    diff = numpy.abs(property[:] - propertyE[:,0])
    mask = propertyE[:,0] != 0.0
    okay = mask*(ratio < tolerance) + ~mask*(diff < tolerance)
    if numpy.sum(okay) != ncells:
      print "Error in values for physical property '%s'." % name
      print "Expected values:",propertyE
      print "Output values:",property
    testcase.assertEqual(numpy.sum(okay), ncells)

  return


# End of file
