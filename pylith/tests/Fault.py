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

## @file pylith/tests/Fault.py
##
## @brief Check fault output from PyLith.

import numpy

def check_vertex_fields(testcase, filename, mesh, fieldNames):
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

  # Check fault information
  tolerance = 1.0e-5

  from spatialdata.units.NondimElasticQuasistatic import NondimElasticQuasistatic
  normalizer = NondimElasticQuasistatic()
  normalizer._configure()

  for name in fieldNames:
    valuesE = testcase.calcFaultField(name, data['vertices'])
    values = data['vertex_fields'][name]

    (nverticesE, dimE) = valuesE.shape
    if 1 == dimE:
      values = values.reshape( (nvertices, dimE) )
    (nvertices, dim) = values.shape

    testcase.assertEqual(nverticesE, nvertices)
    testcase.assertEqual(dimE, dim)

    scale = 1.0
    if name == "traction_change" or name == "traction":
      scale *= normalizer.pressureScale().value

    for i in xrange(dim):
      ratio = numpy.abs(1.0 - values[:,i]/valuesE[:,i])
      diff = numpy.abs(values[:,i] - valuesE[:,i]) / scale
      mask = valuesE[:,i] != 0.0
      okay = mask*(ratio < tolerance) + ~mask*(diff < tolerance)
      if numpy.sum(okay) != nvertices:
        print "Error in component %d of field '%s'." % (i, name)
        print "Expected values:",valuesE
        print "Output values:",values
      testcase.assertEqual(numpy.sum(okay), nvertices)

  return


def check_data(testcase, filename, mesh, fieldNames):
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

  # Check fault information
  tolerance = 1.0e-5

  for name in fieldNames:
    valuesE = testcase.calcFaultInfo(name, data['vertices'])
    values = data['vertex_fields'][name]

    (nverticesE, dim) = valuesE.shape
    values = values.reshape( (nvertices, dim) )
    testcase.assertEqual(nverticesE, nvertices)

    for i in xrange(dim):
      ratio = numpy.abs(1.0 - values[:,i]/valuesE[:,i])
      diff = numpy.abs(values[:,i] - valuesE[:,i])
      mask = valuesE[:,i] != 0.0
      okay = mask*(ratio < tolerance) + ~mask*(diff < tolerance)
      if numpy.sum(okay) != nvertices:
        print "Error in component %d of field '%s'." % (i, name)
        print "Expected values:",valuesE
        print "Output values:",values
      testcase.assertEqual(numpy.sum(okay), nvertices)

  return


# End of file
