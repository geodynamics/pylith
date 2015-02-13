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

## @file pylith/tests/Fault.py
##
## @brief Check fault output from PyLith.

import numpy
import h5py
from spatialdata.units.NondimElasticQuasistatic import NondimElasticQuasistatic

def check_vertex_fields(testcase, filename, mesh, fieldNames):
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

  # Check fault information
  tolerance = 1.0e-5

  normalizer = NondimElasticQuasistatic()
  normalizer._configure()

  for name in fieldNames:
    valuesE = testcase.calcFaultField(name, vertices)
    values = h5['vertex_fields/%s' % name][:]

    (nstepsE, nverticesE, dimE) = valuesE.shape
    (nsteps, nvertices, dim) = values.shape

    testcase.assertEqual(nstepsE, nsteps)
    testcase.assertEqual(nverticesE, nvertices)
    testcase.assertEqual(dimE, dim)

    scale = 1.0
    if name == "traction_change" or name == "traction":
      scale *= normalizer.pressureScale().value

    for istep in xrange(nsteps):
      for idim in xrange(dim):
        ratio = numpy.abs(1.0 - values[istep,:,idim]/valuesE[istep,:,idim])
        diff = numpy.abs(values[istep,:,idim] - valuesE[istep,:,idim]) / scale
        mask = valuesE[istep,:,idim] != 0.0
        okay = mask*(ratio < tolerance) + ~mask*(diff < tolerance)
        if numpy.sum(okay) != nvertices:
          print "Error in component %d of field '%s' for timestep %d." % (idim, name, istep)
          print "Expected values:",valuesE
          print "Output values:",values
          print "Coordinates: ",vertices
        testcase.assertEqual(nvertices, numpy.sum(okay))

  h5.close()
  return


def check_data(testcase, filename, mesh, fieldNames):
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

  # Check fault information
  tolerance = 1.0e-5

  for name in fieldNames:
    valuesE = testcase.calcFaultInfo(name, data['vertices'])
    values = h5['vertex_fields/%s' % name][:]

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
      testcase.assertEqual(nvertices, numpy.sum(okay))
  h5.close()
  return


# End of file
