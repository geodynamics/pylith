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

## @file pylith/tests/StateVariables.py
##
## @brief Check state variables output from PyLith.

import numpy

def check_state_variables(testcase, filename, mesh, stateVarNames):
  """
  Check state variables.
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

  from spatialdata.units.NondimElasticQuasistatic import NondimElasticQuasistatic
  normalizer = NondimElasticQuasistatic()
  normalizer._configure()

  # Check state variables
  tolerance = 1.0e-6

  for name in stateVarNames:
    valuesE = testcase.calcStateVar(name, data['vertices'], data['cells'])
    values = data['cell_fields'][name]

    (ncellsE, dimE) = valuesE.shape
    if 1 == dimE:
      values = values.reshape( (ncells, dimE) )
    (ncells, dim) = values.shape

    testcase.assertEqual(ncellsE, ncells)
    testcase.assertEqual(dimE, dim)

    scale = 1.0
    if name == "stress":
      scale *= normalizer.pressureScale().value

    for i in xrange(dim):
      ratio = numpy.abs(1.0 - values[:,i]/valuesE[:,i])
      diff = numpy.abs(values[:,i] - valuesE[:,i]) / scale
      mask = valuesE[:,i] != 0.0
      okay = mask*(ratio < tolerance) + ~mask*(diff < tolerance)
      if numpy.sum(okay) != ncells:
        print "Error in component %d of state variable '%s'." % (i, name)
        print "Expected values:",valuesE
        print "Output values:",values
      testcase.assertEqual(numpy.sum(okay), ncells)
    
  return


# End of file
