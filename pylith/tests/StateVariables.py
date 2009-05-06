#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# <LicenseText>
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

  # Check state variables
  tolerance = 1.0e-5

  for name in stateVarNames:
    valuesE = testcase.calcStateVar(name, data['vertices'], data['cells'])
    values = data['cell_fields'][name]

    (ncellsE, dim) = valuesE.shape
    values = values.reshape( (ncells, dim) )
    testcase.assertEqual(ncellsE, ncells)

    for i in xrange(dim):
      ratio = numpy.abs(1.0 - values[:,i]/valuesE[:,i])
      diff = numpy.abs(values[:,i] - valuesE[:,i])
      mask = valuesE[:,i] != 0.0
      okay = mask*(ratio < tolerance) + ~mask*(diff < tolerance)
      if numpy.sum(okay) != ncells:
        print "Error in component %d of state variable '%s'." % (i, name)
        print "Expected values:",valuesE
        print "Output values:",values
      testcase.assertEqual(numpy.sum(okay), ncells)
    
  return


# End of file
