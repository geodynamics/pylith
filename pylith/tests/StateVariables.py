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

## @file pylith/tests/StateVariables.py
##
## @brief Check state variables output from PyLith.

import numpy
import h5py

def check_state_variables(testcase, filename, mesh, stateVarNames):
  """
  Check state variables.
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

  # Check state variables
  tolerance = 2.0e-5

  for name in stateVarNames:
    valuesE = testcase.calcStateVar(name, vertices, cells)
    values = h5['cell_fields/%s' % name][:]

    (nstepsE, ncellsE, ncompsE) = valuesE.shape
    (nsteps, ncells, ncomps) = values.shape

    testcase.assertEqual(nstepsE, nsteps)
    testcase.assertEqual(ncellsE, ncells)
    testcase.assertEqual(ncompsE, ncomps)

    if "getValueScale" in dir(testcase):
      scale = testcase.getValueScale(name)
    else:
      if name == "total_strain":
        scale = 1.0
      elif name == "stress" or name == "cauchy_stress":
        scale = 1.0e+9
      else:
        scale = 1.0
      
    for istep in xrange(nsteps):
      for icomp in xrange(ncomps):
        diff = numpy.abs(values[istep,:,icomp] - valuesE[istep,:,icomp]) / scale
        okay = diff < tolerance
        if numpy.sum(okay) != ncells:
          print "Error in component %d of state variable '%s' at time step %d." % (icomp, name, istep)
          print "Expected values:",valuesE
          print "Output values:",values
          print "Expected values (not okay): ",valuesE[istep,~okay,icomp]
          print "Output values (not okay): ",values[istep,~okay,icomp]
          print "Scaled diff (not okay): ",diff[~okay]
          print "Scale",scale

        testcase.assertEqual(ncells, numpy.sum(okay))

  h5.close()
  return


# End of file
