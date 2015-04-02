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

  from spatialdata.units.NondimElasticQuasistatic import NondimElasticQuasistatic
  normalizer = NondimElasticQuasistatic()
  normalizer._configure()

  # Check state variables
  tolerance = 1.0e-6

  for name in stateVarNames:
    valuesE = testcase.calcStateVar(name, vertices, cells)
    values = h5['cell_fields/%s' % name][:]

    (nstepsE, ncellsE, ncompsE) = valuesE.shape
    (nsteps, ncells, ncomps) = values.shape

    testcase.assertEqual(nstepsE, nsteps)
    testcase.assertEqual(ncellsE, ncells)
    testcase.assertEqual(ncompsE, ncomps)

    scale = 1.0
    if name == "stress":
      scale *= normalizer.pressureScale().value

    for istep in xrange(nsteps):
      for icomp in xrange(ncomps):
        ratio = numpy.abs(1.0 - values[istep,:,icomp]/valuesE[istep,:,icomp])
        diff = numpy.abs(values[istep,:,icomp] - valuesE[istep,:,icomp]) / scale
        mask = valuesE[istep,:,icomp] != 0.0
        okay = mask*(ratio < tolerance) + ~mask*(diff < tolerance)
        if numpy.sum(okay) != ncells:
          print "Error in component %d of state variable '%s' at time step %d." % (icomp, name, istep)
          print "Expected values:",valuesE
          print "Output values:",values
        testcase.assertEqual(ncells, numpy.sum(okay))

  h5.close()
  return


# End of file
