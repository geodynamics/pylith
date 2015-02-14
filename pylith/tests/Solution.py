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

## @file pylith/tests/Solution.py
##
## @brief Check displacement solution output from PyLith.

import numpy
import h5py

def check_displacements(testcase, filename, mesh):
  """
  Check displacements.
  """
  h5 = h5py.File(filename, "r", driver="sec2")
  
  # Check vertices
  vertices = h5['geometry/vertices'][:]
  (nvertices, spaceDim) = vertices.shape
  testcase.assertEqual(mesh['nvertices'], nvertices)
  testcase.assertEqual(mesh['spaceDim'], spaceDim)

  # Check displacement solution
  toleranceAbsMask = 0.1
  tolerance = 1.0e-5

  dispE = testcase.calcDisplacements(vertices)
  disp = h5['vertex_fields/displacement'][:]

  (nstepsE, nverticesE, ncompsE) = dispE.shape
  (nsteps, nvertices, ncomps) = disp.shape
  testcase.assertEqual(nstepsE, nsteps)
  testcase.assertEqual(nverticesE, nvertices)
  testcase.assertEqual(ncompsE, ncomps)

  for istep in xrange(nsteps):
    for icomp in xrange(ncomps):

      mask = numpy.abs(dispE[istep,:,icomp]) > toleranceAbsMask
      diff = numpy.abs(disp[istep,:,icomp] - dispE[istep,:,icomp])
      diffR = numpy.abs(1.0 - disp[istep,:,icomp] / dispE[istep,:,icomp])  
      okay = ~mask * (diff < tolerance) + mask * (diffR < tolerance)
      if numpy.sum(okay) != nvertices:
        print "Error in component %d of displacement field at time step %d." % (icomp, istep)
        print "Expected values: ",dispE[istep,:,:]
        print "Output values: ",disp[istep,:,:]
        print "Expected values (not okay): ",dispE[istep,~okay,icomp]
        print "Computed values (not okay): ",disp[istep,~okay,icomp]
        print "Relative diff (not okay): ",diffR[~okay]
        print "Coordinates (not okay): ",vertices[~okay,:]
      testcase.assertEqual(nvertices, numpy.sum(okay))    
    
  h5.close()
  return


# End of file
