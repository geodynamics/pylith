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

## @file pylith/tests/Solution.py
##
## @brief Check displacement solution output from PyLith.

import numpy

def check_displacements(testcase, filename, mesh):
  """
  Check displacements.
  """
  data = testcase.reader.read(filename)
  
  # Check vertices
  (nvertices, spaceDim) = data['vertices'].shape
  testcase.assertEqual(mesh['nvertices'], nvertices)
  testcase.assertEqual(mesh['spaceDim'], spaceDim)

  # Check displacement solution
  tolerance = 1.0e-6

  dispE = testcase.calcDisplacements(data['vertices'])
  disp = data['vertex_fields']['displacement']

  # Check x displacements
  mask = numpy.abs(dispE[:,0]) > tolerance
  diff = numpy.abs(disp[:,0] - dispE[:,0])
  diffR = numpy.abs(1.0 - disp[:,0] / dispE[:,0])  
  okay = ~mask * (diff < tolerance) + mask * (diffR < tolerance)
  if numpy.sum(okay) != nvertices:
    print "Error in x-component of displacement field."
    print "Expected values: ",dispE
    print "Output values: ",disp
  testcase.assertEqual(nvertices, numpy.sum(okay))    
    
  # Check y displacements
  mask = numpy.abs(dispE[:,1]) > tolerance
  diff = numpy.abs(disp[:,1] - dispE[:,1])
  diffR = numpy.abs(1.0 - disp[:,1] / dispE[:,1])  
  okay = ~mask * (diff < tolerance) + mask * (diffR < tolerance)
  if numpy.sum(okay) != nvertices:
    print "Error in y-component of displacement field."
    print "Expected values: ",dispE
    print "Output values: ",disp
  testcase.assertEqual(nvertices, numpy.sum(okay))    

  # Check z displacements
  mask = numpy.abs(dispE[:,2]) > tolerance
  diff = numpy.abs(disp[:,2] - dispE[:,2])
  diffR = numpy.abs(1.0 - disp[:,2] / dispE[:,2])  
  okay = ~mask * (diff < tolerance) + mask * (diffR < tolerance)
  if numpy.sum(okay) != nvertices:
    print "Error in z-component of displacement field."
    print "Expected values: ",dispE
    print "Output values: ",disp
  testcase.assertEqual(nvertices, numpy.sum(okay))    

  return


# End of file
