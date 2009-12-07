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
  tolerance = 1.0e-5

  dispE = testcase.calcDisplacements(data['vertices'])
  disp = data['vertex_fields']['displacement']

  # Check x displacements
  diff = numpy.abs(disp[:,0] - dispE[:,0])
  okay = diff < tolerance
  if numpy.sum(okay) != nvertices:
    "Error in x-component of displacement field."
    print "Expected values: ",dispE
    print "Output values: ",disp
  testcase.assertEqual(nvertices, numpy.sum(okay))    
    
  # Check y displacements
  diff = numpy.abs(disp[:,1] - dispE[:,1])
  okay = diff < tolerance
  if numpy.sum(okay) != nvertices:
    "Error in y-component of displacement field."
    print "Expected values: ",dispE
    print "Output values: ",disp
  testcase.assertEqual(nvertices, numpy.sum(okay))    

  # Check z displacements
  diff = numpy.abs(disp[:,2] - dispE[:,2])
  okay = diff < tolerance
  if numpy.sum(okay) != nvertices:
    "Error in z-component of displacement field."
    print "Expected values: ",dispE
    print "Output values: ",disp
  testcase.assertEqual(nvertices, numpy.sum(okay))    

  return


# End of file
