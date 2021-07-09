// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ======================================================================
//

// ----------------------------------------------------------------------
// List of solution output managers.
%typemap(in) (pylith::meshio::OutputSoln* outputArray[],
        const int numOutputs)
{
  // Check to make sure input is a list.
  if (PyList_Check($input)) {
    const int size = PyList_Size($input);
    $2 = size;
    $1 = (size > 0) ? new pylith::meshio::OutputSoln*[size] : 0;
    for (int i = 0; i < size; i++) {
      PyObject* s = PyList_GetItem($input,i);
      pylith::meshio::OutputSoln* output = 0;
      int err = SWIG_ConvertPtr(s, (void**) &output,
        $descriptor(pylith::meshio::OutputSoln*),
        0);
      if (SWIG_IsOK(err))
  $1[i] = (pylith::meshio::OutputSoln*) output;
      else {
  PyErr_SetString(PyExc_TypeError, "List must contain solution output managers.");
  delete[] $1;
  return NULL;
      } // if
    } // for
  } else {
    PyErr_SetString(PyExc_TypeError, "Expected list of solution output managers.");
    return NULL;
  } // if/else
} // typemap(in) [List of output managers.]

// This cleans up the array we malloc'd before the function call
%typemap(freearg) (pylith::meshio::OutputSoln* outputArray[],
       const int numOutputs) {
  delete[] $1;
}

// End of file
