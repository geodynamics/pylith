// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information. 
// =================================================================================================

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
