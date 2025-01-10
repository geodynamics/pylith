// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information. 
// =================================================================================================
// List of earthquake sources.
%typemap(in) (pylith::faults::KinSrc** ruptures,
	      const int numRuptures)
{
  // Check to make sure input is a list.
  if (PyList_Check($input)) {
    const int size = PyList_Size($input);
    $2 = size;
    $1 = (size > 0) ? new pylith::faults::KinSrc*[size] : 0;
    for (int i = 0; i < size; i++) {
      PyObject* s = PyList_GetItem($input,i);
      pylith::faults::KinSrc** src = 0;
      int err = SWIG_ConvertPtr(s, (void**) &src, 
				$descriptor(pylith::faults::KinSrc*),
				0);
      if (SWIG_IsOK(err))
	$1[i] = (pylith::faults::KinSrc*) src;
      else {
	PyErr_SetString(PyExc_TypeError, 
			"List must contain kinematic earthquake sources.");
	delete[] $1;
	return NULL;
      } // if
    } // for
  } else {
    PyErr_SetString(PyExc_TypeError,
		    "Expected list of kinematic earthquake sources.");
    return NULL;
  } // if/else
} // typemap(in) [List of kinematic earthquake sources.]

// This cleans up the array we malloc'd before the function call
%typemap(freearg) (pylith::faults::KinSrc** ruptures, 
		   const int numRuptures) {
  delete[] $1;
}


// End of file
