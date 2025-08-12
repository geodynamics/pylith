// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information. 
// =================================================================================================

// List of phases.
%typemap(in) (pylith::initializers::InitializePhase* phases[],
	      const size_t numPhases)
{
  // Check to make sure input is a list.
  if (PyList_Check($input)) {
    const size_t size = PyList_Size($input);
    $2 = size;
    $1 = (size > 0) ? new pylith::initializers::InitializePhase*[size] : 0;
    for (size_t i = 0; i < size; i++) {
      PyObject* s = PyList_GetItem($input,i);
      pylith::initializers::InitializePhase* phase = 0;
      int err = SWIG_ConvertPtr(s, (void**) &phase,
				$descriptor(pylith::initializers::InitializePhase*),
				0);
      if (SWIG_IsOK(err))
	$1[i] = (pylith::initializers::InitializePhase*) phase;
      else {
	PyErr_SetString(PyExc_TypeError, "List must contain initialize phases.");
	delete[] $1;
	return NULL;
      } // if
    } // for
  } else {
    PyErr_SetString(PyExc_TypeError, "Expected list of initialize phases.");
    return NULL;
  } // if/else
} // typemap(in) [List of phases.]

// This cleans up the array we malloc'd before the function call
%typemap(freearg) (pylith::initializers::InitializePhase* phases[],
		   const size_t numPhases) {
  delete[] $1;
}


// End of file
