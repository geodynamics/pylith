// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

// ----------------------------------------------------------------------
// List of earthquake sources.
%typemap(in) (pylith::faults::EqKinSrc** sources,
	      const int numSources)
{
  // Check to make sure input is a list.
  if (PyList_Check($input)) {
    const int size = PyList_Size($input);
    $2 = size;
    $1 = (size > 0) ? new pylith::faults::EqKinSrc*[size] : 0;
    for (int i = 0; i < size; i++) {
      PyObject* s = PyList_GetItem($input,i);
      pylith::faults::EqKinSrc** src = 0;
      int err = SWIG_ConvertPtr(s, (void**) &src, 
				$descriptor(pylith::faults::EqKinSrc*),
				0);
      if (SWIG_IsOK(err))
	$1[i] = (pylith::faults::EqKinSrc*) src;
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
%typemap(freearg) (pylith::faults::EqKinSrc** sources, 
		   const int numSources) {
  delete[] $1;
}


// End of file
