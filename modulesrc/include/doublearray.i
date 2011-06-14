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
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

%{
static
int
convert_doublearray(PyObject* input,
		    double* const values,
		    const int size) {
  if (!PySequence_Check(input)) {
    PyErr_SetString(PyExc_TypeError, "Expecting a sequence of floats.");
    return 0;
  } // if
  if (PyObject_Length(input) != size) {
    PyErr_SetString(PyExc_ValueError, "Sequence size mismatch.");
    return 0;
  } // if
  for (int i=0; i < size; i++) {
    PyObject *o = PySequence_GetItem(input,i);
    if (!PyFloat_Check(o)) {
      Py_XDECREF(o);
      PyErr_SetString(PyExc_ValueError,"Expecting a sequence of floats.");
      return 0;
    } // if
    values[i] = PyFloat_AsDouble(o);
    Py_DECREF(o);
  } // for
  return 1;
} // convert_doublearray
%}

// Map a Python sequence of floats into a C double array.
%typemap(in) double [ANY] (double values[$1_dim0]) {
  if (!convert_doublearray($input, values, $1_dim0))
    return NULL;
  $1 = &values[0];
 } // typemap


// End of file
