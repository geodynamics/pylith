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

// Map a Python float scalar to PylistReal.
%typemap(in) PylithReal {
  $1 = PyFloat_AsDouble($input);
 } // typemap

// Map a PylithScalar to a Python float
%typemap(out) PylithReal {
  $result = PyFloat_FromDouble($1);
 } // typemap


%typemap(typecheck,precedence=SWIG_TYPECHECK_FLOAT) PylithReal {
   $1 = PyFloat_Check($input) ? 1 : 0;
}



%{
static
int
convert_scalararray(PyObject* input,
		    PylithReal* const values,
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
} // convert_scalararray
%}

// Map a Python sequence of floats into a C PylithReal array.
%typemap(in) PylithReal [ANY] (PylithReal values[$1_dim0]) {
  if (!convert_scalararray($input, values, $1_dim0))
    return NULL;
  $1 = &values[0];
 } // typemap


// Typemap suite for (PylithReal* IN_ARRAY1, int DIM1)
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY)
  (PylithReal* IN_ARRAY1, int DIM1)
{
  $1 = is_array($input) || PySequence_Check($input);
}
%typemap(in)
  (PylithReal* IN_ARRAY1, int DIM1)
  (PyArrayObject* array=NULL, int is_new_object=0)
{
  if (sizeof(float) == sizeof(PylithReal)) {
    array = obj_to_array_contiguous_allow_conversion($input, NPY_FLOAT, &is_new_object);
  } else if (sizeof(double) == sizeof(PylithReal)) {
    array = obj_to_array_contiguous_allow_conversion($input, NPY_DOUBLE, &is_new_object);
  } else {
    PyErr_Format(PyExc_TypeError, 
		 "Unknown size for PylithReal.  '%ld' given.", 
		 sizeof(PylithReal));
  } // if/else
  npy_intp size[1] = { -1 };
  if (!array || !require_dimensions(array, 1) || !require_size(array, size, 1)) SWIG_fail;
  $1 = (PylithReal*) array_data(array);
  $2 = (int) array_size(array,0);
  
}
%typemap(freearg)
(PylithReal* IN_ARRAY1, int DIM1)
{
  if (is_new_object$argnum && array$argnum)
    { Py_DECREF(array$argnum); }
}


/* Typemap suite for (PylithReal* IN_ARRAY2, int DIM1, int DIM2)
 */
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY)
  (PylithReal* IN_ARRAY2, int DIM1, int DIM2)
{
  $1 = is_array($input) || PySequence_Check($input);
}
%typemap(in)
  (PylithReal* IN_ARRAY2, int DIM1, int DIM2)
  (PyArrayObject* array=NULL, int is_new_object=0)
{
  if (sizeof(float) == sizeof(PylithReal)) {
    array = obj_to_array_contiguous_allow_conversion($input, NPY_FLOAT, &is_new_object);
  } else if (sizeof(double) == sizeof(PylithReal)) {
    array = obj_to_array_contiguous_allow_conversion($input, NPY_DOUBLE, &is_new_object);
  } else {
  } // if/else
  npy_intp size[2] = { -1, -1 };
  if (!array || !require_dimensions(array, 2) || !require_size(array, size, 2)) SWIG_fail;
  $1 = (PylithReal*) array_data(array);
  $2 = (int) array_size(array,0);
  $3 = (int) array_size(array,1);
}
%typemap(freearg)
  (PylithReal* IN_ARRAY2, int DIM1, int DIM2)
{
  if (is_new_object$argnum && array$argnum)
    { Py_DECREF(array$argnum); }
}


// Typemap suite for (PylithReal* IN_ARRAY3, int DIM1, int DIM2, int DIM3)
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY)
  (PylithReal* IN_ARRAY3, int DIM1, int DIM2, int DIM3)
{
  $1 = is_array($input) || PySequence_Check($input);
}
%typemap(in)
  (PylithReal* IN_ARRAY3, int DIM1, int DIM2, int DIM3)
  (PyArrayObject* array=NULL, int is_new_object=0)
{
  if (sizeof(float) == sizeof(PylithReal)) {
    array = obj_to_array_contiguous_allow_conversion($input, NPY_FLOAT, &is_new_object);
  } else if (sizeof(double) == sizeof(PylithReal)) {
    array = obj_to_array_contiguous_allow_conversion($input, NPY_DOUBLE, &is_new_object);
  } else {
    
  } // if/else
  npy_intp size[3] = { -1, -1, -1 };
  if (!array || !require_dimensions(array, 3) || !require_size(array, size, 3)) SWIG_fail;
  $1 = (PylithReal*) array_data(array);
  $2 = (int) array_size(array,0);
  $3 = (int) array_size(array,1);
  $4 = (int) array_size(array,2);
}
%typemap(freearg)
  (PylithReal* IN_ARRAY3, int DIM1, int DIM2, int DIM3)
{
  if (is_new_object$argnum && array$argnum)
    { Py_DECREF(array$argnum); }
}



// End of file
