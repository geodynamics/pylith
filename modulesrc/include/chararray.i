// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information. 
// =================================================================================================
// Treat const char* const* as a special case.
%typemap(in) (const char* const* string_list) {
  // Check to make sure input is a list.
  if (PyList_Check($input)) {
    const int size = PyList_Size($input);
    $1 = (size > 0) ? new char*[size] : NULL;
    for (int i = 0; i < size; i++) {
      PyObject *s = PyList_GetItem($input,i);
      if (PyUnicode_Check(s))
        $1[i] = (char*)PyUnicode_AsUTF8(PyList_GetItem($input,i));
      else {
	PyErr_SetString(PyExc_TypeError, "List must contain strings.");
	delete[] $1;
	return NULL;
      } // else
    } // for
  } else {
    PyErr_SetString(PyExc_TypeError, "Expected list of strings.");
    return NULL;
  } // if/else
} // typemap(in) const char* const*

// This cleans up the char** array we malloc'd before the function call
%typemap(freearg) (const char* const* string_list) {
  delete[] $1;
}


// Treat const char* const* as a special case.
%typemap(in) (const char* const* string_list, const int list_len) {
  // Check to make sure input is a list.
  if (PyList_Check($input)) {
    const int size = PyList_Size($input);
    $2 = size;
    $1 = (size > 0) ? new char*[size] : NULL;
    for (int i = 0; i < size; i++) {
      PyObject *s = PyList_GetItem($input,i);
      if (PyUnicode_Check(s))
        $1[i] = (char*)PyUnicode_AsUTF8(PyList_GetItem($input,i));
      else {
	PyErr_SetString(PyExc_TypeError, "List must contain strings.");
	delete[] $1;
	return NULL;
      } // else
    } // for
  } else {
    PyErr_SetString(PyExc_TypeError, "Expected list of strings.");
    return NULL;
  } // if/else
} // typemap(in) const char* const*

// This cleans up the char** array we malloc'd before the function call
%typemap(freearg) (const char* const* string_list, const int list_len) {
  delete[] $1;
}


// Treat argc & argv as a special case.
%typemap(in) (int argc, char** argv) {
  // Check to make sure input is a list.
  if (PyList_Check($input)) {
    const int size = PyList_Size($input);
    $1 = size;
    $2 = (size > 0) ? new char*[size+1] : NULL;
    for (int i = 0; i < size; i++) {
      PyObject *s = PyList_GetItem($input,i);
      if (PyUnicode_Check(s))
	$2[i] = (char*)PyUnicode_AsUTF8(PyList_GetItem($input,i));
      else {
	PyErr_SetString(PyExc_TypeError, "List must contain strings.");
	delete[] $2;
	return NULL;
      } // else
    } // for
    $2[size] = NULL; // Must terminate argv with null!
  } else {
    PyErr_SetString(PyExc_TypeError, "Expected list of strings.");
    return NULL;
  } // if/else
} // typemap(in) char**

// This cleans up the char** array we malloc'd before the function call
%typemap(freearg) (int argc, char** argv) {
  delete[] $2;
}


// End of file
