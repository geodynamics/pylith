// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ======================================================================
//

// Treat const char* const* as a special case.
%typemap(in) (const char* const* string_list) {
  // Check to make sure input is a list.
  if (PyList_Check($input)) {
    const int size = PyList_Size($input);
    $1 = (size > 0) ? new char*[size] : 0;
    for (int i = 0; i < size; i++) {
      PyObject *s = PyList_GetItem($input,i);
      if (PyString_Check(s))
	$1[i] = PyString_AsString(PyList_GetItem($input,i));
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
} // typemap(in) const char**

// This cleans up the char** array we malloc'd before the function call
%typemap(freearg) (const char** string_list) {
  delete[] $1;
}


// Treat const char* const* as a special case.
%typemap(in) (const char* const* string_list, const int list_len) {
  // Check to make sure input is a list.
  if (PyList_Check($input)) {
    const int size = PyList_Size($input);
    $2 = size;
    $1 = (size > 0) ? new char*[size] : 0;
    for (int i = 0; i < size; i++) {
      PyObject *s = PyList_GetItem($input,i);
      if (PyString_Check(s))
	$1[i] = PyString_AsString(PyList_GetItem($input,i));
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
} // typemap(in) const char**

// This cleans up the char** array we malloc'd before the function call
%typemap(freearg) (const char** string_list, const int list_len) {
  delete[] $1;
}

// Treat argc & argv as a special case.
%typemap(in) (int argc, char** argv) {
  // Check to make sure input is a list.
  if (PyList_Check($input)) {
    const int size = PyList_Size($input);
    $1 = size;
    $2 = (size > 0) ? new char*[size+1] : 0;
    for (int i = 0; i < size; i++) {
      PyObject *s = PyList_GetItem($input,i);
      if (PyString_Check(s))
	$2[i] = PyString_AsString(PyList_GetItem($input,i));
      else {
	PyErr_SetString(PyExc_TypeError, "List must contain strings.");
	delete[] $2;
	return NULL;
      } // else
    } // for
    $2[size] = 0; // Must terminate argv with null!
  } else {
    PyErr_SetString(PyExc_TypeError, "Expected list of strings.");
    return NULL;
  } // if/else
} // typemap(in) char**

// This cleans up the char** array we malloc'd before the function call
%typemap(freearg) (int argc, char** argv) {
  delete[] $2;
}


%typemap(argout) (int *numNames, char ***outNames) {
  int       num   = *$1;
  char    **names = *$2;
  PyObject *l     = PyList_New(num);

  for(Py_ssize_t i = 0; i < (Py_ssize_t) num; ++i) {
    PyList_SetItem(l, i, PyString_FromString(names[i]));
    PetscFree(names[i]);
  }
  PetscFree(names);
  $result = l;
}
%typemap(in,numinputs=0) (int *numNames, char ***outNames)(int tempI, char **tempC) {
    $1 = &tempI;
    $2 = &tempC;
}

// End of file
