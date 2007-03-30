// -*- C++ -*-
// 
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  PyLith by Charles A. Williams, Brad Aagaard, and Matt Knepley
//
//  Copyright (c) 2004-2006 Rensselaer Polytechnic Institute
//
//  Permission is hereby granted, free of charge, to any person obtaining
//  a copy of this software and associated documentation files (the
//  "Software"), to deal in the Software without restriction, including
//  without limitation the rights to use, copy, modify, merge, publish,
//  distribute, sublicense, and/or sell copies of the Software, and to
//  permit persons to whom the Software is furnished to do so, subject to
//  the following conditions:
//
//  The above copyright notice and this permission notice shall be
//  included in all copies or substantial portions of the Software.
//
//  THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
//  EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
//  MERCHANTABILITY,    FITNESS    FOR    A   PARTICULAR    PURPOSE    AND
//  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
//  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
//  OF CONTRACT, TORT OR OTHERWISE,  ARISING FROM, OUT OF OR IN CONNECTION
//  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 

#include <petscsnes.h>
#include <portinfo>
#include "journal/debug.h"

#include <Python.h>

#include "setup.h"
#include "exceptionhandler.h"
#include <stdio.h>
#include <string.h>


// Initialize PETSc

char pypylith3d_petsc_initialize__doc__[] = "";
char pypylith3d_petsc_initialize__name__[] = "PetscInitialize";

PyObject * pypylith3d_petsc_initialize(PyObject *, PyObject *args)
{
  PyObject *argList;
  if (!PyArg_ParseTuple(args, "O:PetscInitialize", &argList)) {
    return 0;
  }

  int argc = PySequence_Length(argList);
  char **argv = (char **) malloc((argc + 1) * sizeof(char *));
  int a;

  for(a = 0; a < argc; a++) {
    char *arg = PyString_AsString(PyList_GetItem(argList, a));

    argv[a] = (char *) malloc((strlen(arg)+1)* sizeof(char));
    strcpy(argv[a], arg);
  }
  argv[argc] = NULL;
  if (PetscInitialize(&argc, &argv, NULL, "Pylith 3D")) {
    PyErr_SetString(PyExc_RuntimeError, "PETSc failed to initialize");
    return 0;
  }
    
  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "Initializing PETSc"
    << journal::endl;

  // Do not free 'argv' -- PETSc saves a reference to it ('PetscGlobalArgs').

  // return
  Py_INCREF(Py_None);
  return Py_None;
}

// Finalize PETSc

char pypylith3d_petsc_finalize__doc__[] = "";
char pypylith3d_petsc_finalize__name__[] = "PetscFinalize";

PyObject * pypylith3d_petsc_finalize(PyObject *, PyObject *)
{
  if (PetscFinalize()) {
    PyErr_SetString(PyExc_RuntimeError, "PETSc failed to finalize");
    return 0;
  }
    
  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "Finalizing PETSc"
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}

// Setup PETSc Logging

char pypylith3d_setup_petsc_logging__doc__[] = "";
char pypylith3d_setup_petsc_logging__name__[] = "setupPETScLogging";

PyObject * pypylith3d_setup_petsc_logging(PyObject *, PyObject *)
{
  PetscInt autoprestrStage, elasticStage, viscousStage;
  PetscEvent iterateEvent;
  PyObject *_pyReturn;

  if (PetscLogStageRegister(&autoprestrStage, "AutoPrestress Solve")) {
    PyErr_SetString(PyExc_RuntimeError, "PETSc failed to create stage");
    return 0;
  }
  if (PetscLogStageRegister(&elasticStage, "Elastic Solve")) {
    PyErr_SetString(PyExc_RuntimeError, "PETSc failed to create stage");
    return 0;
  }
  if (PetscLogStageRegister(&viscousStage, "Viscous Solve")) {
    PyErr_SetString(PyExc_RuntimeError, "PETSc failed to create stage");
    return 0;
  }
  if (PetscLogEventRegister(&iterateEvent, "Iterate", KSP_COOKIE)) {
    PyErr_SetString(PyExc_RuntimeError, "PETSc failed to create event");
    return 0;
  }
    
  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "Setup PETSc logging"
    << journal::endl;

  // return logging handles
  _pyReturn = Py_BuildValue("iiii", autoprestrStage, elasticStage, viscousStage, iterateEvent);
  return _pyReturn;
}

// End of file
