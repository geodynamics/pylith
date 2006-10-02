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
#include "pylith3d_externs.h"
#include <stdio.h>
#include <string.h>


// Initialize PETSc

char pypylith3d_petsc_initialize__doc__[] = "";
char pypylith3d_petsc_initialize__name__[] = "PetscInitialize";

PyObject * pypylith3d_petsc_initialize(PyObject *, PyObject *)
{
  PyObject *sysMod = PyImport_ImportModule("sys");
  PyObject *argList = PyObject_GetAttrString(sysMod, "argv");
  int argc = PySequence_Length(argList);
  char **argv = (char **) malloc(argc * sizeof(char *));
  int a;

  for(a = 0; a < argc; a++) {
    char *arg = PyString_AsString(PyList_GetItem(argList, a));

    argv[a] = (char *) malloc((strlen(arg)+1)* sizeof(char));
    strcpy(argv[a], arg);
  }
  if (PetscInitialize(&argc, &argv, NULL, "Pylith 3D")) {
    PyErr_SetString(PyExc_RuntimeError, "PETSc failed to initialize");
    return 0;
  }
    
  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "Initializing PETSc"
    << journal::endl;

  // cleanup
  for(a = 0; a < argc; a++) {
    free(argv[a]);
  }
  free(argv);
  Py_DECREF(argList);
  Py_DECREF(sysMod);
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

// Initialize material model info

char pypylith3d_matmod_def__doc__[] = "";
char pypylith3d_matmod_def__name__[] = "matmod_def";

PyObject * pypylith3d_matmod_def(PyObject *, PyObject *args)
{
  PyObject* pyPointerToMaterialModelInfo;

  int ok = PyArg_ParseTuple(args, "O:matmod_def",
			    &pyPointerToMaterialModelInfo);

  if (!ok) {
    return 0;
  }

  int* pointerToMaterialModelInfo = (int*) PyCObject_AsVoidPtr(pyPointerToMaterialModelInfo);

  matmod_def_f(pointerToMaterialModelInfo);
    
  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "pointerToMaterialModelInfo:" << pointerToMaterialModelInfo
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}

// Precompute shape function info

char pypylith3d_preshape__doc__[] = "";
char pypylith3d_preshape__name__[] = "preshape";

PyObject * pypylith3d_preshape(PyObject *, PyObject *args)
{
  PyObject* pyPointerToSh;
  PyObject* pyPointerToShj;
  PyObject* pyPointerToGauss;
  int quadratureOrderInt;
  int elementType;
  int numberVolumeElementNodes;
  int numberVolumeElementGaussPoints;

  int ok = PyArg_ParseTuple(args, "OOOiiii:preshape",
			    &pyPointerToSh,
			    &pyPointerToShj,
			    &pyPointerToGauss,
			    &quadratureOrderInt,
			    &elementType,
			    &numberVolumeElementNodes,
			    &numberVolumeElementGaussPoints);

  if (!ok) {
    return 0;
  }

  double* pointerToSh = (double*) PyCObject_AsVoidPtr(pyPointerToSh);
  double* pointerToShj = (double*) PyCObject_AsVoidPtr(pyPointerToShj);
  double* pointerToGauss = (double*) PyCObject_AsVoidPtr(pyPointerToGauss);
  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];

  preshape_f(pointerToSh,
	     pointerToShj,
	     pointerToGauss,
	     &quadratureOrderInt,
	     &elementType,
	     &numberVolumeElementNodes,
	     &numberVolumeElementGaussPoints,
	     &errorcode,
	     errorstring,
	     sizeof(errorstring));

  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }
    
  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "quadratureOrderInt:" << quadratureOrderInt
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}

// Precompute shape function info for element faces

char pypylith3d_preshape2d__doc__[] = "";
char pypylith3d_preshape2d__name__[] = "preshape2d";

PyObject * pypylith3d_preshape2d(PyObject *, PyObject *args)
{
  PyObject* pyPointerToSh2d;
  PyObject* pyPointerToGauss2d;
  int quadratureOrderInt;
  int elementType;
  int numberSurfaceElementNodes;
  int numberSurfaceElementGaussPoints;

  int ok = PyArg_ParseTuple(args, "OOiiii:preshape2d",
			    &pyPointerToSh2d,
			    &pyPointerToGauss2d,
			    &quadratureOrderInt,
			    &elementType,
			    &numberSurfaceElementNodes,
			    &numberSurfaceElementGaussPoints);

  if (!ok) {
    return 0;
  }

  double* pointerToSh2d = (double*) PyCObject_AsVoidPtr(pyPointerToSh2d);
  double* pointerToGauss2d = (double*) PyCObject_AsVoidPtr(pyPointerToGauss2d);
  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];

  preshape2d_f(pointerToSh2d,
	       pointerToGauss2d,
	       &quadratureOrderInt,
	       &elementType,
	       &numberSurfaceElementNodes,
	       &numberSurfaceElementGaussPoints,
	       &errorcode,
	       errorstring,
	       sizeof(errorstring));

  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }
    
  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "quadratureOrderInt:" << quadratureOrderInt
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}
// version
// $Id: setup.cc,v 1.4 2005/06/07 19:39:11 willic3 Exp $

// End of file
