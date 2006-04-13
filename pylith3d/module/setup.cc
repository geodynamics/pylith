// -*- C++ -*-
// 
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
//                               Charles A. Williams
//                        Rensselaer Polytechnic Institute
//                        (C) 2005 All Rights Reserved
// 
//  All worldwide rights reserved.  A license to use, copy, modify and
//  distribute this software for non-commercial research purposes only
//  is hereby granted, provided that this copyright notice and
//  accompanying disclaimer is not modified or removed from the software.
//
//  DISCLAIMER:  The software is distributed "AS IS" without any express
//  or implied warranty, including but not limited to, any implied
//  warranties of merchantability or fitness for a particular purpose
//  or any warranty of non-infringement of any current or pending patent
//  rights.  The authors of the software make no representations about
//  the suitability of this software for any particular purpose.  The
//  entire risk as to the quality and performance of the software is with
//  the user.  Should the software prove defective, the user assumes the
//  cost of all necessary servicing, repair or correction.  In
//  particular, neither Rensselaer Polytechnic Institute, nor the authors
//  of the software are liable for any indirect, special, consequential,
//  or incidental damages related to the software, to the maximum extent
//  the law permits.
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
    
// version
// $Id: setup.cc,v 1.4 2005/06/07 19:39:11 willic3 Exp $

// End of file
