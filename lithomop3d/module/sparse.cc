// -*- C++ -*-
// 
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
//                               Charles A. Williams
//                        Rensselaer Polytechnic Institute
//                        (C) 2004 All Rights Reserved
// 
//  Copyright 2004 Rensselaer Polytechnic Institute.
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

#include <portinfo>
#include "journal/debug.h"

#include <Python.h>

#include "sparse.h"
#include "exceptionhandler.h"
#include "lithomop3d_externs.h"
#include <stdio.h>
#include <string.h>


// Compute number of nonzero entries in stiffness matrix

char pylithomop3d_cmp_stiffsz__doc__[] = "";
char pylithomop3d_cmp_stiffsz__name__[] = "cmp_stiffsz";

PyObject * pylithomop3d_cmp_stiffsz(PyObject *, PyObject *args)
{
  int numberGlobalEquations;
  PyObject* pyPointerToLm;
  PyObject* pyPointerToLmx;
  PyObject* pyPointerToInfiel;
  int connectivitySize;
  int numberElements;
  PyObject* pyPointerToElementTypeInfo;
  int totalNumberSlipperyNodes;

  int ok = PyArg_ParseTuple(args, "iOOOiiOi:cmp_stiffsz",
			    &numberGlobalEquations,
			    &pyPointerToLm,
			    &pyPointerToLmx,
			    &pyPointerToInfiel,
			    &connectivitySize,
			    &numberElements,
			    &pyPointerToElementTypeInfo,
			    &totalNumberSlipperyNodes);

  // printf("Hello from pylithomop3d_cmp_stiffsz!\n");
  // printf("numberGlobalEquations = %d\n", numberGlobalEquations);
  // printf("connectivitySize = %d\n", connectivitySize);
  // printf("numberElements = %d\n", numberElements);
  // printf("workingArraySize = %d\n", workingArraySize);


  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 1024;
  char errorstring[maxsize];
  int* pointerToLm = (int*) PyCObject_AsVoidPtr(pyPointerToLm);
  int* pointerToLmx = (int*) PyCObject_AsVoidPtr(pyPointerToLmx);
  int* pointerToInfiel = (int*) PyCObject_AsVoidPtr(pyPointerToInfiel);
  int* pointerToElementTypeInfo = (int*) PyCObject_AsVoidPtr(pyPointerToElementTypeInfo);
  int workingArraySize = 0;

  // printf("pointerToLm = %d\n", pointerToLm);
  // printf("pointerToLmx = %d\n", pointerToLmx);
  // printf("pointerToInfiel = %d\n", pointerToInfiel);
  // printf("pointerToElementTypeInfo = %d\n", pointerToElementTypeInfo);
  // printf("connectivitySize = %d\n", connectivitySize);
  // printf("numberElements = %d\n", numberElements);
  // printf("workingArraySize = %d\n", workingArraySize);
  // printf("totalNumberSlipperyNodes = %d\n", totalNumberSlipperyNodes);

  cmp_stiffsz_f(&numberGlobalEquations,
	   	pointerToLm,
	   	pointerToLmx,
	   	pointerToInfiel,
	   	&connectivitySize,
	   	&numberElements,
	   	pointerToElementTypeInfo,
	   	&workingArraySize,
	   	&totalNumberSlipperyNodes,
	   	&errorcode,
	   	errorstring,
	   	strlen(errorstring));

  // printf("workingArraySize = %d\n", workingArraySize);

  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }
		  
  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "workingArraySize:" << workingArraySize
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue("i", workingArraySize);
}



// Create linked list of nonzero row and column entries in the stiffness matrix

char pylithomop3d_lnklst__doc__[] = "";
char pylithomop3d_lnklst__name__[] = "lnklst";

PyObject * pylithomop3d_lnklst(PyObject *, PyObject *args)
{
  int numberGlobalEquations;
  PyObject* pyPointerToLm;
  PyObject* pyPointerToLmx;
  PyObject* pyPointerToInfiel;
  int connectivitySize;
  int numberElements;
  PyObject* pyPointerToElementTypeInfo;
  PyObject* pyPointerToIndx;
  PyObject* pyPointerToLink;
  PyObject* pyPointerToNbrs;
  int workingArraySize;
  int totalNumberSlipperyNodes;

  int ok = PyArg_ParseTuple(args, "iOOOiiOOOOii:lnklst",
			    &numberGlobalEquations,
			    &pyPointerToLm,
			    &pyPointerToLmx,
			    &pyPointerToInfiel,
			    &connectivitySize,
			    &numberElements,
			    &pyPointerToElementTypeInfo,
			    &pyPointerToIndx,
			    &pyPointerToLink,
			    &pyPointerToNbrs,
			    &workingArraySize,
			    &totalNumberSlipperyNodes);

  // printf("Hello from pylithomop3d_lnklst!\n");
  // printf("numberGlobalEquations = %d\n", numberGlobalEquations);
  // printf("connectivitySize = %d\n", connectivitySize);
  // printf("numberElements = %d\n", numberElements);
  // printf("workingArraySize = %d\n", workingArraySize);


  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 1024;
  char errorstring[maxsize];
  int* pointerToLm = (int*) PyCObject_AsVoidPtr(pyPointerToLm);
  int* pointerToLmx = (int*) PyCObject_AsVoidPtr(pyPointerToLmx);
  int* pointerToInfiel = (int*) PyCObject_AsVoidPtr(pyPointerToInfiel);
  int* pointerToElementTypeInfo = (int*) PyCObject_AsVoidPtr(pyPointerToElementTypeInfo);
  int* pointerToIndx = (int*) PyCObject_AsVoidPtr(pyPointerToIndx);
  int* pointerToLink = (int*) PyCObject_AsVoidPtr(pyPointerToLink);
  int* pointerToNbrs = (int*) PyCObject_AsVoidPtr(pyPointerToNbrs);
  int stiffnessMatrixSize =0;
  int stiffnessOffDiagonalSize =0;

  // printf("pointerToLm = %d\n", pointerToLm);
  // printf("pointerToLmx = %d\n", pointerToLmx);
  // printf("pointerToInfiel = %d\n", pointerToInfiel);
  // printf("pointerToElementTypeInfo = %d\n", pointerToElementTypeInfo);
  // printf("pointerToIndx = %d\n", pointerToIndx);
  // printf("pointerToLink = %d\n", pointerToLink);
  // printf("pointerToNbrs = %d\n", pointerToNbrs);
  // printf("connectivitySize = %d\n", connectivitySize);
  // printf("numberElements = %d\n", numberElements);
  // printf("stiffnessMatrixSize = %d\n", stiffnessMatrixSize);
  // printf("stiffnessOffDiagonalSize = %d\n", stiffnessOffDiagonalSize);
  // printf("totalNumberSlipperyNodes = %d\n", totalNumberSlipperyNodes);

  lnklst_f(&numberGlobalEquations,
	   pointerToLm,
	   pointerToLmx,
	   pointerToInfiel,
	   &connectivitySize,
	   &numberElements,
	   pointerToElementTypeInfo,
	   pointerToIndx,
	   pointerToLink,
	   pointerToNbrs,
	   &workingArraySize,
	   &stiffnessOffDiagonalSize,
	   &stiffnessMatrixSize,
	   &totalNumberSlipperyNodes,
	   &errorcode,
	   errorstring,
	   strlen(errorstring));

  // printf("stiffnessOffDiagonalSize = %d\n", stiffnessOffDiagonalSize);
  // printf("stiffnessMatrixSize = %d\n", stiffnessMatrixSize);

  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }
		  
  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "stiffnessMatrixSize:" << stiffnessMatrixSize
    << journal::endl;

  // return
  return Py_BuildValue("ii", stiffnessMatrixSize,
   		  stiffnessOffDiagonalSize);
}

// Create a PETSc Mat
#include <petscmat.h>
char pylithomop3d_createPETScMat__doc__[] = "";
char pylithomop3d_createPETScMat__name__[] = "createPETScMat";

PyObject * pylithomop3d_createPETScMat(PyObject *, PyObject *args)
{
  PyObject *pyA;
  Mat A;
  int size;

  int ok = PyArg_ParseTuple(args, "i:createPETScMat", &size);
  if (!ok) {
    return 0;
  }

  if (MatCreate(PETSC_COMM_WORLD, size, size, PETSC_DETERMINE, PETSC_DETERMINE, &A)) {
    PyErr_SetString(PyExc_RuntimeError, "Could not create PETSc Mat");
    return 0;
  }

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "Created PETSc Mat:" << size
    << journal::endl;

  // return Py_None;
  pyA = PyCObject_FromVoidPtr(A, NULL);
  return Py_BuildValue("N", pyA);
}

// Destroy a PETSc Mat

char pylithomop3d_destroyPETScMat__doc__[] = "";
char pylithomop3d_destroyPETScMat__name__[] = "destroyPETScMat";

PyObject * pylithomop3d_destroyPETScMat(PyObject *, PyObject *args)
{
  PyObject *pyA;
  Mat A;
  int size;

  int ok = PyArg_ParseTuple(args, "O:destroyPETScMat", &pyA);
  if (!ok) {
    return 0;
  }

  A = (Mat) PyCObject_AsVoidPtr(pyA);
  if (MatDestroy(A)) {
    PyErr_SetString(PyExc_RuntimeError, "Could not destroy PETSc Mat");
    return 0;
  }

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "Destroyed PETSc Mat"
    << journal::endl;

  // return Py_None;
  Py_INCREF(Py_None);
  return Py_None;
}

// Transform linked list into index array for modified sparse row format

char pylithomop3d_makemsr__doc__[] = "";
char pylithomop3d_makemsr__name__[] = "makemsr";

PyObject * pylithomop3d_makemsr(PyObject *, PyObject *args)
{
  PyObject* pyA;
  PyObject* pyPointerToIndx;
  PyObject* pyPointerToLink;
  PyObject* pyPointerToNbrs;
  int numberGlobalEquations;
  int stiffnessMatrixSize;
  int workingArraySize;

  int ok = PyArg_ParseTuple(args, "OOOOiii:makemsr",
                            &pyA,
			    &pyPointerToIndx,
			    &pyPointerToLink,
			    &pyPointerToNbrs,
			    &numberGlobalEquations,
			    &stiffnessMatrixSize,
			    &workingArraySize);

  // printf("Hello from pylithomop3d_makemsr!\n");
  // printf("numberGlobalEquations = %d\n", numberGlobalEquations);
  // printf("stiffnessMatrixSize = %d\n", stiffnessMatrixSize);
  // printf("workingArraySize = %d\n", workingArraySize);

  if (!ok) {
    return 0;
  }
  Mat A = (Mat) PyCObject_AsVoidPtr(pyA);
  int* pointerToIndx = (int*) PyCObject_AsVoidPtr(pyPointerToIndx);
  int* pointerToLink = (int*) PyCObject_AsVoidPtr(pyPointerToLink);
  int* pointerToNbrs = (int*) PyCObject_AsVoidPtr(pyPointerToNbrs);
  int minimumNonzeroTermsPerRow = 0;
  int maximumNonzeroTermsPerRow = 0;
  double averageNonzeroTermsPerRow = 0.0;

  makemsr_f(&A,
	    pointerToIndx,
	    pointerToLink,
	    pointerToNbrs,
	    &numberGlobalEquations,
	    &stiffnessMatrixSize,
	    &workingArraySize,
	    &minimumNonzeroTermsPerRow,
	    &maximumNonzeroTermsPerRow,
	    &averageNonzeroTermsPerRow);

  // printf("workingArraySize = %d\n", workingArraySize);
  // printf("minimumNonzeroTermsPerRow = %d\n", minimumNonzeroTermsPerRow);
  // printf("maximumNonzeroTermsPerRow = %d\n", maximumNonzeroTermsPerRow);
  // printf("averageNonzeroTermsPerRow = %g\n", averageNonzeroTermsPerRow);
		  
  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "workingArraySize:" << workingArraySize
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  // return Py_None;
  return Py_BuildValue("iid",
		       minimumNonzeroTermsPerRow,
		       maximumNonzeroTermsPerRow,
		       averageNonzeroTermsPerRow);
}


// version
// $Id: sparse.cc,v 1.7 2005/03/11 04:07:42 knepley Exp $

// End of file
