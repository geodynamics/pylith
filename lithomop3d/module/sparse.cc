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

#include <petscmat.h>
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
  int numberVolumeElements;
  int totalNumberSlipperyNodes;
  int numberVolumeElementNodes;

  int ok = PyArg_ParseTuple(args, "iOOiii:cmp_stiffsz",
			    &numberGlobalEquations,
			    &pyPointerToLm,
			    &pyPointerToLmx,
			    &numberVolumeElements,
			    &totalNumberSlipperyNodes,
			    &numberVolumeElementNodes);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int* pointerToLm = (int*) PyCObject_AsVoidPtr(pyPointerToLm);
  int* pointerToLmx = (int*) PyCObject_AsVoidPtr(pyPointerToLmx);
  int workingArraySize = 0;

  cmp_stiffsz_f(&numberGlobalEquations,
	   	pointerToLm,
	   	pointerToLmx,
	   	&numberVolumeElements,
	   	&workingArraySize,
	   	&totalNumberSlipperyNodes,
	   	&numberVolumeElementNodes,
	   	&errorcode,
	   	errorstring,
	   	sizeof(errorstring));

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
  int numberVolumeElements;
  int numberVolumeElementNodes;
  int numberVolumeElementEquations;
  PyObject* pyPointerToIndx;
  PyObject* pyPointerToLink;
  PyObject* pyPointerToNbrs;
  int workingArraySize;
  int totalNumberSlipperyNodes;

  int ok = PyArg_ParseTuple(args, "iOOiiiOOOii:lnklst",
			    &numberGlobalEquations,
			    &pyPointerToLm,
			    &pyPointerToLmx,
			    &numberVolumeElements,
			    &numberVolumeElementNodes,
			    &numberVolumeElementEquations,
			    &pyPointerToIndx,
			    &pyPointerToLink,
			    &pyPointerToNbrs,
			    &workingArraySize,
			    &totalNumberSlipperyNodes);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int* pointerToLm = (int*) PyCObject_AsVoidPtr(pyPointerToLm);
  int* pointerToLmx = (int*) PyCObject_AsVoidPtr(pyPointerToLmx);
  int* pointerToIndx = (int*) PyCObject_AsVoidPtr(pyPointerToIndx);
  int* pointerToLink = (int*) PyCObject_AsVoidPtr(pyPointerToLink);
  int* pointerToNbrs = (int*) PyCObject_AsVoidPtr(pyPointerToNbrs);
  int stiffnessMatrixSize =0;
  int stiffnessOffDiagonalSize =0;

  lnklst_f(&numberGlobalEquations,
	   pointerToLm,
	   pointerToLmx,
	   &numberVolumeElements,
	   &numberVolumeElementNodes,
	   &numberVolumeElementEquations,
	   pointerToIndx,
	   pointerToLink,
	   pointerToNbrs,
	   &workingArraySize,
	   &stiffnessOffDiagonalSize,
	   &stiffnessMatrixSize,
	   &totalNumberSlipperyNodes,
	   &errorcode,
	   errorstring,
	   sizeof(errorstring));

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


// Localize id array for reference by element

char pylithomop3d_local__doc__[] = "";
char pylithomop3d_local__name__[] = "local";

PyObject * pylithomop3d_local(PyObject *, PyObject *args)
{
  PyObject* pyPointerToId;
  int numberNodes;
  PyObject* pyPointerToIens;
  PyObject* pyPointerToLm;
  int numberVolumeElements;
  int numberVolumeElementNodes;

  int ok = PyArg_ParseTuple(args, "OiOOii:local",
			    &pyPointerToId,
			    &numberNodes,
			    &pyPointerToIens,
			    &pyPointerToLm,
			    &numberVolumeElements,
			    &numberVolumeElementNodes);

  if (!ok) {
    return 0;
  }

  int* pointerToId = (int*) PyCObject_AsVoidPtr(pyPointerToId);
  int* pointerToIens = (int*) PyCObject_AsVoidPtr(pyPointerToIens);
  int* pointerToLm = (int*) PyCObject_AsVoidPtr(pyPointerToLm);

  local_f(pointerToId,
	  &numberNodes,
	  pointerToIens,
	  pointerToLm,
	  &numberVolumeElements,
	  &numberVolumeElementNodes);

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberVolumeElements:" << numberVolumeElements
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}
    

// Localize nfault array for reference by element

char pylithomop3d_localf__doc__[] = "";
char pylithomop3d_localf__name__[] = "localf";

PyObject * pylithomop3d_localf(PyObject *, PyObject *args)
{
  PyObject* pyPointerToIens;
  PyObject* pyPointerToLmf;
  int numberVolumeElements;
  PyObject* pyPointerToNfault;
  int numberSplitNodeEntries;
  int numberVolumeElementNodes;

  int ok = PyArg_ParseTuple(args, "OOiOii:localf",
			    &pyPointerToIens,
			    &pyPointerToLmf,
			    &numberVolumeElements,
			    &pyPointerToNfault,
			    &numberSplitNodeEntries,
			    &numberVolumeElementNodes);

  if (!ok) {
    return 0;
  }

  int* pointerToIens = (int*) PyCObject_AsVoidPtr(pyPointerToIens);
  int* pointerToLmf = (int*) PyCObject_AsVoidPtr(pyPointerToLmf);
  int* pointerToNfault = (int*) PyCObject_AsVoidPtr(pyPointerToNfault);

  localf_f(pointerToIens,
	   pointerToLmf,
	   &numberVolumeElements,
	   pointerToNfault,
	   &numberSplitNodeEntries,
	   &numberVolumeElementNodes);
		  
  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberSplitNodeEntries:" << numberSplitNodeEntries
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}
    

// Localize idx array for reference by element

char pylithomop3d_localx__doc__[] = "";
char pylithomop3d_localx__name__[] = "localx";

PyObject * pylithomop3d_localx(PyObject *, PyObject *args)
{
  PyObject* pyPointerToIdx;
  int numberNodes;
  PyObject* pyPointerToIens;
  PyObject* pyPointerToLmx;
  int numberVolumeElements;
  PyObject* pyPointerToNslip;
  int numberSlipperyNodeEntries;
  int numberVolumeElementNodes;

  int ok = PyArg_ParseTuple(args, "OiOOiOii:localx",
			    &pyPointerToIdx,
  			    &numberNodes,
  			    &pyPointerToIens,
  			    &pyPointerToLmx,
  			    &numberVolumeElements,
  			    &pyPointerToNslip,
			    &numberSlipperyNodeEntries,
  			    &numberVolumeElementNodes);

  if (!ok) {
    return 0;
  }

  int* pointerToIdx = (int*) PyCObject_AsVoidPtr(pyPointerToIdx);
  int* pointerToIens = (int*) PyCObject_AsVoidPtr(pyPointerToIens);
  int* pointerToLmx = (int*) PyCObject_AsVoidPtr(pyPointerToLmx);
  int* pointerToNslip = (int*) PyCObject_AsVoidPtr(pyPointerToNslip);

  localx_f(pointerToIdx,
	   &numberNodes,
	   pointerToIens,
	   pointerToLmx,
	   &numberVolumeElements,
	   pointerToNslip,
	   &numberSlipperyNodeEntries,
	   &numberVolumeElementNodes);

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberVolumeElements:" << numberVolumeElements
    << journal::endl;

  // return
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
// $Id: sparse.cc,v 1.18 2005/06/07 19:39:11 willic3 Exp $

// End of file
