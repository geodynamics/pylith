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

#include "numbering.h"
#include "exceptionhandler.h"
#include "lithomop3d_externs.h"
#include <stdio.h>
#include <string.h>


// Create Id and Idx arrays, and determine the number of equations.

char pylithomop3d_create_id__doc__[] = "";
char pylithomop3d_create_id__name__[] = "create_id";

PyObject * pylithomop3d_create_id(PyObject *, PyObject *args)
{
  PyObject* pyPointerToId;
  PyObject* pyPointerToIdx;
  PyObject* pyPointerToIbond;
  PyObject* pyPointerToNslip;
  PyObject* pyPointerToIdslp;
  int numberSlipperyNodeEntries;
  int numberNodes;
  int totalNumberSlipperyNodes;

  int ok = PyArg_ParseTuple(args, "OOOOOiii:create_id",
			    &pyPointerToId,
			    &pyPointerToIdx,
			    &pyPointerToIbond,
			    &pyPointerToNslip,
			    &pyPointerToIdslp,
			    &numberSlipperyNodeEntries,
			    &numberNodes,
			    &totalNumberSlipperyNodes);

  if (!ok) {
    return 0;
  }

  int* pointerToId = (int*) PyCObject_AsVoidPtr(pyPointerToId);
  int* pointerToIdx = (int*) PyCObject_AsVoidPtr(pyPointerToIdx);
  int* pointerToIbond = (int*) PyCObject_AsVoidPtr(pyPointerToIbond);
  int* pointerToNslip = (int*) PyCObject_AsVoidPtr(pyPointerToNslip);
  int* pointerToIdslp = (int*) PyCObject_AsVoidPtr(pyPointerToIdslp);
  int numberGlobalEquations = 0;

  create_id_f(pointerToId,
	  pointerToIdx,
	  pointerToIbond,
	  pointerToNslip,
	  pointerToIdslp,
	  &numberSlipperyNodeEntries,
	  &numberNodes,
	  &totalNumberSlipperyNodes,
	  &numberGlobalEquations);

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberGlobalEquations:" << numberGlobalEquations
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue("i", numberGlobalEquations);
}


// Form id array for split nodes

char pylithomop3d_id_split__doc__[] = "";
char pylithomop3d_id_split__name__[] = "id_split";

PyObject * pylithomop3d_id_split(PyObject *, PyObject *args)
{
  PyObject* pyPointerToNfault;
  PyObject* pyPointerToIdftn;
  int numberNodes;
  int numberSplitNodeEntries;
  int totalNumberSplitNodes;
  int f77PlotOutput;
  int plotOutputInt;
  char* plotOutputFile;

  int ok = PyArg_ParseTuple(args, "OOiiiiis:id_split",
			    &pyPointerToNfault,
			    &pyPointerToIdftn,
			    &numberNodes,
			    &numberSplitNodeEntries,
			    &totalNumberSplitNodes,
			    &f77PlotOutput,
			    &plotOutputInt,
			    &plotOutputFile);

  if (!ok) {
    return 0;
  }

  int* pointerToNfault = (int*) PyCObject_AsVoidPtr(pyPointerToNfault);
  int* pointerToIdftn = (int*) PyCObject_AsVoidPtr(pyPointerToIdftn);

  id_split_f(pointerToNfault,
	     pointerToIdftn,
	     &numberNodes,
	     &numberSplitNodeEntries,
	     &totalNumberSplitNodes,
	     &f77PlotOutput,
	     &plotOutputInt,
	     plotOutputFile,
	     strlen(plotOutputFile));

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberSplitNodeEntries:" << numberSplitNodeEntries
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}
    

// Localize id array for reference by element

char pylithomop3d_local__doc__[] = "";
char pylithomop3d_local__name__[] = "local";

PyObject * pylithomop3d_local(PyObject *, PyObject *args)
{
  PyObject* pyPointerToId;
  int numberNodes;
  PyObject* pyPointerToIen;
  PyObject* pyPointerToLm;
  PyObject* pyPointerToInfiel;
  int connectivitySize;
  int numberElements;
  PyObject* pyPointerToElementTypeInfo;

  int ok = PyArg_ParseTuple(args, "OiOOOiiO:local",
			    &pyPointerToId,
			    &numberNodes,
			    &pyPointerToIen,
			    &pyPointerToLm,
			    &pyPointerToInfiel,
			    &connectivitySize,
			    &numberElements,
			    &pyPointerToElementTypeInfo);

  if (!ok) {
    return 0;
  }

  int* pointerToId = (int*) PyCObject_AsVoidPtr(pyPointerToId);
  int* pointerToIen = (int*) PyCObject_AsVoidPtr(pyPointerToIen);
  int* pointerToLm = (int*) PyCObject_AsVoidPtr(pyPointerToLm);
  int* pointerToInfiel = (int*) PyCObject_AsVoidPtr(pyPointerToInfiel);
  int* pointerToElementTypeInfo = (int*) PyCObject_AsVoidPtr(pyPointerToElementTypeInfo);

  local_f(pointerToId,
	  &numberNodes,
	  pointerToIen,
	  pointerToLm,
	  pointerToInfiel,
	  &connectivitySize,
	  &numberElements,
	  pointerToElementTypeInfo);

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "connectivitySize:" << connectivitySize
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
  PyObject* pyPointerToIen;
  PyObject* pyPointerToLmf;
  PyObject* pyPointerToInfiel;
  int connectivitySize;
  int numberElements;
  PyObject* pyPointerToElementTypeInfo;
  PyObject* pyPointerToNfault;
  int numberSplitNodeEntries;

  int ok = PyArg_ParseTuple(args, "OOOiiOOi:localf",
			    &pyPointerToIen,
			    &pyPointerToLmf,
			    &pyPointerToInfiel,
			    &connectivitySize,
			    &numberElements,
			    &pyPointerToElementTypeInfo,
			    &pyPointerToNfault,
			    &numberSplitNodeEntries);

  if (!ok) {
    return 0;
  }

  int* pointerToIen = (int*) PyCObject_AsVoidPtr(pyPointerToIen);
  int* pointerToLmf = (int*) PyCObject_AsVoidPtr(pyPointerToLmf);
  int* pointerToInfiel = (int*) PyCObject_AsVoidPtr(pyPointerToInfiel);
  int* pointerToElementTypeInfo = (int*) PyCObject_AsVoidPtr(pyPointerToElementTypeInfo);
  int* pointerToNfault = (int*) PyCObject_AsVoidPtr(pyPointerToNfault);

  localf_f(pointerToIen,
	   pointerToLmf,
	   pointerToInfiel,
	   &connectivitySize,
	   &numberElements,
	   pointerToElementTypeInfo,
	   pointerToNfault,
	   &numberSplitNodeEntries);
		  
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
  PyObject* pyPointerToIen;
  PyObject* pyPointerToLmx;
  PyObject* pyPointerToInfiel;
  int connectivitySize;
  int numberElements;
  PyObject* pyPointerToElementTypeInfo;
  PyObject* pyPointerToNslip;
  int numberSlipperyNodeEntries;

  int ok = PyArg_ParseTuple(args, "OiOOOiiOOi:localx",
			    &pyPointerToIdx,
  			    &numberNodes,
  			    &pyPointerToIen,
  			    &pyPointerToLmx,
  			    &pyPointerToInfiel,
  			    &connectivitySize,
  			    &numberElements,
  			    &pyPointerToElementTypeInfo,
  			    &pyPointerToNslip,
			    &numberSlipperyNodeEntries);

  if (!ok) {
    return 0;
  }

  int* pointerToIdx = (int*) PyCObject_AsVoidPtr(pyPointerToIdx);
  int* pointerToIen = (int*) PyCObject_AsVoidPtr(pyPointerToIen);
  int* pointerToLmx = (int*) PyCObject_AsVoidPtr(pyPointerToLmx);
  int* pointerToInfiel = (int*) PyCObject_AsVoidPtr(pyPointerToInfiel);
  int* pointerToElementTypeInfo = (int*) PyCObject_AsVoidPtr(pyPointerToElementTypeInfo);
  int* pointerToNslip = (int*) PyCObject_AsVoidPtr(pyPointerToNslip);

  localx_f(pointerToIdx,
	   &numberNodes,
	   pointerToIen,
	   pointerToLmx,
	   pointerToInfiel,
	   &connectivitySize,
	   &numberElements,
	   pointerToElementTypeInfo,
	   pointerToNslip,
	   &numberSlipperyNodeEntries);

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberElements:" << numberElements
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


// Find closest fault neighbors for slippery nodes

char pylithomop3d_nfind__doc__[] = "";
char pylithomop3d_nfind__name__[] = "nfind";

PyObject * pylithomop3d_nfind(PyObject *, PyObject *args)
{
  PyObject* pyPointerToX;
  PyObject* pyPointerToXtmp;
  PyObject* pyPointerToIdslp;
  PyObject* pyPointerToIpslp;
  PyObject* pyPointerToItmp;
  PyObject* pyPointerToItmp1;
  PyObject* pyPointerToItmp2;
  PyObject* pyPointerToNslip;
  int numberSlipperyNodeEntries;
  int totalNumberSlipperyNodes;
  int numberNodes;

  int ok = PyArg_ParseTuple(args, "OOOOOOOOiii:nfind",
			    &pyPointerToX,
			    &pyPointerToXtmp,
			    &pyPointerToIdslp,
			    &pyPointerToIpslp,
			    &pyPointerToItmp,
			    &pyPointerToItmp1,
			    &pyPointerToItmp2,
			    &pyPointerToNslip,
			    &numberSlipperyNodeEntries,
			    &totalNumberSlipperyNodes,
			    &numberNodes);

  if (!ok) {
    return 0;
  }

  double* pointerToX = (double*) PyCObject_AsVoidPtr(pyPointerToX);
  double* pointerToXtmp = (double*) PyCObject_AsVoidPtr(pyPointerToXtmp);
  int* pointerToIdslp = (int*) PyCObject_AsVoidPtr(pyPointerToIdslp);
  int* pointerToIpslp = (int*) PyCObject_AsVoidPtr(pyPointerToIpslp);
  int* pointerToItmp = (int*) PyCObject_AsVoidPtr(pyPointerToItmp);
  int* pointerToItmp1 = (int*) PyCObject_AsVoidPtr(pyPointerToItmp1);
  int* pointerToItmp2 = (int*) PyCObject_AsVoidPtr(pyPointerToItmp2);
  int* pointerToNslip = (int*) PyCObject_AsVoidPtr(pyPointerToNslip);

  nfind_f(pointerToX,
	  pointerToXtmp,
	  pointerToIdslp,
	  pointerToIpslp,
	  pointerToItmp,
	  pointerToItmp1,
	  pointerToItmp2,
	  pointerToNslip,
	  &numberSlipperyNodeEntries,
	  &totalNumberSlipperyNodes,
	  &numberNodes);

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "pointerToX:" << pointerToX
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}

// version
// $Id: numbering.cc,v 1.1 2005/03/12 02:03:18 willic3 Exp $

// End of file
