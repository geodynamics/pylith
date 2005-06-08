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

#include "numbering.h"
#include "exceptionhandler.h"
#include "lithomop3d_externs.h"
#include <stdio.h>
#include <string.h>



// Assign equation numbers to Iwink and Iwinkx arrays, and compact
// Wink and Winkx arrays to correspond to assigned BC.

char pylithomop3d_assign_wink__doc__[] = "";
char pylithomop3d_assign_wink__name__[] = "assign_wink";

PyObject * pylithomop3d_assign_wink(PyObject *, PyObject *args)
{
  PyObject* pyPointerToWinkdef;
  PyObject* pyPointerToWink;
  PyObject* pyPointerToIwinkdef;
  PyObject* pyPointerToIwinkid;
  PyObject* pyPointerToIwink;
  PyObject* pyPointerToId;
  int numberNodes;
  int numberWinklerForces;
  int numberWinklerEntries;

  int ok = PyArg_ParseTuple(args, "OOOOOOiii:assign_wink",
			    &pyPointerToWinkdef,
			    &pyPointerToWink,
			    &pyPointerToIwinkdef,
			    &pyPointerToIwinkid,
			    &pyPointerToIwink,
			    &pyPointerToId,
			    &numberNodes,
			    &numberWinklerForces,
			    &numberWinklerEntries);

  if (!ok) {
    return 0;
  }

  double* pointerToWinkdef = (double*) PyCObject_AsVoidPtr(pyPointerToWinkdef);
  double* pointerToWink = (double*) PyCObject_AsVoidPtr(pyPointerToWink);
  int* pointerToIwinkdef = (int*) PyCObject_AsVoidPtr(pyPointerToIwinkdef);
  int* pointerToIwinkid = (int*) PyCObject_AsVoidPtr(pyPointerToIwinkid);
  int* pointerToIwink = (int*) PyCObject_AsVoidPtr(pyPointerToIwink);
  int* pointerToId = (int*) PyCObject_AsVoidPtr(pyPointerToId);

  assign_wink_f(pointerToWinkdef,
		pointerToWink,
		pointerToIwinkdef,
		pointerToIwinkid,
		pointerToIwink,
		pointerToId,
		&numberNodes,
		&numberWinklerForces,
		&numberWinklerEntries);

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberWinklerForces:" << numberWinklerForces
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}
    

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

  int ok = PyArg_ParseTuple(args, "OOiii:id_split",
			    &pyPointerToNfault,
			    &pyPointerToIdftn,
			    &numberNodes,
			    &numberSplitNodeEntries,
			    &totalNumberSplitNodes);

  if (!ok) {
    return 0;
  }

  int* pointerToNfault = (int*) PyCObject_AsVoidPtr(pyPointerToNfault);
  int* pointerToIdftn = (int*) PyCObject_AsVoidPtr(pyPointerToIdftn);

  id_split_f(pointerToNfault,
	     pointerToIdftn,
	     &numberNodes,
	     &numberSplitNodeEntries,
	     &totalNumberSplitNodes);

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberSplitNodeEntries:" << numberSplitNodeEntries
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
// $Id: numbering.cc,v 1.4 2005/06/07 19:39:11 willic3 Exp $

// End of file
