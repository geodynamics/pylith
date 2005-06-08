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

#include "sorting.h"
#include "exceptionhandler.h"
#include "lithomop3d_externs.h"
#include <stdio.h>
#include <string.h>


// Sort elements into element families

char pylithomop3d_sort_elements__doc__[] = "";
char pylithomop3d_sort_elements__name__[] = "sort_elements";

PyObject * pylithomop3d_sort_elements(PyObject *, PyObject *args)
{
  PyObject* pyPointerToIen;
  PyObject* pyPointerToMat;
  PyObject* pyPointerToMaterialModelInfo;
  PyObject* pyPointerToVolumeElementFamilyList;
  PyObject* pyPointerToIvfamily;
  PyObject* pyPointerToIens;
  PyObject* pyPointerToIvftmp;
  PyObject* pyPointerToIndxiel;
  int numberVolumeElementNodes;
  int numberVolumeElementGaussPoints;
  int maxNumberVolumeElementFamilies;
  int numberVolumeElementFamilies;
  int prestressFlag;
  int numberVolumeElements;
  int numberNodes;

  int ok = PyArg_ParseTuple(args, "OOOOOOOOiiiiiii:sort_elements",
			    &pyPointerToIen,
			    &pyPointerToMat,
			    &pyPointerToMaterialModelInfo,
			    &pyPointerToVolumeElementFamilyList,
			    &pyPointerToIvfamily,
			    &pyPointerToIens,
			    &pyPointerToIvftmp,
			    &pyPointerToIndxiel,
			    &numberVolumeElementNodes,
			    &numberVolumeElementGaussPoints,
			    &maxNumberVolumeElementFamilies,
			    &numberVolumeElementFamilies,
			    &prestressFlag,
			    &numberVolumeElements,
			    &numberNodes);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int* pointerToIen = (int*) PyCObject_AsVoidPtr(pyPointerToIen);
  int* pointerToMat = (int*) PyCObject_AsVoidPtr(pyPointerToMat);
  int* pointerToMaterialModelInfo = (int*) PyCObject_AsVoidPtr(pyPointerToMaterialModelInfo);
  int* pointerToVolumeElementFamilyList = (int*) PyCObject_AsVoidPtr(pyPointerToVolumeElementFamilyList);
  int* pointerToIvfamily = (int*) PyCObject_AsVoidPtr(pyPointerToIvfamily);
  int* pointerToIens = (int*) PyCObject_AsVoidPtr(pyPointerToIens);
  int* pointerToIvftmp = (int*) PyCObject_AsVoidPtr(pyPointerToIvftmp);
  int* pointerToIndxiel = (int*) PyCObject_AsVoidPtr(pyPointerToIndxiel);
  int stateSize = 0;
  int state0Size = 0;
  int propertySize = 0;

  sort_elements_f(pointerToIen,
		  pointerToMat,
		  pointerToMaterialModelInfo,
		  pointerToVolumeElementFamilyList,
		  pointerToIvfamily,
		  pointerToIens,
		  pointerToIvftmp,
		  pointerToIndxiel,
		  &numberVolumeElementNodes,
		  &numberVolumeElementGaussPoints,
		  &maxNumberVolumeElementFamilies,
		  &numberVolumeElementFamilies,
		  &prestressFlag,
		  &numberVolumeElements,
		  &numberNodes,
		  &stateSize,
		  &state0Size,
		  &propertySize,
		  &errorcode,
		  errorstring,
		  strlen(errorstring));

  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "stateSize:" << stateSize
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue("iii",
		       stateSize,
		       state0Size,
		       propertySize);
}

// Sort slippery nodes according to reordered elements

char pylithomop3d_sort_slip_nodes__doc__[] = "";
char pylithomop3d_sort_slip_nodes__name__[] = "sort_slip_nodes";

PyObject * pylithomop3d_sort_slip_nodes(PyObject *, PyObject *args)
{
  PyObject* pyPointerToNslip;
  PyObject* pyPointerToIndxiel;
  int numberSlipperyNodeEntries;
  int numberVolumeElements;

  int ok = PyArg_ParseTuple(args, "OOii:sort_slip_nodes",
			    &pyPointerToNslip,
			    &pyPointerToIndxiel,
			    &numberSlipperyNodeEntries,
			    &numberVolumeElements);

  if (!ok) {
    return 0;
  }

  int* pointerToNslip = (int*) PyCObject_AsVoidPtr(pyPointerToNslip);
  int* pointerToIndxiel = (int*) PyCObject_AsVoidPtr(pyPointerToIndxiel);

  sort_slip_nodes_f(pointerToNslip,
		    pointerToIndxiel,
		    &numberSlipperyNodeEntries,
		    &numberVolumeElements);

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberSlipperyNodeEntries:" << numberSlipperyNodeEntries
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}
    


// Sort split nodes according to reordered elements

char pylithomop3d_sort_split_nodes__doc__[] = "";
char pylithomop3d_sort_split_nodes__name__[] = "sort_split_nodes";

PyObject * pylithomop3d_sort_split_nodes(PyObject *, PyObject *args)
{
  PyObject* pyPointerToNfault;
  PyObject* pyPointerToIndxiel;
  int numberSplitNodeEntries;
  int numberVolumeElements;

  int ok = PyArg_ParseTuple(args, "OOii:sort_split_nodes",
			    &pyPointerToNfault,
			    &pyPointerToIndxiel,
			    &numberSplitNodeEntries,
			    &numberVolumeElements);

  if (!ok) {
    return 0;
  }

  int* pointerToNfault = (int*) PyCObject_AsVoidPtr(pyPointerToNfault);
  int* pointerToIndxiel = (int*) PyCObject_AsVoidPtr(pyPointerToIndxiel);

  sort_split_nodes_f(pointerToNfault,
		    pointerToIndxiel,
		    &numberSplitNodeEntries,
		    &numberVolumeElements);

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberSplitNodeEntries:" << numberSplitNodeEntries
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}

// version
// $Id: sorting.cc,v 1.2 2005/06/07 19:39:11 willic3 Exp $

// End of file
