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

#include "setup.h"
#include "exceptionhandler.h"
#include "lithomop3d_externs.h"
#include <stdio.h>
#include <string.h>


// Initialize material model info

char pylithomop3d_matmod_def__doc__[] = "";
char pylithomop3d_matmod_def__name__[] = "matmod_def";

PyObject * pylithomop3d_matmod_def(PyObject *, PyObject *args)
{
  PyObject* pyPointerToMaterialModelInfo;
  PyObject* pyPointerToMaterialModelStates;

  int ok = PyArg_ParseTuple(args, "OO:matmod_def",
			    &pyPointerToMaterialModelInfo,
			    &pyPointerToMaterialModelStates);

  if (!ok) {
    return 0;
  }

  int* pointerToMaterialModelInfo = (int*) PyCObject_AsVoidPtr(pyPointerToMaterialModelInfo);
  int* pointerToMaterialModelStates = (int*) PyCObject_AsVoidPtr(pyPointerToMaterialModelStates);

  matmod_def_f(pointerToMaterialModelInfo,
	       pointerToMaterialModelStates);
    
  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "pointerToMaterialModelInfo:" << pointerToMaterialModelInfo
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}

// Precompute shape function info

char pylithomop3d_preshape__doc__[] = "";
char pylithomop3d_preshape__name__[] = "preshape";

PyObject * pylithomop3d_preshape(PyObject *, PyObject *args)
{
  PyObject* pyPointerToSh;
  PyObject* pyPointerToShj;
  PyObject* pyPointerToGauss;
  PyObject* pyPointerToElementTypeInfo;
  int quadratureOrderInt;

  int ok = PyArg_ParseTuple(args, "OOOOi:preshape",
			    &pyPointerToSh,
			    &pyPointerToShj,
			    &pyPointerToGauss,
			    &pyPointerToElementTypeInfo,
			    &quadratureOrderInt);

  if (!ok) {
    return 0;
  }

  double* pointerToSh = (double*) PyCObject_AsVoidPtr(pyPointerToSh);
  double* pointerToShj = (double*) PyCObject_AsVoidPtr(pyPointerToShj);
  double* pointerToGauss = (double*) PyCObject_AsVoidPtr(pyPointerToGauss);
  int* pointerToElementTypeInfo = (int*) PyCObject_AsVoidPtr(pyPointerToElementTypeInfo);

  preshape_f(pointerToSh,
	     pointerToShj,
	     pointerToGauss,
	     pointerToElementTypeInfo,
	     &quadratureOrderInt);
    
  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "quadratureOrderInt:" << quadratureOrderInt
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}
    
// version
// $Id: setup.cc,v 1.1 2004/07/19 19:07:25 willic3 Exp $

// End of file
