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
//
#include <portinfo>
#include <Python.h>

#include "array.h"


// allocateInt
char pylithomop3d_allocateInt__doc__[] = "";
char pylithomop3d_allocateInt__name__[] = "allocateInt";

PyObject * pylithomop3d_allocateInt(PyObject *, PyObject * args)
{

  int len;
  //  int i;

  int ok = PyArg_ParseTuple(args, "i:allocateInt",&len);
  
  if (!ok) {
    return 0;
  }
  
  int * p = (int *)malloc(sizeof(int)*len);

  // Experimental portion to initialize all values to zero.
  //  for (i = 0; i < len; ++i)
      //p[i] = 0;

  // Experimental section to increment the reference count for the
  // array pointer.
  PyObject *arrayPointer;
  arrayPointer = PyCObject_FromVoidPtr(p,clearp);
    //Py_XINCREF(arrayPointer);

  // return
  return arrayPointer;

}


// allocateDouble
char pylithomop3d_allocateDouble__doc__[] = "";
char pylithomop3d_allocateDouble__name__[] = "allocateDouble";

PyObject * pylithomop3d_allocateDouble(PyObject *, PyObject * args)
{

  int len;
  // int i;

  int ok = PyArg_ParseTuple(args, "i:allocateDouble",&len);
  
  if (!ok) {
    return 0;
  }
  
  double * p = (double *)malloc(sizeof(double)*len);

  // Experimental portion to initialize all values to zero.
    //for (i = 0; i < len; ++i)
      //p[i] = 0.0;

  // Experimental section to increment the reference count for the
  // array pointer.
  PyObject *arrayPointer;
  arrayPointer = PyCObject_FromVoidPtr(p,clearp);
    //Py_XINCREF(arrayPointer);

  // return
  return arrayPointer;

}


// intListToArray
char pylithomop3d_intListToArray__doc__[] = "";
char pylithomop3d_intListToArray__name__[] = "intListToArray";

PyObject * pylithomop3d_intListToArray(PyObject *, PyObject * args)
{
  PyObject *listobj;
  PyObject *listitem;
  int len,i;

  int ok = PyArg_ParseTuple(args, "O!:intListToArray",&PyList_Type,&listobj);

  if (!ok) {
    return 0;
  }

  len = PyList_Size(listobj);
  int * p = (int *)malloc(sizeof(int)*len);
  for(i=0;i<len;i++) {
    listitem = PyList_GetItem(listobj,i);
    if(!PyInt_Check(listitem)) {  //error if list item is not int
      return 0;			  //may want better error handling
    }
    p[i] = (int)PyInt_AsLong(listitem);
  }

  // Experimental section to increment the reference count for the
  // array pointer.
  PyObject *arrayPointer;
  arrayPointer = PyCObject_FromVoidPtr(p,clearp);
    //Py_XINCREF(arrayPointer);

  // return
  return arrayPointer;

}


// doubleListToArray
char pylithomop3d_doubleListToArray__doc__[] = "";
char pylithomop3d_doubleListToArray__name__[] = "doubleListToArray";

PyObject * pylithomop3d_doubleListToArray(PyObject *, PyObject * args)
{
  PyObject *listobj;
  PyObject *listitem;
  int len,i;

  int ok = PyArg_ParseTuple(args, "O!:doubleListToArray",&PyList_Type,&listobj);

  if (!ok) {
    return 0;
  }

  len = PyList_Size(listobj);
  double * p = (double *)malloc(sizeof(double)*len);
  for(i=0;i<len;i++) {
    listitem = PyList_GetItem(listobj,i);
    if(!PyFloat_Check(listitem)) {  //error if list item is not float
      return 0;			  //may want better error handling
    }
    p[i] = PyFloat_AsDouble(listitem);
  }

  // Experimental section to increment the reference count for the
  // array pointer.
  PyObject *arrayPointer;
  arrayPointer = PyCObject_FromVoidPtr(p,clearp);
    //Py_XINCREF(arrayPointer);

  // return
  return arrayPointer;

}


// clearp 
// deconstructor
void clearp(void *p)
{
  free(p);
  return;
}


// version
// $Id: array.cc,v 1.1 2004/04/14 21:24:47 willic3 Exp $

// End of file
