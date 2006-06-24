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

#include <petscmat.h>
#include <portinfo>
#include "journal/debug.h"

#include <Python.h>

#include "sorting.h"
#include "exceptionhandler.h"
#include "pylith3d_externs.h"
#include <stdio.h>
#include <string.h>


// Sort elements into element families

char pypylith3d_sort_elements__doc__[] = "";
char pypylith3d_sort_elements__name__[] = "sort_elements";

PyObject * pypylith3d_sort_elements(PyObject *, PyObject *args)
{
  PyObject* pyPointerToIen;
  PyObject* pyPointerToMat;
  PyObject* pyPointerToMaterialModelInfo;
  PyObject* pyPointerToVolumeElementFamilyList;
  PyObject* pyPointerToIvfamily;
  PyObject* pyPointerToIens;
  PyObject* pyPointerToIvftmp;
  PyObject* pyPointerToIndxiel;
  PyObject* pyPointerToIelindx;
  int numberVolumeElementNodes;
  int numberVolumeElementGaussPoints;
  int maxNumberVolumeElementFamilies;
  int numberVolumeElementFamilies;
  int prestressFlag;
  int numberVolumeElements;
  int numberNodes;

  int ok = PyArg_ParseTuple(args, "OOOOOOOOOiiiiiii:sort_elements",
			    &pyPointerToIen,
			    &pyPointerToMat,
			    &pyPointerToMaterialModelInfo,
			    &pyPointerToVolumeElementFamilyList,
			    &pyPointerToIvfamily,
			    &pyPointerToIens,
			    &pyPointerToIvftmp,
			    &pyPointerToIndxiel,
			    &pyPointerToIelindx,
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
  int* pointerToIelindx = (int*) PyCObject_AsVoidPtr(pyPointerToIelindx);
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
		  pointerToIelindx,
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
		  sizeof(errorstring));

  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
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

char pypylith3d_sort_slip_nodes__doc__[] = "";
char pypylith3d_sort_slip_nodes__name__[] = "sort_slip_nodes";

PyObject * pypylith3d_sort_slip_nodes(PyObject *, PyObject *args)
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

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numberSlipperyNodeEntries:" << numberSlipperyNodeEntries
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}
    


// Sort split nodes according to reordered elements

char pypylith3d_sort_split_nodes__doc__[] = "";
char pypylith3d_sort_split_nodes__name__[] = "sort_split_nodes";

PyObject * pypylith3d_sort_split_nodes(PyObject *, PyObject *args)
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

  journal::debug_t debug("pylith3d");
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
