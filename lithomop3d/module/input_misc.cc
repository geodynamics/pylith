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

#include "input_misc.h"
#include "exceptionhandler.h"
#include "lithomop3d_externs.h"
#include <stdio.h>
#include <string.h>


// Write out element info

char pylithomop3d_write_element_info__doc__[] = "";
char pylithomop3d_write_element_info__name__[] = "write_element_info";

PyObject * pylithomop3d_write_element_info(PyObject *, PyObject *args)
{
  int numberElements;
  int quadratureOrderInt;
  int prestressAutoComputeInt;
  int prestressAutoChangeElasticPropsInt;
  double prestressAutoComputePoisson;
  double prestressAutoComputeYoungs;
  int f77AsciiOutput;
  int asciiOutputInt;
  char* asciiOutputFile;

  int ok = PyArg_ParseTuple(args, "iiiiddiis:write_element_info",
			    &numberElements,
			    &quadratureOrderInt,
			    &prestressAutoComputeInt,
			    &prestressAutoChangeElasticPropsInt,
			    &prestressAutoComputePoisson,
			    &prestressAutoComputeYoungs,
			    &f77AsciiOutput,
			    &asciiOutputInt,
			    &asciiOutputFile);

  if (!ok) {
    return 0;
  }

  write_element_info_f(&numberElements,
		       &quadratureOrderInt,
		       &prestressAutoComputeInt,
		       &prestressAutoChangeElasticPropsInt,
		       &prestressAutoComputePoisson,
		       &prestressAutoComputeYoungs,
		       &f77AsciiOutput,
		       &asciiOutputInt,
		       asciiOutputFile,strlen(asciiOutputFile));

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberElements:" << numberElements
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


// Write out global info

char pylithomop3d_write_global_info__doc__[] = "";
char pylithomop3d_write_global_info__name__[] = "write_global_info";

PyObject * pylithomop3d_write_global_info(PyObject *, PyObject *args)
{
  char* title;
  int asciiOutputInt;
  int plotOutputInt;
  int numberNodes;
  int analysisTypeInt;
  int debuggingOutput;
  int f77AsciiOutput;
  int f77PlotOutput;
  char* asciiOutputFile;
  char* plotOutputFile;

  int ok = PyArg_ParseTuple(args, "siiiiiiiss:write_global_info",
			    &title,
			    &asciiOutputInt,
			    &plotOutputInt,
			    &numberNodes,
			    &analysisTypeInt,
			    &debuggingOutput,
			    &f77AsciiOutput,
			    &f77PlotOutput,
			    &asciiOutputFile,
			    &plotOutputFile);

  if (!ok) {
    return 0;
  }

  write_global_info_f(title,
		      &asciiOutputInt,
		      &plotOutputInt,
		      &numberNodes,
		      &analysisTypeInt,
		      &debuggingOutput,
		      &f77AsciiOutput,
		      &f77PlotOutput,
		      asciiOutputFile,
		      plotOutputFile,
		      strlen(title),
		      strlen(asciiOutputFile),
		      strlen(plotOutputFile));

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "debuggingOutput:" << debuggingOutput
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


// Write out material property info

char pylithomop3d_write_props__doc__[] = "";
char pylithomop3d_write_props__name__[] = "write_props";

PyObject * pylithomop3d_write_props(PyObject *, PyObject *args)
{
  PyObject* pyPointerToListArrayPropertyList;
  PyObject* pyPointerToListArrayGrav;
  PyObject* pyPointerToMaterialInfo;
  PyObject* pyPointerToMaterialModelInfo;
  int numberMaterials;
  int propertyListSize;
  int asciiOutputInt;
  int plotOutputInt;
  int f77AsciiOutput;
  int f77PlotOutput;
  char* asciiOutputFile;
  char* plotOutputFile;

  int ok = PyArg_ParseTuple(args, "OOOOiiiiiiss:write_props",
			    &pyPointerToListArrayPropertyList,
			    &pyPointerToListArrayGrav,
			    &pyPointerToMaterialInfo,
			    &pyPointerToMaterialModelInfo,
			    &numberMaterials,
			    &propertyListSize,
			    &asciiOutputInt,
			    &plotOutputInt,
			    &f77AsciiOutput,
			    &f77PlotOutput,
			    &asciiOutputFile,
			    &plotOutputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 1024;
  char errorstring[maxsize];
  double* pointerToListArrayPropertyList = (double*) PyCObject_AsVoidPtr(pyPointerToListArrayPropertyList);
  double* pointerToListArrayGrav = (double*) PyCObject_AsVoidPtr(pyPointerToListArrayGrav);
  int* pointerToMaterialInfo = (int*) PyCObject_AsVoidPtr(pyPointerToMaterialInfo);
  int* pointerToMaterialModelInfo = (int*) PyCObject_AsVoidPtr(pyPointerToMaterialModelInfo);

  write_props_f(pointerToListArrayPropertyList,
		pointerToListArrayGrav,
		pointerToMaterialInfo,
		pointerToMaterialModelInfo,
		&numberMaterials,
		&propertyListSize,
		&asciiOutputInt,
		&plotOutputInt,
		&f77AsciiOutput,
		&f77PlotOutput,
		asciiOutputFile,
		plotOutputFile,
		&errorcode,
		errorstring,
		strlen(asciiOutputFile),
		strlen(plotOutputFile),
		strlen(errorstring));

  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberMaterials:" << numberMaterials
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}

// Write out sparse matrix info

char pylithomop3d_write_sparse_info__doc__[] = "";
char pylithomop3d_write_sparse_info__name__[] = "write_sparse_info";

PyObject * pylithomop3d_write_sparse_info(PyObject *, PyObject *args)
{
  int numberGlobalEquations;
  int stiffnessMatrixSize;
  int minimumNonzeroTermsPerRow;
  int maximumNonzeroTermsPerRow;
  double averageNonzeroTermsPerRow;
  int asciiOutputInt;
  int f77AsciiOutput;
  char* asciiOutputFile;

  int ok = PyArg_ParseTuple(args, "iiiidiis:write_sparse_info",
			    &numberGlobalEquations,
			    &stiffnessMatrixSize,
			    &minimumNonzeroTermsPerRow,
			    &maximumNonzeroTermsPerRow,
			    &averageNonzeroTermsPerRow,
			    &asciiOutputInt,
			    &f77AsciiOutput,
			    &asciiOutputFile);

  //debug printf("Hello from pylithomop3d_write_sparse_info!\n");
  //debug printf("numberGlobalEquations = %d\n", numberGlobalEquations);
  //debug printf("stiffnessMatrixSize = %d\n", stiffnessMatrixSize);
  //debug printf("minimumNonzeroTermsPerRow = %d\n", minimumNonzeroTermsPerRow);
  //debug printf("maximumNonzeroTermsPerRow = %d\n", maximumNonzeroTermsPerRow);
  //debug printf("averageNonzeroTermsPerRow = %g\n", averageNonzeroTermsPerRow);
  //debug printf("asciiOutputInt = %d\n", asciiOutputInt);
  //debug printf("f77AsciiOutput = %d\n", f77AsciiOutput);
  //debug printf("asciiOutputFile = %s\n", asciiOutputFile);
  if (!ok) {
    return 0;
  }

  write_sparse_info_f(&numberGlobalEquations,
		      &stiffnessMatrixSize,
		      &minimumNonzeroTermsPerRow,
		      &maximumNonzeroTermsPerRow,
		      &averageNonzeroTermsPerRow,
		      &asciiOutputInt,
		      &f77AsciiOutput,
		      asciiOutputFile,strlen(asciiOutputFile));
		  
  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "minimumNonzeroTermsPerRow:" << minimumNonzeroTermsPerRow
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


// Write out stress integration parameters

char pylithomop3d_write_strscomp__doc__[] = "";
char pylithomop3d_write_strscomp__name__[] = "write_strscomp";

PyObject * pylithomop3d_write_strscomp(PyObject *, PyObject *args)
{
  double stressTolerance;
  double minimumStrainPerturbation;
  double initialStrainPerturbation;
  int f77AsciiOutput;
  int asciiOutputInt;
  char* asciiOutputFile;

  int ok = PyArg_ParseTuple(args, "dddiis:write_strscomp",
			    &stressTolerance,
			    &minimumStrainPerturbation,
			    &initialStrainPerturbation,
			    &f77AsciiOutput,
			    &asciiOutputInt,
			    &asciiOutputFile);

  if (!ok) {
    return 0;
  }

  write_strscomp_f(&stressTolerance,
		   &minimumStrainPerturbation,
		   &initialStrainPerturbation,
		   &f77AsciiOutput,
		   &asciiOutputInt,
		   asciiOutputFile,strlen(asciiOutputFile));

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "initialStrainPerturbation:" << initialStrainPerturbation
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


// Write out subiteration convergence parameters

char pylithomop3d_write_subiter__doc__[] = "";
char pylithomop3d_write_subiter__name__[] = "write_subiter";

PyObject * pylithomop3d_write_subiter(PyObject *, PyObject *args)
{
  PyObject* pyPointerToListArrayRmult;
  PyObject* pyPointerToListArrayRmin;
  int preconditionerTypeInt;
  int maxPcgIterations;
  int f77AsciiOutput;
  int asciiOutputInt;
  char* asciiOutputFile;

  int ok = PyArg_ParseTuple(args, "OOiiiis:write_subiter",
			    &pyPointerToListArrayRmult,
			    &pyPointerToListArrayRmin,
			    &preconditionerTypeInt,
			    &maxPcgIterations,
			    &f77AsciiOutput,
			    &asciiOutputInt,
			    &asciiOutputFile);

  if (!ok) {
    return 0;
  }

  double* pointerToListArrayRmult = (double*) PyCObject_AsVoidPtr(pyPointerToListArrayRmult);
  double* pointerToListArrayRmin = (double*) PyCObject_AsVoidPtr(pyPointerToListArrayRmin);

  write_subiter_f(pointerToListArrayRmult,
		  pointerToListArrayRmin,
		  &preconditionerTypeInt,
		  &maxPcgIterations,
		  &f77AsciiOutput,
		  &asciiOutputInt,
		  asciiOutputFile,strlen(asciiOutputFile));
		  
  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "maxPcgIterations:" << maxPcgIterations
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


// Write mesh info to UCD file

char pylithomop3d_write_ucd_mesh__doc__[] = "";
char pylithomop3d_write_ucd_mesh__name__[] = "write_ucd_mesh";

PyObject * pylithomop3d_write_ucd_mesh(PyObject *, PyObject *args)
{
  PyObject* pyPointerToX;
  int numberNodes;
  PyObject* pyPointerToIen;
  PyObject* pyPointerToInfiel;
  int numberElements;
  int connectivitySize;
  PyObject* pyPointerToSh;
  PyObject* pyPointerToElementTypeInfo;
  PyObject* pyPointerToIstatout;
  int f77UcdOutput;
  char* ucdOutputRoot;

  int ok = PyArg_ParseTuple(args, "OiOOiiOOOis:write_ucd_mesh",
			    &pyPointerToX,
			    &numberNodes,
			    &pyPointerToIen,
			    &pyPointerToInfiel,
			    &numberElements,
			    &connectivitySize,
			    &pyPointerToSh,
			    &pyPointerToElementTypeInfo,
			    &pyPointerToIstatout,
			    &f77UcdOutput,
			    &ucdOutputRoot);

  if (!ok) {
    return 0;
  }

  double* pointerToX = (double*) PyCObject_AsVoidPtr(pyPointerToX);
  int* pointerToIen = (int*) PyCObject_AsVoidPtr(pyPointerToIen);
  int* pointerToInfiel = (int*) PyCObject_AsVoidPtr(pyPointerToInfiel);
  double* pointerToSh = (double*) PyCObject_AsVoidPtr(pyPointerToSh);
  int* pointerToElementTypeInfo = (int*) PyCObject_AsVoidPtr(pyPointerToElementTypeInfo);
  int* pointerToIstatout = (int*) PyCObject_AsVoidPtr(pyPointerToIstatout);

  write_ucd_mesh_f(pointerToX,
		   &numberNodes,
		   pointerToIen,
		   pointerToInfiel,
		   &numberElements,
		   &connectivitySize,
		   pointerToSh,
		   pointerToElementTypeInfo,
		   pointerToIstatout,
		   &f77UcdOutput,
		   ucdOutputRoot,strlen(ucdOutputRoot));
		  
  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberElements:" << numberElements
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}



// version
// $Id: input_misc.cc,v 1.7 2005/03/12 02:03:18 willic3 Exp $

// End of file
