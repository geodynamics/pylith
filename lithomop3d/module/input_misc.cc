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
  int numberVolumeElements;
  int quadratureOrderInt;
  int prestressAutoComputeInt;
  int prestressAutoChangeElasticPropsInt;
  double prestressAutoComputePoisson;
  double prestressAutoComputeYoungs;
  int f77AsciiOutput;
  int asciiOutputInt;
  char* asciiOutputFile;

  int ok = PyArg_ParseTuple(args, "iiiiddiis:write_element_info",
			    &numberVolumeElements,
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

  write_element_info_f(&numberVolumeElements,
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
    << "numberVolumeElements:" << numberVolumeElements
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
  int debuggingOutputInt;
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
			    &debuggingOutputInt,
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
		      &debuggingOutputInt,
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
    << "debuggingOutputInt:" << debuggingOutputInt
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
  PyObject* pyPointerToIvfamily;
  PyObject* pyPointerToMaterialModelInfo;
  int numberVolumeElementFamilies;
  int propertySize;
  int asciiOutputInt;
  int plotOutputInt;
  int f77AsciiOutput;
  int f77PlotOutput;
  char* asciiOutputFile;
  char* plotOutputFile;

  int ok = PyArg_ParseTuple(args, "OOOOiiiiiiss:write_props",
			    &pyPointerToListArrayPropertyList,
			    &pyPointerToListArrayGrav,
			    &pyPointerToIvfamily,
			    &pyPointerToMaterialModelInfo,
			    &numberVolumeElementFamilies,
			    &propertySize,
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
  const int maxsize = 4096;
  char errorstring[maxsize];
  double* pointerToListArrayPropertyList = (double*) PyCObject_AsVoidPtr(pyPointerToListArrayPropertyList);
  double* pointerToListArrayGrav = (double*) PyCObject_AsVoidPtr(pyPointerToListArrayGrav);
  int* pointerToIvfamily = (int*) PyCObject_AsVoidPtr(pyPointerToIvfamily);
  int* pointerToMaterialModelInfo = (int*) PyCObject_AsVoidPtr(pyPointerToMaterialModelInfo);

  write_props_f(pointerToListArrayPropertyList,
		pointerToListArrayGrav,
		pointerToIvfamily,
		pointerToMaterialModelInfo,
		&numberVolumeElementFamilies,
		&propertySize,
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
    << "numberVolumeElementFamilies:" << numberVolumeElementFamilies
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
  int usePreviousDisplacementFlag;
  int f77AsciiOutput;
  int asciiOutputInt;
  char* asciiOutputFile;

  int ok = PyArg_ParseTuple(args, "iiis:write_subiter",
			    &usePreviousDisplacementFlag,
			    &f77AsciiOutput,
			    &asciiOutputInt,
			    &asciiOutputFile);

  if (!ok) {
    return 0;
  }


  write_subiter_f(&usePreviousDisplacementFlag,
		  &f77AsciiOutput,
		  &asciiOutputInt,
		  asciiOutputFile,strlen(asciiOutputFile));
		  
  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "usePreviousDisplacementFlag:" << usePreviousDisplacementFlag
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
  PyObject* pyPointerToIvfamily;
  int numberVolumeElements;
  int numberVolumeElementFamilies;
  PyObject* pyPointerToSh;
  int numberVolumeElementNodes;
  int numberVolumeElementGaussPoints;
  int volumeElementType;
  PyObject* pyPointerToIstatout;
  PyObject* pyPointerToNstatout;
  int f77UcdOutput;
  char* ucdOutputRoot;

  int ok = PyArg_ParseTuple(args, "OiOOiiOiiiOOis:write_ucd_mesh",
			    &pyPointerToX,
			    &numberNodes,
			    &pyPointerToIen,
			    &pyPointerToIvfamily,
			    &numberVolumeElements,
			    &numberVolumeElementFamilies,
			    &pyPointerToSh,
			    &numberVolumeElementNodes,
			    &numberVolumeElementGaussPoints,
			    &volumeElementType,
			    &pyPointerToIstatout,
			    &pyPointerToNstatout,
			    &f77UcdOutput,
			    &ucdOutputRoot);

  if (!ok) {
    return 0;
  }

  double* pointerToX = (double*) PyCObject_AsVoidPtr(pyPointerToX);
  int* pointerToIen = (int*) PyCObject_AsVoidPtr(pyPointerToIen);
  int* pointerToIvfamily = (int*) PyCObject_AsVoidPtr(pyPointerToIvfamily);
  double* pointerToSh = (double*) PyCObject_AsVoidPtr(pyPointerToSh);
  int* pointerToIstatout = (int*) PyCObject_AsVoidPtr(pyPointerToIstatout);
  int* pointerToNstatout = (int*) PyCObject_AsVoidPtr(pyPointerToNstatout);

  write_ucd_mesh_f(pointerToX,
		   &numberNodes,
		   pointerToIen,
		   pointerToIvfamily,
		   &numberVolumeElements,
		   &numberVolumeElementFamilies,
		   pointerToSh,
		   &numberVolumeElementNodes,
		   &numberVolumeElementGaussPoints,
		   &volumeElementType,
		   pointerToIstatout,
		   pointerToNstatout,
		   &f77UcdOutput,
		   ucdOutputRoot,strlen(ucdOutputRoot));
		  
  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberVolumeElements:" << numberVolumeElements
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}



// version
// $Id: input_misc.cc,v 1.8 2005/03/31 23:27:57 willic3 Exp $

// End of file
