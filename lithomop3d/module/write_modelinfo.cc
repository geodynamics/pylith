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

#include "write_modelinfo.h"
#include "exceptionhandler.h"
#include "lithomop3d_externs.h"
#include <stdio.h>
#include <string.h>


// Write boundary conditions

char pylithomop3d_write_bc__doc__[] = "";
char pylithomop3d_write_bc__name__[] = "write_bc";

PyObject * pylithomop3d_write_bc(PyObject *, PyObject *args)
{
  PyObject* pyPointerToBond;
  PyObject* pyPointerToIbond;
  int numberNodes;
  int f77AsciiOutput;
  int asciiOutputInt;
  char* asciiOutputFile;

  int ok = PyArg_ParseTuple(args, "OOiiis:write_bc",
			    &pyPointerToBond,
			    &pyPointerToIbond,
			    &numberNodes,
			    &f77AsciiOutput,
			    &asciiOutputInt,
			    &asciiOutputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  double* pointerToBond = (double*) PyCObject_AsVoidPtr(pyPointerToBond);
  int* pointerToIbond = (int*) PyCObject_AsVoidPtr(pyPointerToIbond);

  write_bc_f(pointerToBond,
	     pointerToIbond,
	     &numberNodes,
	     &f77AsciiOutput,
	     &asciiOutputInt,
	     asciiOutputFile,
	     &errorcode,
	     errorstring,
	     strlen(asciiOutputFile),
	     strlen(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberNodes:" << numberNodes
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


// Write connectivities

char pylithomop3d_write_connect__doc__[] = "";
char pylithomop3d_write_connect__name__[] = "write_connect";

PyObject * pylithomop3d_write_connect(PyObject *, PyObject *args)
{
  PyObject* pyPointerToIens;
  PyObject* pyPointerToIvfamily;
  PyObject* pyPointerToIndxiel;
  int numberVolumeElementNodes;
  int numberVolumeElementGaussPoints;
  int numberVolumeElements;
  int volumeElementType;
  int numberVolumeElementFamilies;
  int f77AsciiOutput;
  int f77PlotOutput;
  int asciiOutputInt;
  int plotOutputInt;
  char* asciiOutputFile;
  char* plotOutputFile;

  int ok = PyArg_ParseTuple(args, "OOOiiiiiiiiiss:write_connect",
                            &pyPointerToIens,
                            &pyPointerToIvfamily,
                            &pyPointerToIndxiel,
                            &numberVolumeElementNodes,
                            &numberVolumeElementGaussPoints,
                            &numberVolumeElements,
                            &volumeElementType,
                            &numberVolumeElementFamilies,
                            &f77AsciiOutput,
                            &f77PlotOutput,
                            &asciiOutputInt,
                            &plotOutputInt,
                            &asciiOutputFile,
                            &plotOutputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int* pointerToIens = (int*) PyCObject_AsVoidPtr(pyPointerToIens);
  int* pointerToIvfamily = (int*) PyCObject_AsVoidPtr(pyPointerToIvfamily);
  int* pointerToIndxiel = (int*) PyCObject_AsVoidPtr(pyPointerToIndxiel);

  write_connect_f(pointerToIens,
		  pointerToIvfamily,
		  pointerToIndxiel,
		  &numberVolumeElementNodes,
		  &numberVolumeElementGaussPoints,
		  &numberVolumeElements,
		  &volumeElementType,
		  &numberVolumeElementFamilies,
		  &f77AsciiOutput,
		  &f77PlotOutput,
		  &asciiOutputInt,
		  &plotOutputInt,
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
    << "numberVolumeElements:" << numberVolumeElements
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


// Write coordinates

char pylithomop3d_write_coords__doc__[] = "";
char pylithomop3d_write_coords__name__[] = "write_coords";

PyObject * pylithomop3d_write_coords(PyObject *, PyObject *args)
{
  PyObject* pyPointerToX;
  int numberNodes;
  int f77AsciiOutput;
  int f77PlotOutput;
  int asciiOutputInt;
  int plotOutputInt;
  char* asciiOutputFile;
  char* plotOutputFile;

  int ok = PyArg_ParseTuple(args, "Oiiiiiss:write_coords",
			    &pyPointerToX,
			    &numberNodes,
			    &f77AsciiOutput,
			    &f77PlotOutput,
			    &asciiOutputInt,
			    &plotOutputInt,
			    &asciiOutputFile,
			    &plotOutputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  double* pointerToX = (double*) PyCObject_AsVoidPtr(pyPointerToX);

  write_coords_f(pointerToX,
		 &numberNodes,
		 &f77AsciiOutput,
		 &f77PlotOutput,
		 &asciiOutputInt,
		 &plotOutputInt,
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
    << "numberNodes:" << numberNodes
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


// Write differential forces

char pylithomop3d_write_diff__doc__[] = "";
char pylithomop3d_write_diff__name__[] = "write_diff";

PyObject * pylithomop3d_write_diff(PyObject *, PyObject *args)
{
  PyObject* pyPointerToDiforc;
  PyObject* pyPointerToNslip;
  PyObject* pyPointerToIdhist;
  int numberSlipperyNodeEntries;
  int numberDifferentialForceEntries;
  int numberNodes;
  int f77AsciiOutput;
  int asciiOutputInt;
  char* asciiOutputFile;

  int ok = PyArg_ParseTuple(args, "OOOiiiiis:write_diff",
			    &pyPointerToDiforc,
			    &pyPointerToNslip,
			    &pyPointerToIdhist,
			    &numberSlipperyNodeEntries,
			    &numberDifferentialForceEntries,
			    &numberNodes,
			    &f77AsciiOutput,
			    &asciiOutputInt,
			    &asciiOutputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  double* pointerToDiforc = (double*) PyCObject_AsVoidPtr(pyPointerToDiforc);
  int* pointerToNslip = (int*) PyCObject_AsVoidPtr(pyPointerToNslip);
  int* pointerToIdhist = (int*) PyCObject_AsVoidPtr(pyPointerToIdhist);

  write_diff_f(pointerToDiforc,
	       pointerToNslip,
	       pointerToIdhist,
	       &numberSlipperyNodeEntries,
	       &numberDifferentialForceEntries,
	       &numberNodes,
	       &f77AsciiOutput,
	       &asciiOutputInt,
	       asciiOutputFile,
	       &errorcode,
	       errorstring,
	       strlen(asciiOutputFile),
	       strlen(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberDifferentialForceEntries:" << numberDifferentialForceEntries
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}

// Write out element info

char pylithomop3d_write_element_info__doc__[] = "";
char pylithomop3d_write_element_info__name__[] = "write_element_info";

PyObject * pylithomop3d_write_element_info(PyObject *, PyObject *args)
{
  int numberVolumeElements;
  int numberVolumeElementNodes;
  int numberVolumeElementGaussPoints;
  int volumeElementType;
  int quadratureOrderInt;
  int prestressAutoComputeInt;
  int prestressAutoChangeElasticPropsInt;
  double prestressAutoComputePoisson;
  double prestressAutoComputeYoungs;
  int f77AsciiOutput;
  int asciiOutputInt;
  char* asciiOutputFile;

  int ok = PyArg_ParseTuple(args, "iiiiiiiddiis:write_element_info",
			    &numberVolumeElements,
			    &numberVolumeElementNodes,
			    &numberVolumeElementGaussPoints,
			    &volumeElementType,
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
		       &numberVolumeElementNodes,
		       &numberVolumeElementGaussPoints,
		       &volumeElementType,
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


// Write time steps at which full output is desired

char pylithomop3d_write_fuldat__doc__[] = "";
char pylithomop3d_write_fuldat__name__[] = "write_fuldat";

PyObject * pylithomop3d_write_fuldat(PyObject *, PyObject *args)
{
  PyObject* pyPointerToIprint;
  int numberFullOutputs;
  int analysisTypeInt;
  int numberCycles;
  int totalNumberTimeSteps;
  int f77AsciiOutput;
  int f77PlotOutput;
  int asciiOutputInt;
  int plotOutputInt;
  char* asciiOutputFile;
  char* plotOutputFile;

  int ok = PyArg_ParseTuple(args, "Oiiiiiiiiss:write_fuldat",
			    &pyPointerToIprint,
			    &numberFullOutputs,
			    &analysisTypeInt,
			    &numberCycles,
			    &totalNumberTimeSteps,
			    &f77AsciiOutput,
			    &f77PlotOutput,
			    &asciiOutputInt,
			    &plotOutputInt,
			    &asciiOutputFile,
			    &plotOutputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int* pointerToIprint = (int*) PyCObject_AsVoidPtr(pyPointerToIprint);

  write_fuldat_f(pointerToIprint,
		 &numberFullOutputs,
		 &analysisTypeInt,
		 &numberCycles,
		 &totalNumberTimeSteps,
		 &f77AsciiOutput,
		 &f77PlotOutput,
		 &asciiOutputInt,
		 &plotOutputInt,
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
    << "numberFullOutputs:" << numberFullOutputs
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


// Write load histories

char pylithomop3d_write_hist__doc__[] = "";
char pylithomop3d_write_hist__name__[] = "write_hist";

PyObject * pylithomop3d_write_hist(PyObject *, PyObject *args)
{
  PyObject* pyPointerToHistry;
  PyObject* pyPointerToTimes;
  int numberLoadHistories;
  int totalNumberTimeSteps;
  int f77AsciiOutput;
  int asciiOutputInt;
  char* asciiOutputFile;

  int ok = PyArg_ParseTuple(args, "OOiiiis:write_hist",
			    &pyPointerToHistry,
			    &pyPointerToTimes,
			    &numberLoadHistories,
			    &totalNumberTimeSteps,
			    &f77AsciiOutput,
			    &asciiOutputInt,
			    &asciiOutputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  double* pointerToHistry = (double*) PyCObject_AsVoidPtr(pyPointerToHistry);
  double* pointerToTimes = (double*) PyCObject_AsVoidPtr(pyPointerToTimes);

  write_hist_f(pointerToHistry,
	       pointerToTimes,
	       &numberLoadHistories,
	       &totalNumberTimeSteps,
	       &f77AsciiOutput,
	       &asciiOutputInt,
	       asciiOutputFile,
	       &errorcode,
	       errorstring,
	       strlen(asciiOutputFile),
	       strlen(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberLoadHistories:" << numberLoadHistories
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


// Write local coordinate rotations

char pylithomop3d_write_skew__doc__[] = "";
char pylithomop3d_write_skew__name__[] = "write_skew";

PyObject * pylithomop3d_write_skew(PyObject *, PyObject *args)
{
  PyObject* pyPointerToSkew;
  int numberRotationEntries;
  int autoRotateSlipperyNodesInt;
  int numberNodes;
  int f77AsciiOutput;
  int asciiOutputInt;
  char* asciiOutputFile;

  int ok = PyArg_ParseTuple(args, "Oiiiiis:write_skew",
			    &pyPointerToSkew,
			    &numberRotationEntries,
			    &autoRotateSlipperyNodesInt,
			    &numberNodes,
			    &f77AsciiOutput,
			    &asciiOutputInt,
			    &asciiOutputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  double* pointerToSkew = (double*) PyCObject_AsVoidPtr(pyPointerToSkew);

  write_skew_f(pointerToSkew,
	       &numberRotationEntries,
	       &autoRotateSlipperyNodesInt,
	       &numberNodes,
	       &f77AsciiOutput,
	       &asciiOutputInt,
	       asciiOutputFile,
	       &errorcode,
	       errorstring,
	       strlen(asciiOutputFile),
	       strlen(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberRotationEntries:" << numberRotationEntries
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


// Write slippery node entries

char pylithomop3d_write_slip__doc__[] = "";
char pylithomop3d_write_slip__name__[] = "write_slip";

PyObject * pylithomop3d_write_slip(PyObject *, PyObject *args)
{
  PyObject* pyPointerToNslip;
  int numberSlipperyNodeEntries;
  int totalNumberSlipperyNodes;
  int f77AsciiOutput;
  int f77PlotOutput;
  int asciiOutputInt;
  int plotOutputInt;
  char* asciiOutputFile;
  char* plotOutputFile;

  int ok = PyArg_ParseTuple(args, "Oiiiiiiss:write_slip",
			    &pyPointerToNslip,
			    &numberSlipperyNodeEntries,
			    &totalNumberSlipperyNodes,
			    &f77AsciiOutput,
			    &f77PlotOutput,
			    &asciiOutputInt,
			    &plotOutputInt,
			    &asciiOutputFile,
			    &plotOutputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int* pointerToNslip = (int*) PyCObject_AsVoidPtr(pyPointerToNslip);

  write_slip_f(pointerToNslip,
	       &numberSlipperyNodeEntries,
	       &totalNumberSlipperyNodes,
	       &f77AsciiOutput,
	       &f77PlotOutput,
	       &asciiOutputInt,
	       &plotOutputInt,
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
    << "totalNumberSlipperyNodes:" << totalNumberSlipperyNodes
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


// Write split node entries

char pylithomop3d_write_split__doc__[] = "";
char pylithomop3d_write_split__name__[] = "write_split";

PyObject * pylithomop3d_write_split(PyObject *, PyObject *args)
{
  PyObject* pyPointerToFault;
  PyObject* pyPointerToNfault;
  int numberSplitNodeEntries;
  int f77AsciiOutput;
  int f77PlotOutput;
  int asciiOutputInt;
  int plotOutputInt;
  char* asciiOutputFile;
  char* plotOutputFile;

  int ok = PyArg_ParseTuple(args, "OOiiiiiss:write_split",
			    &pyPointerToFault,
			    &pyPointerToNfault,
			    &numberSplitNodeEntries,
			    &f77AsciiOutput,
			    &f77PlotOutput,
			    &asciiOutputInt,
			    &plotOutputInt,
			    &asciiOutputFile,
			    &plotOutputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  double* pointerToFault = (double*) PyCObject_AsVoidPtr(pyPointerToFault);
  int* pointerToNfault = (int*) PyCObject_AsVoidPtr(pyPointerToNfault);

  write_split_f(pointerToFault,
		pointerToNfault,
		&numberSplitNodeEntries,
		&f77AsciiOutput,
		&f77PlotOutput,
		&asciiOutputInt,
		&plotOutputInt,
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
    << "numberSplitNodeEntries:" << numberSplitNodeEntries
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


// Write split node entries for plot output

char pylithomop3d_write_split_plot__doc__[] = "";
char pylithomop3d_write_split_plot__name__[] = "write_split_plot";

PyObject * pylithomop3d_write_split_plot(PyObject *, PyObject *args)
{
  PyObject* pyPointerToIdftn;
  int totalNumberSplitNodes;
  int f77PlotOutput;
  int plotOutputInt;
  char* plotOutputFile;

  int ok = PyArg_ParseTuple(args, "Oiiis:write_split_plot",
			    &pyPointerToIdftn,
			    &totalNumberSplitNodes,
			    &f77PlotOutput,
			    &plotOutputInt,
			    &plotOutputFile);

  if (!ok) {
    return 0;
  }

  int* pointerToIdftn = (int*) PyCObject_AsVoidPtr(pyPointerToIdftn);

  write_split_plot_f(pointerToIdftn,
		     &totalNumberSplitNodes,
		     &f77PlotOutput,
		     &plotOutputInt,
		     plotOutputFile,
		     strlen(plotOutputFile));
    

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "totalNumberSplitNodes:" << totalNumberSplitNodes
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


// Write state output information

char pylithomop3d_write_stateout__doc__[] = "";
char pylithomop3d_write_stateout__name__[] = "write_stateout";

PyObject * pylithomop3d_write_stateout(PyObject *, PyObject *args)
{
  PyObject* pyPointerToIstatout;
  PyObject* pyPointerToNstatout;
  int f77AsciiOutput;
  int f77PlotOutput;
  int asciiOutputInt;
  int plotOutputInt;
  char* asciiOutputFile;
  char* plotOutputFile;

  int ok = PyArg_ParseTuple(args, "OOiiiiss:write_stateout",
			    &pyPointerToIstatout,
			    &pyPointerToNstatout,
			    &f77AsciiOutput,
			    &f77PlotOutput,
			    &asciiOutputInt,
			    &plotOutputInt,
			    &asciiOutputFile,
			    &plotOutputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int* pointerToIstatout = (int*) PyCObject_AsVoidPtr(pyPointerToIstatout);
  int* pointerToNstatout = (int*) PyCObject_AsVoidPtr(pyPointerToNstatout);

  write_stateout_f(pointerToIstatout,
		   pointerToNstatout,
		   &f77AsciiOutput,
		   &f77PlotOutput,
		   &asciiOutputInt,
		   &plotOutputInt,
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
    << "pointerToIstatout:" << pointerToIstatout
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


// Write time step data

char pylithomop3d_write_timdat__doc__[] = "";
char pylithomop3d_write_timdat__name__[] = "write_timdat";

PyObject * pylithomop3d_write_timdat(PyObject *, PyObject *args)
{
  PyObject* pyPointerToDelt;
  PyObject* pyPointerToAlfa;
  PyObject* pyPointerToUtol;
  PyObject* pyPointerToFtol;
  PyObject* pyPointerToEtol;
  PyObject* pyPointerToTimes;
  PyObject* pyPointerToMaxstp;
  PyObject* pyPointerToMaxit;
  PyObject* pyPointerToNtdinit;
  PyObject* pyPointerToLgdef;
  PyObject* pyPointerToItmax;
  int numberTimeStepGroups;
  int totalNumberTimeSteps;
  int f77AsciiOutput;
  int asciiOutputInt;
  char* asciiOutputFile;

  int ok = PyArg_ParseTuple(args, "OOOOOOOOOOOiiiis:write_timdat",
			    &pyPointerToDelt,
			    &pyPointerToAlfa,
			    &pyPointerToUtol,
			    &pyPointerToFtol,
			    &pyPointerToEtol,
			    &pyPointerToTimes,
			    &pyPointerToMaxstp,
			    &pyPointerToMaxit,
			    &pyPointerToNtdinit,
			    &pyPointerToLgdef,
			    &pyPointerToItmax,
			    &numberTimeStepGroups,
			    &totalNumberTimeSteps,
			    &f77AsciiOutput,
			    &asciiOutputInt,
			    &asciiOutputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  double* pointerToDelt = (double*) PyCObject_AsVoidPtr(pyPointerToDelt);
  double* pointerToAlfa = (double*) PyCObject_AsVoidPtr(pyPointerToAlfa);
  double* pointerToUtol = (double*) PyCObject_AsVoidPtr(pyPointerToUtol);
  double* pointerToFtol = (double*) PyCObject_AsVoidPtr(pyPointerToFtol);
  double* pointerToEtol = (double*) PyCObject_AsVoidPtr(pyPointerToEtol);
  double* pointerToTimes = (double*) PyCObject_AsVoidPtr(pyPointerToTimes);
  int* pointerToMaxstp = (int*) PyCObject_AsVoidPtr(pyPointerToMaxstp);
  int* pointerToMaxit = (int*) PyCObject_AsVoidPtr(pyPointerToMaxit);
  int* pointerToNtdinit = (int*) PyCObject_AsVoidPtr(pyPointerToNtdinit);
  int* pointerToLgdef = (int*) PyCObject_AsVoidPtr(pyPointerToLgdef);
  int* pointerToItmax = (int*) PyCObject_AsVoidPtr(pyPointerToItmax);

  write_timdat_f(pointerToDelt,
		 pointerToAlfa,
		 pointerToUtol,
		 pointerToFtol,
		 pointerToEtol,
		 pointerToTimes,
		 pointerToMaxstp,
		 pointerToMaxit,
		 pointerToNtdinit,
		 pointerToLgdef,
		 pointerToItmax,
		 &numberTimeStepGroups,
		 &totalNumberTimeSteps,
		 &f77AsciiOutput,
		 &asciiOutputInt,
		 asciiOutputFile,
		 &errorcode,
		 errorstring,
		 strlen(asciiOutputFile),
		 strlen(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "totalNumberTimeSteps:" << totalNumberTimeSteps
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
  PyObject* pyPointerToIens;
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
  int ucdOutputInt;
  char* ucdOutputRoot;

  int ok = PyArg_ParseTuple(args, "OiOOiiOiiiOOiis:write_ucd_mesh",
			    &pyPointerToX,
			    &numberNodes,
			    &pyPointerToIens,
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
			    &ucdOutputInt,
			    &ucdOutputRoot);

  if (!ok) {
    return 0;
  }

  double* pointerToX = (double*) PyCObject_AsVoidPtr(pyPointerToX);
  int* pointerToIens = (int*) PyCObject_AsVoidPtr(pyPointerToIens);
  int* pointerToIvfamily = (int*) PyCObject_AsVoidPtr(pyPointerToIvfamily);
  double* pointerToSh = (double*) PyCObject_AsVoidPtr(pyPointerToSh);
  int* pointerToIstatout = (int*) PyCObject_AsVoidPtr(pyPointerToIstatout);
  int* pointerToNstatout = (int*) PyCObject_AsVoidPtr(pyPointerToNstatout);

  write_ucd_mesh_f(pointerToX,
		   &numberNodes,
		   pointerToIens,
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
		   &ucdOutputInt,
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


// Write winkler BC

char pylithomop3d_write_wink__doc__[] = "";
char pylithomop3d_write_wink__name__[] = "write_wink";

PyObject * pylithomop3d_write_wink(PyObject *, PyObject *args)
{
  PyObject* pyPointerToWinkdef;
  PyObject* pyPointerToIwinkdef;
  PyObject* pyPointerToIwinkid;
  int numberWinklerEntries;
  int f77AsciiOutput;
  int asciiOutputInt;
  char* asciiOutputFile;

  int ok = PyArg_ParseTuple(args, "OOOiiis:write_wink",
			    &pyPointerToWinkdef,
			    &pyPointerToIwinkdef,
			    &pyPointerToIwinkid,
			    &numberWinklerEntries,
			    &f77AsciiOutput,
			    &asciiOutputInt,
			    &asciiOutputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  double* pointerToWinkdef = (double*) PyCObject_AsVoidPtr(pyPointerToWinkdef);
  int* pointerToIwinkdef = (int*) PyCObject_AsVoidPtr(pyPointerToIwinkdef);
  int* pointerToIwinkid = (int*) PyCObject_AsVoidPtr(pyPointerToIwinkid);

  write_wink_f(pointerToWinkdef,
	       pointerToIwinkdef,
	       pointerToIwinkid,
	       &numberWinklerEntries,
	       &f77AsciiOutput,
	       &asciiOutputInt,
	       asciiOutputFile,
	       &errorcode,
	       errorstring,
	       strlen(asciiOutputFile),
	       strlen(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberWinklerEntries:" << numberWinklerEntries
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


// Write slippery winkler BC

char pylithomop3d_write_winkx__doc__[] = "";
char pylithomop3d_write_winkx__name__[] = "write_winkx";

PyObject * pylithomop3d_write_winkx(PyObject *, PyObject *args)
{
  PyObject* pyPointerToWinkxdef;
  PyObject* pyPointerToIwinkxdef;
  PyObject* pyPointerToIwinkxid;
  int numberSlipperyWinklerEntries;
  int f77AsciiOutput;
  int asciiOutputInt;
  char* asciiOutputFile;

  int ok = PyArg_ParseTuple(args, "OOOiiis:write_winkx",
			    &pyPointerToWinkxdef,
			    &pyPointerToIwinkxdef,
			    &pyPointerToIwinkxid,
			    &numberSlipperyWinklerEntries,
			    &f77AsciiOutput,
			    &asciiOutputInt,
			    &asciiOutputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  double* pointerToWinkxdef = (double*) PyCObject_AsVoidPtr(pyPointerToWinkxdef);
  int* pointerToIwinkxdef = (int*) PyCObject_AsVoidPtr(pyPointerToIwinkxdef);
  int* pointerToIwinkxid = (int*) PyCObject_AsVoidPtr(pyPointerToIwinkxid);

  write_winkx_f(pointerToWinkxdef,
		pointerToIwinkxdef,
		pointerToIwinkxid,
		&numberSlipperyWinklerEntries,
		&f77AsciiOutput,
		&asciiOutputInt,
		asciiOutputFile,
		&errorcode,
		errorstring,
		strlen(asciiOutputFile),
		strlen(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberSlipperyWinklerEntries:" << numberSlipperyWinklerEntries
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}

// version
// $Id: write_modelinfo.cc,v 1.5 2005/06/24 20:17:43 willic3 Exp $

// End of file
