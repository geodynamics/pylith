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

#include "parser.h"
#include "exceptionhandler.h"
#include "lithomop3d_externs.h"
#include <stdio.h>
#include <string.h>


// Read boundary conditions

char pylithomop3d_read_bc__doc__[] = "";
char pylithomop3d_read_bc__name__[] = "read_bc";

PyObject * pylithomop3d_read_bc(PyObject *, PyObject *args)
{
  PyObject* pyPointerToBond;
  double displacementScaleFactor;
  double velocityScaleFactor;
  double forceScaleFactor;
  PyObject* pyPointerToIbond;
  PyObject* pyPointerToId;
  int numberNodes;
  int numberBcEntries;
  int f77FileInput;
  int f77AsciiOutput;
  int asciiOutputInt;
  char* bcInputFile;
  char* asciiOutputFile;

  int ok = PyArg_ParseTuple(args, "OdddOOiiiiiss:read_bc",
			    &pyPointerToBond,
			    &displacementScaleFactor,
			    &velocityScaleFactor,
			    &forceScaleFactor,
			    &pyPointerToIbond,
			    &pyPointerToId,
			    &numberNodes,
			    &numberBcEntries,
			    &f77FileInput,
			    &f77AsciiOutput,
			    &asciiOutputInt,
			    &bcInputFile,
			    &asciiOutputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 1024;
  char errorstring[maxsize];
  double* pointerToBond = (double*) PyCObject_AsVoidPtr(pyPointerToBond);
  int* pointerToIbond = (int*) PyCObject_AsVoidPtr(pyPointerToIbond);
  int* pointerToId = (int*) PyCObject_AsVoidPtr(pyPointerToId);
  int numberGlobalEquations = 0;
  int numberConcForces = 0;

  read_bc_f(pointerToBond,
	    &displacementScaleFactor,
	    &velocityScaleFactor,
	    &forceScaleFactor,
	    pointerToIbond,
	    pointerToId,
	    &numberNodes,
	    &numberBcEntries,
	    &numberGlobalEquations,
	    &numberConcForces,
	    &f77FileInput,
	    &f77AsciiOutput,
	    &asciiOutputInt,
	    bcInputFile,
	    asciiOutputFile,
	    &errorcode,
	    errorstring,
	    strlen(bcInputFile),
	    strlen(asciiOutputFile),
	    strlen(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberBcEntries:" << numberBcEntries
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue("ii", numberGlobalEquations,
		  numberConcForces);
}


// Read connectivities

char pylithomop3d_read_connect__doc__[] = "";
char pylithomop3d_read_connect__name__[] = "read_connect";

PyObject * pylithomop3d_read_connect(PyObject *, PyObject *args)
{
  PyObject* pyPointerToListArrayNumberElementNodesBase;
  PyObject* pyPointerToElementTypeInfo;
  PyObject* pyPointerToMaterialInfo;
  PyObject* pyPointerToMaterialModelInfo;
  PyObject* pyPointerToIen;
  PyObject* pyPointerToInfiel;
  PyObject* pyPointerToIndmat;
  PyObject* pyPointerToImgrp;
  int connectivitySize;
  int prestressFlag;
  int numberElements;
  int numberNodes;
  int numberMaterials;
  int f77FileInput;
  int f77AsciiOutput;
  int f77PlotOutput;
  int asciiOutputInt;
  int plotOutputInt;
  char* connectivityInputFile;
  char* asciiOutputFile;
  char* plotOutputFile;

  int ok = PyArg_ParseTuple(args, "OOOOOOOOiiiiiiiiiisss:read_connect",
			    &pyPointerToListArrayNumberElementNodesBase,
			    &pyPointerToElementTypeInfo,
			    &pyPointerToMaterialInfo,
			    &pyPointerToMaterialModelInfo,
			    &pyPointerToIen,
			    &pyPointerToInfiel,
			    &pyPointerToIndmat,
			    &pyPointerToImgrp,
			    &connectivitySize,
			    &prestressFlag,
			    &numberElements,
			    &numberNodes,
			    &numberMaterials,
			    &f77FileInput,
			    &f77AsciiOutput,
			    &f77PlotOutput,
			    &asciiOutputInt,
			    &plotOutputInt,
			    &connectivityInputFile,
			    &asciiOutputFile,
			    &plotOutputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 1024;
  char errorstring[maxsize];
  int* pointerToListArrayNumberElementNodesBase = (int*) PyCObject_AsVoidPtr(pyPointerToListArrayNumberElementNodesBase);
  int* pointerToElementTypeInfo = (int*) PyCObject_AsVoidPtr(pyPointerToElementTypeInfo);
  int* pointerToMaterialInfo = (int*) PyCObject_AsVoidPtr(pyPointerToMaterialInfo);
  int* pointerToMaterialModelInfo = (int*) PyCObject_AsVoidPtr(pyPointerToMaterialModelInfo);
  int* pointerToIen = (int*) PyCObject_AsVoidPtr(pyPointerToIen);
  int* pointerToInfiel = (int*) PyCObject_AsVoidPtr(pyPointerToInfiel);
  int* pointerToIndmat = (int*) PyCObject_AsVoidPtr(pyPointerToIndmat);
  int* pointerToImgrp = (int*) PyCObject_AsVoidPtr(pyPointerToImgrp);
  int stateSize = 0;
  int dmatSize = 0;
  int totalVolumeGaussPoints = 0;
  int state0Size = 0;

  read_connect_f(pointerToListArrayNumberElementNodesBase,
		 pointerToElementTypeInfo,
		 pointerToMaterialInfo,
		 pointerToMaterialModelInfo,
		 pointerToIen,
		 pointerToInfiel,
		 pointerToIndmat,
		 pointerToImgrp,
		 &connectivitySize,
		 &prestressFlag,
		 &numberElements,
		 &totalVolumeGaussPoints,
		 &numberNodes,
		 &numberMaterials,
		 &stateSize,
		 &state0Size,
		 &dmatSize,
		 &f77FileInput,
		 &f77AsciiOutput,
		 &f77PlotOutput,
		 &asciiOutputInt,
		 &plotOutputInt,
		 connectivityInputFile,
		 asciiOutputFile,
		 plotOutputFile,
		 &errorcode,
		 errorstring,
		 strlen(connectivityInputFile),
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
  return Py_BuildValue("iiii",
		  	stateSize,
		       	dmatSize,
		       	totalVolumeGaussPoints,
		       	state0Size);
}


// Read coordinates

char pylithomop3d_read_coords__doc__[] = "";
char pylithomop3d_read_coords__name__[] = "read_coords";

PyObject * pylithomop3d_read_coords(PyObject *, PyObject *args)
{
  PyObject* pyPointerToX;
  double coordinateScaleFactor;
  int numberNodes;
  int f77FileInput;
  int f77AsciiOutput;
  int f77PlotOutput;
  int asciiOutputInt;
  int plotOutputInt;
  char* coordinateInputFile;
  char* asciiOutputFile;
  char* plotOutputFile;

  int ok = PyArg_ParseTuple(args, "Odiiiiiisss:read_coords",
			    &pyPointerToX,
			    &coordinateScaleFactor,
			    &numberNodes,
			    &f77FileInput,
			    &f77AsciiOutput,
			    &f77PlotOutput,
			    &asciiOutputInt,
			    &plotOutputInt,
			    &coordinateInputFile,
			    &asciiOutputFile,
			    &plotOutputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 1024;
  char errorstring[maxsize];
  double* pointerToX = (double*) PyCObject_AsVoidPtr(pyPointerToX);

  read_coords_f(pointerToX,
		&coordinateScaleFactor,
		&numberNodes,
		&f77FileInput,
		&f77AsciiOutput,
		&f77PlotOutput,
		&asciiOutputInt,
		&plotOutputInt,
		coordinateInputFile,
		asciiOutputFile,
		plotOutputFile,
		&errorcode,
		errorstring,
		strlen(coordinateInputFile),
		strlen(asciiOutputFile),
		strlen(plotOutputFile),
		strlen(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "coordinateScaleFactor:" << coordinateScaleFactor
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


// Read differential forces

char pylithomop3d_read_diff__doc__[] = "";
char pylithomop3d_read_diff__name__[] = "read_diff";

PyObject * pylithomop3d_read_diff(PyObject *, PyObject *args)
{
  PyObject* pyPointerToDiforc;
  PyObject* pyPointerToNslip;
  PyObject* pyPointerToIdhist;
  int numberSlipperyNodeEntries;
  int numberDifferentialForceEntries;
  int numberNodes;
  int f77FileInput;
  int f77AsciiOutput;
  int asciiOutputInt;
  char* differentialForceInputFile;
  char* asciiOutputFile;

  int ok = PyArg_ParseTuple(args, "OOOiiiiiiss:read_diff",
			    &pyPointerToDiforc,
			    &pyPointerToNslip,
			    &pyPointerToIdhist,
			    &numberSlipperyNodeEntries,
			    &numberDifferentialForceEntries,
			    &numberNodes,
			    &f77FileInput,
			    &f77AsciiOutput,
			    &asciiOutputInt,
			    &differentialForceInputFile,
			    &asciiOutputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 1024;
  char errorstring[maxsize];
  double* pointerToDiforc = (double*) PyCObject_AsVoidPtr(pyPointerToDiforc);
  int* pointerToNslip = (int*) PyCObject_AsVoidPtr(pyPointerToNslip);
  int* pointerToIdhist = (int*) PyCObject_AsVoidPtr(pyPointerToIdhist);

  read_diff_f(pointerToDiforc,
	      pointerToNslip,
	      pointerToIdhist,
	      &numberSlipperyNodeEntries,
	      &numberDifferentialForceEntries,
	      &numberNodes,
	      &f77FileInput,
	      &f77AsciiOutput,
	      &asciiOutputInt,
	      differentialForceInputFile,
	      asciiOutputFile,
	      &errorcode,
	      errorstring,
	      strlen(differentialForceInputFile),
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


// Read time steps at which full output is desired

char pylithomop3d_read_fuldat__doc__[] = "";
char pylithomop3d_read_fuldat__name__[] = "read_fuldat";

PyObject * pylithomop3d_read_fuldat(PyObject *, PyObject *args)
{
  PyObject* pyPointerToIprint;
  int numberFullOutputs;
  int analysisTypeInt;
  int numberCycles;
  int totalNumberTimeSteps;
  int f77FileInput;
  int f77AsciiOutput;
  int f77PlotOutput;
  int asciiOutputInt;
  int plotOutputInt;
  char* fullOutputInputFile;
  char* asciiOutputFile;
  char* plotOutputFile;

  int ok = PyArg_ParseTuple(args, "Oiiiiiiiiisss:read_fuldat",
			    &pyPointerToIprint,
			    &numberFullOutputs,
			    &analysisTypeInt,
			    &numberCycles,
			    &totalNumberTimeSteps,
			    &f77FileInput,
			    &f77AsciiOutput,
			    &f77PlotOutput,
			    &asciiOutputInt,
			    &plotOutputInt,
			    &fullOutputInputFile,
			    &asciiOutputFile,
			    &plotOutputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 1024;
  char errorstring[maxsize];
  int* pointerToIprint = (int*) PyCObject_AsVoidPtr(pyPointerToIprint);

  read_fuldat_f(pointerToIprint,
		&numberFullOutputs,
		&analysisTypeInt,
		&numberCycles,
		&totalNumberTimeSteps,
		&f77FileInput,
		&f77AsciiOutput,
		&f77PlotOutput,
		&asciiOutputInt,
		&plotOutputInt,
		fullOutputInputFile,
		asciiOutputFile,
		plotOutputFile,
		&errorcode,
		errorstring,
		strlen(fullOutputInputFile),
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


// Read load histories

char pylithomop3d_read_hist__doc__[] = "";
char pylithomop3d_read_hist__name__[] = "read_hist";

PyObject * pylithomop3d_read_hist(PyObject *, PyObject *args)
{
  PyObject* pyPointerToHistry;
  PyObject* pyPointerToTimes;
  int numberLoadHistories;
  int totalNumberTimeSteps;
  int f77FileInput;
  int f77AsciiOutput;
  int asciiOutputInt;
  char* loadHistoryInputFile;
  char* asciiOutputFile;

  int ok = PyArg_ParseTuple(args, "OOiiiiiss:read_hist",
			    &pyPointerToHistry,
			    &pyPointerToTimes,
			    &numberLoadHistories,
			    &totalNumberTimeSteps,
			    &f77FileInput,
			    &f77AsciiOutput,
			    &asciiOutputInt,
			    &loadHistoryInputFile,
			    &asciiOutputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 1024;
  char errorstring[maxsize];
  double* pointerToHistry = (double*) PyCObject_AsVoidPtr(pyPointerToHistry);
  double* pointerToTimes = (double*) PyCObject_AsVoidPtr(pyPointerToTimes);

  read_hist_f(pointerToHistry,
	      pointerToTimes,
	      &numberLoadHistories,
	      &totalNumberTimeSteps,
	      &f77FileInput,
	      &f77AsciiOutput,
	      &asciiOutputInt,
	      loadHistoryInputFile,
	      asciiOutputFile,
	      &errorcode,
	      errorstring,
	      strlen(loadHistoryInputFile),
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


// Read material histories

char pylithomop3d_read_mathist__doc__[] = "";
char pylithomop3d_read_mathist__name__[] = "read_mathist";

PyObject * pylithomop3d_read_mathist(PyObject *, PyObject *args)
{
  PyObject* pyPointerToMhist;
  PyObject* pyPointerToMaterialInfo;
  int numberMaterials;
  int propertyListSize;
  int numberLoadHistories;
  int f77FileInput;
  int f77AsciiOutput;
  int f77PlotOutput;
  int asciiOutputInt;
  int plotOutputInt;
  char* materialHistoryInputFile;
  char* asciiOutputFile;
  char* plotOutputFile;

  int ok = PyArg_ParseTuple(args, "OOiiiiiiiisss:read_mathist",
			    &pyPointerToMhist,
  			    &pyPointerToMaterialInfo,
  			    &numberMaterials,
  			    &propertyListSize,
  			    &numberLoadHistories,
  			    &f77FileInput,
  			    &f77AsciiOutput,
  			    &f77PlotOutput,
  			    &asciiOutputInt,
  			    &plotOutputInt,
  			    &materialHistoryInputFile,
  			    &asciiOutputFile,
			    &plotOutputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 1024;
  char errorstring[maxsize];
  int* pointerToMhist = (int*) PyCObject_AsVoidPtr(pyPointerToMhist);
  int* pointerToMaterialInfo = (int*) PyCObject_AsVoidPtr(pyPointerToMaterialInfo);

  read_mathist_f(pointerToMhist,
		 pointerToMaterialInfo,
		 &numberMaterials,
		 &propertyListSize,
		 &numberLoadHistories,
		 &f77FileInput,
		 &f77AsciiOutput,
		 &f77PlotOutput,
		 &asciiOutputInt,
		 &plotOutputInt,
		 materialHistoryInputFile,
		 asciiOutputFile,
		 plotOutputFile,
		 &errorcode,
		 errorstring,
		 strlen(materialHistoryInputFile),
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

// Read element prestresses

// char pylithomop3d_read_prestr__doc__[] = "";
// char pylithomop3d_read_prestr__name__[] = "read_prestr";

// PyObject * pylithomop3d_read_prestr(PyObject *, PyObject *args)
// {
  // PyObject* pyPointerToStn;
  // PyObject* pyPointerToSt0;
  // PyObject* pyPointerToListArrayPrscal;
  // int numberStressComponents;
  // int numberGaussPoints;
  // int numberPrestressGaussPoints;
  // int numberElements;
  // int numberPrestressEntries;
  // int prestressAutoComputeInt;
  // int asciiOutputInt;
  // int f77FileInput;
  // int f77AsciiOutput;
  // char* prestressInputFile;
  // char* asciiOutputFile;

  // int ok = PyArg_ParseTuple(args, "OOOiiiiiiiiiss:read_prestr",
			    // &pyPointerToStn,
			    // &pyPointerToSt0,
			    // &pyPointerToListArrayPrscal,
			    // &numberStressComponents,
			    // &numberGaussPoints,
			    // &numberPrestressGaussPoints,
			    // &numberElements,
			    // &numberPrestressEntries,
			    // &prestressAutoComputeInt,
			    // &asciiOutputInt,
			    // &f77FileInput,
			    // &f77AsciiOutput,
			    // &prestressInputFile,
			    // &asciiOutputFile);

  // if (!ok) {
    // return 0;
  // }

  // int errorcode = 0;
  // const int maxsize = 1024;
  // char errorstring[maxsize];
  // double* pointerToStn = (double*) PyCObject_AsVoidPtr(pyPointerToStn);
  // double* pointerToSt0 = (double*) PyCObject_AsVoidPtr(pyPointerToSt0);
  // double* pointerToListArrayPrscal = (double*) PyCObject_AsVoidPtr(pyPointerToListArrayPrscal);

  // read_prestr_f(pointerToStn,
		// pointerToSt0,
		// pointerToListArrayPrscal,
		// &numberStressComponents,
		// &numberGaussPoints,
		// &numberPrestressGaussPoints,
		// &numberElements,
		// &numberPrestressEntries,
		// &prestressAutoComputeInt,
		// &asciiOutputInt,
		// &f77FileInput,
		// &f77AsciiOutput,
		// &errorcode,
		// prestressInputFile,
		// asciiOutputFile,
		// strlen(prestressInputFile),
		// strlen(asciiOutputFile));
    
// if(0 != exceptionhandler(errorcode, errorstring)) {
    // return 0;
  // }

  // journal::debug_t debug("lithomop3d");
  // debug
    // << journal::at(__HERE__)
    // << "numberPrestressEntries:" << numberPrestressEntries
    // << journal::endl;

  // return
  // Py_INCREF(Py_None);
  // return Py_None;
// }


// Read local coordinate rotations

char pylithomop3d_read_skew__doc__[] = "";
char pylithomop3d_read_skew__name__[] = "read_skew";

PyObject * pylithomop3d_read_skew(PyObject *, PyObject *args)
{
  PyObject* pyPointerToSkew;
  double rotationScaleFactor;
  int numberRotationEntries;
  int numberNodes;
  int autoRotateSlipperyNodesInt;
  int f77FileInput;
  int f77AsciiOutput;
  int asciiOutputInt;
  char* rotationInputFile;
  char* asciiOutputFile;

  int ok = PyArg_ParseTuple(args, "Odiiiiiiss:read_skew",
			    &pyPointerToSkew,
			    &rotationScaleFactor,
			    &numberRotationEntries,
			    &numberNodes,
			    &autoRotateSlipperyNodesInt,
			    &f77FileInput,
			    &f77AsciiOutput,
			    &asciiOutputInt,
			    &rotationInputFile,
			    &asciiOutputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 1024;
  char errorstring[maxsize];
  double* pointerToSkew = (double*) PyCObject_AsVoidPtr(pyPointerToSkew);

  read_skew_f(pointerToSkew,
	      &rotationScaleFactor,
	      &numberRotationEntries,
	      &numberNodes,
	      &autoRotateSlipperyNodesInt,
	      &f77FileInput,
	      &f77AsciiOutput,
	      &asciiOutputInt,
	      rotationInputFile,
	      asciiOutputFile,
	      &errorcode,
	      errorstring,
	      strlen(rotationInputFile),
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


// Read slippery node entries

char pylithomop3d_read_slip__doc__[] = "";
char pylithomop3d_read_slip__name__[] = "read_slip";

PyObject * pylithomop3d_read_slip(PyObject *, PyObject *args)
{
  PyObject* pyPointerToNslip;
  int numberSlipperyNodeEntries;
  int numberNodes;
  int autoRotateSlipperyNodesInt;
  int f77FileInput;
  int f77AsciiOutput;
  int f77PlotOutput;
  int asciiOutputInt;
  int plotOutputInt;
  char* slipperyNodeInputFile;
  char* asciiOutputFile;
  char* plotOutputFile;

  int ok = PyArg_ParseTuple(args, "Oiiiiiiiisss:read_slip",
			    &pyPointerToNslip,
			    &numberSlipperyNodeEntries,
			    &numberNodes,
			    &autoRotateSlipperyNodesInt,
			    &f77FileInput,
			    &f77AsciiOutput,
			    &f77PlotOutput,
			    &asciiOutputInt,
			    &plotOutputInt,
			    &slipperyNodeInputFile,
			    &asciiOutputFile,
			    &plotOutputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 1024;
  char errorstring[maxsize];
  int* pointerToNslip = (int*) PyCObject_AsVoidPtr(pyPointerToNslip);
  int totalNumberSlipperyNodes = 0;

  read_slip_f(pointerToNslip,
	      &numberSlipperyNodeEntries,
	      &totalNumberSlipperyNodes,
	      &numberNodes,
	      &autoRotateSlipperyNodesInt,
	      &f77FileInput,
	      &f77AsciiOutput,
	      &f77PlotOutput,
	      &asciiOutputInt,
	      &plotOutputInt,
	      slipperyNodeInputFile,
	      asciiOutputFile,
	      plotOutputFile,
	      &errorcode,
	      errorstring,
	      strlen(slipperyNodeInputFile),
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
  return Py_BuildValue("i", totalNumberSlipperyNodes);
}


// Read split node entries

char pylithomop3d_read_split__doc__[] = "";
char pylithomop3d_read_split__name__[] = "read_split";

PyObject * pylithomop3d_read_split(PyObject *, PyObject *args)
{
  PyObject* pyPointerToFault;
  PyObject* pyPointerToNfault;
  int numberSplitNodeEntries;
  int numberNodes;
  int f77FileInput;
  int f77AsciiOutput;
  int f77PlotOutput;
  int asciiOutputInt;
  int plotOutputInt;
  char* splitNodeInputFile;
  char* asciiOutputFile;
  char* plotOutputFile;

  int ok = PyArg_ParseTuple(args, "OOiiiiiiisss:read_split",
			    &pyPointerToFault,
			    &pyPointerToNfault,
			    &numberSplitNodeEntries,
			    &numberNodes,
			    &f77FileInput,
			    &f77AsciiOutput,
			    &f77PlotOutput,
			    &asciiOutputInt,
			    &plotOutputInt,
			    &splitNodeInputFile,
			    &asciiOutputFile,
			    &plotOutputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 1024;
  char errorstring[maxsize];
  double* pointerToFault = (double*) PyCObject_AsVoidPtr(pyPointerToFault);
  int* pointerToNfault = (int*) PyCObject_AsVoidPtr(pyPointerToNfault);
  int totalNumberSplitNodes = 0;

  read_split_f(pointerToFault,
	       pointerToNfault,
	       &numberSplitNodeEntries,
	       &totalNumberSplitNodes,
	       &numberNodes,
	       &f77FileInput,
	       &f77AsciiOutput,
	       &f77PlotOutput,
	       &asciiOutputInt,
	       &plotOutputInt,
	       splitNodeInputFile,
	       asciiOutputFile,
	       plotOutputFile,
	       &errorcode,
	       errorstring,
	       strlen(splitNodeInputFile),
	       strlen(asciiOutputFile),
	       strlen(plotOutputFile),
	       strlen(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "totalNumberSplitNodes:" << totalNumberSplitNodes
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue("i", totalNumberSplitNodes);
}


// Read state output information

char pylithomop3d_read_stateout__doc__[] = "";
char pylithomop3d_read_stateout__name__[] = "read_stateout";

PyObject * pylithomop3d_read_stateout(PyObject *, PyObject *args)
{
  PyObject* pyPointerToIstatout;
  int f77FileInput;
  int f77AsciiOutput;
  int f77PlotOutput;
  int asciiOutputInt;
  int plotOutputInt;
  char* stateVariableInputFile;
  char* asciiOutputFile;
  char* plotOutputFile;

  int ok = PyArg_ParseTuple(args, "Oiiiiisss:read_stateout",
			    &pyPointerToIstatout,
			    &f77FileInput,
			    &f77AsciiOutput,
			    &f77PlotOutput,
			    &asciiOutputInt,
			    &plotOutputInt,
			    &stateVariableInputFile,
			    &asciiOutputFile,
			    &plotOutputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 1024;
  char errorstring[maxsize];
  int* pointerToIstatout = (int*) PyCObject_AsVoidPtr(pyPointerToIstatout);

  read_stateout_f(pointerToIstatout,
	       &f77FileInput,
	       &f77AsciiOutput,
	       &f77PlotOutput,
	       &asciiOutputInt,
	       &plotOutputInt,
	       stateVariableInputFile,
	       asciiOutputFile,
	       plotOutputFile,
	       &errorcode,
	       errorstring,
	       strlen(stateVariableInputFile),
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


// Read time step data

char pylithomop3d_read_timdat__doc__[] = "";
char pylithomop3d_read_timdat__name__[] = "read_timdat";

PyObject * pylithomop3d_read_timdat(PyObject *, PyObject *args)
{
  PyObject* pyPointerToDelt;
  PyObject* pyPointerToAlfa;
  PyObject* pyPointerToUtol;
  PyObject* pyPointerToFtol;
  PyObject* pyPointerToEtol;
  PyObject* pyPointerToTimes;
  double timeScaleFactor;
  PyObject* pyPointerToMaxstp;
  PyObject* pyPointerToMaxit;
  PyObject* pyPointerToNtdinit;
  PyObject* pyPointerToLgdef;
  PyObject* pyPointerToItmax;
  int numberTimeStepGroups;
  int totalNumberTimeSteps;
  int f77FileInput;
  int f77AsciiOutput;
  int asciiOutputInt;
  char* timeStepInputFile;
  char* asciiOutputFile;

  int ok = PyArg_ParseTuple(args, "OOOOOOdOOOOOiiiiiss:read_timdat",
			    &pyPointerToDelt,
			    &pyPointerToAlfa,
			    &pyPointerToUtol,
			    &pyPointerToFtol,
			    &pyPointerToEtol,
			    &pyPointerToTimes,
			    &timeScaleFactor,
			    &pyPointerToMaxstp,
			    &pyPointerToMaxit,
			    &pyPointerToNtdinit,
			    &pyPointerToLgdef,
			    &pyPointerToItmax,
			    &numberTimeStepGroups,
			    &totalNumberTimeSteps,
			    &f77FileInput,
			    &f77AsciiOutput,
			    &asciiOutputInt,
			    &timeStepInputFile,
			    &asciiOutputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 1024;
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

  read_timdat_f(pointerToDelt,
		pointerToAlfa,
		pointerToUtol,
		pointerToFtol,
		pointerToEtol,
		pointerToTimes,
		&timeScaleFactor,
		pointerToMaxstp,
		pointerToMaxit,
		pointerToNtdinit,
		pointerToLgdef,
		pointerToItmax,
		&numberTimeStepGroups,
		&totalNumberTimeSteps,
		&f77FileInput,
		&f77AsciiOutput,
		&asciiOutputInt,
		timeStepInputFile,
		asciiOutputFile,
		&errorcode,
		errorstring,
		strlen(timeStepInputFile),
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


// Read traction BC

// char pylithomop3d_read_traction__doc__[] = "";
// char pylithomop3d_read_traction__name__[] = "read_traction";

// PyObject * pylithomop3d_read_traction(PyObject *, PyObject *args)
// {
  // PyObject* pyPointerToPres;
  // PyObject* pyPointerToPdir;
  // double tractionBcScaleFactor;
  // PyObject* pyPointerToIelno;
  // PyObject* pyPointerToIside;
  // PyObject* pyPointerToIhistry;
  // int numberTractionBc;
  // int numberElementNodes;
  // int numberTractionDirections;
  // int f77FileInput;
  // int f77AsciiOutput;
  // int asciiOutputInt;
  // char* tractionInputFile;
  // char* asciiOutputFile;

  // int ok = PyArg_ParseTuple(args, "OOdOOOiiiiiiss:read_traction",
			    // &pyPointerToPres,
			    // &pyPointerToPdir,
			    // &tractionBcScaleFactor,
			    // &pyPointerToIelno,
			    // &pyPointerToIside,
			    // &pyPointerToIhistry,
			    // &numberTractionBc,
			    // &numberElementNodes,
			    // &numberTractionDirections,
			    // &f77FileInput,
			    // &f77AsciiOutput,
			    // &asciiOutputInt,
			    // &tractionInputFile,
			    // &asciiOutputFile);

  // if (!ok) {
    // return 0;
  // }

  // int errorcode = 0;
  // const int maxsize = 1024;
  // char errorstring[maxsize];
  // double* pointerToPres = (double*) PyCObject_AsVoidPtr(pyPointerToPres);
  // double* pointerToPdir = (double*) PyCObject_AsVoidPtr(pyPointerToPdir);
  // int* pointerToIelno = (int*) PyCObject_AsVoidPtr(pyPointerToIelno);
  // int* pointerToIside = (int*) PyCObject_AsVoidPtr(pyPointerToIside);
  // int* pointerToIhistry = (int*) PyCObject_AsVoidPtr(pyPointerToIhistry);

  // read_traction_f(pointerToPres,
		  // pointerToPdir,
		  // &tractionBcScaleFactor,
		  // pointerToIelno,
		  // pointerToIside,
		  // pointerToIhistry,
		  // &numberTractionBc,
		  // &numberElementNodes,
		  // &numberTractionDirections,
		  // &f77FileInput,
		  // &f77AsciiOutput,
		  // &asciiOutputInt,
		  // &errorcode,
		  // tractionInputFile,
		  // asciiOutputFile,
		  // strlen(tractionInputFile),
		  // strlen(asciiOutputFile));
    
// if(0 != exceptionhandler(errorcode, errorstring)) {
    // return 0;
  // }

  // journal::debug_t debug("lithomop3d");
  // debug
    // << journal::at(__HERE__)
    // << "numberTractionBc:" << numberTractionBc
    // << journal::endl;

  // return
  // Py_INCREF(Py_None);
  // return Py_None;
// }


// Read winkler BC

char pylithomop3d_read_wink__doc__[] = "";
char pylithomop3d_read_wink__name__[] = "read_wink";

PyObject * pylithomop3d_read_wink(PyObject *, PyObject *args)
{
  PyObject* pyPointerToWink;
  PyObject* pyPointerToListArrayWscal;
  PyObject* pyPointerToIwink;
  PyObject* pyPointerToId;
  int numberNodes;
  int numberWinklerForces;
  int numberWinklerEntries;
  int f77FileInput;
  int f77AsciiOutput;
  int asciiOutputInt;
  char* winklerInputFile;
  char* asciiOutputFile;

  int ok = PyArg_ParseTuple(args, "OOOOiiiiiiss:read_wink",
			    &pyPointerToWink,
			    &pyPointerToListArrayWscal,
			    &pyPointerToIwink,
			    &pyPointerToId,
			    &numberNodes,
			    &numberWinklerForces,
			    &numberWinklerEntries,
			    &f77FileInput,
			    &f77AsciiOutput,
			    &asciiOutputInt,
			    &winklerInputFile,
			    &asciiOutputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 1024;
  char errorstring[maxsize];
  double* pointerToWink = (double*) PyCObject_AsVoidPtr(pyPointerToWink);
  double* pointerToListArrayWscal = (double*) PyCObject_AsVoidPtr(pyPointerToListArrayWscal);
  int* pointerToIwink = (int*) PyCObject_AsVoidPtr(pyPointerToIwink);
  int* pointerToId = (int*) PyCObject_AsVoidPtr(pyPointerToId);

  read_wink_f(pointerToWink,
	      pointerToListArrayWscal,
	      pointerToIwink,
	      pointerToId,
	      &numberNodes,
	      &numberWinklerForces,
	      &numberWinklerEntries,
	      &f77FileInput,
	      &f77AsciiOutput,
	      &asciiOutputInt,
	      winklerInputFile,
	      asciiOutputFile,
	      &errorcode,
	      errorstring,
	      strlen(winklerInputFile),
	      strlen(asciiOutputFile),
	      strlen(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberWinklerForces:" << numberWinklerForces
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


// Read winkler BC for slippery nodes

char pylithomop3d_read_winkx__doc__[] = "";
char pylithomop3d_read_winkx__name__[] = "read_winkx";

PyObject * pylithomop3d_read_winkx(PyObject *, PyObject *args)
{
  PyObject* pyPointerToWinkx;
  PyObject* pyPointerToListArrayWxscal;
  PyObject* pyPointerToIwinkx;
  PyObject* pyPointerToIdx;
  int numberNodes;
  int numberSlipperyWinklerForces;
  int numberSlipperyWinklerEntries;
  int f77FileInput;
  int f77AsciiOutput;
  int asciiOutputInt;
  char* slipperyWinklerInputFile;
  char* asciiOutputFile;

  int ok = PyArg_ParseTuple(args, "OOOOiiiiiiss:read_winkx",
			    &pyPointerToWinkx,
			    &pyPointerToListArrayWxscal,
			    &pyPointerToIwinkx,
			    &pyPointerToIdx,
			    &numberNodes,
			    &numberSlipperyWinklerForces,
			    &numberSlipperyWinklerEntries,
			    &f77FileInput,
			    &f77AsciiOutput,
			    &asciiOutputInt,
			    &slipperyWinklerInputFile,
			    &asciiOutputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 1024;
  char errorstring[maxsize];
  double* pointerToWinkx = (double*) PyCObject_AsVoidPtr(pyPointerToWinkx);
  double* pointerToListArrayWxscal = (double*) PyCObject_AsVoidPtr(pyPointerToListArrayWxscal);
  int* pointerToIwinkx = (int*) PyCObject_AsVoidPtr(pyPointerToIwinkx);
  int* pointerToIdx = (int*) PyCObject_AsVoidPtr(pyPointerToIdx);

  read_winkx_f(pointerToWinkx,
	       pointerToListArrayWxscal,
	       pointerToIwinkx,
	       pointerToIdx,
	       &numberNodes,
	       &numberSlipperyWinklerForces,
	       &numberSlipperyWinklerEntries,
	       &f77FileInput,
	       &f77AsciiOutput,
	       &asciiOutputInt,
	       slipperyWinklerInputFile,
	       asciiOutputFile,
	       &errorcode,
	       errorstring,
	       strlen(slipperyWinklerInputFile),
	       strlen(asciiOutputFile),
	       strlen(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberSlipperyWinklerForces:" << numberSlipperyWinklerForces
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}
    
// version
// $Id: parser.cc,v 1.6 2005/02/24 00:41:01 willic3 Exp $

// End of file
