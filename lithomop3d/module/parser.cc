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
#include "lithomop3d_externs.h"
// #include <stdio.h>


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
  int numberDegreesFreedom;
  int numberNodes;
  int numberBcEntries;
  int f77FileInput;
  int f77AsciiOutput;
  int asciiOutputInt;
  char* bcInputFile;
  char* asciiOutputFile;

  int ok = PyArg_ParseTuple(args, "OdddOOiiiiiiss:read_bc",
			    &pyPointerToBond,
			    &displacementScaleFactor,
			    &velocityScaleFactor,
			    &forceScaleFactor,
			    &pyPointerToIbond,
			    &pyPointerToId,
			    &numberDegreesFreedom,
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

  double* pointerToBond = (double*) PyCObject_AsVoidPtr(pyPointerToBond);
  int* pointerToIbond = (int*) PyCObject_AsVoidPtr(pyPointerToIbond);
  int* pointerToId = (int*) PyCObject_AsVoidPtr(pyPointerToId);
  int errorCode = 0;
  int numberGlobalEquations = 0;

  read_bc_f(pointerToBond,
	    &displacementScaleFactor,
	    &velocityScaleFactor,
	    &forceScaleFactor,
	    pointerToIbond,
	    pointerToId,
	    &numberDegreesFreedom,
	    &numberNodes,
	    &numberBcEntries,
	    &numberGlobalEquations,
	    &f77FileInput,
	    &f77AsciiOutput,
	    &asciiOutputInt,
	    &errorCode,
	    bcInputFile,
	    asciiOutputFile,
	    strlen(bcInputFile),
	    strlen(asciiOutputFile));
    
  switch(errorCode)
    {
    case 0:
      break;
    case 1:
      PyErr_SetString(PyExc_IOError,
		      "Parser: Error opening BC input file:");
      break;
    case 2:
      PyErr_SetString(PyExc_ValueError,
		      "Parser: All units not specified for BC:");
      break;
    case 3:
      PyErr_SetString(PyExc_IOError,
		      "Parser: Error reading BC input file:");
      break;
    case 4:
      PyErr_SetString(PyExc_ValueError,
		      "Parser: BC assigned for nonexistent node:");
      break;
    default:
      PyErr_SetString(PyExc_ValueError,
		      "Parser: Unknown error reading BC input file:");
    }

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberBcEntries:" << numberBcEntries
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue("i", numberGlobalEquations);
}


// Read connectivities

char pylithomop3d_read_connect__doc__[] = "";
char pylithomop3d_read_connect__name__[] = "read_connect";

PyObject * pylithomop3d_read_connect(PyObject *, PyObject *args)
{
  PyObject* pyPointerToIen;
  PyObject* pyPointerToMat;
  PyObject* pyPointerToInfin;
  int numberElementNodes;
  int numberElements;
  int numberNodes;
  int numberMaterialTypes;
  int numberGaussPoints;
  int f77FileInput;
  int f77AsciiOutput;
  int f77PlotOutput;
  int asciiOutputInt;
  int plotOutputInt;
  char* connectivityInputFile;
  char* asciiOutputFile;
  char* plotOutputFile;

  int ok = PyArg_ParseTuple(args, "OOOiiiiiiiiiisss:read_connect",
			    &pyPointerToIen,
			    &pyPointerToMat,
			    &pyPointerToInfin,
			    &numberElementNodes,
			    &numberElements,
			    &numberNodes,
			    &numberMaterialTypes,
			    &numberGaussPoints,
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

  int* pointerToIen = (int*) PyCObject_AsVoidPtr(pyPointerToIen);
  int* pointerToMat = (int*) PyCObject_AsVoidPtr(pyPointerToMat);
  int* pointerToInfin = (int*) PyCObject_AsVoidPtr(pyPointerToInfin);
  int errorCode = 0;

  read_connect_f(pointerToIen,
		 pointerToMat,
		 pointerToInfin,
		 &numberElementNodes,
		 &numberElements,
		 &numberNodes,
		 &numberMaterialTypes,
		 &numberGaussPoints,
		 &f77FileInput,
		 &f77AsciiOutput,
		 &f77PlotOutput,
		 &asciiOutputInt,
		 &plotOutputInt,
		 &errorCode,
		 connectivityInputFile,
		 asciiOutputFile,
		 plotOutputFile,
		 strlen(connectivityInputFile),
		 strlen(asciiOutputFile),
		 strlen(plotOutputFile));
    
  switch(errorCode)
    {
    case 0:
      break;
    case 1:
      PyErr_SetString(PyExc_IOError,
		      "Parser: Error opening connectivity input file:");
      break;
      //    case 2:
              //PyErr_SetString(PyExc_ValueError,
	                  //"Connectivity units not specified:");
              //break;
    case 3:
      PyErr_SetString(PyExc_IOError,
		      "Parser: Error reading connectivity input file:");
      break;
    case 4:
      PyErr_SetString(PyExc_ValueError,
		      "Parser: Bad material number:");
      break;
    case 5:
      PyErr_SetString(PyExc_ValueError,
		      "Parser: Bad connectivity:");
      break;
    default:
      PyErr_SetString(PyExc_ValueError,
		      "Parser: Unknown error reading connectivity input file:");
    }

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberMaterialTypes:" << numberMaterialTypes
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


// Read coordinates

char pylithomop3d_read_coords__doc__[] = "";
char pylithomop3d_read_coords__name__[] = "read_coords";

PyObject * pylithomop3d_read_coords(PyObject *, PyObject *args)
{
  PyObject* pyPointerToX;
  double coordinateScaleFactor;
  int numberSpaceDimensions;
  int numberNodes;
  int f77FileInput;
  int f77AsciiOutput;
  int f77PlotOutput;
  int asciiOutputInt;
  int plotOutputInt;
  char* coordinateInputFile;
  char* asciiOutputFile;
  char* plotOutputFile;

  int ok = PyArg_ParseTuple(args, "Odiiiiiiisss:read_coords",
			    &pyPointerToX,
			    &coordinateScaleFactor,
			    &numberSpaceDimensions,
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

  double* pointerToX = (double*) PyCObject_AsVoidPtr(pyPointerToX);
  int errorCode = 0;

  read_coords_f(pointerToX,
		&coordinateScaleFactor,
		&numberSpaceDimensions,
		&numberNodes,
		&f77FileInput,
		&f77AsciiOutput,
		&f77PlotOutput,
		&asciiOutputInt,
		&plotOutputInt,
		&errorCode,
		coordinateInputFile,
		asciiOutputFile,
		plotOutputFile,
		strlen(coordinateInputFile),
		strlen(asciiOutputFile),
		strlen(plotOutputFile));
    
  switch(errorCode)
    {
    case 0:
      break;
    case 1:
      PyErr_SetString(PyExc_IOError,
		      "Parser: Error opening coordinate input file:");
      break;
    case 2:
      PyErr_SetString(PyExc_ValueError,
		      "Parser: Coordinate units not specified:");
      break;
    case 3:
      PyErr_SetString(PyExc_IOError,
		      "Parser: Error reading coordinate input file:");
      break;
    case 4:
      PyErr_SetString(PyExc_ValueError,
		      "Parser: Wrong number of nodes read:");
      break;
    default:
      PyErr_SetString(PyExc_ValueError,
		      "Parser: Unknown error reading coordinate input file:");
    }

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberSpaceDimensions:" << numberSpaceDimensions
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
  int numberDegreesFreedom;
  int numberNodes;
  int numberSlipDimensions;
  int f77FileInput;
  int f77AsciiOutput;
  int asciiOutputInt;
  char* differentialForceInputFile;
  char* asciiOutputFile;

  int ok = PyArg_ParseTuple(args, "OOOiiiiiiiiss:read_diff",
			    &pyPointerToDiforc,
			    &pyPointerToNslip,
			    &pyPointerToIdhist,
			    &numberSlipperyNodeEntries,
			    &numberDifferentialForceEntries,
			    &numberDegreesFreedom,
			    &numberNodes,
			    &numberSlipDimensions,
			    &f77FileInput,
			    &f77AsciiOutput,
			    &asciiOutputInt,
			    &differentialForceInputFile,
			    &asciiOutputFile);

  if (!ok) {
    return 0;
  }

  double* pointerToDiforc = (double*) PyCObject_AsVoidPtr(pyPointerToDiforc);
  int* pointerToNslip = (int*) PyCObject_AsVoidPtr(pyPointerToNslip);
  int* pointerToIdhist = (int*) PyCObject_AsVoidPtr(pyPointerToIdhist);
  int errorCode = 0;

  read_diff_f(pointerToDiforc,
	      pointerToNslip,
	      pointerToIdhist,
	      &numberSlipperyNodeEntries,
	      &numberDifferentialForceEntries,
	      &numberDegreesFreedom,
	      &numberNodes,
	      &numberSlipDimensions,
	      &f77FileInput,
	      &f77AsciiOutput,
	      &asciiOutputInt,
	      &errorCode,
	      differentialForceInputFile,
	      asciiOutputFile,
	      strlen(differentialForceInputFile),
	      strlen(asciiOutputFile));
    
  switch(errorCode)
    {
    case 0:
      break;
      //    case 1:
            //PyErr_SetString(PyExc_IOError,
		            //"Error opening coordinate input file:");
            //break;
          //case 2:
            //PyErr_SetString(PyExc_ValueError,
		            //"Coordinate units not specified:");
            //break;
    case 3:
      PyErr_SetString(PyExc_IOError,
		      "Parser: Error reading differential force input file:");
      break;
    case 4:
      PyErr_SetString(PyExc_ValueError,
		      "Parser: Differential force applied to non-slippery node:");
      break;
    default:
      PyErr_SetString(PyExc_ValueError,
		      "Parser: Unknown error reading differential force input file:");
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

  int* pointerToIprint = (int*) PyCObject_AsVoidPtr(pyPointerToIprint);
  int errorCode = 0;

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
		&errorCode,
		fullOutputInputFile,
		asciiOutputFile,
		plotOutputFile,
		strlen(fullOutputInputFile),
		strlen(asciiOutputFile),
		strlen(plotOutputFile));
    
  switch(errorCode)
    {
    case 0:
      break;
    case 1:
      PyErr_SetString(PyExc_IOError,
		      "Parser: Error opening file specifying time steps for output:");
      break;
      //case 2:
            //PyErr_SetString(PyExc_ValueError,
		            //"Coordinate units not specified:");
            //break;
    case 3:
      PyErr_SetString(PyExc_IOError,
		      "Parser: Error reading file specifying time steps for output:");
      break;
      //case 4:
            //PyErr_SetString(PyExc_ValueError,
		            //"Wrong number of nodes read:");
            //break;
    default:
      PyErr_SetString(PyExc_ValueError,
		      "Parser: Unknown error reading time steps to output:");
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
  PyObject* pyPointerToHist;
  PyObject* pyPointerToTimes;
  int numberLoadHistories;
  int numberTimeSteps;
  int f77FileInput;
  int f77AsciiOutput;
  int asciiOutputInt;
  char* loadHistoryInputFile;
  char* asciiOutputFile;

  int ok = PyArg_ParseTuple(args, "OOiiiiiss:read_hist",
			    &pyPointerToHist,
			    &pyPointerToTimes,
			    &numberLoadHistories,
			    &numberTimeSteps,
			    &f77FileInput,
			    &f77AsciiOutput,
			    &asciiOutputInt,
			    &loadHistoryInputFile,
			    &asciiOutputFile);

  if (!ok) {
    return 0;
  }

  double* pointerToHist = (double*) PyCObject_AsVoidPtr(pyPointerToHist);
  double* pointerToTimes = (double*) PyCObject_AsVoidPtr(pyPointerToTimes);
  int errorCode = 0;

  read_hist_f(pointerToHist,
	      pointerToTimes,
	      &numberLoadHistories,
	      &numberTimeSteps,
	      &f77FileInput,
	      &f77AsciiOutput,
	      &asciiOutputInt,
	      &errorCode,
	      loadHistoryInputFile,
	      asciiOutputFile,
	      strlen(loadHistoryInputFile),
	      strlen(asciiOutputFile));
    
  switch(errorCode)
    {
    case 0:
      break;
      //case 1:
            //PyErr_SetString(PyExc_IOError,
		            //"Error opening load history input file:");
            //break;
      //case 2:
            //PyErr_SetString(PyExc_ValueError,
		            //"Load history units not specified:");
      //break;
    case 3:
      PyErr_SetString(PyExc_IOError,
		      "Parser: Error reading load history input file:");
      break;
    case 4:
      PyErr_SetString(PyExc_ValueError,
		      "Parser: Load history times are out of order:");
      break;
    default:
      PyErr_SetString(PyExc_ValueError,
		      "Parser: Unknown error reading load history input file:");
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


// Read element prestresses

char pylithomop3d_read_prestr__doc__[] = "";
char pylithomop3d_read_prestr__name__[] = "read_prestr";

PyObject * pylithomop3d_read_prestr(PyObject *, PyObject *args)
{
  PyObject* pyPointerToStn;
  PyObject* pyPointerToSt0;
  PyObject* pyPointerToListArrayPrscal;
  int numberStressComponents;
  int numberGaussPoints;
  int numberPrestressGaussPoints;
  int numberElements;
  int numberPrestressEntries;
  int prestressAutoComputeInt;
  int asciiOutputInt;
  int f77FileInput;
  int f77AsciiOutput;
  char* prestressInputFile;
  char* asciiOutputFile;

  int ok = PyArg_ParseTuple(args, "OOOiiiiiiiiiss:read_prestr",
			    &pyPointerToStn,
			    &pyPointerToSt0,
			    &pyPointerToListArrayPrscal,
			    &numberStressComponents,
			    &numberGaussPoints,
			    &numberPrestressGaussPoints,
			    &numberElements,
			    &numberPrestressEntries,
			    &prestressAutoComputeInt,
			    &asciiOutputInt,
			    &f77FileInput,
			    &f77AsciiOutput,
			    &prestressInputFile,
			    &asciiOutputFile);

  if (!ok) {
    return 0;
  }

  double* pointerToStn = (double*) PyCObject_AsVoidPtr(pyPointerToStn);
  double* pointerToSt0 = (double*) PyCObject_AsVoidPtr(pyPointerToSt0);
  double* pointerToListArrayPrscal = (double*) PyCObject_AsVoidPtr(pyPointerToListArrayPrscal);
  int errorCode = 0;

  read_prestr_f(pointerToStn,
		pointerToSt0,
		pointerToListArrayPrscal,
		&numberStressComponents,
		&numberGaussPoints,
		&numberPrestressGaussPoints,
		&numberElements,
		&numberPrestressEntries,
		&prestressAutoComputeInt,
		&asciiOutputInt,
		&f77FileInput,
		&f77AsciiOutput,
		&errorCode,
		prestressInputFile,
		asciiOutputFile,
		strlen(prestressInputFile),
		strlen(asciiOutputFile));
    
  switch(errorCode)
    {
    case 0:
      break;
      //case 1:
            //PyErr_SetString(PyExc_IOError,
		            //"Error opening coordinate input file:");
            //break;
          //case 2:
            //PyErr_SetString(PyExc_ValueError,
		            //"Coordinate units not specified:");
            //break;
    case 3:
      PyErr_SetString(PyExc_IOError,
		      "Parser: Error reading prestress input file:");
      break;
      //case 4:
            //PyErr_SetString(PyExc_ValueError,
		            //"Wrong number of nodes read:");
            //break;
    default:
      PyErr_SetString(PyExc_ValueError,
		      "Parser: Unknown error reading prestress input file:");
    }

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberPrestressEntries:" << numberPrestressEntries
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


// Read material properties

char pylithomop3d_read_prop__doc__[] = "";
char pylithomop3d_read_prop__name__[] = "read_prop";

PyObject * pylithomop3d_read_prop(PyObject *, PyObject *args)
{
  PyObject* pyPointerToProp;
  PyObject* pyPointerToListArrayGrav;
  double densityScaleFactor;
  double youngScaleFactor;
  double viscosityCoefficientScaleFactor;
  double cohesionScaleFactor;
  int numberDegreesFreedom;
  int numberMaterialProperties;
  int numberMaterialTypes;
  int asciiOutputInt;
  int plotOutputInt;
  int f77FileInput;
  int f77AsciiOutput;
  int f77PlotOutput;
  int viscousFlagInt;
  int plasticFlagInt;
  char* materialPropertiesInputFile;
  char* asciiOutputFile;
  char* plotOutputFile;

  int ok = PyArg_ParseTuple(args, "OOddddiiiiiiiiiisss:read_prop",
			    &pyPointerToProp,
			    &pyPointerToListArrayGrav,
			    &densityScaleFactor,
			    &youngScaleFactor,
			    &viscosityCoefficientScaleFactor,
			    &cohesionScaleFactor,
			    &numberDegreesFreedom,
			    &numberMaterialProperties,
			    &numberMaterialTypes,
			    &asciiOutputInt,
			    &plotOutputInt,
			    &f77FileInput,
			    &f77AsciiOutput,
			    &f77PlotOutput,
			    &viscousFlagInt,
			    &plasticFlagInt,
			    &materialPropertiesInputFile,
			    &asciiOutputFile,
			    &plotOutputFile);

  if (!ok) {
    return 0;
  }

  double* pointerToProp = (double*) PyCObject_AsVoidPtr(pyPointerToProp);
  double* pointerToListArrayGrav = (double*) PyCObject_AsVoidPtr(pyPointerToListArrayGrav);
  int errorCode = 0;

  read_prop_f(pointerToProp,
	      pointerToListArrayGrav,
	      &densityScaleFactor,
	      &youngScaleFactor,
	      &viscosityCoefficientScaleFactor,
	      &cohesionScaleFactor,
	      &numberDegreesFreedom,
	      &numberMaterialProperties,
	      &numberMaterialTypes,
	      &asciiOutputInt,
	      &plotOutputInt,
	      &f77FileInput,
	      &f77AsciiOutput,
	      &f77PlotOutput,
	      &errorCode,
	      &viscousFlagInt,
	      &plasticFlagInt,
	      materialPropertiesInputFile,
	      asciiOutputFile,
	      plotOutputFile,
	      strlen(materialPropertiesInputFile),
	      strlen(asciiOutputFile),
	      strlen(plotOutputFile));
    
  switch(errorCode)
    {
    case 0:
      break;
    case 1:
      PyErr_SetString(PyExc_IOError,
		      "Parser: Error opening material properties input file:");
      break;
    case 2:
      PyErr_SetString(PyExc_ValueError,
		      "Parser: All material property units not specified:");
      break;
    case 3:
      PyErr_SetString(PyExc_IOError,
		      "Parser: Error reading material properties input file:");
      break;
      //case 4:
            //PyErr_SetString(PyExc_ValueError,
		            //"Wrong number of nodes read:");
            //break;
    default:
      PyErr_SetString(PyExc_ValueError,
		      "Parser: Unknown error reading material properties input file:");
    }

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberMaterialTypes:" << numberMaterialTypes
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


// Read local coordinate rotations

char pylithomop3d_read_skew__doc__[] = "";
char pylithomop3d_read_skew__name__[] = "read_skew";

PyObject * pylithomop3d_read_skew(PyObject *, PyObject *args)
{
  PyObject* pyPointerToSkew;
  double rotationScaleFactor;
  int numberRotationEntries;
  int numberSkewDimensions;
  int numberNodes;
  int autoRotateSlipperyNodesInt;
  int f77FileInput;
  int f77AsciiOutput;
  int asciiOutputInt;
  char* rotationInputFile;
  char* asciiOutputFile;

  int ok = PyArg_ParseTuple(args, "Odiiiiiiiss:read_skew",
			    &pyPointerToSkew,
			    &rotationScaleFactor,
			    &numberRotationEntries,
			    &numberSkewDimensions,
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

  double* pointerToSkew = (double*) PyCObject_AsVoidPtr(pyPointerToSkew);
  int errorCode = 0;

  read_skew_f(pointerToSkew,
	      &rotationScaleFactor,
	      &numberRotationEntries,
	      &numberSkewDimensions,
	      &numberNodes,
	      &autoRotateSlipperyNodesInt,
	      &f77FileInput,
	      &f77AsciiOutput,
	      &asciiOutputInt,
	      &errorCode,
	      rotationInputFile,
	      asciiOutputFile,
	      strlen(rotationInputFile),
	      strlen(asciiOutputFile));
    
  switch(errorCode)
    {
    case 0:
      break;
      //case 1:
            //PyErr_SetString(PyExc_IOError,
		            //"Error opening coordinate input file:");
            //break;
    case 2:
      PyErr_SetString(PyExc_ValueError,
		      "Parser: Rotation units not specified:");
      break;
    case 3:
      PyErr_SetString(PyExc_IOError,
		      "Parser: Error reading rotation input file:");
      break;
      //case 4:
            //PyErr_SetString(PyExc_ValueError,
		            //"Wrong number of nodes read:");
            //break;
    default:
      PyErr_SetString(PyExc_ValueError,
		      "Parser: Unknown error reading rotation input file:");
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
  int numberDegreesFreedom;
  int numberNodes;
  int autoRotateSlipperyNodesInt;
  int numberSlipDimensions;
  int f77FileInput;
  int f77AsciiOutput;
  int f77PlotOutput;
  int asciiOutputInt;
  int plotOutputInt;
  char* slipperyNodeInputFile;
  char* asciiOutputFile;
  char* plotOutputFile;

  int ok = PyArg_ParseTuple(args, "Oiiiiiiiiiisss:read_slip",
			    &pyPointerToNslip,
			    &numberSlipperyNodeEntries,
			    &numberDegreesFreedom,
			    &numberNodes,
			    &autoRotateSlipperyNodesInt,
			    &numberSlipDimensions,
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

  int* pointerToNslip = (int*) PyCObject_AsVoidPtr(pyPointerToNslip);
  int totalNumberSlipperyNodes = 0;
  int errorCode = 0;

  read_slip_f(pointerToNslip,
	      &numberSlipperyNodeEntries,
	      &totalNumberSlipperyNodes,
	      &numberDegreesFreedom,
	      &numberNodes,
	      &autoRotateSlipperyNodesInt,
	      &numberSlipDimensions,
	      &f77FileInput,
	      &f77AsciiOutput,
	      &f77PlotOutput,
	      &asciiOutputInt,
	      &plotOutputInt,
	      &errorCode,
	      slipperyNodeInputFile,
	      asciiOutputFile,
	      plotOutputFile,
	      strlen(slipperyNodeInputFile),
	      strlen(asciiOutputFile),
	      strlen(plotOutputFile));
    
  switch(errorCode)
    {
    case 0:
      break;
      //case 1:
            //PyErr_SetString(PyExc_IOError,
		            //"Error opening coordinate input file:");
            //break;
          //case 2:
            //PyErr_SetString(PyExc_ValueError,
		            //"Coordinate units not specified:");
            //break;
    case 3:
      PyErr_SetString(PyExc_IOError,
		      "Parser: Error reading slippery nodes input file:");
      break;
      //case 4:
            //PyErr_SetString(PyExc_ValueError,
		            //"Wrong number of nodes read:");
            //break;
    default:
      PyErr_SetString(PyExc_ValueError,
		      "Parser: Unknown error reading slippery nodes input file:");
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
  int numberDegreesFreedom;
  int numberNodes;
  int f77FileInput;
  int f77AsciiOutput;
  int f77PlotOutput;
  int asciiOutputInt;
  int plotOutputInt;
  char* splitNodeInputFile;
  char* asciiOutputFile;
  char* plotOutputFile;

  int ok = PyArg_ParseTuple(args, "OOiiiiiiiisss:read_split",
			    &pyPointerToFault,
			    &pyPointerToNfault,
			    &numberSplitNodeEntries,
			    &numberDegreesFreedom,
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

  double* pointerToFault = (double*) PyCObject_AsVoidPtr(pyPointerToFault);
  int* pointerToNfault = (int*) PyCObject_AsVoidPtr(pyPointerToNfault);
  int totalNumberSplitNodes = 0;
  int errorCode = 0;

  read_split_f(pointerToFault,
	       pointerToNfault,
	       &numberSplitNodeEntries,
	       &totalNumberSplitNodes,
	       &numberDegreesFreedom,
	       &numberNodes,
	       &f77FileInput,
	       &f77AsciiOutput,
	       &f77PlotOutput,
	       &asciiOutputInt,
	       &plotOutputInt,
	       &errorCode,
	       splitNodeInputFile,
	       asciiOutputFile,
	       plotOutputFile,
	       strlen(splitNodeInputFile),
	       strlen(asciiOutputFile),
	       strlen(plotOutputFile));
    
  switch(errorCode)
    {
    case 0:
      break;
      //case 1:
            //PyErr_SetString(PyExc_IOError,
		            //"Error opening coordinate input file:");
            //break;
          //case 2:
            //PyErr_SetString(PyExc_ValueError,
		            //"Coordinate units not specified:");
            //break;
    case 3:
      PyErr_SetString(PyExc_IOError,
		      "Parser: Error reading split nodes input file:");
      break;
      //case 4:
            //PyErr_SetString(PyExc_ValueError,
		            //"Wrong number of nodes read:");
            //break;
    default:
      PyErr_SetString(PyExc_ValueError,
		      "Parser: Unknown error reading split nodes input file:");
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
  PyObject* pyPointerToMaxitc;
  PyObject* pyPointerToLgdef;
  PyObject* pyPointerToIbbar;
  PyObject* pyPointerToItmax;
  int numberTimeStepGroups;
  int totalNumberTimeSteps;
  int f77FileInput;
  int f77AsciiOutput;
  int asciiOutputInt;
  char* timeStepInputFile;
  char* asciiOutputFile;

  int ok = PyArg_ParseTuple(args, "OOOOOOdOOOOOOiiiiiss:read_timdat",
			    &pyPointerToDelt,
			    &pyPointerToAlfa,
			    &pyPointerToUtol,
			    &pyPointerToFtol,
			    &pyPointerToEtol,
			    &pyPointerToTimes,
			    &timeScaleFactor,
			    &pyPointerToMaxstp,
			    &pyPointerToMaxit,
			    &pyPointerToMaxitc,
			    &pyPointerToLgdef,
			    &pyPointerToIbbar,
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

  double* pointerToDelt = (double*) PyCObject_AsVoidPtr(pyPointerToDelt);
  double* pointerToAlfa = (double*) PyCObject_AsVoidPtr(pyPointerToAlfa);
  double* pointerToUtol = (double*) PyCObject_AsVoidPtr(pyPointerToUtol);
  double* pointerToFtol = (double*) PyCObject_AsVoidPtr(pyPointerToFtol);
  double* pointerToEtol = (double*) PyCObject_AsVoidPtr(pyPointerToEtol);
  double* pointerToTimes = (double*) PyCObject_AsVoidPtr(pyPointerToTimes);
  int* pointerToMaxstp = (int*) PyCObject_AsVoidPtr(pyPointerToMaxstp);
  int* pointerToMaxit = (int*) PyCObject_AsVoidPtr(pyPointerToMaxit);
  int* pointerToMaxitc = (int*) PyCObject_AsVoidPtr(pyPointerToMaxitc);
  int* pointerToLgdef = (int*) PyCObject_AsVoidPtr(pyPointerToLgdef);
  int* pointerToIbbar = (int*) PyCObject_AsVoidPtr(pyPointerToIbbar);
  int* pointerToItmax = (int*) PyCObject_AsVoidPtr(pyPointerToItmax);
  int errorCode = 0;

  read_timdat_f(pointerToDelt,
		pointerToAlfa,
		pointerToUtol,
		pointerToFtol,
		pointerToEtol,
		pointerToTimes,
		&timeScaleFactor,
		pointerToMaxstp,
		pointerToMaxit,
		pointerToMaxitc,
		pointerToLgdef,
		pointerToIbbar,
		pointerToItmax,
		&numberTimeStepGroups,
		&totalNumberTimeSteps,
		&f77FileInput,
		&f77AsciiOutput,
		&asciiOutputInt,
		&errorCode,
		timeStepInputFile,
		asciiOutputFile,
		strlen(timeStepInputFile),
		strlen(asciiOutputFile));
    
  switch(errorCode)
    {
    case 0:
      break;
    case 1:
      PyErr_SetString(PyExc_IOError,
		      "Parser: Error opening time step input file:");
      break;
    case 2:
      PyErr_SetString(PyExc_ValueError,
		      "Parser: Time units not specified:");
      break;
    case 3:
      PyErr_SetString(PyExc_IOError,
		      "Parser: Error reading time step input file:");
      break;
      //case 4:
            //PyErr_SetString(PyExc_ValueError,
		            //"Wrong number of nodes read:");
            //break;
    default:
      PyErr_SetString(PyExc_ValueError,
		      "Parser: Unknown error reading time step input file:");
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

char pylithomop3d_read_traction__doc__[] = "";
char pylithomop3d_read_traction__name__[] = "read_traction";

PyObject * pylithomop3d_read_traction(PyObject *, PyObject *args)
{
  PyObject* pyPointerToPres;
  PyObject* pyPointerToPdir;
  double tractionBcScaleFactor;
  PyObject* pyPointerToIelno;
  PyObject* pyPointerToIside;
  PyObject* pyPointerToIhistry;
  int numberTractionBc;
  int numberElementNodes;
  int numberTractionDirections;
  int f77FileInput;
  int f77AsciiOutput;
  int asciiOutputInt;
  char* tractionInputFile;
  char* asciiOutputFile;

  int ok = PyArg_ParseTuple(args, "OOdOOOiiiiiiss:read_traction",
			    &pyPointerToPres,
			    &pyPointerToPdir,
			    &tractionBcScaleFactor,
			    &pyPointerToIelno,
			    &pyPointerToIside,
			    &pyPointerToIhistry,
			    &numberTractionBc,
			    &numberElementNodes,
			    &numberTractionDirections,
			    &f77FileInput,
			    &f77AsciiOutput,
			    &asciiOutputInt,
			    &tractionInputFile,
			    &asciiOutputFile);

  if (!ok) {
    return 0;
  }

  double* pointerToPres = (double*) PyCObject_AsVoidPtr(pyPointerToPres);
  double* pointerToPdir = (double*) PyCObject_AsVoidPtr(pyPointerToPdir);
  int* pointerToIelno = (int*) PyCObject_AsVoidPtr(pyPointerToIelno);
  int* pointerToIside = (int*) PyCObject_AsVoidPtr(pyPointerToIside);
  int* pointerToIhistry = (int*) PyCObject_AsVoidPtr(pyPointerToIhistry);
  int errorCode = 0;

  read_traction_f(pointerToPres,
		  pointerToPdir,
		  &tractionBcScaleFactor,
		  pointerToIelno,
		  pointerToIside,
		  pointerToIhistry,
		  &numberTractionBc,
		  &numberElementNodes,
		  &numberTractionDirections,
		  &f77FileInput,
		  &f77AsciiOutput,
		  &asciiOutputInt,
		  &errorCode,
		  tractionInputFile,
		  asciiOutputFile,
		  strlen(tractionInputFile),
		  strlen(asciiOutputFile));
    
  switch(errorCode)
    {
    case 0:
      break;
    case 1:
      PyErr_SetString(PyExc_IOError,
		      "Parser: Error opening traction input file:");
      break;
    case 2:
      PyErr_SetString(PyExc_ValueError,
		      "Parser: Traction units not specified:");
      break;
    case 3:
      PyErr_SetString(PyExc_IOError,
		      "Parser: Error reading traction input file:");
      break;
      //case 4:
            //PyErr_SetString(PyExc_ValueError,
		            //"Wrong number of nodes read:");
            //break;
    default:
      PyErr_SetString(PyExc_ValueError,
		      "Parser: Unknown error reading traction input file:");
    }

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberTractionBc:" << numberTractionBc
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


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
  int numberDegreesFreedom;
  int numberWinklerForces;
  int numberWinklerEntries;
  int f77FileInput;
  int f77AsciiOutput;
  int asciiOutputInt;
  char* winklerInputFile;
  char* asciiOutputFile;

  int ok = PyArg_ParseTuple(args, "OOOOiiiiiiiss:read_wink",
			    &pyPointerToWink,
			    &pyPointerToListArrayWscal,
			    &pyPointerToIwink,
			    &pyPointerToId,
			    &numberNodes,
			    &numberDegreesFreedom,
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

  double* pointerToWink = (double*) PyCObject_AsVoidPtr(pyPointerToWink);
  double* pointerToListArrayWscal = (double*) PyCObject_AsVoidPtr(pyPointerToListArrayWscal);
  int* pointerToIwink = (int*) PyCObject_AsVoidPtr(pyPointerToIwink);
  int* pointerToId = (int*) PyCObject_AsVoidPtr(pyPointerToId);
  int errorCode = 0;

  read_wink_f(pointerToWink,
	      pointerToListArrayWscal,
	      pointerToIwink,
	      pointerToId,
	      &numberNodes,
	      &numberDegreesFreedom,
	      &numberWinklerForces,
	      &numberWinklerEntries,
	      &f77FileInput,
	      &f77AsciiOutput,
	      &asciiOutputInt,
	      &errorCode,
	      winklerInputFile,
	      asciiOutputFile,
	      strlen(winklerInputFile),
	      strlen(asciiOutputFile));
    
  switch(errorCode)
    {
    case 0:
      break;
    case 1:
      PyErr_SetString(PyExc_IOError,
		      "Parser: Error opening winkler input file:");
      break;
      //case 2:
            //PyErr_SetString(PyExc_ValueError,
		            //"Coordinate units not specified:");
            //break;
    case 3:
      PyErr_SetString(PyExc_IOError,
		      "Parser: Error reading winkler input file:");
      break;
      //case 4:
            //PyErr_SetString(PyExc_ValueError,
		            //"Wrong number of nodes read:");
            //break;
    default:
      PyErr_SetString(PyExc_ValueError,
		      "Parser: Unknown error reading winkler input file:");
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
  int numberDegreesFreedom;
  int numberSlipperyWinklerForces;
  int numberSlipperyWinklerEntries;
  int f77FileInput;
  int f77AsciiOutput;
  int asciiOutputInt;
  char* slipperyWinklerInputFile;
  char* asciiOutputFile;

  int ok = PyArg_ParseTuple(args, "OOOOiiiiiiiss:read_winkx",
			    &pyPointerToWinkx,
			    &pyPointerToListArrayWxscal,
			    &pyPointerToIwinkx,
			    &pyPointerToIdx,
			    &numberNodes,
			    &numberDegreesFreedom,
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

  double* pointerToWinkx = (double*) PyCObject_AsVoidPtr(pyPointerToWinkx);
  double* pointerToListArrayWxscal = (double*) PyCObject_AsVoidPtr(pyPointerToListArrayWxscal);
  int* pointerToIwinkx = (int*) PyCObject_AsVoidPtr(pyPointerToIwinkx);
  int* pointerToIdx = (int*) PyCObject_AsVoidPtr(pyPointerToIdx);
  int errorCode = 0;

  read_winkx_f(pointerToWinkx,
	       pointerToListArrayWxscal,
	       pointerToIwinkx,
	       pointerToIdx,
	       &numberNodes,
	       &numberDegreesFreedom,
	       &numberSlipperyWinklerForces,
	       &numberSlipperyWinklerEntries,
	       &f77FileInput,
	       &f77AsciiOutput,
	       &asciiOutputInt,
	       &errorCode,
	       slipperyWinklerInputFile,
	       asciiOutputFile,
	       strlen(slipperyWinklerInputFile),
	       strlen(asciiOutputFile));
    
  switch(errorCode)
    {
    case 0:
      break;
    case 1:
      PyErr_SetString(PyExc_IOError,
		      "Parser: Error opening winkler input file for slippery nodes:");
      break;
      //case 2:
            //PyErr_SetString(PyExc_ValueError,
		            //"Coordinate units not specified:");
            //break;
    case 3:
      PyErr_SetString(PyExc_IOError,
		      "Parser: Error reading winkler input file for slippery nodes:");
      break;
      //case 4:
            //PyErr_SetString(PyExc_ValueError,
		            //"Wrong number of nodes read:");
            //break;
    default:
      PyErr_SetString(PyExc_ValueError,
		      "Parser: Unknown error reading winkler input file for slippery nodes:");
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
// $Id: parser.cc,v 1.1 2004/04/14 21:24:47 willic3 Exp $

// End of file
