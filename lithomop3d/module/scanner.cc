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

#include "scanner.h"
#include "lithomop3d_externs.h"
// #include <stdio.h>


// Scan boundary conditions

char pylithomop3d_scan_bc__doc__[] = "";
char pylithomop3d_scan_bc__name__[] = "scan_bc";

PyObject * pylithomop3d_scan_bc(PyObject *, PyObject *args)
{
  int numberDegreesFreedom;
  int f77FileInput;
  char* displacementUnits;
  char* velocityUnits;
  char* forceUnits;
  char* bcInputFile;

  int ok = PyArg_ParseTuple(args, "iissss:scan_bc",
			    &numberDegreesFreedom,
			    &f77FileInput,
			    &displacementUnits,
			    &velocityUnits,
			    &forceUnits,
			    &bcInputFile);

  if (!ok) {
    return 0;
  }

  int errorCode = 0;
  int numberBcEntries = 0;

  scan_bc_f(&numberDegreesFreedom,
	    &numberBcEntries,
	    &f77FileInput,
	    &errorCode,
	    displacementUnits,
	    velocityUnits,
	    forceUnits,
	    bcInputFile,
	    strlen(displacementUnits),
	    strlen(velocityUnits),
	    strlen(forceUnits),
	    strlen(bcInputFile));
    
  switch(errorCode)
    {
    case 0:
      break;
    case 1:
      PyErr_SetString(PyExc_IOError,
		      "Scanner: Error opening BC input file:");
      break;
    case 2:
      PyErr_SetString(PyExc_ValueError,
		      "Scanner: All units not specified for BC:");
      break;
    case 3:
      PyErr_SetString(PyExc_IOError,
		      "Scanner: Error reading BC input file:");
      break;
      //case 4:
            //PyErr_SetString(PyExc_ValueError,
		            //"BC assigned for nonexistent node:");
            //break;
    default:
      PyErr_SetString(PyExc_ValueError,
		      "Scanner: Unknown error reading BC input file:");
    }

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberBcEntries:" << numberBcEntries
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue("i", numberBcEntries);
}


// Scan connectivities

char pylithomop3d_scan_connect__doc__[] = "";
char pylithomop3d_scan_connect__name__[] = "scan_connect";

PyObject * pylithomop3d_scan_connect(PyObject *, PyObject *args)
{
  int numberElementNodes;
  int f77FileInput;
  char* connectivityInputFile;

  int ok = PyArg_ParseTuple(args, "iis:scan_connect",
			    &numberElementNodes,
			    &f77FileInput,
			    &connectivityInputFile);

  if (!ok) {
    return 0;
  }

  int errorCode = 0;
  int numberElements = 0;

  scan_connect_f(&numberElementNodes,
		 &numberElements,
		 &f77FileInput,
		 &errorCode,
		 connectivityInputFile,strlen(connectivityInputFile));
    
  switch(errorCode)
    {
    case 0:
      break;
    case 1:
      PyErr_SetString(PyExc_IOError,
		      "Scanner: Error opening connectivity input file:");
      break;
      //    case 2:
              //PyErr_SetString(PyExc_ValueError,
	                  //"Connectivity units not specified:");
              //break;
    case 3:
      PyErr_SetString(PyExc_IOError,
		      "Scanner: Error reading connectivity input file:");
      break;
      //case 4:
            //PyErr_SetString(PyExc_ValueError,
		            //"Bad material number:");
            //break;
          //case 5:
            //PyErr_SetString(PyExc_ValueError,
		            //"Bad connectivity:");
            //break;
    default:
      PyErr_SetString(PyExc_ValueError,
		      "Scanner: Unknown error reading connectivity input file:");
    }

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberElements:" << numberElements
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue("i", numberElements);
}


// Read coordinates

char pylithomop3d_scan_coords__doc__[] = "";
char pylithomop3d_scan_coords__name__[] = "scan_coords";

PyObject * pylithomop3d_scan_coords(PyObject *, PyObject *args)
{
  int numberSpaceDimensions;
  int f77FileInput;
  char *coordinateUnits;
  char *coordinateInputFile;

  int ok = PyArg_ParseTuple(args, "iiss:scan_coords",
			    &numberSpaceDimensions,
			    &f77FileInput,
			    &coordinateUnits,
			    &coordinateInputFile);

  if (!ok) {
    return 0;
  }

  int errorCode = 0;
  int numberNodes = 0;

  scan_coords_f(&numberSpaceDimensions,
		&numberNodes,
		&f77FileInput,
		&errorCode,
		coordinateUnits,
		coordinateInputFile,
		strlen(coordinateUnits),
		strlen(coordinateInputFile));
    
  switch(errorCode)
    {
    case 0:
      break;
    case 1:
      PyErr_SetString(PyExc_IOError,
		      "Scanner: Error opening coordinate input file:");
      break;
    case 2:
      PyErr_SetString(PyExc_ValueError,
		      "Scanner: Coordinate units not specified:");
      break;
    case 3:
      PyErr_SetString(PyExc_IOError,
		      "Scanner: Error reading coordinate input file:");
      break;
      //case 4:
            //PyErr_SetString(PyExc_ValueError,
		            //"Wrong number of nodes read:");
            //break;
    default:
      PyErr_SetString(PyExc_ValueError,
		      "Scanner: Unknown error reading coordinate input file:");
    }

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberSpaceDimensions:" << numberSpaceDimensions
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue("i", numberNodes);
}


// Read differential forces

char pylithomop3d_scan_diff__doc__[] = "";
char pylithomop3d_scan_diff__name__[] = "scan_diff";

PyObject * pylithomop3d_scan_diff(PyObject *, PyObject *args)
{
  int numberDegreesFreedom;
  int numberSlipperyNodeEntries;
  int f77FileInput;
  char* differentialForceInputFile;

  int ok = PyArg_ParseTuple(args, "iiis:scan_diff",
			    &numberDegreesFreedom,
			    &numberSlipperyNodeEntries,
			    &f77FileInput,
			    &differentialForceInputFile);

  if (!ok) {
    return 0;
  }

  int errorCode = 0;
  int numberDifferentialForceEntries = 0;

  scan_diff_f(&numberDegreesFreedom,
	      &numberSlipperyNodeEntries,
	      &numberDifferentialForceEntries,
	      &f77FileInput,
	      &errorCode,
	      differentialForceInputFile,strlen(differentialForceInputFile));
    
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
		      "Scanner: Error reading differential force input file:");
      break;
      //case 4:
            //PyErr_SetString(PyExc_ValueError,
		            //"Differential force applied to non-slippery node:");
            //break;
    default:
      PyErr_SetString(PyExc_ValueError,
		      "Scanner: Unknown error reading differential force input file:");
    }

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberDifferentialForceEntries:" << numberDifferentialForceEntries
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue("i",numberDifferentialForceEntries);
}


// Read time steps at which full output is desired

char pylithomop3d_scan_fuldat__doc__[] = "";
char pylithomop3d_scan_fuldat__name__[] = "scan_fuldat";

PyObject * pylithomop3d_scan_fuldat(PyObject *, PyObject *args)
{
  int analysisTypeInt;
  int totalNumberTimeSteps;
  int f77FileInput;
  char* fullOutputInputFile;

  int ok = PyArg_ParseTuple(args, "iiis:scan_fuldat",
			    &analysisTypeInt,
			    &totalNumberTimeSteps,
			    &f77FileInput,
			    &fullOutputInputFile);

  if (!ok) {
    return 0;
  }

  int errorCode = 0;
  int numberFullOutputs = 0;

  scan_fuldat_f(&analysisTypeInt,
		&totalNumberTimeSteps,
		&numberFullOutputs,
		&f77FileInput,
		&errorCode,
		fullOutputInputFile,strlen(fullOutputInputFile));
    
  switch(errorCode)
    {
    case 0:
      break;
    case 1:
      PyErr_SetString(PyExc_IOError,
		      "Scanner: Error opening file specifying time steps for output:");
      break;
      //case 2:
            //PyErr_SetString(PyExc_ValueError,
		            //"Coordinate units not specified:");
            //break;
    case 3:
      PyErr_SetString(PyExc_IOError,
		      "Scanner: Error reading file specifying time steps for output:");
      break;
      //case 4:
            //PyErr_SetString(PyExc_ValueError,
		            //"Wrong number of nodes read:");
            //break;
    default:
      PyErr_SetString(PyExc_ValueError,
		      "Scanner: Unknown error reading time steps to output:");
    }

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberFullOutputs:" << numberFullOutputs
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue("i",numberFullOutputs);
}


// Read load histories

char pylithomop3d_scan_hist__doc__[] = "";
char pylithomop3d_scan_hist__name__[] = "scan_hist";

PyObject * pylithomop3d_scan_hist(PyObject *, PyObject *args)
{
  int f77FileInput;
  char* loadHistoryInputFile;

  int ok = PyArg_ParseTuple(args, "is:scan_hist",
			    &f77FileInput,
			    &loadHistoryInputFile);

  if (!ok) {
    return 0;
  }

  int errorCode = 0;
  int numberLoadHistories = 0;

  scan_hist_f(&numberLoadHistories,
	      &f77FileInput,
	      &errorCode,
	      loadHistoryInputFile,strlen(loadHistoryInputFile));
    
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
		      "Scanner: Error reading load history input file:");
      break;
    //case 4:
      //PyErr_SetString(PyExc_ValueError,
		      //"Load history times are out of order:");
      //break;
    default:
      PyErr_SetString(PyExc_ValueError,
		      "Scanner: Unknown error reading load history input file:");
    }

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberLoadHistories:" << numberLoadHistories
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue("i",numberLoadHistories);
}


// Read element prestresses

char pylithomop3d_scan_prestr__doc__[] = "";
char pylithomop3d_scan_prestr__name__[] = "scan_prestr";

PyObject * pylithomop3d_scan_prestr(PyObject *, PyObject *args)
{
  int numberStressComponents;
  int numberPrestressGaussPoints;
  int numberElements;
  int prestressAutoComputeInt;
  int f77FileInput;
  char* prestressInputFile;

  int ok = PyArg_ParseTuple(args, "iiiiis:scan_prestr",
			    &numberStressComponents,
			    &numberPrestressGaussPoints,
			    &numberElements,
			    &prestressAutoComputeInt,
			    &f77FileInput,
			    &prestressInputFile);

  if (!ok) {
    return 0;
  }

  int errorCode = 0;
  int numberPrestressEntries = 0;

  scan_prestr_f(&numberStressComponents,
		&numberPrestressGaussPoints,
		&numberPrestressEntries,
		&numberElements,
		&prestressAutoComputeInt,
		&f77FileInput,
		&errorCode,
		prestressInputFile,strlen(prestressInputFile));
    
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
		      "Scanner: Error reading prestress input file:");
      break;
      //case 4:
            //PyErr_SetString(PyExc_ValueError,
		            //"Wrong number of nodes read:");
            //break;
    default:
      PyErr_SetString(PyExc_ValueError,
		      "Scanner: Unknown error reading prestress input file:");
    }

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberPrestressEntries:" << numberPrestressEntries
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue("i",numberPrestressEntries);
}


// Read material properties

char pylithomop3d_scan_prop__doc__[] = "";
char pylithomop3d_scan_prop__name__[] = "scan_prop";

PyObject * pylithomop3d_scan_prop(PyObject *, PyObject *args)
{
  int numberMaterialProperties;
  int f77FileInput;
  char* densityUnits;
  char* youngUnits;
  char* viscosityCoefficientUnits;
  char* cohesionUnits;
  char* materialPropertiesInputFile;

  int ok = PyArg_ParseTuple(args, "iisssss:scan_prop",
			    &numberMaterialProperties,
			    &f77FileInput,
			    &densityUnits,
			    &youngUnits,
			    &viscosityCoefficientUnits,
			    &cohesionUnits,
			    &materialPropertiesInputFile);

  if (!ok) {
    return 0;
  }

  int errorCode = 0;
  int numberMaterialTypes = 0;
  int viscousFlagInt = 0;
  int plasticFlagInt = 0;
  int materialHistoryFlagInt = 0;

  scan_prop_f(&numberMaterialProperties,
	      &numberMaterialTypes,
	      &f77FileInput,
	      &errorCode,
	      densityUnits,
	      youngUnits,
	      viscosityCoefficientUnits,
	      cohesionUnits,
	      &viscousFlagInt,
	      &plasticFlagInt,
	      &materialHistoryFlagInt,
	      materialPropertiesInputFile,
	      strlen(densityUnits),
	      strlen(youngUnits),
	      strlen(viscosityCoefficientUnits),
	      strlen(cohesionUnits),
	      strlen(materialPropertiesInputFile));
    
  switch(errorCode)
    {
    case 0:
      break;
    case 1:
      PyErr_SetString(PyExc_IOError,
		      "Scanner: Error opening material properties input file:");
      break;
    case 2:
      PyErr_SetString(PyExc_ValueError,
		      "Scanner: All material property units not specified:");
      break;
    case 3:
      PyErr_SetString(PyExc_IOError,
		      "Scanner: Error reading material properties input file:");
      break;
      //case 4:
            //PyErr_SetString(PyExc_ValueError,
		            //"Wrong number of nodes read:");
            //break;
    default:
      PyErr_SetString(PyExc_ValueError,
		      "Scanner: Unknown error reading material properties input file:");
    }

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberMaterialTypes:" << numberMaterialTypes
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue("iiii",numberMaterialTypes,
		       viscousFlagInt,
		       plasticFlagInt,
		       materialHistoryFlagInt);
}


// Read local coordinate rotations

char pylithomop3d_scan_skew__doc__[] = "";
char pylithomop3d_scan_skew__name__[] = "scan_skew";

PyObject * pylithomop3d_scan_skew(PyObject *, PyObject *args)
{
  int numberSkewDimensions;
  int f77FileInput;
  char* rotationUnits;
  char* rotationInputFile;

  int ok = PyArg_ParseTuple(args, "iiss:scan_skew",
			    &numberSkewDimensions,
			    &f77FileInput,
			    &rotationUnits,
			    &rotationInputFile);

  if (!ok) {
    return 0;
  }

  int errorCode = 0;
  int numberRotationEntries = 0;

  scan_skew_f(&numberSkewDimensions,
	      &numberRotationEntries,
	      &f77FileInput,
	      &errorCode,
	      rotationUnits,
	      rotationInputFile,
	      strlen(rotationUnits),
	      strlen(rotationInputFile));
    
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
		      "Scanner: Rotation units not specified:");
      break;
    case 3:
      PyErr_SetString(PyExc_IOError,
		      "Scanner: Error reading rotation input file:");
      break;
      //case 4:
            //PyErr_SetString(PyExc_ValueError,
		            //"Wrong number of nodes read:");
            //break;
    default:
      PyErr_SetString(PyExc_ValueError,
		      "Scanner: Unknown error reading rotation input file:");
    }

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberRotationEntries:" << numberRotationEntries
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue("i",numberRotationEntries);
}


// Read slippery node entries

char pylithomop3d_scan_slip__doc__[] = "";
char pylithomop3d_scan_slip__name__[] = "scan_slip";

PyObject * pylithomop3d_scan_slip(PyObject *, PyObject *args)
{
  int numberDegreesFreedom;
  int f77FileInput;
  char* slipperyNodeInputFile;

  int ok = PyArg_ParseTuple(args, "iis:scan_slip",
			    &numberDegreesFreedom,
			    &f77FileInput,
			    &slipperyNodeInputFile);

  if (!ok) {
    return 0;
  }

  int errorCode = 0;
  int numberSlipperyNodeEntries = 0;

  scan_slip_f(&numberDegreesFreedom,
	      &numberSlipperyNodeEntries,
	      &f77FileInput,
	      &errorCode,
	      slipperyNodeInputFile,strlen(slipperyNodeInputFile));
    
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
		      "Scanner: Error reading slippery nodes input file:");
      break;
      //case 4:
            //PyErr_SetString(PyExc_ValueError,
		            //"Wrong number of nodes read:");
            //break;
    default:
      PyErr_SetString(PyExc_ValueError,
		      "Scanner: Unknown error reading slippery nodes input file:");
    }

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberSlipperyNodeEntries:" << numberSlipperyNodeEntries
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue("i",numberSlipperyNodeEntries);
}


// Read split node entries

char pylithomop3d_scan_split__doc__[] = "";
char pylithomop3d_scan_split__name__[] = "scan_split";

PyObject * pylithomop3d_scan_split(PyObject *, PyObject *args)
{
  int numberDegreesFreedom;
  int f77FileInput;
  char* splitNodeInputFile;

  int ok = PyArg_ParseTuple(args, "iis:scan_split",
			    &numberDegreesFreedom,
			    &f77FileInput,
			    &splitNodeInputFile);

  if (!ok) {
    return 0;
  }

  int errorCode = 0;
  int numberSplitNodeEntries = 0;

  scan_split_f(&numberDegreesFreedom,
	       &numberSplitNodeEntries,
	       &f77FileInput,
	       &errorCode,
	       splitNodeInputFile,strlen(splitNodeInputFile));
    
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
		      "Scanner: Error reading split nodes input file:");
      break;
      //case 4:
            //PyErr_SetString(PyExc_ValueError,
		            //"Wrong number of nodes read:");
            //break;
    default:
      PyErr_SetString(PyExc_ValueError,
		      "Scanner: Unknown error reading split nodes input file:");
    }

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberSplitNodeEntries:" << numberSplitNodeEntries
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue("i",numberSplitNodeEntries);
}


// Read time step data

char pylithomop3d_scan_timdat__doc__[] = "";
char pylithomop3d_scan_timdat__name__[] = "scan_timdat";

PyObject * pylithomop3d_scan_timdat(PyObject *, PyObject *args)
{
  int f77FileInput;
  char* timeUnits;
  char* timeStepInputFile;

  int ok = PyArg_ParseTuple(args, "iss:scan_timdat",
			    &f77FileInput,
			    &timeUnits,
			    &timeStepInputFile);

  if (!ok) {
    return 0;
  }

  int errorCode = 0;
  int numberTimeStepGroups = 0;
  int totalNumberTimeSteps = 0;

  scan_timdat_f(&totalNumberTimeSteps,
		&numberTimeStepGroups,
		&f77FileInput,
		&errorCode,
		timeUnits,
		timeStepInputFile,
		strlen(timeUnits),
		strlen(timeStepInputFile));
    
  switch(errorCode)
    {
    case 0:
      break;
    case 1:
      PyErr_SetString(PyExc_IOError,
		      "Scanner: Error opening time step input file:");
      break;
    case 2:
      PyErr_SetString(PyExc_ValueError,
		      "Scanner: Time units not specified:");
      break;
    case 3:
      PyErr_SetString(PyExc_IOError,
		      "Scanner: Error reading time step input file:");
      break;
      //case 4:
            //PyErr_SetString(PyExc_ValueError,
		            //"Wrong number of nodes read:");
            //break;
    default:
      PyErr_SetString(PyExc_ValueError,
		      "Scanner: Unknown error reading time step input file:");
    }

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "totalNumberTimeSteps:" << totalNumberTimeSteps
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue("ii",numberTimeStepGroups, totalNumberTimeSteps);
}


// Read traction BC

char pylithomop3d_scan_traction__doc__[] = "";
char pylithomop3d_scan_traction__name__[] = "scan_traction";

PyObject * pylithomop3d_scan_traction(PyObject *, PyObject *args)
{
  int numberElementNodes;
  int numberTractionDirections;
  char* tractionBcUnits;
  int f77FileInput;
  char* tractionInputFile;

  int ok = PyArg_ParseTuple(args, "iisis:scan_traction",
			    &numberElementNodes,
			    &numberTractionDirections,
			    &tractionBcUnits,
			    &f77FileInput,
			    &tractionInputFile);

  if (!ok) {
    return 0;
  }

  int errorCode = 0;
  int numberTractionBc = 0;

  scan_traction_f(&numberElementNodes,
		  &numberTractionDirections,
		  &numberTractionBc,
		  &f77FileInput,
		  &errorCode,
		  tractionBcUnits,
		  tractionInputFile,
		  strlen(tractionBcUnits),
		  strlen(tractionInputFile));
    
  switch(errorCode)
    {
    case 0:
      break;
      //case 1:
            //PyErr_SetString(PyExc_IOError,
		            //"Error opening traction input file:");
            //break;
    case 2:
      PyErr_SetString(PyExc_ValueError,
		      "Scanner: Traction units not specified:");
      break;
    case 3:
      PyErr_SetString(PyExc_IOError,
		      "Scanner: Error reading traction input file:");
      break;
      //case 4:
            //PyErr_SetString(PyExc_ValueError,
		            //"Wrong number of nodes read:");
            //break;
    default:
      PyErr_SetString(PyExc_ValueError,
		      "Scanner: Unknown error reading traction input file:");
    }

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberTractionBc:" << numberTractionBc
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue("i",numberTractionBc);
}


// Read winkler BC

char pylithomop3d_scan_wink__doc__[] = "";
char pylithomop3d_scan_wink__name__[] = "scan_wink";

PyObject * pylithomop3d_scan_wink(PyObject *, PyObject *args)
{
  int numberDegreesFreedom;
  int numberWinklerForces;
  int f77FileInput;
  char* winklerInputFile;

  int ok = PyArg_ParseTuple(args, "iiis:scan_wink",
			    &numberDegreesFreedom,
			    &numberWinklerForces,
			    &f77FileInput,
			    &winklerInputFile);

  if (!ok) {
    return 0;
  }

  int errorCode = 0;
  int numberWinklerEntries = 0;

  scan_wink_f(&numberDegreesFreedom,
	      &numberWinklerEntries,
	      &numberWinklerForces,
	      &f77FileInput,
	      &errorCode,
	      winklerInputFile,strlen(winklerInputFile));
    
  switch(errorCode)
    {
    case 0:
      break;
      //case 1:
            //PyErr_SetString(PyExc_IOError,
		            //"Error opening winkler input file:");
            //break;
      //case 2:
            //PyErr_SetString(PyExc_ValueError,
		            //"Coordinate units not specified:");
            //break;
    case 3:
      PyErr_SetString(PyExc_IOError,
		      "Scanner: Error reading winkler input file:");
      break;
      //case 4:
            //PyErr_SetString(PyExc_ValueError,
		            //"Wrong number of nodes read:");
            //break;
    default:
      PyErr_SetString(PyExc_ValueError,
		      "Scanner: Unknown error reading winkler input file:");
    }

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberWinklerForces:" << numberWinklerForces
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue("i",numberWinklerEntries);
}


// Read winkler BC for slippery nodes

char pylithomop3d_scan_winkx__doc__[] = "";
char pylithomop3d_scan_winkx__name__[] = "scan_winkx";

PyObject * pylithomop3d_scan_winkx(PyObject *, PyObject *args)
{
  int numberDegreesFreedom;
  int numberSlipperyWinklerForces;
  int numberSlipperyNodeEntries;
  int f77FileInput;
  char* slipperyWinklerInputFile;

  int ok = PyArg_ParseTuple(args, "iiiis:scan_winkx",
			    &numberDegreesFreedom,
			    &numberSlipperyWinklerForces,
			    &numberSlipperyNodeEntries,
			    &f77FileInput,
			    &slipperyWinklerInputFile);

  if (!ok) {
    return 0;
  }

  int errorCode = 0;
  int numberSlipperyWinklerEntries = 0;

  scan_winkx_f(&numberDegreesFreedom,
	       &numberSlipperyNodeEntries,
	       &numberSlipperyWinklerEntries,
	       &numberSlipperyWinklerForces,
	       &f77FileInput,
	       &errorCode,
	       slipperyWinklerInputFile,strlen(slipperyWinklerInputFile));
    
  switch(errorCode)
    {
    case 0:
      break;
      //case 1:
            //PyErr_SetString(PyExc_IOError,
		            //"Error opening winkler input file for slippery nodes:");
            //break;
      //case 2:
            //PyErr_SetString(PyExc_ValueError,
		            //"Coordinate units not specified:");
            //break;
    case 3:
      PyErr_SetString(PyExc_IOError,
		      "Scanner: Error reading winkler input file for slippery nodes:");
      break;
      //case 4:
            //PyErr_SetString(PyExc_ValueError,
		            //"Wrong number of nodes read:");
            //break;
    default:
      PyErr_SetString(PyExc_ValueError,
		      "Scanner: Unknown error reading winkler input file for slippery nodes:");
    }

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberSlipperyWinklerForces:" << numberSlipperyWinklerForces
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue("i",numberSlipperyWinklerEntries);
}
    
// version
// $Id: scanner.cc,v 1.1 2004/04/14 21:24:47 willic3 Exp $

// End of file
