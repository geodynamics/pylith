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
#include <Python.h>
#include <stdio.h>

#include <string.h> // USES strncpy()

#include "exceptionhandler.h"

char*
safestrcat(char* dest, const char* src, const int destsize)
{
  // Perform a safe string concatenation. If the destination string is
  // too small to append all of the src string, only the portion of
  // the src string that will fit is appended.

  // maximum number of chars we can add to dest string from src string
  // 1 accounts for terminating '\0' not included in strlen
  const int maxadd = destsize - (1+strlen(dest));
  if (strlen(src) < maxadd)
    // we can safely add all chars from src to dest
    strcat(dest, src);
  else
    // we can only add maxadd chars from src to dest, mesg will be truncated
    strncat(dest, src, maxadd);
  return dest;
}

int
exceptionhandler(const int errorcode, 
		 const char* errorstring)
{
  // Generate python exceptions from error codes

  const int maxsize = 4096;
  char errormsg[maxsize];
  // copy at most 4095 chars (account for terminating '\0')
  strncpy(errormsg, errorstring, maxsize-1);

  switch(errorcode)
    {
    case 0:
      break;
    case 1:
      PyErr_SetString(PyExc_IOError, 
		      safestrcat(errormsg,
				 ": Error opening file for reading.",
				 maxsize));
      break;
    case 2:
      PyErr_SetString(PyExc_IOError, 
		      safestrcat(errormsg,
				 ": Error opening file for writing.",
				 maxsize));
      break;
    case 3:
      PyErr_SetString(PyExc_IOError, 
		      safestrcat(errormsg,
				 ": Error reading from file.",
				 maxsize));
      break;
    case 4:
      PyErr_SetString(PyExc_IOError, 
		      safestrcat(errormsg,
				 ": Error writing to file.",
				 maxsize));
      break;
    case 5:
      PyErr_SetString(PyExc_IOError, 
		      safestrcat(errormsg,
				 ": Invalid or missing units specification.",
				 maxsize));
      break;
    case 100:
      PyErr_SetString(PyExc_ValueError, 
		      safestrcat(errormsg,
				 ": Attempt to use undefined history.",
				 maxsize));
      break;
    case 101:
      PyErr_SetString(PyExc_ValueError, 
		      safestrcat(errormsg,
				 ": Attempt to use undefined material model.",
				 maxsize));
      break;
    case 102:
      PyErr_SetString(PyExc_ValueError, 
		      safestrcat(errormsg,
				 ": Invalid number of state variables.",
				 maxsize));
      break;
    case 103:
      PyErr_SetString(PyExc_ValueError, 
		      safestrcat(errormsg,
				 ": Zero or negative diagonal in stiffness matrix.",
				 maxsize));
      break;
    case 104:
      PyErr_SetString(PyExc_ValueError, 
		      safestrcat(errormsg,
				 ": BC assigned for nonexistent node.",
				 maxsize));
      break;
    case 105:
      PyErr_SetString(PyExc_ValueError, 
		      safestrcat(errormsg,
				 ": Wrong number of nodes read.",
				 maxsize));
      break;
    case 106:
      PyErr_SetString(PyExc_ValueError, 
		      safestrcat(errormsg,
				 ": Invalid element type specified.",
				 maxsize));
      break;
    case 107:
      PyErr_SetString(PyExc_ValueError, 
		      safestrcat(errormsg,
				 ": Invalid material type specified.",
				 maxsize));
      break;
    case 108:
      PyErr_SetString(PyExc_ValueError, 
		      safestrcat(errormsg,
				 ": Undefined node used in connectivity.",
				 maxsize));
      break;
    case 109:
      PyErr_SetString(PyExc_ValueError, 
		      safestrcat(errormsg,
				 ": Differential force applied to non-slippery node.",
				 maxsize));
      break;
    case 110:
      PyErr_SetString(PyExc_ValueError, 
		      safestrcat(errormsg,
				 ": Matrix not positive-definite.",
				 maxsize));
      break;
    case 111:
      PyErr_SetString(PyExc_ValueError, 
		      safestrcat(errormsg,
				 ": Load history times are out of order.",
				 maxsize));
      break;
    case 112:
      PyErr_SetString(PyExc_ValueError, 
		      safestrcat(errormsg,
				 ": Not enough points to define a plane for slippery nodes.",
				 maxsize));
      break;
    case 113:
      PyErr_SetString(PyExc_ValueError, 
		      safestrcat(errormsg,
				 ": Zero or negative jacobian for element.",
				 maxsize));
      break;
    case 300:
      PyErr_SetString(PyExc_MemoryError, 
		      safestrcat(errormsg,
				 ": Insufficient memory assigned.",
				 maxsize));
      break;
    default:
	PyErr_SetString(PyExc_ValueError, 
			safestrcat(errormsg,
				   ": Unknown error.",
				   maxsize));
	break;
    }
  return errorcode;
}
    
// version
// $Id: exceptionhandler.cc,v 1.3 2005/03/31 23:27:57 willic3 Exp $

// End of file
