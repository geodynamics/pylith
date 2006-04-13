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

#include "exceptionhandler.h"

int
exceptionhandler(const int errorcode, 
		 const char* errorstring)
{
  // Generate python exceptions from error codes

  const int maxsize = 4096;
  char errormsg[maxsize];

  if (errorcode == 0) {
    return 0;
  }

  // copy string to 'errormsg', trim spaces, and null-terminate
  char *dest = errormsg + maxsize;
  const char *src = errorstring + maxsize;
  while (dest > errormsg && (*--dest = *--src) == ' ')
    ;
  dest[1] = '\0';
  while (dest > errormsg)
    *--dest = *--src;

  PyObject *exception = NULL;
  const char *format = NULL;

  switch(errorcode)
    {
    // IOErrors
    case 1: exception = PyExc_IOError; format = "%s: Error opening file for reading."; break;
    case 2: exception = PyExc_IOError; format = "%s: Error opening file for writing."; break;
    case 3: exception = PyExc_IOError; format = "%s: Error reading from file."; break;
    case 4: exception = PyExc_IOError; format = "%s: Error writing to file."; break;
    // ValueErrors
    case   5: exception = PyExc_ValueError; format = "%s: Invalid or missing units specification."; break;
    case 100: exception = PyExc_ValueError; format = "%s: Attempt to use undefined history."; break;
    case 101: exception = PyExc_ValueError; format = "%s: Attempt to use undefined material model."; break;
    case 102: exception = PyExc_ValueError; format = "%s: Invalid number of state variables."; break;
    case 103: exception = PyExc_ValueError; format = "%s: Zero or negative diagonal in stiffness matrix."; break;
    case 104: exception = PyExc_ValueError; format = "%s: BC assigned for nonexistent node."; break;
    case 105: exception = PyExc_ValueError; format = "%s: Wrong number of nodes read."; break;
    case 106: exception = PyExc_ValueError; format = "%s: Invalid element type specified."; break;
    case 107: exception = PyExc_ValueError; format = "%s: Invalid material type specified."; break;
    case 108: exception = PyExc_ValueError; format = "%s: Undefined node used in connectivity."; break;
    case 109: exception = PyExc_ValueError; format = "%s: Differential force applied to non-slippery node."; break;
    case 110: exception = PyExc_ValueError; format = "%s: Matrix not positive-definite."; break;
    case 111: exception = PyExc_ValueError; format = "%s: Load history times are out of order."; break;
    case 112: exception = PyExc_ValueError; format = "%s: Not enough points to define a plane for slippery nodes."; break;
    case 113: exception = PyExc_ValueError; format = "%s: Zero or negative jacobian for element."; break;
    // MemoryErrors
    case 300: exception = PyExc_MemoryError; format = "%s: Insufficient memory assigned."; break;
    // binary I/O
    case 400: exception = PyExc_RuntimeError; format = "%s: Don't know how to open a binary file."; break;
    case 402: exception = PyExc_RuntimeError; format = "%s: Error opening binary file for writing."; break;
    case 403: exception = PyExc_RuntimeError; format = "%s: Error reading from binary file."; break;
    case 404: exception = PyExc_RuntimeError; format = "%s: Error writing to binary file."; break;
    case 405: exception = PyExc_RuntimeError; format = "%s: Unexpected 'inquire' result for binary file."; break;
    case 406: exception = PyExc_RuntimeError; format = "%s: Binary data read does not match data written."; break;
    // ???
    default: exception = PyExc_RuntimeError; format = "%s: Unknown error."; break;
    }
  
  PyErr_Format(exception, format, errormsg);
  
  return errorcode;
}
    
// version
// $Id: exceptionhandler.cc,v 1.3 2005/03/31 23:27:57 willic3 Exp $

// End of file
