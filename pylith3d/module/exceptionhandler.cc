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
    case 114: exception = PyExc_ValueError; format = "%s: Initial bracketing values are identical."; break;
    case 115: exception = PyExc_ValueError; format = "%s: Root is not initially bracketed."; break;
    // MemoryErrors
    case 300: exception = PyExc_MemoryError; format = "%s: Insufficient memory assigned."; break;
    // binary I/O
    case 400: exception = PyExc_RuntimeError; format = "%s: Don't know how to open a binary file."; break;
    case 402: exception = PyExc_RuntimeError; format = "%s: Error opening binary file for writing."; break;
    case 403: exception = PyExc_RuntimeError; format = "%s: Error reading from binary file."; break;
    case 404: exception = PyExc_RuntimeError; format = "%s: Error writing to binary file."; break;
    case 405: exception = PyExc_RuntimeError; format = "%s: Unexpected 'inquire' result for binary file."; break;
    case 406: exception = PyExc_RuntimeError; format = "%s: Binary data read does not match data written."; break;
    case 410: exception = PyExc_RuntimeError; format = "%s: Bracketing values not found."; break;
    case 411: exception = PyExc_RuntimeError; format = "%s: Maximum iterations exceeded."; break;
    // ???
    default: exception = PyExc_RuntimeError; format = "%s: Unknown error."; break;
    }
  
  PyErr_Format(exception, format, errormsg);
  
  return errorcode;
}
    
// version
// $Id: exceptionhandler.cc,v 1.3 2005/03/31 23:27:57 willic3 Exp $

// End of file
