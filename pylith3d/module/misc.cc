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

#include "misc.h"
#include "exceptionhandler.h"
#include "pylith3d_externs.h"


// try_binio

char pypylith3d_try_binio__doc__[] = "";
char pypylith3d_try_binio__name__[] = "try_binio";

PyObject * pypylith3d_try_binio(PyObject *, PyObject *args)
{
  int unit;

  int ok = PyArg_ParseTuple(args, "i:try_binio",
			    &unit);
  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];

  try_binio_f(&unit,
	      &errorcode,
	      errorstring,
	      sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  // return
  Py_INCREF(Py_None);
  return Py_None;
}
    

// copyright

char pypylith3d_copyright__doc__[] = "";
char pypylith3d_copyright__name__[] = "copyright";

static char pypylith3d_copyright_note[] = 
    "pylith3d python module: Copyright (c) 2005 Rensselaer Polytechnic Institute";


PyObject * pypylith3d_copyright(PyObject *, PyObject *)
{
    return Py_BuildValue("s", pypylith3d_copyright_note);
}
    
// version
// $Id: misc.cc,v 1.2 2005/03/31 23:27:57 willic3 Exp $

// End of file
