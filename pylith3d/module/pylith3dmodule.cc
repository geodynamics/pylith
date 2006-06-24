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

#include "pylith3dmodule.h"
#include "exceptions.h"
#include "bindings.h"


char pypylith3d_module__doc__[] = "";

void pypylith3d_init(const char *name)
{
  // create the module and add the functions
  PyObject * m = Py_InitModule4(
				(char *)name, pypylith3d_methods,
				pypylith3d_module__doc__, 0, PYTHON_API_VERSION);

  // get its dictionary
  PyObject * d = PyModule_GetDict(m);

  // check for errors
  if (PyErr_Occurred()) {
    Py_FatalError("can't initialize module pylith3d");
  }

  // install the module exceptions
  pypylith3d_runtimeError = PyErr_NewException("pylith3d.runtime", 0, 0);
  PyDict_SetItemString(d, "RuntimeException", pypylith3d_runtimeError);

  return;
}

// Initialization function for the module (*must* be called initpylith3d)
extern "C"
void
initpylith3d()
{
  pypylith3d_init("pylith3d");
  return;
}

// version
// $Id: pylith3dmodule.cc,v 1.2 2005/03/31 23:27:57 willic3 Exp $

// End of file
