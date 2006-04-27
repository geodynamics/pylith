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
