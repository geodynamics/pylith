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
