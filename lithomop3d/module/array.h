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

#if !defined(pylithomop3d_array_h)
#define pylithomop3d_array_h

// allocateInt
extern char pylithomop3d_allocateInt__doc__[];
extern char pylithomop3d_allocateInt__name__[];

extern "C"
PyObject * pylithomop3d_allocateInt(PyObject *, PyObject *);

// allocateDouble
extern char pylithomop3d_allocateDouble__doc__[];
extern char pylithomop3d_allocateDouble__name__[];

extern "C"
PyObject * pylithomop3d_allocateDouble(PyObject *, PyObject *);

// intListToArray
extern char pylithomop3d_intListToArray__doc__[];
extern char pylithomop3d_intListToArray__name__[];

extern "C"
PyObject * pylithomop3d_intListToArray(PyObject *, PyObject *);

// doubleListToArray
extern char pylithomop3d_doubleListToArray__doc__[];
extern char pylithomop3d_doubleListToArray__name__[];

extern "C"
PyObject * pylithomop3d_doubleListToArray(PyObject *, PyObject *);

// clear p
void clearp(void *);

#endif


// version
// $Id: array.h,v 1.1 2004/04/14 21:24:47 willic3 Exp $

// End of file
