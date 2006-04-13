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

#if !defined(pypylith3d_setup_h)
#define pypylith3d_setup_h

// Initialize PETSc
extern char pypylith3d_petsc_initialize__doc__[];
extern char pypylith3d_petsc_initialize__name__[];
extern "C"
PyObject * pypylith3d_petsc_initialize(PyObject *, PyObject *);

// Finalize PETSc
extern char pypylith3d_petsc_finalize__doc__[];
extern char pypylith3d_petsc_finalize__name__[];
extern "C"
PyObject * pypylith3d_petsc_finalize(PyObject *, PyObject *);

// Setup PETSc Logging
extern char pypylith3d_setup_petsc_logging__doc__[];
extern char pypylith3d_setup_petsc_logging__name__[];
extern "C"
PyObject * pypylith3d_setup_petsc_logging(PyObject *, PyObject *);

// Initialize material model info
extern char pypylith3d_matmod_def__name__[];
extern char pypylith3d_matmod_def__doc__[];
extern "C"
PyObject * pypylith3d_matmod_def(PyObject *, PyObject *);

// Precompute shape function info
extern char pypylith3d_preshape__name__[];
extern char pypylith3d_preshape__doc__[];
extern "C"
PyObject * pypylith3d_preshape(PyObject *, PyObject *);

#endif

// version
// $Id: setup.h,v 1.3 2005/03/31 23:27:58 willic3 Exp $

// End of file
