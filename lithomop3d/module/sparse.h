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

#if !defined(pylithomop3d_sparse_h)
#define pylithomop3d_sparse_h

// compute maximum number of nonzero entries in stiffness matrix
extern char pylithomop3d_cmp_stiffsz__name__[];
extern char pylithomop3d_cmp_stiffsz__doc__[];
extern "C"
PyObject * pylithomop3d_cmp_stiffsz(PyObject *, PyObject *);

// create linked list for sparse matrix
extern char pylithomop3d_lnklst__name__[];
extern char pylithomop3d_lnklst__doc__[];
extern "C"
PyObject * pylithomop3d_lnklst(PyObject *, PyObject *);

// create a PETSc Mat
extern char pylithomop3d_createPETScMat__name__[];
extern char pylithomop3d_createPETScMat__doc__[];
extern "C"
PyObject * pylithomop3d_createPETScMat(PyObject *, PyObject *);

// create a PETSc Mat
extern char pylithomop3d_destroyPETScMat__name__[];
extern char pylithomop3d_destroyPETScMat__doc__[];
extern "C"
PyObject * pylithomop3d_destroyPETScMat(PyObject *, PyObject *);

// create sparse matrix in modified sparse row format
extern char pylithomop3d_makemsr__name__[];
extern char pylithomop3d_makemsr__doc__[];
extern "C"
PyObject * pylithomop3d_makemsr(PyObject *, PyObject *);

#endif

// version
// $Id: sparse.h,v 1.3 2005/03/10 01:10:39 knepley Exp $

// End of file
