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

#if !defined(pylithomop3d_numbering_h)
#define pylithomop3d_numbering_h

// create id array
extern char pylithomop3d_create_id__name__[];
extern char pylithomop3d_create_id__doc__[];
extern "C"
PyObject * pylithomop3d_create_id(PyObject *, PyObject *);

// form id array for split nodes
extern char pylithomop3d_id_split__name__[];
extern char pylithomop3d_id_split__doc__[];
extern "C"
PyObject * pylithomop3d_id_split(PyObject *, PyObject *);

// localize id array
extern char pylithomop3d_local__name__[];
extern char pylithomop3d_local__doc__[];
extern "C"
PyObject * pylithomop3d_local(PyObject *, PyObject *);

// localize id array for split nodes
extern char pylithomop3d_localf__name__[];
extern char pylithomop3d_localf__doc__[];
extern "C"
PyObject * pylithomop3d_localf(PyObject *, PyObject *);

// localize id array for slippery nodes
extern char pylithomop3d_localx__name__[];
extern char pylithomop3d_localx__doc__[];
extern "C"
PyObject * pylithomop3d_localx(PyObject *, PyObject *);

// find closest fault neignbors for slippery nodes
extern char pylithomop3d_nfind__name__[];
extern char pylithomop3d_nfind__doc__[];
extern "C"
PyObject * pylithomop3d_nfind(PyObject *, PyObject *);

#endif

// version
// $Id: numbering.h,v 1.2 2005/03/31 23:27:58 willic3 Exp $

// End of file
