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

#if !defined(pylithomop3d_input_misc_h)
#define pylithomop3d_input_misc_h

// write out element and prestress information
extern char pylithomop3d_write_element_info__name__[];
extern char pylithomop3d_write_element_info__doc__[];
extern "C"
PyObject * pylithomop3d_write_element_info(PyObject *, PyObject *);

// write out global control information
extern char pylithomop3d_write_global_info__name__[];
extern char pylithomop3d_write_global_info__doc__[];
extern "C"
PyObject * pylithomop3d_write_global_info(PyObject *, PyObject *);

// write out material property info
extern char pylithomop3d_write_props__name__[];
extern char pylithomop3d_write_props__doc__[];
extern "C"
PyObject * pylithomop3d_write_props(PyObject *, PyObject *);

// write out sparse matrix information
extern char pylithomop3d_write_sparse_info__name__[];
extern char pylithomop3d_write_sparse_info__doc__[];
extern "C"
PyObject * pylithomop3d_write_sparse_info(PyObject *, PyObject *);

// write out stress computation information
extern char pylithomop3d_write_strscomp__name__[];
extern char pylithomop3d_write_strscomp__doc__[];
extern "C"
PyObject * pylithomop3d_write_strscomp(PyObject *, PyObject *);

// write out subiteration convergence information
extern char pylithomop3d_write_subiter__name__[];
extern char pylithomop3d_write_subiter__doc__[];
extern "C"
PyObject * pylithomop3d_write_subiter(PyObject *, PyObject *);

// write mesh info to UCD file
extern char pylithomop3d_write_ucd_mesh__name__[];
extern char pylithomop3d_write_ucd_mesh__doc__[];
extern "C"
PyObject * pylithomop3d_write_ucd_mesh(PyObject *, PyObject *);

#endif

// version
// $Id: input_misc.h,v 1.5 2005/03/31 23:27:58 willic3 Exp $

// End of file
