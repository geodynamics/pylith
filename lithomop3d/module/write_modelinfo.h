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

#if !defined(pylithomop3d_write_modelinfo_h)
#define pylithomop3d_write_modelinfo_h

// write boundary conditions
extern char pylithomop3d_write_bc__name__[];
extern char pylithomop3d_write_bc__doc__[];
extern "C"
PyObject * pylithomop3d_write_bc(PyObject *, PyObject *);

// write element node array
extern char pylithomop3d_write_connect__name__[];
extern char pylithomop3d_write_connect__doc__[];
extern "C"
PyObject * pylithomop3d_write_connect(PyObject *, PyObject *);

// write nodal coordinates
extern char pylithomop3d_write_coords__name__[];
extern char pylithomop3d_write_coords__doc__[];
extern "C"
PyObject * pylithomop3d_write_coords(PyObject *, PyObject *);

// write differential forces
extern char pylithomop3d_write_diff__name__[];
extern char pylithomop3d_write_diff__doc__[];
extern "C"
PyObject * pylithomop3d_write_diff(PyObject *, PyObject *);

// write out element and prestress information
extern char pylithomop3d_write_element_info__name__[];
extern char pylithomop3d_write_element_info__doc__[];
extern "C"
PyObject * pylithomop3d_write_element_info(PyObject *, PyObject *);

// write time steps at which full output is desired
extern char pylithomop3d_write_fuldat__name__[];
extern char pylithomop3d_write_fuldat__doc__[];
extern "C"
PyObject * pylithomop3d_write_fuldat(PyObject *, PyObject *);

// write out global control information
extern char pylithomop3d_write_global_info__name__[];
extern char pylithomop3d_write_global_info__doc__[];
extern "C"
PyObject * pylithomop3d_write_global_info(PyObject *, PyObject *);

// write load histories
extern char pylithomop3d_write_hist__name__[];
extern char pylithomop3d_write_hist__doc__[];
extern "C"
PyObject * pylithomop3d_write_hist(PyObject *, PyObject *);

// write out material property info
extern char pylithomop3d_write_props__name__[];
extern char pylithomop3d_write_props__doc__[];
extern "C"
PyObject * pylithomop3d_write_props(PyObject *, PyObject *);

// write local coordinate rotations
extern char pylithomop3d_write_skew__name__[];
extern char pylithomop3d_write_skew__doc__[];
extern "C"
PyObject * pylithomop3d_write_skew(PyObject *, PyObject *);

// write slippery node entries
extern char pylithomop3d_write_slip__name__[];
extern char pylithomop3d_write_slip__doc__[];
extern "C"
PyObject * pylithomop3d_write_slip(PyObject *, PyObject *);

// write out sparse matrix information
extern char pylithomop3d_write_sparse_info__name__[];
extern char pylithomop3d_write_sparse_info__doc__[];
extern "C"
PyObject * pylithomop3d_write_sparse_info(PyObject *, PyObject *);

// write split node entries
extern char pylithomop3d_write_split__name__[];
extern char pylithomop3d_write_split__doc__[];
extern "C"
PyObject * pylithomop3d_write_split(PyObject *, PyObject *);

// write split node entries for plot output
extern char pylithomop3d_write_split_plot__name__[];
extern char pylithomop3d_write_split_plot__doc__[];
extern "C"
PyObject * pylithomop3d_write_split_plot(PyObject *, PyObject *);

// write state output information
extern char pylithomop3d_write_stateout__name__[];
extern char pylithomop3d_write_stateout__doc__[];
extern "C"
PyObject * pylithomop3d_write_stateout(PyObject *, PyObject *);

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

// write time step data
extern char pylithomop3d_write_timdat__name__[];
extern char pylithomop3d_write_timdat__doc__[];
extern "C"
PyObject * pylithomop3d_write_timdat(PyObject *, PyObject *);

// write mesh info to UCD file
extern char pylithomop3d_write_ucd_mesh__name__[];
extern char pylithomop3d_write_ucd_mesh__doc__[];
extern "C"
PyObject * pylithomop3d_write_ucd_mesh(PyObject *, PyObject *);

// write winkler BC
extern char pylithomop3d_write_wink__name__[];
extern char pylithomop3d_write_wink__doc__[];
extern "C"
PyObject * pylithomop3d_write_wink(PyObject *, PyObject *);

// write slippery winkler BC
extern char pylithomop3d_write_winkx__name__[];
extern char pylithomop3d_write_winkx__doc__[];
extern "C"
PyObject * pylithomop3d_write_winkx(PyObject *, PyObject *);

#endif

// version
// $Id: write_modelinfo.h,v 1.1 2005/04/21 00:04:36 willic3 Exp $

// End of file
