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

#if !defined(pylithomop3d_scanner_h)
#define pylithomop3d_scanner_h

// scan boundary condition file
extern char pylithomop3d_scan_bc__name__[];
extern char pylithomop3d_scan_bc__doc__[];
extern "C"
PyObject * pylithomop3d_scan_bc(PyObject *, PyObject *);

// scan connectivity file
extern char pylithomop3d_scan_connect__name__[];
extern char pylithomop3d_scan_connect__doc__[];
extern "C"
PyObject * pylithomop3d_scan_connect(PyObject *, PyObject *);

// scan coordinates file
extern char pylithomop3d_scan_coords__name__[];
extern char pylithomop3d_scan_coords__doc__[];
extern "C"
PyObject * pylithomop3d_scan_coords(PyObject *, PyObject *);

// scan differential forces file
extern char pylithomop3d_scan_diff__name__[];
extern char pylithomop3d_scan_diff__doc__[];
extern "C"
PyObject * pylithomop3d_scan_diff(PyObject *, PyObject *);

// scan time step output info file
extern char pylithomop3d_scan_fuldat__name__[];
extern char pylithomop3d_scan_fuldat__doc__[];
extern "C"
PyObject * pylithomop3d_scan_fuldat(PyObject *, PyObject *);

// scan time history definition file
extern char pylithomop3d_scan_hist__name__[];
extern char pylithomop3d_scan_hist__doc__[];
extern "C"
PyObject * pylithomop3d_scan_hist(PyObject *, PyObject *);

// scan prestress file
// extern char pylithomop3d_scan_prestr__name__[];
// extern char pylithomop3d_scan_prestr__doc__[];
// extern "C"
// PyObject * pylithomop3d_scan_prestr(PyObject *, PyObject *);

// scan local coordinate rotations file
extern char pylithomop3d_scan_skew__name__[];
extern char pylithomop3d_scan_skew__doc__[];
extern "C"
PyObject * pylithomop3d_scan_skew(PyObject *, PyObject *);

// scan slippery node definitions file
extern char pylithomop3d_scan_slip__name__[];
extern char pylithomop3d_scan_slip__doc__[];
extern "C"
PyObject * pylithomop3d_scan_slip(PyObject *, PyObject *);

// scan split node definitions file
extern char pylithomop3d_scan_split__name__[];
extern char pylithomop3d_scan_split__doc__[];
extern "C"
PyObject * pylithomop3d_scan_split(PyObject *, PyObject *);

// scan time step group info file
extern char pylithomop3d_scan_timdat__name__[];
extern char pylithomop3d_scan_timdat__doc__[];
extern "C"
PyObject * pylithomop3d_scan_timdat(PyObject *, PyObject *);

// scan traction boundary conditions file
// extern char pylithomop3d_scan_traction__name__[];
// extern char pylithomop3d_scan_traction__doc__[];
// extern "C"
// PyObject * pylithomop3d_scan_traction(PyObject *, PyObject *);

// scan winkler forces info file
extern char pylithomop3d_scan_wink__name__[];
extern char pylithomop3d_scan_wink__doc__[];
extern "C"
PyObject * pylithomop3d_scan_wink(PyObject *, PyObject *);

// scan winkler forces info file for slippery nodes
extern char pylithomop3d_scan_winkx__name__[];
extern char pylithomop3d_scan_winkx__doc__[];
extern "C"
PyObject * pylithomop3d_scan_winkx(PyObject *, PyObject *);

#endif

// version
// $Id: scanner.h,v 1.3 2005/03/31 23:27:58 willic3 Exp $

// End of file
