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

#if !defined(pylithomop3d_parser_h)
#define pylithomop3d_parser_h

// read boundary conditions
extern char pylithomop3d_read_bc__name__[];
extern char pylithomop3d_read_bc__doc__[];
extern "C"
PyObject * pylithomop3d_read_bc(PyObject *, PyObject *);

// read connectivities
extern char pylithomop3d_read_connect__name__[];
extern char pylithomop3d_read_connect__doc__[];
extern "C"
PyObject * pylithomop3d_read_connect(PyObject *, PyObject *);

// read coordinates
extern char pylithomop3d_read_coords__name__[];
extern char pylithomop3d_read_coords__doc__[];
extern "C"
PyObject * pylithomop3d_read_coords(PyObject *, PyObject *);

// read differential forces
extern char pylithomop3d_read_diff__name__[];
extern char pylithomop3d_read_diff__doc__[];
extern "C"
PyObject * pylithomop3d_read_diff(PyObject *, PyObject *);

// read time steps at which full output is desired
extern char pylithomop3d_read_fuldat__name__[];
extern char pylithomop3d_read_fuldat__doc__[];
extern "C"
PyObject * pylithomop3d_read_fuldat(PyObject *, PyObject *);

// read time history info
extern char pylithomop3d_read_hist__name__[];
extern char pylithomop3d_read_hist__doc__[];
extern "C"
PyObject * pylithomop3d_read_hist(PyObject *, PyObject *);

// read material history info
extern char pylithomop3d_read_mathist__name__[];
extern char pylithomop3d_read_mathist__doc__[];
extern "C"
PyObject * pylithomop3d_read_mathist(PyObject *, PyObject *);

// read prestresses
// extern char pylithomop3d_read_prestr__name__[];
// extern char pylithomop3d_read_prestr__doc__[];
// extern "C"
// PyObject * pylithomop3d_read_prestr(PyObject *, PyObject *);

// read local coordinate rotations
extern char pylithomop3d_read_skew__name__[];
extern char pylithomop3d_read_skew__doc__[];
extern "C"
PyObject * pylithomop3d_read_skew(PyObject *, PyObject *);

// read slippery node definitions
extern char pylithomop3d_read_slip__name__[];
extern char pylithomop3d_read_slip__doc__[];
extern "C"
PyObject * pylithomop3d_read_slip(PyObject *, PyObject *);

// read split node definitions
extern char pylithomop3d_read_split__name__[];
extern char pylithomop3d_read_split__doc__[];
extern "C"
PyObject * pylithomop3d_read_split(PyObject *, PyObject *);

// read state output information
extern char pylithomop3d_read_stateout__name__[];
extern char pylithomop3d_read_stateout__doc__[];
extern "C"
PyObject * pylithomop3d_read_stateout(PyObject *, PyObject *);

// read time step group informations
extern char pylithomop3d_read_timdat__name__[];
extern char pylithomop3d_read_timdat__doc__[];
extern "C"
PyObject * pylithomop3d_read_timdat(PyObject *, PyObject *);

// read traction boundary conditions
// extern char pylithomop3d_read_traction__name__[];
// extern char pylithomop3d_read_traction__doc__[];
// extern "C"
// PyObject * pylithomop3d_read_traction(PyObject *, PyObject *);

// read winkler force information
extern char pylithomop3d_read_wink__name__[];
extern char pylithomop3d_read_wink__doc__[];
extern "C"
PyObject * pylithomop3d_read_wink(PyObject *, PyObject *);

// read winkler force information for slippery nodes
extern char pylithomop3d_read_winkx__name__[];
extern char pylithomop3d_read_winkx__doc__[];
extern "C"
PyObject * pylithomop3d_read_winkx(PyObject *, PyObject *);

#endif

// version
// $Id: parser.h,v 1.2 2004/07/20 00:44:48 willic3 Exp $

// End of file
