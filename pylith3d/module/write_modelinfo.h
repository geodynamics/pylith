// -*- C++ -*-
// 
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  PyLith by Charles A. Williams, Brad Aagaard, and Matt Knepley
//
//  Copyright (c) 2004-2006 Rensselaer Polytechnic Institute
//
//  Permission is hereby granted, free of charge, to any person obtaining
//  a copy of this software and associated documentation files (the
//  "Software"), to deal in the Software without restriction, including
//  without limitation the rights to use, copy, modify, merge, publish,
//  distribute, sublicense, and/or sell copies of the Software, and to
//  permit persons to whom the Software is furnished to do so, subject to
//  the following conditions:
//
//  The above copyright notice and this permission notice shall be
//  included in all copies or substantial portions of the Software.
//
//  THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
//  EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
//  MERCHANTABILITY,    FITNESS    FOR    A   PARTICULAR    PURPOSE    AND
//  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
//  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
//  OF CONTRACT, TORT OR OTHERWISE,  ARISING FROM, OUT OF OR IN CONNECTION
//  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 

#if !defined(pypylith3d_write_modelinfo_h)
#define pypylith3d_write_modelinfo_h

// write boundary conditions
extern char pypylith3d_write_bc__name__[];
extern char pypylith3d_write_bc__doc__[];
extern "C"
PyObject * pypylith3d_write_bc(PyObject *, PyObject *);

// write element node array
extern char pypylith3d_write_connect__name__[];
extern char pypylith3d_write_connect__doc__[];
extern "C"
PyObject * pypylith3d_write_connect(PyObject *, PyObject *);

// write nodal coordinates
extern char pypylith3d_write_coords__name__[];
extern char pypylith3d_write_coords__doc__[];
extern "C"
PyObject * pypylith3d_write_coords(PyObject *, PyObject *);

// write differential forces
extern char pypylith3d_write_diff__name__[];
extern char pypylith3d_write_diff__doc__[];
extern "C"
PyObject * pypylith3d_write_diff(PyObject *, PyObject *);

// write out element and prestress information
extern char pypylith3d_write_element_info__name__[];
extern char pypylith3d_write_element_info__doc__[];
extern "C"
PyObject * pypylith3d_write_element_info(PyObject *, PyObject *);

// write time steps at which full output is desired
extern char pypylith3d_write_fuldat__name__[];
extern char pypylith3d_write_fuldat__doc__[];
extern "C"
PyObject * pypylith3d_write_fuldat(PyObject *, PyObject *);

// write out global control information
extern char pypylith3d_write_global_info__name__[];
extern char pypylith3d_write_global_info__doc__[];
extern "C"
PyObject * pypylith3d_write_global_info(PyObject *, PyObject *);

// write load histories
extern char pypylith3d_write_hist__name__[];
extern char pypylith3d_write_hist__doc__[];
extern "C"
PyObject * pypylith3d_write_hist(PyObject *, PyObject *);

// write out material property info
extern char pypylith3d_write_props__name__[];
extern char pypylith3d_write_props__doc__[];
extern "C"
PyObject * pypylith3d_write_props(PyObject *, PyObject *);

// write local coordinate rotations
extern char pypylith3d_write_skew__name__[];
extern char pypylith3d_write_skew__doc__[];
extern "C"
PyObject * pypylith3d_write_skew(PyObject *, PyObject *);

// write slippery node entries
extern char pypylith3d_write_slip__name__[];
extern char pypylith3d_write_slip__doc__[];
extern "C"
PyObject * pypylith3d_write_slip(PyObject *, PyObject *);

// write out sparse matrix information
extern char pypylith3d_write_sparse_info__name__[];
extern char pypylith3d_write_sparse_info__doc__[];
extern "C"
PyObject * pypylith3d_write_sparse_info(PyObject *, PyObject *);

// write split node entries
extern char pypylith3d_write_split__name__[];
extern char pypylith3d_write_split__doc__[];
extern "C"
PyObject * pypylith3d_write_split(PyObject *, PyObject *);

// write split node entries for plot output
extern char pypylith3d_write_split_plot__name__[];
extern char pypylith3d_write_split_plot__doc__[];
extern "C"
PyObject * pypylith3d_write_split_plot(PyObject *, PyObject *);

// write state output information
extern char pypylith3d_write_stateout__name__[];
extern char pypylith3d_write_stateout__doc__[];
extern "C"
PyObject * pypylith3d_write_stateout(PyObject *, PyObject *);

// write out stress computation information
extern char pypylith3d_write_strscomp__name__[];
extern char pypylith3d_write_strscomp__doc__[];
extern "C"
PyObject * pypylith3d_write_strscomp(PyObject *, PyObject *);

// write out subiteration convergence information
extern char pypylith3d_write_subiter__name__[];
extern char pypylith3d_write_subiter__doc__[];
extern "C"
PyObject * pypylith3d_write_subiter(PyObject *, PyObject *);

// write time step data
extern char pypylith3d_write_timdat__name__[];
extern char pypylith3d_write_timdat__doc__[];
extern "C"
PyObject * pypylith3d_write_timdat(PyObject *, PyObject *);

// write mesh info to UCD file
extern char pypylith3d_write_ucd_mesh__name__[];
extern char pypylith3d_write_ucd_mesh__doc__[];
extern "C"
PyObject * pypylith3d_write_ucd_mesh(PyObject *, PyObject *);

// write winkler BC
extern char pypylith3d_write_wink__name__[];
extern char pypylith3d_write_wink__doc__[];
extern "C"
PyObject * pypylith3d_write_wink(PyObject *, PyObject *);

// write slippery winkler BC
extern char pypylith3d_write_winkx__name__[];
extern char pypylith3d_write_winkx__doc__[];
extern "C"
PyObject * pypylith3d_write_winkx(PyObject *, PyObject *);

#endif

// version
// $Id: write_modelinfo.h,v 1.1 2005/04/21 00:04:36 willic3 Exp $

// End of file
