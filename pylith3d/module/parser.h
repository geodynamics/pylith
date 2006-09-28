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

#if !defined(pypylith3d_parser_h)
#define pypylith3d_parser_h

// read boundary conditions
extern char pypylith3d_read_bc__name__[];
extern char pypylith3d_read_bc__doc__[];
extern "C"
PyObject * pypylith3d_read_bc(PyObject *, PyObject *);

// read connectivities
extern char pypylith3d_read_connect__name__[];
extern char pypylith3d_read_connect__doc__[];
extern "C"
PyObject * pypylith3d_read_connect(PyObject *, PyObject *);

// read coordinates
extern char pypylith3d_read_coords__name__[];
extern char pypylith3d_read_coords__doc__[];
extern "C"
PyObject * pypylith3d_read_coords(PyObject *, PyObject *);

// read differential forces
extern char pypylith3d_read_diff__name__[];
extern char pypylith3d_read_diff__doc__[];
extern "C"
PyObject * pypylith3d_read_diff(PyObject *, PyObject *);

// read time steps at which full output is desired
extern char pypylith3d_read_fuldat__name__[];
extern char pypylith3d_read_fuldat__doc__[];
extern "C"
PyObject * pypylith3d_read_fuldat(PyObject *, PyObject *);

// read time history info
extern char pypylith3d_read_hist__name__[];
extern char pypylith3d_read_hist__doc__[];
extern "C"
PyObject * pypylith3d_read_hist(PyObject *, PyObject *);

// read material history info
// these are not being used for now
// extern char pypylith3d_read_mathist__name__[];
// extern char pypylith3d_read_mathist__doc__[];
// extern "C"
// PyObject * pypylith3d_read_mathist(PyObject *, PyObject *);

// read prestresses
// extern char pypylith3d_read_prestr__name__[];
// extern char pypylith3d_read_prestr__doc__[];
// extern "C"
// PyObject * pypylith3d_read_prestr(PyObject *, PyObject *);

// read local coordinate rotations
extern char pypylith3d_read_skew__name__[];
extern char pypylith3d_read_skew__doc__[];
extern "C"
PyObject * pypylith3d_read_skew(PyObject *, PyObject *);

// read slippery node definitions
extern char pypylith3d_read_slip__name__[];
extern char pypylith3d_read_slip__doc__[];
extern "C"
PyObject * pypylith3d_read_slip(PyObject *, PyObject *);

// read split node definitions
extern char pypylith3d_read_split__name__[];
extern char pypylith3d_read_split__doc__[];
extern "C"
PyObject * pypylith3d_read_split(PyObject *, PyObject *);

// read state output information
extern char pypylith3d_read_stateout__name__[];
extern char pypylith3d_read_stateout__doc__[];
extern "C"
PyObject * pypylith3d_read_stateout(PyObject *, PyObject *);

// read time step group informations
extern char pypylith3d_read_timdat__name__[];
extern char pypylith3d_read_timdat__doc__[];
extern "C"
PyObject * pypylith3d_read_timdat(PyObject *, PyObject *);

// read traction boundary conditions
extern char pypylith3d_read_tractions__name__[];
extern char pypylith3d_read_tractions__doc__[];
extern "C"
PyObject * pypylith3d_read_tractions(PyObject *, PyObject *);

// read winkler force information
extern char pypylith3d_read_wink__name__[];
extern char pypylith3d_read_wink__doc__[];
extern "C"
PyObject * pypylith3d_read_wink(PyObject *, PyObject *);

#endif

// version
// $Id: parser.h,v 1.3 2005/03/31 23:27:58 willic3 Exp $

// End of file
