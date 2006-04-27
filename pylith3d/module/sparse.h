// -*- C++ -*-
// 
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  PyLith by Charles A. Williams
//  Copyright (c) 2003-2006 Rensselaer Polytechnic Institute
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

#if !defined(pypylith3d_sparse_h)
#define pypylith3d_sparse_h

// compute maximum number of nonzero entries in stiffness matrix
extern char pypylith3d_cmp_stiffsz__name__[];
extern char pypylith3d_cmp_stiffsz__doc__[];
extern "C"
PyObject * pypylith3d_cmp_stiffsz(PyObject *, PyObject *);

// create linked list for sparse matrix
extern char pypylith3d_lnklst__name__[];
extern char pypylith3d_lnklst__doc__[];
extern "C"
PyObject * pypylith3d_lnklst(PyObject *, PyObject *);

// localize id array for reference by element
extern char pypylith3d_local__name__[];
extern char pypylith3d_local__doc__[];
extern "C"
PyObject * pypylith3d_local(PyObject *, PyObject *);

// localize nfault array for reference by element
extern char pypylith3d_localf__name__[];
extern char pypylith3d_localf__doc__[];
extern "C"
PyObject * pypylith3d_localf(PyObject *, PyObject *);

// localize idx array for reference by element
extern char pypylith3d_localx__name__[];
extern char pypylith3d_localx__doc__[];
extern "C"
PyObject * pypylith3d_localx(PyObject *, PyObject *);

// create sparse matrix in modified sparse row format
extern char pypylith3d_makemsr__name__[];
extern char pypylith3d_makemsr__doc__[];
extern "C"
PyObject * pypylith3d_makemsr(PyObject *, PyObject *);

#endif

// version
// $Id: sparse.h,v 1.5 2005/04/21 23:21:00 willic3 Exp $

// End of file
