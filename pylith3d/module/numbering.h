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

#if !defined(pypylith3d_numbering_h)
#define pypylith3d_numbering_h

// assign equation numbers for Winkler BC
extern char pypylith3d_assign_wink__name__[];
extern char pypylith3d_assign_wink__doc__[];
extern "C"
PyObject * pypylith3d_assign_wink(PyObject *, PyObject *);

// create id array
extern char pypylith3d_create_id__name__[];
extern char pypylith3d_create_id__doc__[];
extern "C"
PyObject * pypylith3d_create_id(PyObject *, PyObject *);

// form id array for split nodes
extern char pypylith3d_id_split__name__[];
extern char pypylith3d_id_split__doc__[];
extern "C"
PyObject * pypylith3d_id_split(PyObject *, PyObject *);

// find closest fault neignbors for slippery nodes
extern char pypylith3d_nfind__name__[];
extern char pypylith3d_nfind__doc__[];
extern "C"
PyObject * pypylith3d_nfind(PyObject *, PyObject *);

#endif

// version
// $Id: numbering.h,v 1.3 2005/04/21 23:17:42 willic3 Exp $

// End of file
