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

#if !defined(pylithomop3d_sorting_h)
#define pylithomop3d_sorting_h

// sort elements into element families
extern char pylithomop3d_sort_elements__name__[];
extern char pylithomop3d_sort_elements__doc__[];
extern "C"
PyObject * pylithomop3d_sort_elements(PyObject *, PyObject *);

// Sort slippery nodes according to reordered elements
extern char pylithomop3d_sort_slip__name__[];
extern char pylithomop3d_sort_slip__doc__[];
extern "C"
PyObject * pylithomop3d_sort_slip(PyObject *, PyObject *);

// Sort split nodes according to reordered elements
extern char pylithomop3d_sort_split__name__[];
extern char pylithomop3d_sort_split__doc__[];
extern "C"
PyObject * pylithomop3d_sort_split(PyObject *, PyObject *);

#endif

// version
// $Id: sorting.h,v 1.1 2005/04/21 01:05:24 willic3 Exp $

// End of file
