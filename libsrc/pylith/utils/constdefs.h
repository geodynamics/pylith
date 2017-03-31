// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/utils/constdefs.h
 *
 * @brief Macro definitions for PyLith.
 */

#if !defined(pylith_utils_constdefs_h)
#define pylith_utils_constdefs_h

#include "types.hh" // HASA PylithScalar

#include "limits.h" // USES INT_MAX, INT_MIN, LONG_MAX, LONG_MIN

namespace pylith {
    static const double PYLITH_MAXDOUBLE = 1.0e+99;
    static const float PYLITH_MAXFLOAT = 1.0e+30;
    static const PylithInt PYLITH_MAXINT = sizeof(PylithInt) == sizeof(long) ? LONG_MAX : INT_MAX;
    static const PylithInt PYLITH_MININT = sizeof(PylithInt) == sizeof(long) ? LONG_MIN : INT_MIN;
    static const PylithScalar PYLITH_MAXSCALAR = (sizeof(PylithScalar) == sizeof(double)) ? PYLITH_MAXDOUBLE : PYLITH_MAXFLOAT;
}


#endif // pylith_utils_constdefs_h


// End of file
