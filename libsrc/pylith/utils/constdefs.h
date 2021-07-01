// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
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

namespace pylith {
    static const double PYLITH_MAXDOUBLE = 1.0e+99;
    static const float PYLITH_MAXFLOAT = 1.0e+30;
    static const PylithInt PYLITH_MAXINT = PETSC_MAX_INT;
    static const PylithInt PYLITH_MININT = PETSC_MIN_INT;;
    static const PylithScalar PYLITH_MAXSCALAR = (sizeof(PylithScalar) == sizeof(double)) ? PYLITH_MAXDOUBLE : PYLITH_MAXFLOAT;
}


#endif // pylith_utils_constdefs_h


// End of file
