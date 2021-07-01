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
 * @file libsrc/utils/macrodefs.h
 *
 * @brief Macro definitions for PyLith.
 */

#if !defined(pylith_utils_macrodefs_h)
#define pylith_utils_macrodefs_h

#if !defined(CALL_MEMBER_FN)
#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember))
#endif

#endif // pylith_utils_macro_defs_h


// End of file
