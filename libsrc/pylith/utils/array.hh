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
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/utils/array.hh
 *
 * @brief Header file for PyLith array objects.
 *
 * Since the arrays are really C++ STL objects, we simply include the
 * STL header files.
 */

#if !defined(pylith_utils_array_hh)
#define pylith_utils_array_hh

#include "types.hh"
#include "arrayfwd.hh"

#include "sievetypes.hh" // ensure we include petscsys.h BEFORE valarray to prevent clash over isinf() and isnan().

#include <vector>
#include <valarray>

namespace pylith {
  typedef std::valarray<PylithScalar> scalar_array;
} // namespace pylith

#endif // pylith_utils_array_hh

// End of file
