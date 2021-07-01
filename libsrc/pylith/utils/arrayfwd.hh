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
 * @file libsrc/utils/arrayfwd.hh
 *
 * @brief Forward declarations for PyLith array objects.
 *
 * These are generally just typenames for C++ STL objects.
 *
 * For simple types (i.e., int and PylithScalar) std::valarray provides some
 * features that std::vector does not have, such as operating on the
 * whole array at once.
 */

#if !defined(pylith_utils_arrayfwd_hh)
#define pylith_utils_arrayfwd_hh

#include "types.hh" // USES PylithScalar, PylithReal

#include <string> // USES std::string
#include <vector> // USES std::vector
#include <valarray> // USES std::valarray

/// Aliases
namespace pylith {
    /// Alias for std::vector<int>
    typedef std::vector<int, std::allocator<int> > int_vector;

    /// Alias for std::vector<double>
    typedef std::vector<double, std::allocator<double> > double_vector;

    /// Alias for std::vector<std::string>
    typedef std::vector<std::string, std::allocator<std::string> > string_vector;

    /// Alias for std::valarray<float>
    typedef std::valarray<float> float_array;

    /// Alias for std::valarray<double>
    typedef std::valarray<double> double_array;

    /// Alias for std::valarray<PylithInt>
    typedef std::valarray<PylithInt> int_array;

    /// Alias for std::valarray<char>
    typedef std::valarray<char> char_array;

    /// Alias for std::valarray<PylithReal>
    typedef std::valarray<PylithReal> real_array;

    /// Alias for std::valarray<PylithScalar>
    typedef std::valarray<PylithScalar> scalar_array;

} // pylith

#endif // pylith_utils_arrayfwd_hh

// End of file
