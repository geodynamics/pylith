// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "pylith/utils/utilsfwd.hh" // forward declarations

#include "array.hh" // USES scalar_array

// TestArray ------------------------------------------------------------
/** @brief C++ object for testing array values.
 *
 * This object is used in unit testing of SWIG interfaces where the
 * C++ object has an accessor returning a std::valarray. The TestArray
 * methods provide the ability to compare the array returned by the
 * accessor against the expected values, which are supplied via a
 * pointer and a size (number of values).
 */
class pylith::utils::TestArray { // TestArray
    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    /** Check to make sure array of values match expected values.
     *
     * @param valuesE Array of expected values.
     * @param nvalues Array size.
     * @param values Array of values to check.
     */
    static
    bool check(const PylithScalar* valuesE,
               const int nvalues,
               const scalar_array& values);

}; // EventLogger

// End of file
