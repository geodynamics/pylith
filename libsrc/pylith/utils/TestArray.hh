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
 * @file libsrc/utils/TestArray.hh
 *
 * @brief C++ object for testing array values.
 */

#if !defined(pylith_utils_testarray_hh)
#define pylith_utils_testarray_hh

// Include directives ---------------------------------------------------
#include "utilsfwd.hh" // forward declarations

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
class pylith::utils::TestArray
{ // TestArray

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /** Check to make sure array of values match expected values.
   *
   * @param valuesE Array of expected values.
   * @param nvalues Array size.
   * @param values Array of values to check.
   */
  static
  bool
  check(const PylithScalar* valuesE,
	const int nvalues,
	const scalar_array& values);

}; // EventLogger

#endif // pylith_utils_testarray_hh


// End of file 
