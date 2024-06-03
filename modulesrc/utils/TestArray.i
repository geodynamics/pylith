// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information. 
// =================================================================================================

/**
 * @file modulesrc/utils/TestArray.i
 *
 * @brief Python interface to C++ TestArray object.
 */

%apply(PylithScalar* IN_ARRAY1, int DIM1) {
  (const PylithScalar* valuesE,
   const int nvalues)
    };
%inline %{
  /** Check to make sure array of values match expected values.
   *
   * @param valuesE Array of expected values.
   * @param nvalues Array size.
   * @param values Array of values to check.
   */
  bool
  TestArray_checkScalar(const PylithScalar* valuesE,
			const int nvalues,
			const pylith::scalar_array& values) {
    return pylith::utils::TestArray::check(valuesE, nvalues, values);
  } // check(PylithScalar)
%} // inline
%clear(const PylithScalar* valuesE, const int nvalues);


// End of file 
