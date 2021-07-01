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
