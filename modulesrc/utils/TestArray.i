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
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file modulesrc/utils/TestArray.i
 *
 * @brief Python interface to C++ TestArray object.
 */

%apply(double* IN_ARRAY1, int DIM1) {
  (const double* valuesE,
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
  TestArray_checkDouble(const double* valuesE,
			const int nvalues,
			const pylith::double_array& values) {
    pylith::utils::TestArray::check(valuesE, nvalues, values);
  } // check(double)
%} // inline
%clear(const double* valuesE, const int nvalues);


// End of file 
