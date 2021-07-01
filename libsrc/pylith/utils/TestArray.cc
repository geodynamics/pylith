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

#include <portinfo>

#include "TestArray.hh" // implementation of class methods

#include "array.hh" // USES scalar_array

#include <iostream> // USES std::cerr
#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Check to make sure array of values match expected values.
bool
pylith::utils::TestArray::check(const PylithScalar* valuesE,
				const int nvalues,
				const scalar_array& values)
{ // check(PylithScalar)
  assert( (0 == nvalues && 0 == valuesE) ||
	  (0 < nvalues && 0 != valuesE) );

  if (size_t(nvalues) != values.size()) {
    std::cerr << "Array size mismatch, expected: " << nvalues
	      << " actual: " << values.size() << std::endl;
    return false;
  } // if

  const PylithScalar tolerance = 1.0e-06;
  bool okay = true;
  for (int i=0; i < nvalues; ++i) {
    okay = true;
    if (0.0 != valuesE[i]) {
      if (fabs(1.0 - values[i]/valuesE[i]) > tolerance)
	okay = false;
    } else if (fabs(values[i] - valuesE[i]) > tolerance)
      okay = false;

    if (!okay) {
      std::cerr << "Mismatch in array at index " << i << ", expected: "
		<< valuesE[i] << ", actual: " << values[i] << std::endl;
      return false;
    } // if
  } // for

  return true;
} // check(PylithScalar)


// End of file 
