// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ======================================================================
//

#include "TestArray.hh" // implementation of class methods

#include "array.hh" // USES double_array

#include <iostream> // USES std::cerr
#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Check to make sure array of values match expected values.
bool
pylith::utils::TestArray::check(const double* valuesE,
				const int nvalues,
				const double_array& values)
{ // check(double)
  assert( (0 == nvalues && 0 == valuesE) ||
	  (0 < nvalues && 0 != valuesE) );

  if (nvalues != values.size()) {
    std::cerr << "Array size mismatch, expected: " << nvalues
	      << " actual: " << values.size() << std::endl;
    return false;
  } // if

  const double tolerance = 1.0e-06;
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
} // check(double)


// End of file 
