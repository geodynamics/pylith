// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "EffectiveStress.hh" // implementation of object methods

#include "petsc.h" // USES PetscLogFlops

#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
// Compute effective stress, given an initial guess, a vector of parameters,
// and functions to compute the effective stress function and it's derivative.
// The first entry in effStressParams should be a stress scaling factor
// (e.g., mu) to provide a reasonable initial guess in the case where the
// actual initial guess is zero.
double
pylith::materials::EffectiveStress::getEffStress(
					const double effStressInitialGuess,
					const double* effStressParams,
					***effStressFunc,
					***effStressFuncDFunc)
{ // getEffStress
  // Check parameters
  if (effStressInitialGuess < 0.0)
    throw std::runtime_error("Effective stress initial guess must be >= 0.");

  // Bracket the root.
  double x1 = 0.0;
  double x2 = 0.0;
  if (effStressInitialGuess > 0.0) {
    x1 = effStressInitialGuess - 0.5 * effStressInitialGuess;
    x2 = effStressInitialGuess + 0.5 * effStressInitialGuess;
  } else {
    x1 = effStressParams[0] - 0.5 * effStressParams[0];
    x2 = effStressParams[0] + 0.5 * effStressParams[0];
  } // else

  PetscLogFlops(4);
  _bracketEffStress(x1, x2, effStressParams, effStressFunc);

  // Find effective stress using Newton's method with bisection.
  const double effStress = _findEffStress(x1,
					  x2,
					  effStressParams,
					  effStressFunc,
					  effStressFuncDFunc);

  return effStress;

} // getEffStress

// ----------------------------------------------------------------------
// Bracket effective stress.
void
pylith::materials::EffectiveStress::_bracketEffStress(
					double* const x1,
					double* const x2,
					const double* effStressParams,
					??? effStressFunc)
{ // _bracketEffStress
  // Arbitrary number of iterations to bracket the root
  const int maxIterations = 50;

  // Arbitrary factor by which to increase the brackets.
  const double bracketFactor = 1.6;

  double funcValue1 = effStressFunc(x1, effStressParams);
  double funcValue2 = effStressFunc(x2, effStressParams);

  int iteration = 0;
  bool bracketed = false;

  while ((iteration < maxIterations) && (bracketed == false)) {
    if ((funcValue1 * funcValue2) < 0.0) {
      bracketed = true;
    } else {
      if (std::abs(f1) < std::abs(f2)) {
	x1 += bracketFactor * (x1 - x2);
	x1 = std::max(x1, 0.0);
	funcValue1 = effStressFunc(x1, effStressParams);
      } else {
	x2 += bracketFactor * (x1 - x2);
	x2 = std::max(x2, 0.0);
	funcValue2 = effStressFunc(x2, effStressParams);
      } // else
    } // else
    ++iteration;
  } // while

  PetscLogFlops(5 * iteration);
  if (bracketed == false)
    throw std::runtime_error("Unable to bracket effective stress.");

} // _bracketEffStress

// ----------------------------------------------------------------------
// Find root using Newton's method with bisection.
void
pylith::materials::EffectiveStress::_findEffStress(
					double* const x1,
					double* const x2,
					const double* effStressParams,
					??? effStressFunc,
					??? effStressFuncDFunc)
{ // _findEffStress
  // Arbitrary number of iterations to find the root
  const int maxIterations = 100;

  // Desired accuracy for root. This is a bit arbitrary for now.
  const double accuracy = 1.0e-10;

  /// Determine if root has already been found, or if root is not bracketed.
  // Otherwise, organize search so that effStressFunc(xLow) is less than zero.
  double funcValueLow = effStressFunc(x1, effStressParams);
  double funcValueHigh = effStressFunc(x2, effStressParams);
  if (funcValueLow * funcValueHigh > 0.0)
    throw std::runtime_error("Effective stress is not bracketed.");

  double effStress = 0.0;
  double xLow = 0.0;
  double xHigh = 0.0;

  if (std::abs(funcValueLow) < accuracy) {
    effStress = x1;
    return effStress;
  } else if (std::abs(funcValueHigh) < accuracy) {
    effStress = x2;
    return effStress;
  } else if (funcValueLow < 0.0) {
    xLow = x1;
    xHigh = x2;
  } else {
    xHigh = x1;
    xLow = x2;
  }

  effStress = 0.5 * (x1 + x2);
  double dxPrevious = std::abs(x2 - x1);
  double dx = dxPrevious;
  double funcValue = 0.0;
  double funcDeriv = 0.0;
  effStressFuncDFunc(effStress, effStressParams, funcValue, funcDeriv);
  int iteration = 0;
  bool converged = false;

  // *******  finish fixing through here ***********
  while ((iteration < maxIterations) && (converged == false)) {
    if ((funcValue1 * funcValue2) < 0.0) {
      bracketed = true;
    } else {
      if (std::abs(f1) < std::abs(f2)) {
	x1 += bracketFactor * (x1 - x2);
	x1 = std::max(x1, 0.0);
	funcValue1 = effStressFunc(x1, effStressParams);
      } else {
	x2 += bracketFactor * (x1 - x2);
	x2 = std::max(x2, 0.0);
	funcValue2 = effStressFunc(x2, effStressParams);
      } // else
    } // else
  } // while

  if (bracketed == false)
    throw std::runtime_error("Unable to bracket effective stress.");

} // _bracketEffStress


// End of file 
