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

#include <cmath> // USES fabs()
#include <cassert> // USES assert()
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
				 const double stressScale,
				 const EffStressStruct& effStressParams,
				 effStressFunc_fn_type effStressFunc,
				 effStressFuncDFunc_fn_type effStressFuncDFunc)
{ // getEffStress
  // Check parameters
  assert(effStressInitialGuess >= 0.0);

  // Bracket the root.
  double x1 = 0.0;
  double x2 = 0.0;
  if (effStressInitialGuess > 0.0) {
    x1 = effStressInitialGuess - 0.5 * effStressInitialGuess;
    x2 = effStressInitialGuess + 0.5 * effStressInitialGuess;
  } else {
    x1 = stressScale - 0.5 * stressScale;
    x2 = stressScale + 0.5 * stressScale;
  } // else

  PetscLogFlops(4);
  _bracketEffStress(&x1, &x2, effStressParams, effStressFunc);

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
				     double* px1,
				     double* px2,
				     const EffStressStruct& effStressParams,
				     effStressFunc_fn_type effStressFunc)
{ // _bracketEffStress
  // Arbitrary number of iterations to bracket the root
  const int maxIterations = 50;

  // Arbitrary factor by which to increase the brackets.
  const double bracketFactor = 1.6;
  double x1 = *px1;
  double x2 = *px2;

  double funcValue1 = effStressFunc(x1, effStressParams);
  double funcValue2 = effStressFunc(x2, effStressParams);

  int iteration = 0;
  bool bracketed = false;
  while (iteration < maxIterations) {
    if ((funcValue1 * funcValue2) < 0.0) {
      bracketed = true;
      break;
    } // if

    if (fabs(funcValue1) < fabs(funcValue2)) {
      x1 += bracketFactor * (x1 - x2);
      x1 = std::max(x1, 0.0);
      funcValue1 = effStressFunc(x1, effStressParams);
    } else {
      x2 += bracketFactor * (x1 - x2);
      x2 = std::max(x2, 0.0);
      funcValue2 = effStressFunc(x2, effStressParams);
    } // else
    ++iteration;
  } // while

  *px1 = x1;
  *px2 = x2;

  PetscLogFlops(5 * iteration);
  if (!bracketed)
    throw std::runtime_error("Unable to bracket effective stress.");
} // _bracketEffStress

// ----------------------------------------------------------------------
// Find root using Newton's method with bisection.
double
pylith::materials::EffectiveStress::_findEffStress(
				 const double x1,
				 const double x2,
				 const EffStressStruct& effStressParams,
				 effStressFunc_fn_type effStressFunc,
				 effStressFuncDFunc_fn_type effStressFuncDFunc)
{ // _findEffStress
  // Arbitrary number of iterations to find the root
  const int maxIterations = 100;

  // Desired accuracy for root. This is a bit arbitrary for now.
  const double accuracy = 1.0e-6;

  /// Determine if root has already been found, or if root is not bracketed.
  // Otherwise, organize search so that effStressFunc(xLow) is less than zero.
  double funcValueLow = effStressFunc(x1, effStressParams);
  double funcValueHigh = effStressFunc(x2, effStressParams);
  assert(funcValueLow * funcValueHigh <= 0.0);

  double effStress = 0.0;
  double xLow = 0.0;
  double xHigh = 0.0;
  bool converged = false;

  if (fabs(funcValueLow) < accuracy) {
    effStress = x1;
    converged = true;
    return effStress;
  } else if (fabs(funcValueHigh) < accuracy) {
    effStress = x2;
    converged = true;
    return effStress;
  } else if (funcValueLow < 0.0) {
    xLow = x1;
    xHigh = x2;
  } else {
    xHigh = x1;
    xLow = x2;
  }

  effStress = 0.5 * (x1 + x2);
  double dxPrevious = fabs(x2 - x1);
  double dx = dxPrevious;
  double funcValue = 0.0;
  double funcDeriv = 0.0;
  double funcXHigh = 0.0;
  double funcXLow = 0.0;
  effStressFuncDFunc(effStress, effStressParams, &funcValue, &funcDeriv);
  int iteration = 0;

  while (iteration < maxIterations) {
    funcXHigh = (effStress - xHigh) * funcDeriv - funcValue;
    funcXLow = (effStress - xLow) * funcDeriv - funcValue;
    // Use bisection if solution goes out of bounds or is not converging
    // fast enough.
    if ( (funcXHigh * funcXLow >= 0.0) ||
	 (fabs(2.0 * funcValue) > fabs(dxPrevious * funcDeriv))) {
      dxPrevious = dx;
      dx = 0.5 * (xHigh - xLow);
      effStress = xLow + dx;
    } else {
      dxPrevious = dx;
      dx = funcValue/funcDeriv;
      effStress = effStress - dx;
    } // else
    if (fabs(dx) < accuracy) {
      converged = true;
      break;
    } // if
    effStressFuncDFunc(effStress, effStressParams, &funcValue, &funcDeriv);
    if (funcValue < 0.0) {
      xLow = effStress;
    } else {
      xHigh = effStress;
    } // else
    ++iteration;
  } // while

  if (converged == false)
    throw std::runtime_error("Cannot find root of effective stress function.");

  PetscLogFlops(5 + 15 * iteration);
  return effStress;

  return effStress;
} // _findEffStress


// End of file 
