// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "ViscoelasticMaxwell.hh" // implementation of object methods

#include "petsc.h" // USES PetscLogFlops

#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
// Compute viscous strain parameter for a linear Maxwell model.
double
pylith::materials::ViscoelasticMaxwell::viscousStrainParam(const double dt,
							   const double maxwellTime)
{ // viscousStrainParam
  // Check parameters
  if (maxwellTime <= 0.0)
    throw std::runtime_error("Maxwell time must be greater than 0.");

  // Define cutoff values
  const double timeFrac = 1.0e-10;

  // Compute viscous strain parameter.  The ratio of dt and
  // maxwellTime should never approach timeFrac for any reasonable
  // computation, but I have put in alternative solutions just in
  // case.

  double dq = 0.0;

  // Use series expansion if dt is very small, since default solution
  // blows up otherwise.

  if (dt < timeFrac*maxwellTime) {
    double fSign = 1.0;
    double factorial = 1.0;
    double fraction = 1.0;
    dq = 1.0;

    const int numTerms = 5;
    for (int iTerm=2; iTerm <= numTerms; ++iTerm) {
      factorial *= iTerm;
      fSign *= -1.0;
      fraction *= dt / maxwellTime;
      dq += fSign * fraction / factorial;
    } // for
    PetscLogFlops(8*(numTerms-1));
  } else if (maxwellTime < timeFrac*dt) {
    // Throw away exponential term if maxwellTime is very small.
    dq = maxwellTime / dt;
    PetscLogFlops(1);
  } else{
    // Default solution.
    dq = maxwellTime*(1.0-exp(-dt/maxwellTime))/dt;
    PetscLogFlops(6);
  } // else

  return dq;
} // viscousStrainParam
  

// End of file 
