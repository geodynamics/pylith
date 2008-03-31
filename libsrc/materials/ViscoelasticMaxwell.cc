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

#include "ViscoelasticMaxwell.hh" // implementation of object methods

#include "petsc.h" // USES PetscLogFlops

#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
// Compute viscous strain parameter for a linear Maxwell model.
double
pylith::materials::ViscoelasticMaxwell::computeVisStrain(const double dt,
							 const double maxwelltime)
{ // check parameters and define cutoff values
  if (maxwelltime <= 0.0)
    throw std::runtime_error("Maxwell time must be > 0.");
  const double timeFrac = 1.0e-10;
  const int numTerms = 5;

  // Compute viscous strain parameter.
  // The ratio of dt and maxwelltime should never approach timeFrac for any
  // reasonable computation, but I have put in alternative solutions just in
  // case.
  double dq = 0.0;
  // Use series expansion if dt is very small, since default solution blows
  // up otherwise.
  if(dt < timeFrac*maxwelltime) {
    double fSign = 1.0;
    double factorial = 1.0;
    double fraction = 1.0;
    dq = 1.0;
    for (int iTerm=2; iTerm <= numTerms; ++iTerm) {
      factorial *= iTerm;
      fSign *= -1.0;
      fraction *= dt/maxwelltime;
      dq += fSign*fraction/factorial;
    } // for
    PetscLogFlops(8*(numTerms-1));
  // Throw away exponential term if maxwelltime is very small.
  } else if (maxwelltime < timeFrac*dt) {
    dq = maxwelltime/dt;
    PetscLogFlops(1);
  // Default solution.
  } else{
    dq = maxwelltime*(1.0-exp(-dt/maxwelltime))/dt;
    PetscLogFlops(6);
  } // else

  return dq;
} // computeVisStrain
  

// End of file 
