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

#include <portinfo>

#include "Solver.hh" // implementation of class methods

#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Constructor
pylith::problems::Solver::Solver(void) :
  _formulation(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::problems::Solver::~Solver(void)
{ // destructor
  _formulation = 0; // Handle only, do not manage memory.
} // destructor

// ----------------------------------------------------------------------
// Initialize solver.
void
pylith::problems::Solver::initialize(const topology::SolutionFields& fields,
				     const topology::Jacobian& jacobian,
				     Formulation* const formulation)
{ // initialize
  assert(0 != formulation);

  _formulation = formulation;
} // initialize


// End of file
