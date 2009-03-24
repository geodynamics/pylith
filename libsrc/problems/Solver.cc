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

// ----------------------------------------------------------------------
// Constructor
pylith::problems::Solver::Solver(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::problems::Solver::~Solver(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Initialize solver.
void
pylith::problems::Solver::initialize(topology::SolutionFields* fields,
				     Formulation* const formulation)
{ // initialize
  _formulation = formulation;
} // initialize


// End of file
