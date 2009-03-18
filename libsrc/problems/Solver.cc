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

#include "pylith/topology/SolutionFields.hh" // USES SolutionFields

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
pylith::problems::Solver::initialize(topology::SolutionFields* fields)
{ // initialize
  assert(0 != fields);

  topology::Field<topology::Mesh>& solution = fields->solution();
  solution.createVector(); // Move this to use of SolutionFields and copyLayout?
  solution.createScatter();

  topology::Field<topology::Mesh>& residual = fields->get("residual");
  residual.createVector();
  solution.createScatter(); // :TODO: eliminate duplication
} // initialize


// End of file
