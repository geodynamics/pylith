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

#include <petscksp.h> // USES PetscKSP

// ----------------------------------------------------------------------
// Constructor
pylith::problems::Solver::Solver(void) :
  _ksp(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::problems::Solver::~Solver(void)
{ // destructor
  if (0 != _ksp) {
    PetscErrorCode err = KSPDestroy(_ksp); _ksp = 0;
    CHECK_PETSC_ERROR(err);
  } // if
} // destructor

// ----------------------------------------------------------------------
// Initialize solver.
void
pylith::problems::Solver::initialize(topology::SolutionFields* fields)
{ // initialize
  assert(0 != fields);

  fields->createScatter();

  topology::Field<topology::Mesh>& solution = fields->solution();
  solution.createVector();

  topology::Field<topology::Mesh>& residual = fields->get("residual");
  residual.createVector();

  PetscErrorCode err = 0;

  if (0 != _ksp) {
    err = KSPDestroy(_ksp); _ksp = 0;
    CHECK_PETSC_ERROR(err);
  } // if    
  err = KSPCreate(fields->mesh().comm(), &_ksp); CHECK_PETSC_ERROR(err);

  err = KSPSetFromOptions(_ksp); CHECK_PETSC_ERROR(err);
} // initialize


// End of file
