// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "Discretization.hh" // implementation of class methods

#include "pylith/utils/error.hh" // USES PYLITH_CHECK_ERROR

#include "petscdm.h" // USES PetscDM
#include "petscds.h" // USES PetscDS


// ----------------------------------------------------------------------
// Default constructor.
pylith::topology::Discretization::Discretization() :
  _basisOrder(1),
  _quadOrder(-1), // use basis order
  _basisContinuity(true)
{ // constructor
} // constructor


// ----------------------------------------------------------------------
// Destructor.
pylith::topology::Discretization::~Discretization(void)
{ // destructor
} // destructor


// ----------------------------------------------------------------------
// Set order for basis functions.
void
pylith::topology::Discretization::basisOrder(const int value)
{ // basisOrder
  PYLITH_METHOD_BEGIN;

  if (value < 0) {
    std::ostringstream msg;
    msg << "Basis order (" << value << ") must be positive.";
    throw std::runtime_error(msg.str());
  } // if

  _basisOrder = value;

  PYLITH_METHOD_END;
} // basisOrder


// ----------------------------------------------------------------------
// Get order of basis functions.
int
pylith::topology::Discretization::basisOrder(void) const
{ // basisOrder
  return _basisOrder;
} // basisOrder


// ----------------------------------------------------------------------
// Set basis continuity flag.
void
pylith::topology::Discretization::isBasisContinuous(const bool value)
{ // isBasisContinuous
  _basisContinuity = true;
} // isBasisContinuous


// ----------------------------------------------------------------------
// Get basis continuity flag.
int
pylith::topology::Discretization::isBasisContinuous(void) const
{ // isBasisContinuous
  return _basisContinuity;
} // isBasisContinuous


// ----------------------------------------------------------------------
// Set order of quadrature scheme.
void
pylith::topology::Discretization::quadratureOrder(const int value)
{ // quadratureOrder
  _quadOrder = value;
} // quadratureOrder


// ----------------------------------------------------------------------
// Get order of quadrature scheme.
int
pylith::topology::Discretization::quadratureOrder(void) const
{ // quadratureOrder
  return _quadOrder;
} // quadratureOrder


// ----------------------------------------------------------------------
// Create PetscFE object for discretization.
PetscFE
pylith::topology::Discretization::createFE(const PetscDM dm,
					    const int numComponents)
{ // createFE
  PYLITH_METHOD_BEGIN;


  PetscErrorCode err;

  // Get spatial dimension of mesh.
  int dim = 0;
  err = DMGetDimension(dm, &dim);PYLITH_CHECK_ERROR(err);

  // Determine type of cell in mesh.
  PetscBool isSimplex = PETSC_FALSE;
  PetscInt closureSize, vStart, vEnd;
  PetscInt* closure = NULL;
  err = DMPlexGetDepthStratum(dm, 0, &vStart, &vEnd);PYLITH_CHECK_ERROR(err);
  err = DMPlexGetTransitiveClosure(dm, 0, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
  PetscInt numVertices = 0;
  for (PetscInt c = 0; c < closureSize*2; c+=2) {
    if ((closure[c] >= vStart) && (closure[c] < vEnd)) {
      ++numVertices;
    } // if
  } // for
  if (numVertices == dim+1) {
    isSimplex = PETSC_TRUE;
  } // if
  err = DMPlexRestoreTransitiveClosure(dm, 0, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);

  // Create space 
  PetscSpace space = NULL;
  err = PetscSpaceCreate(PetscObjectComm((PetscObject) dm), &space);PYLITH_CHECK_ERROR(err);
  err = PetscSpacePolynomialSetNumVariables(space, dim);PYLITH_CHECK_ERROR(err);
  err = PetscSpaceSetOrder(space, _basisOrder);
  err = PetscSpaceSetUp(space);PYLITH_CHECK_ERROR(err);

  // Create dual space
  PetscDualSpace dualspace = NULL;
  PetscDM dmCell = NULL;
  err = PetscDualSpaceCreate(PetscObjectComm((PetscObject) dm), &dualspace);PYLITH_CHECK_ERROR(err);
  err = PetscDualSpaceCreateReferenceCell(dualspace, dim, isSimplex, &dmCell);PYLITH_CHECK_ERROR(err);
  err = PetscDualSpaceSetDM(dualspace, dmCell);PYLITH_CHECK_ERROR(err);
  err = DMDestroy(&dmCell);PYLITH_CHECK_ERROR(err);
  err = PetscDualSpaceSetOrder(dualspace, _basisOrder);PYLITH_CHECK_ERROR(err);
  err = PetscDualSpaceSetFromOptions(dualspace);PYLITH_CHECK_ERROR(err);
  err = PetscDualSpaceSetUp(dualspace);PYLITH_CHECK_ERROR(err);

  // Create element
  PetscFE fe = NULL;
  err = PetscFECreate(PetscObjectComm((PetscObject) dm), &fe);PYLITH_CHECK_ERROR(err);
  err = PetscFESetBasisSpace(fe, space);PYLITH_CHECK_ERROR(err);
  err = PetscFESetDualSpace(fe, dualspace);PYLITH_CHECK_ERROR(err);
  err = PetscFESetNumComponents(fe, numComponents);PYLITH_CHECK_ERROR(err);
  err = PetscFESetUp(fe);PYLITH_CHECK_ERROR(err);
  err = PetscSpaceDestroy(&space);PYLITH_CHECK_ERROR(err);
  err = PetscDualSpaceDestroy(&dualspace);PYLITH_CHECK_ERROR(err);

  // Create quadrature
  PetscQuadrature quadrature = NULL;
  const int qorder = PetscMax(_quadOrder > 0 ? _quadOrder : _basisOrder, 1);
  if (isSimplex) {
    err = PetscDTGaussJacobiQuadrature(dim, qorder, -1.0, 1.0, &quadrature);PYLITH_CHECK_ERROR(err);
  } else {
    err = PetscDTGaussTensorQuadrature(dim, qorder, -1.0, 1.0, &quadrature);PYLITH_CHECK_ERROR(err);
  }
  err = PetscFESetQuadrature(fe, quadrature);PYLITH_CHECK_ERROR(err);
  err = PetscQuadratureDestroy(&quadrature);PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_RETURN(fe);
} // createFE


// End of file 
