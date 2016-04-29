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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "MaterialNew.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // HOLDSA Field
#include "pylith/topology/FieldQuery.hh" // HOLDSA FieldQuery
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields
#include "pylith/topology/Jacobian.hh" // USES Jacobian
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/topology/Stratum.hh" // USES StratumIS
#include "pylith/utils/array.hh" // USES scalar_array, std::vector

#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <petscds.h> // USES PetscDS

#include <strings.h> // USES strcasecmp()
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

extern "C" PetscErrorCode DMPlexComputeResidual_Internal(DM dm, PetscInt cStart, PetscInt cEnd, PetscReal time, Vec locX, Vec locX_t, Vec locF, void *user);
extern "C" PetscErrorCode DMPlexComputeJacobian_Internal(DM dm, PetscInt cStart, PetscInt cEnd, PetscReal t, PetscReal X_tShift, Vec X, Vec X_t, Mat Jac, Mat JacP,void *user);
extern "C" PetscErrorCode DMPlexComputeJacobianAction_Internal(DM dm, PetscInt cStart, PetscInt cEnd, PetscReal t, PetscReal X_tShift, Vec X, Vec X_t, Vec Y, Vec z, void *user);


// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::MaterialNew::MaterialNew(const int dimension) :
  _materialIS(NULL),
  _dimension(dimension),
  _id(0),
  _label("")
{ // constructor
  const topology::FieldBase::DiscretizeInfo defaultInfo = {-1, -1, true};
  _auxFieldsFEInfo["default"] = defaultInfo;
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::MaterialNew::~MaterialNew(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::materials::MaterialNew::deallocate(void)
{ // deallocate
  PYLITH_METHOD_BEGIN;

  IntegratorPointwise::deallocate();
  delete _materialIS; _materialIS = NULL;

  PYLITH_METHOD_END;
} // deallocate
  
// ----------------------------------------------------------------------
// Get physical property parameters and initial state (if used) from database.
void
pylith::materials::MaterialNew::initialize(const pylith::topology::Field& solution)
{ // initialize
  PYLITH_METHOD_BEGIN;

  // Get cells associated with material
  const pylith::topology::Mesh& mesh = solution.mesh();
  PetscDM dmMesh = mesh.dmMesh();assert(dmMesh);
  topology::CoordsVisitor::optimizeClosure(dmMesh);
  
  const bool includeOnlyCells = true;
  delete _materialIS; _materialIS = new topology::StratumIS(dmMesh, "material-id", _id, includeOnlyCells);assert(_materialIS);

  delete _auxFields; _auxFields = new topology::Field(mesh);assert(_auxFields);
  delete _auxFieldsQuery; _auxFieldsQuery = new topology::FieldQuery(*_auxFields);assert(_auxFieldsQuery);
  _auxFields->label("auxiliary fields");
  _auxFieldsSetup();
  _auxFields->subfieldsSetup();
  _auxFields->allocate();
  _auxFields->zeroAll();

  if (_auxFieldsDB) {
    assert(_normalizer);
    _auxFieldsQuery->openDB(_auxFieldsDB, _normalizer->lengthScale());
    _auxFieldsQuery->queryDB();
    _auxFieldsQuery->closeDB(_auxFieldsDB);
  } else { // else
    assert(0);
    throw std::logic_error("Unknown case for setting up auxiliary fields.");
  } // if/else
  _auxFields->complete();

  PYLITH_METHOD_END;
} // initialize


// ----------------------------------------------------------------------
// Compute RHS residual for G(t,s).
void
pylith::materials::MaterialNew::computeRHSResidual(pylith::topology::Field* residual,
						   const PylithReal t,
						   const PylithReal dt,
						   const pylith::topology::Field& solution)
{ // computeRHSResidual
  PYLITH_METHOD_BEGIN;

  _setFEKernelsRHSResidual(solution);
  pylith::topology::Field solutionDotNull(solution.mesh()); // :KLUDGE: fake field to satisfy general interface to _computeResidual()
  _computeResidual(residual, t, dt, solution, solutionDotNull);

  PYLITH_METHOD_END;
} // computeRHSResidual


// ----------------------------------------------------------------------
// Compute RHS Jacobian for G(t,s).
void
pylith::materials::MaterialNew::computeRHSJacobian(pylith::topology::Jacobian* jacobian,
						   pylith::topology::Jacobian* preconditioner,
						   const PylithReal t,
						   const PylithReal dt,
						   const pylith::topology::Field& solution)
{ // computeRHSJacobian
  PYLITH_METHOD_BEGIN;

  _setFEKernelsRHSJacobian(solution);
  pylith::topology::Field solutionDotNull(solution.mesh()); // :KLUDGE: fake field to satisfy general interface to _computeResidual()
  const PylithReal tshift = 0.0; // No dependence on time derivative of solution in RHS, so shift isn't applicable.
  _computeJacobian(jacobian, preconditioner, t, dt, tshift, solution, solutionDotNull);
  _needNewRHSJacobian = false;

  PYLITH_METHOD_END;
} // computeRHSJacobian


// ----------------------------------------------------------------------
// Compute LHS residual for F(t,s,\dot{s}).
void
pylith::materials::MaterialNew::computeLHSResidual(pylith::topology::Field* residual,
						   const PylithReal t,
						   const PylithReal dt,
						   const pylith::topology::Field& solution,
						   const pylith::topology::Field& solutionDot)
{ // computeLHSResidual
  PYLITH_METHOD_BEGIN;

  _setFEKernelsLHSResidual(solution);
  _computeResidual(residual, t, dt, solution, solutionDot);

  PYLITH_METHOD_END;
} // computeLHSResidual

// ----------------------------------------------------------------------
// Compute LHS Jacobian for F(t,s,\dot{s}).
void
pylith::materials::MaterialNew::computeLHSJacobianImplicit(pylith::topology::Jacobian* jacobian,
							   pylith::topology::Jacobian* preconditioner,
							   const PylithReal t,
							   const PylithReal dt,
							   const PylithReal tshift,
							   const pylith::topology::Field& solution,
							   const pylith::topology::Field& solutionDot)
{ // computeLHSJacobianImplicit
  PYLITH_METHOD_BEGIN;

  _setFEKernelsLHSJacobianImplicit(solution);
  _computeJacobian(jacobian, preconditioner, t, dt, tshift, solution, solutionDot);
  _needNewLHSJacobian = false;

  PYLITH_METHOD_END;
} // computeLHSJacobianImplicit


// ----------------------------------------------------------------------
// Compute LHS Jacobian for F(t,s,\dot{s}).
void
pylith::materials::MaterialNew::computeLHSJacobianInverseExplicit(pylith::topology::Field* jacobian,
								  const PylithReal t,
								  const PylithReal dt,
								  const PylithReal tshift,
								  const pylith::topology::Field& solution,
								  const pylith::topology::Field& solutionDot)
{ // computeLHSJacobianInverseExplicit
  PYLITH_METHOD_BEGIN;

  _setFEKernelsLHSJacobianExplicit(solution);
  _computeJacobianInverseLumped(jacobian, t, dt, tshift, solution, solutionDot);
  _needNewLHSJacobian = false;

  PYLITH_METHOD_END;
} // computeLHSJacobianInverseExplicit


// ----------------------------------------------------------------------
// Update state variables as needed.
void
pylith::materials::MaterialNew::updateStateVars(const pylith::topology::Field& solution)
{ // updateStateVars
  PYLITH_METHOD_BEGIN;

  throw std::logic_error("MaterialNew::updateStateVars() not implemented");

  PYLITH_METHOD_END;
} // updateStateVars

// ----------------------------------------------------------------------
// Compute residual using current kernels.
void
pylith::materials::MaterialNew::_computeResidual(pylith::topology::Field* residual,
						 const PylithReal t,
						 const PylithReal dt,
						 const pylith::topology::Field& solution,
						 const pylith::topology::Field& solutionDot)
{ // _computeResidual
  PYLITH_METHOD_BEGIN;

  assert(residual);
  assert(_auxFields);

  PetscDS prob = NULL;
  PetscInt cStart = 0, cEnd = 0;
  PetscErrorCode err;

  PetscDM dmMesh = solution.dmMesh();
  PetscDM dmAux = _auxFields->dmMesh();
  PetscDMLabel dmLabel;

  // Pointwise function have been set in DS
  err = DMGetDS(dmMesh, &prob);PYLITH_CHECK_ERROR(err);

  // Get auxiliary data
  err = PetscObjectCompose((PetscObject) dmMesh, "dmAux", (PetscObject) dmAux);PYLITH_CHECK_ERROR(err);
  err = PetscObjectCompose((PetscObject) dmMesh, "A", (PetscObject) _auxFields->localVector());PYLITH_CHECK_ERROR(err);

  // Compute the local residual
  err = DMGetLabel(dmMesh, "material-id", &dmLabel);PYLITH_CHECK_ERROR(err);
  err = DMLabelGetStratumBounds(dmLabel, id(), &cStart, &cEnd);PYLITH_CHECK_ERROR(err);
  err = DMPlexComputeResidual_Internal(dmMesh, cStart, cEnd, PETSC_MIN_REAL, solution.localVector(), solutionDot.localVector(), residual->localVector(), NULL);PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_END;
} // _computeResidual
  

// ----------------------------------------------------------------------
// Compute Jacobian using current kernels.
void
pylith::materials::MaterialNew::_computeJacobian(pylith::topology::Jacobian* jacobian,
						 pylith::topology::Jacobian* preconditioner,
						 const PylithReal t,
						 const PylithReal dt,
						 const PylithReal tshift,
						 const pylith::topology::Field& solution,
						 const pylith::topology::Field& solutionDot)
{ // _computeJacobian
  PYLITH_METHOD_BEGIN;

  assert(jacobian);
  assert(preconditioner);

  PetscDS prob = NULL;
  PetscInt cStart = 0, cEnd = 0;
  PetscErrorCode err;

  const PetscMat jacobianMat = jacobian->matrix();assert(jacobianMat);
  const PetscMat precondMat = preconditioner->matrix();assert(precondMat);
  PetscDM dmMesh = solution.dmMesh();
  PetscDM dmAux = _auxFields->dmMesh();
  PetscDMLabel dmLabel;

  // Pointwise function have been set in DS
  err = DMGetDS(dmMesh, &prob);PYLITH_CHECK_ERROR(err);

  // Get auxiliary data
  err = PetscObjectCompose((PetscObject) dmMesh, "dmAux", (PetscObject) dmAux);PYLITH_CHECK_ERROR(err);
  err = PetscObjectCompose((PetscObject) dmMesh, "A", (PetscObject) auxFields().localVector());PYLITH_CHECK_ERROR(err);

  // Compute the local Jacobian
  err = DMGetLabel(dmMesh, "material-id", &dmLabel);PYLITH_CHECK_ERROR(err);
  err = DMLabelGetStratumBounds(dmLabel, id(), &cStart, &cEnd);PYLITH_CHECK_ERROR(err);
  err = DMPlexComputeJacobian_Internal(dmMesh, cStart, cEnd, t, tshift, solution.localVector(), solutionDot.localVector(), jacobianMat, precondMat, NULL);PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_END;
} // _computeJacobian


// ----------------------------------------------------------------------
// Compute Jacobian using current kernels.
void
pylith::materials::MaterialNew::_computeJacobianInverseLumped(pylith::topology::Field* jacobian,
							      const PylithReal t,
							      const PylithReal dt,
							      const PylithReal tshift,
							      const pylith::topology::Field& solution,
							      const pylith::topology::Field& solutionDot)
{ // _computeJacobianInverseLumped
  PYLITH_METHOD_BEGIN;

  assert(jacobian);

  PetscDS prob = NULL;
  PetscInt cStart = 0, cEnd = 0;
  PetscErrorCode err;

  PetscDM dmMesh = solution.dmMesh();
  PetscDM dmAux = _auxFields->dmMesh();
  PetscDMLabel dmLabel;

  // Pointwise function have been set in DS
  err = DMGetDS(dmMesh, &prob);PYLITH_CHECK_ERROR(err);

  // Get auxiliary data
  err = PetscObjectCompose((PetscObject) dmMesh, "dmAux", (PetscObject) dmAux);PYLITH_CHECK_ERROR(err);
  err = PetscObjectCompose((PetscObject) dmMesh, "A", (PetscObject) auxFields().localVector());PYLITH_CHECK_ERROR(err);

  PetscVec vecRowSum = NULL;
  err = DMGetGlobalVector(dmMesh, &vecRowSum);PYLITH_CHECK_ERROR(err);
  err = VecSet(vecRowSum, 1.0);PYLITH_CHECK_ERROR(err);

  // Compute the local Jacobian action
  err = DMGetLabel(dmMesh, "material-id", &dmLabel);PYLITH_CHECK_ERROR(err);
  err = DMLabelGetStratumBounds(dmLabel, id(), &cStart, &cEnd);PYLITH_CHECK_ERROR(err);
#if 0 // NOT YET IMPLEMENTED IN petsc-dev knepley/pylith
  err = DMPlexComputeJacobianAction_Internal(dmMesh, cStart, cEnd, t, tshift, vecRowSum, NULL, vecRowSum, jacobian->localVector(), NULL);PYLITH_CHECK_ERROR(err);

  // Compute the Jacobian inverse.
  err = VecReciprocal(jacobian->localVector());PYLITH_CHECK_ERROR(err);
#endif

  PYLITH_METHOD_END;
} // _computeJacobianInverseLumped


// End of file 
