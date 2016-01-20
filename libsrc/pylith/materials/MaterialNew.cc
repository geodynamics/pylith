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
  _computeResidual(residual, t, dt, solution);

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
  _computeJacobian(jacobian, preconditioner, t, dt, solution);

  PYLITH_METHOD_END;
} // computeRHSJacobian


// ----------------------------------------------------------------------
// Compute LHS residual for F(t,s,\dot{s}).
void
pylith::materials::MaterialNew::computeLHSResidual(pylith::topology::Field* residual,
						   const PylithReal t,
						   const PylithReal dt,
						   const pylith::topology::Field& solution)
{ // computeLHSResidual
  PYLITH_METHOD_BEGIN;

  _setFEKernelsLHSResidual(solution);
  _computeResidual(residual, t, dt, solution);

  PYLITH_METHOD_END;
} // computeLHSResidual

// ----------------------------------------------------------------------
// Compute LHS Jacobian for F(t,s,\dot{s}).
void
pylith::materials::MaterialNew::computeLHSJacobianImplicit(pylith::topology::Jacobian* jacobian,
							   pylith::topology::Jacobian* preconditioner,
							   const PylithReal t,
							   const PylithReal dt,
							   const pylith::topology::Field& solution)
{ // computeLHSJacobianImplicit
  PYLITH_METHOD_BEGIN;

  _setFEKernelsLHSJacobianImplicit(solution);
  _computeJacobian(jacobian, preconditioner, t, dt, solution);

  PYLITH_METHOD_END;
} // computeLHSJacobianImplicit


// ----------------------------------------------------------------------
// Compute LHS Jacobian for F(t,s,\dot{s}).
void
pylith::materials::MaterialNew::computeLHSJacobianExplicit(pylith::topology::Jacobian* jacobian,
							   pylith::topology::Jacobian* preconditioner,
							   const PylithReal t,
							   const PylithReal dt,
							   const pylith::topology::Field& solution)
{ // computeLHSJacobianExplicit
  PYLITH_METHOD_BEGIN;

  _setFEKernelsLHSJacobianExplicit(solution);
  _computeJacobian(jacobian, preconditioner, t, dt, solution);

  PYLITH_METHOD_END;
} // computeLHSJacobianExplicit


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
						 const pylith::topology::Field& solution)
{ // _computeResidual
  PYLITH_METHOD_BEGIN;

  assert(residual);
  assert(_auxFields);

  PetscDS prob = NULL;
  PetscInt cStart = 0, cEnd = 0;
  PetscErrorCode err;

  PetscDM dmMesh = solution.dmMesh();
  PetscDM dmAux = _auxFields->dmMesh();
  PetscDMLabel label;

  // Pointwise function have been set in DS
  err = DMGetDS(dmMesh, &prob);PYLITH_CHECK_ERROR(err);

  // Get auxiliary data
  err = PetscObjectCompose((PetscObject) dmMesh, "dmAux", (PetscObject) dmAux);PYLITH_CHECK_ERROR(err);
  err = PetscObjectCompose((PetscObject) dmMesh, "A", (PetscObject) _auxFields->localVector());PYLITH_CHECK_ERROR(err);

  // Compute the local residual
  err = DMGetLabel(dmMesh, "material-id", &label);PYLITH_CHECK_ERROR(err);
  err = DMLabelGetStratumBounds(label, id(), &cStart, &cEnd);PYLITH_CHECK_ERROR(err);
  err = DMPlexComputeResidual_Internal(dmMesh, cStart, cEnd, PETSC_MIN_REAL, solution.localVector(), NULL, residual->localVector(), NULL);PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_END;
} // _computeResidual
  

// ----------------------------------------------------------------------
// Compute Jacobian using current kernels.
void
pylith::materials::MaterialNew::_computeJacobian(pylith::topology::Jacobian* jacobian,
						 pylith::topology::Jacobian* preconditioner,
						 const PylithReal t,
						 const PylithReal dt,
						 const pylith::topology::Field& solution)
{ // _computeJacobian
  PYLITH_METHOD_BEGIN;

  assert(_logger);
  assert(jacobian);

  PetscDS prob = NULL;
  PetscErrorCode err;

  const PetscMat jacobianMat = jacobian->matrix();assert(jacobianMat);
  const PetscMat precondMat = preconditioner->matrix();assert(precondMat);
  PetscDM dmMesh = solution.dmMesh();
  PetscDM dmAux = _auxFields->dmMesh();

  // Pointwise function have been set in DS
  err = DMGetDS(dmMesh, &prob);PYLITH_CHECK_ERROR(err);

  // Get auxiliary data
  err = PetscObjectCompose((PetscObject) dmMesh, "dmAux", (PetscObject) dmAux);PYLITH_CHECK_ERROR(err);
  err = PetscObjectCompose((PetscObject) dmMesh, "A", (PetscObject) auxFields().localVector());PYLITH_CHECK_ERROR(err);

  // Compute the local Jacobian
  err = DMPlexSNESComputeJacobianFEM(dmMesh, solution.localVector(), jacobianMat, precondMat, NULL);PYLITH_CHECK_ERROR(err);

  _needNewJacobian = false;

  PYLITH_METHOD_END;
} // _computeJacobian


// End of file 
