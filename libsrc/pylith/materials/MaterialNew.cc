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

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::MaterialNew::MaterialNew(const int dimension) :
  _materialIS(NULL),
  _auxFieldsDB(NULL),
  _auxFieldsQuery(NULL),
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

  delete _materialIS; _materialIS = NULL;
  delete _auxFieldsQuery; _auxFieldsQuery = NULL;

  _auxFieldsDB = 0; // :TODO: Use shared pointer.

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

  // Set finite-element kernels
  _setFEKernels(solution);

  PYLITH_METHOD_END;
} // initialize


// ----------------------------------------------------------------------
// Set discretization information for subfield.
void
pylith::materials::MaterialNew::discretization(const char* name,
					       const pylith::topology::FieldBase::DiscretizeInfo& feInfo)
{ // discretization
  _auxFieldsFEInfo[name] = feInfo;
} // discretization


// ----------------------------------------------------------------------
// Get discretization information for subfield.
const pylith::topology::FieldBase::DiscretizeInfo& 
pylith::materials::MaterialNew::discretization(const char* name) const
{ // discretization
  PYLITH_METHOD_BEGIN;

  discretizations_type::const_iterator iter = _auxFieldsFEInfo.find(name);
  if (iter != _auxFieldsFEInfo.end()) {
    PYLITH_METHOD_RETURN(iter->second);
  } else { // not found so try default
    iter = _auxFieldsFEInfo.find("default");
    if (iter == _auxFieldsFEInfo.end()) {
      throw std::logic_error("Default discretization not set for auxiliary fields.");
    } // if
  } // if/else

  PYLITH_METHOD_RETURN(iter->second); // default
} // discretization


extern "C" PetscErrorCode DMPlexComputeResidual_Internal(DM dm, PetscInt cStart, PetscInt cEnd, PetscReal time, Vec locX, Vec locX_t, Vec locF, void *user);
// ----------------------------------------------------------------------
void
pylith::materials::MaterialNew::integrateResidual(const topology::Field& residual,
							  const PylithScalar t,
							  topology::SolutionFields* const fields)
{ // integrateResidual
  PYLITH_METHOD_BEGIN;

  assert(_logger);
  assert(fields);

  assert(fields);
  assert(_auxFields);

  PetscDS prob = NULL;
  PetscVec dispTpdtVec = NULL;
  PetscInt cStart = 0, cEnd = 0;
  PetscErrorCode err;

  PetscDM dmMesh = fields->get("dispIncr(t->t+dt)").dmMesh();
  PetscDM dmAux = _auxFields->dmMesh();
  PetscDMLabel label;

  // Pointwise function have been set in DS
  err = DMGetDS(dmMesh, &prob);PYLITH_CHECK_ERROR(err);

  // Create full solution
  err = VecDuplicate(residual.localVector(), &dispTpdtVec);PYLITH_CHECK_ERROR(err);
  err = VecWAXPY(dispTpdtVec, 1.0, fields->get("disp(t)").localVector(), fields->get("dispIncr(t->t+dt)").localVector());PYLITH_CHECK_ERROR(err);

  // Get auxiliary data
  err = PetscObjectCompose((PetscObject) dmMesh, "dmAux", (PetscObject) dmAux);PYLITH_CHECK_ERROR(err);
  err = PetscObjectCompose((PetscObject) dmMesh, "A", (PetscObject) auxFields().localVector());PYLITH_CHECK_ERROR(err);

  // Compute the local residual
  err = DMPlexGetLabel(dmMesh, "material-id", &label);PYLITH_CHECK_ERROR(err);
  err = DMLabelGetStratumBounds(label, id(), &cStart, &cEnd);PYLITH_CHECK_ERROR(err);
  err = DMPlexComputeResidual_Internal(dmMesh, cStart, cEnd, PETSC_MIN_REAL, dispTpdtVec, NULL, residual.localVector(), NULL);PYLITH_CHECK_ERROR(err);
  err = VecDestroy(&dispTpdtVec);PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_END;
} // integrateResidual

// ----------------------------------------------------------------------
// Compute stiffness matrix.
void
pylith::materials::MaterialNew::integrateJacobian(topology::Jacobian* jacobian,
							  const PylithScalar t,
							  topology::SolutionFields* fields)
{ // integrateJacobian
  PYLITH_METHOD_BEGIN;

  assert(_logger);
  assert(jacobian);
  assert(fields);

  PetscDS prob = NULL;
  PetscVec dispTpdtVec = NULL;
  PetscErrorCode err;

  const PetscMat jacobianMat = jacobian->matrix();assert(jacobianMat);
  const topology::Field& solnField = fields->get("dispIncr(t->t+dt)");
  PetscDM dmMesh = solnField.dmMesh();
  PetscDM dmAux = _auxFields->dmMesh();

  // Pointwise function have been set in DS
  err = DMGetDS(dmMesh, &prob);PYLITH_CHECK_ERROR(err);

  // Create full solution
  err = VecDuplicate(fields->get("disp(t)").localVector(), &dispTpdtVec);PYLITH_CHECK_ERROR(err);
  err = VecWAXPY(dispTpdtVec, 1.0, fields->get("disp(t)").localVector(), solnField.localVector());PYLITH_CHECK_ERROR(err);

  // Get auxiliary data
  err = PetscObjectCompose((PetscObject) dmMesh, "dmAux", (PetscObject) dmAux);PYLITH_CHECK_ERROR(err);
  err = PetscObjectCompose((PetscObject) dmMesh, "A", (PetscObject) auxFields().localVector());PYLITH_CHECK_ERROR(err);

  // Compute the local Jacobian
  err = DMPlexSNESComputeJacobianFEM(dmMesh, dispTpdtVec, jacobianMat, jacobianMat, NULL);PYLITH_CHECK_ERROR(err);
  err = VecDestroy(&dispTpdtVec);PYLITH_CHECK_ERROR(err);

  _needNewJacobian = false;

  PYLITH_METHOD_END;
} // integrateJacobian

// End of file 
