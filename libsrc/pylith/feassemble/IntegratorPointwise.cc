// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, Rice University
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

#include "IntegratorPointwise.hh" // implementation of class methods

#include "Quadrature.hh" // USES Quadrature
#include "CellGeometry.hh" // USES CellGeometry

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/Jacobian.hh" // USES Jacobian
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional
#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField

#include "pylith/utils/array.hh" // USES scalar_array
#include "pylith/utils/EventLogger.hh" // USES EventLogger

#include <strings.h> // USES strcasecmp()
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <iostream> // USES std::cerr

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::IntegratorPointwise::IntegratorPointwise(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::IntegratorPointwise::~IntegratorPointwise(void)
{ // destructor
  deallocate();
} // destructor
  
// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::feassemble::IntegratorPointwise::deallocate(void)
{ // deallocate
  PYLITH_METHOD_BEGIN;

  Integrator::deallocate();

  PYLITH_METHOD_END;
} // deallocate

// ----------------------------------------------------------------------
// Return auxiliary fields for this problem
const pylith::topology::Field&
pylith::feassemble::IntegratorPointwise::auxFields() const
{ // auxFields
  PYLITH_METHOD_BEGIN;

  PYLITH_METHOD_RETURN(*_fieldsAux);
} // auxFields

// ----------------------------------------------------------------------
// Determine whether we need to recompute the Jacobian.
bool
pylith::feassemble::IntegratorPointwise::needNewJacobian(void)
{ // needNewJacobian
  PYLITH_METHOD_BEGIN;

  if (!_needNewJacobian)
    _needNewJacobian = true;
  PYLITH_METHOD_RETURN(_needNewJacobian);
} // needNewJacobian

// ----------------------------------------------------------------------
// Initialize integrator.
void
pylith::feassemble::IntegratorPointwise::initialize(const topology::Mesh& mesh)
{ // initialize
  PYLITH_METHOD_BEGIN;

  assert(_quadrature);

  _initializeLogger();

  // Compute geometry for quadrature operations.
  _quadrature->initializeGeometry();

  // Optimize coordinate retrieval in closure
  PetscDM dmMesh = mesh.dmMesh();assert(dmMesh);
  topology::CoordsVisitor::optimizeClosure(dmMesh);

  _isJacobianSymmetric = false;

  // Set up gravity field database for querying
  if (_gravityField) {
    const int spaceDim = _quadrature->spaceDim();
    _gravityField->open();
    if (1 == spaceDim){
      const char* queryNames[] = { "x"};
      _gravityField->queryVals(queryNames, spaceDim);
    } else if (2 == spaceDim){
      const char* queryNames[] = { "x", "y"};
      _gravityField->queryVals(queryNames, spaceDim);
    } else if (3 == spaceDim){
      const char* queryNames[] = { "x", "y", "z"};
      _gravityField->queryVals(queryNames, spaceDim);
    } else {
      std::cerr << "Bad spatial dimension '" << spaceDim << "'." << std::endl;
      assert(0);
      throw std::logic_error("Bad spatial dimension in IntegratorPointwise.");
    } // else
  } // if

  PYLITH_METHOD_END;
} // initialize
  
// ----------------------------------------------------------------------
// Set time step for advancing from time t to time t+dt.
void
pylith::feassemble::IntegratorPointwise::timeStep(const PylithScalar dt)
{ // timeStep
  PYLITH_METHOD_BEGIN;

  _dt = dt;

  PYLITH_METHOD_END;
} // timeStep

// ----------------------------------------------------------------------
// Get stable time step for advancing from time t to time t+dt.
PylithScalar
pylith::feassemble::IntegratorPointwise::stableTimeStep(const topology::Mesh& mesh) const
{ // stableTimeStep
  PYLITH_METHOD_BEGIN;

  PYLITH_METHOD_RETURN(0.0);
} // stableTimeStep

// ----------------------------------------------------------------------
// Update state variables as needed.
void
pylith::feassemble::IntegratorPointwise::updateStateVars(const PylithScalar t,
							  topology::SolutionFields* const fields)
{ // updateState
  PYLITH_METHOD_BEGIN;

  assert(_quadrature);
  assert(fields);

#if 0
  // No need to update state vars if material doesn't have any.
  if (!_material->hasStateVars())
    PYLITH_METHOD_END;

  // Get cell information that doesn't depend on particular cell
  const int cellDim = _quadrature->cellDim();
  const int numQuadPts = _quadrature->numQuadPts();
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int numCorners = _quadrature->refGeometry().numCorners();
  const int tensorSize = _material->tensorSize();
  totalStrain_fn_type calcTotalStrainFn;
  if (2 == cellDim) {
    calcTotalStrainFn = &pylith::feassemble::IntegratorPointwise::_calcTotalStrain2D;
  } else if (3 == cellDim) {
    calcTotalStrainFn = &pylith::feassemble::IntegratorPointwise::_calcTotalStrain3D;
  } else {
      std::cerr << "Bad cell dimension '" << cellDim << "'." << std::endl;
      assert(0);
      throw std::logic_error("Bad cell dimension in IntegratorPointwise::updateStateVars().");
  } // else

  // Allocate arrays for cell data.
  scalar_array strainCell(numQuadPts*tensorSize);
  strainCell = 0.0;

  // Get cell information
  PetscDM dmMesh = fields->mesh().dmMesh();assert(dmMesh);
  assert(_materialIS);
  const PetscInt* cells = _materialIS->points();
  const PetscInt numCells = _materialIS->size();

  // Setup visitors.
  scalar_array dispCell(numBasis*spaceDim);
  topology::VecVisitorMesh dispVisitor(fields->get("disp(t)"), "displacement");
  dispVisitor.optimizeClosure();

  scalar_array coordsCell(numCorners*spaceDim);
  topology::CoordsVisitor coordsVisitor(dmMesh);

  // Loop over cells
  for(PetscInt c = 0; c < numCells; ++c) {
    const PetscInt cell = cells[c];
    // Retrieve geometry information for current cell
    coordsVisitor.getClosure(&coordsCell, cell);
    _quadrature->computeGeometry(&coordsCell[0], coordsCell.size(), cell);
    const scalar_array& basisDeriv = _quadrature->basisDeriv();

    // Get physical properties and state variables for cell.
    _material->retrievePropsAndVars(cell);

    // Restrict input fields to cell
    dispVisitor.getClosure(&dispCell, cell);

    // Compute strains
    calcTotalStrainFn(&strainCell, basisDeriv, &dispCell[0], numBasis, spaceDim, numQuadPts);

    // Update material state
    _material->updateStateVars(strainCell, cell);
  } // for
#endif

  PYLITH_METHOD_END;
} // updateStateVars

// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::feassemble::IntegratorPointwise::verifyConfiguration(const topology::Mesh& mesh) const
{ // verifyConfiguration
  PYLITH_METHOD_BEGIN;

  assert(_quadrature);

  const int dimension = mesh.dimension();

  // check compatibility of mesh and quadrature scheme
  if (_quadrature->cellDim() != dimension) {
    std::ostringstream msg;
    msg << "Quadrature is incompatible with cells'.\n"
	<< "Dimension of mesh: " << dimension
	<< ", dimension of quadrature: " << _quadrature->cellDim()
	<< ".";
    throw std::runtime_error(msg.str());
  } // if

#if 0
  PetscDM dmMesh = mesh.dmMesh();assert(dmMesh);
  const bool includeOnlyCells = true;
  topology::StratumIS materialIS(dmMesh, "material-id", _material->id(), includeOnlyCells);
  const PetscInt* cells = NULL /*materialIS.points()*/;
  const PetscInt numCells = 0 /*materialIS.size()*/;

  const int numCorners = _quadrature->refGeometry().numCorners();
  PetscInt vStart, vEnd;
  PetscErrorCode err;

  err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);PYLITH_CHECK_ERROR(err);
  for(PetscInt c = 0; c < numCells; ++c) {
    const PetscInt cell = cells[c];
    PetscInt cellNumCorners = 0, closureSize, *closure = NULL;

    err = DMPlexGetTransitiveClosure(dmMesh, cell, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
    for (PetscInt cl = 0; cl < closureSize*2; cl += 2) {
      if ((closure[cl] >= vStart) && (closure[cl] < vEnd)) ++cellNumCorners;
    }
    err = DMPlexRestoreTransitiveClosure(dmMesh, cell, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
    if (numCorners != cellNumCorners) {
      std::ostringstream msg;
      msg << "Quadrature is incompatible with cell in material '"
	  << _material->label() << "'.\n"
	  << "Cell " << cell << " has " << cellNumCorners
	  << " vertices but quadrature reference cell has "
	  << numCorners << " vertices.";
      throw std::runtime_error(msg.str());
    } // if
  } // for
#endif

  PYLITH_METHOD_END;
} // verifyConfiguration

// ----------------------------------------------------------------------
// Initialize logger.
void
pylith::feassemble::IntegratorPointwise::_initializeLogger(void)
{ // initializeLogger
  PYLITH_METHOD_BEGIN;

  delete _logger; _logger = new utils::EventLogger;
  assert(_logger);
  _logger->className("PointwiseIntegrator");
  _logger->initialize();
  _logger->registerEvent("ElIR setup");
  _logger->registerEvent("ElIR geometry");
  _logger->registerEvent("ElIR compute");
  _logger->registerEvent("ElIR restrict");
  _logger->registerEvent("ElIR stateVars");
  _logger->registerEvent("ElIR update");
 
  _logger->registerEvent("ElIJ setup");
  _logger->registerEvent("ElIJ geometry");
  _logger->registerEvent("ElIJ compute");
  _logger->registerEvent("ElIJ restrict");
  _logger->registerEvent("ElIJ stateVars");
  _logger->registerEvent("ElIJ update");

  PYLITH_METHOD_END;
} // initializeLogger

// ----------------------------------------------------------------------
void
pylith::feassemble::IntegratorPointwise::integrateResidual(const topology::Field& residual,
							  const PylithScalar t,
							  topology::SolutionFields* const fields)
{ // integrateResidual
  PYLITH_METHOD_BEGIN;

  assert(_logger);
  assert(fields);

  PetscDM        dmMesh = _dm; /* MUST match dm in 'fields' */
  PetscDS        prob;
  Vec            dispTpdtVec;
  PetscErrorCode err;

  // Pointwise function have been set in DS
  // :TODO: call _setFEKernels from initialize()
  err = DMGetDS(dmMesh, &prob);PYLITH_CHECK_ERROR(err);
  // Create full solution
  err = VecDuplicate(residual.localVector(), &dispTpdtVec);PYLITH_CHECK_ERROR(err);
  err = VecWAXPY(dispTpdtVec, 1.0, fields->get("disp(t)").localVector(), fields->get("dispIncr(t->t+dt)").localVector());PYLITH_CHECK_ERROR(err);
  // Get auxiliary data
  err = PetscObjectCompose((PetscObject) dmMesh, "dmAux", (PetscObject) _dmAux);PYLITH_CHECK_ERROR(err);
  err = PetscObjectCompose((PetscObject) dmMesh, "A", (PetscObject) auxFields().localVector());PYLITH_CHECK_ERROR(err);
  // Compute the local residual
  err = DMPlexSNESComputeResidualFEM(dmMesh, dispTpdtVec, residual.localVector(), NULL);PYLITH_CHECK_ERROR(err);
  err = VecDestroy(&dispTpdtVec);PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_END;
} // integrateResidual

// ----------------------------------------------------------------------
// Compute stiffness matrix.
void
pylith::feassemble::IntegratorPointwise::integrateJacobian(topology::Jacobian* jacobian,
							  const PylithScalar t,
							  topology::SolutionFields* fields)
{ // integrateJacobian
  PYLITH_METHOD_BEGIN;

  assert(_logger);
  assert(jacobian);
  assert(fields);

  PetscDM        dmMesh = _dm; /* MUST match dm in 'fields' */
  PetscDS        prob;
  const PetscMat jacobianMat = jacobian->matrix();assert(jacobianMat);
  Vec            dispTpdtVec;
  PetscErrorCode err;

  // Pointwise function have been set in DS
  // :TODO: call _setFEKernels from initialize()
  err = DMGetDS(dmMesh, &prob);PYLITH_CHECK_ERROR(err);
  // Create full solution
  err = VecDuplicate(fields->get("disp(t)").localVector(), &dispTpdtVec);PYLITH_CHECK_ERROR(err);
  err = VecWAXPY(dispTpdtVec, 1.0, fields->get("disp(t)").localVector(), fields->get("dispIncr(t->t+dt)").localVector());PYLITH_CHECK_ERROR(err);
  // Get auxiliary data
  err = PetscObjectCompose((PetscObject) dmMesh, "dmAux", (PetscObject) _dmAux);PYLITH_CHECK_ERROR(err);
  err = PetscObjectCompose((PetscObject) dmMesh, "A", (PetscObject) auxFields().localVector());PYLITH_CHECK_ERROR(err);
  // Compute the local Jacobian
  err = DMPlexSNESComputeJacobianFEM(dmMesh, dispTpdtVec, jacobianMat, jacobianMat, NULL);PYLITH_CHECK_ERROR(err);
  err = VecDestroy(&dispTpdtVec);PYLITH_CHECK_ERROR(err);

  _needNewJacobian = false;

  PYLITH_METHOD_END;
} // integrateJacobian

// End of file 
