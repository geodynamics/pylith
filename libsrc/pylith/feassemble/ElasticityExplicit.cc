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

#include "ElasticityExplicit.hh" // implementation of class methods

#include "Quadrature.hh" // USES Quadrature
#include "CellGeometry.hh" // USES CellGeometry

#include "pylith/materials/ElasticMaterial.hh" // USES ElasticMaterial
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields
#include "pylith/topology/Jacobian.hh" // USES Jacobian
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor

#include "pylith/utils/array.hh" // USES scalar_array
#include "pylith/utils/EventLogger.hh" // USES EventLogger
#include "pylith/utils/macrodefs.h" // USES CALL_MEMBER_FN

#include "pylith/utils/error.h" // USES PYLITH_CHECK_ERROR

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField
#include "spatialdata/units/Nondimensional.hh" // USES Nondimendional

#include "petscmat.h" // USES PetscMat

#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

//#define DETAILED_EVENT_LOGGING

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::ElasticityExplicit::ElasticityExplicit(void) :
  _dtm1(-1.0),
  _normViscosity(0.1)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::ElasticityExplicit::~ElasticityExplicit(void)
{ // destructor
  deallocate();
} // destructor
  
// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::feassemble::ElasticityExplicit::deallocate(void)
{ // deallocate
  PYLITH_METHOD_BEGIN;

  IntegratorElasticity::deallocate();

  PYLITH_METHOD_END;
} // deallocate
  
// ----------------------------------------------------------------------
// Set time step for advancing from time t to time t+dt.
void
pylith::feassemble::ElasticityExplicit::timeStep(const PylithScalar dt)
{ // timeStep
  PYLITH_METHOD_BEGIN;
  
  if (_dt != -1.0) {
    _dtm1 = _dt;
  } else {
    _dtm1 = dt;
  } // if/else
  _dt = dt;
  assert(_dt == _dtm1); // For now, don't allow variable time step
  if (_material) {
    _material->timeStep(_dt);
  } // if

  PYLITH_METHOD_END;
} // timeStep

// ----------------------------------------------------------------------
// Get stable time step for advancing from time t to time t+dt.
PylithScalar
pylith::feassemble::ElasticityExplicit::stableTimeStep(const topology::Mesh& mesh) const
{ // stableTimeStep
  PYLITH_METHOD_BEGIN;
  
  assert(_material);
  PYLITH_METHOD_RETURN(_material->stableTimeStepExplicit(mesh, _quadrature));
} // stableTimeStep

// ----------------------------------------------------------------------
// Set normalized viscosity for numerical damping.
void
pylith::feassemble::ElasticityExplicit::normViscosity(const PylithScalar viscosity)
{ // normViscosity
  PYLITH_METHOD_BEGIN;
  
  if (viscosity < 0.0) {
    std::ostringstream msg;
    msg << "Normalized viscosity (" << viscosity << ") must be nonnegative.";
    throw std::runtime_error(msg.str());
  } // if

  _normViscosity = viscosity;

  PYLITH_METHOD_END;
} // normViscosity

// ----------------------------------------------------------------------
// Integrate constributions to residual term (r) for operator.
void
pylith::feassemble::ElasticityExplicit::integrateResidual(const topology::Field& residual,
							  const PylithScalar t,
							  topology::SolutionFields* const fields)
{ // integrateResidual
  PYLITH_METHOD_BEGIN;
  
  /// Member prototype for _elasticityResidualXD()
  typedef void (pylith::feassemble::ElasticityExplicit::*elasticityResidual_fn_type)
    (const scalar_array&);

  assert(_quadrature);
  assert(_material);
  assert(_logger);
  assert(fields);

  const int setupEvent = _logger->eventId("ElIR setup");
  const int computeEvent = _logger->eventId("ElIR compute");
#if defined(DETAILED_EVENT_LOGGING)
  const int geometryEvent = _logger->eventId("ElIR geometry");
  const int restrictEvent = _logger->eventId("ElIR restrict");
  const int stateVarsEvent = _logger->eventId("ElIR stateVars");
  const int stressEvent = _logger->eventId("ElIR stress");
  const int updateEvent = _logger->eventId("ElIR update");
#endif

  _logger->eventBegin(setupEvent);

  // Get cell geometry information that doesn't depend on cell
  const int numQuadPts = _quadrature->numQuadPts();
  const scalar_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == size_t(numQuadPts));
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int cellDim = _quadrature->cellDim();
  const int tensorSize = _material->tensorSize();
  /** :TODO:
   *
   * If cellDim and spaceDim are different, we need to transform
   * displacements into cellDim, compute action, and transform result
   * back into spaceDim. We get this information from the Jacobian and
   * inverse of the Jacobian.
   */
  if (cellDim != spaceDim)
    throw std::logic_error("Integration for cells with spatial dimensions "
         "different than the spatial dimension of the "
         "domain not implemented yet.");

  // Set variables dependent on dimension of cell
  totalStrain_fn_type calcTotalStrainFn;
  elasticityResidual_fn_type elasticityResidualFn;
  if (2 == cellDim) {
    elasticityResidualFn =
      &pylith::feassemble::ElasticityExplicit::_elasticityResidual2D;
    calcTotalStrainFn =
      &pylith::feassemble::IntegratorElasticity::_calcTotalStrain2D;
  } else if (3 == cellDim) {
    elasticityResidualFn =
      &pylith::feassemble::ElasticityExplicit::_elasticityResidual3D;
    calcTotalStrainFn =
      &pylith::feassemble::IntegratorElasticity::_calcTotalStrain3D;
  } else {
    assert(0);
    throw std::runtime_error("Error unknown cell dimension.");
  } // if/else

  // Allocate vectors for cell values.
  scalar_array strainCell(numQuadPts*tensorSize);
  strainCell = 0.0;
  scalar_array gravVec(spaceDim);
  scalar_array quadPtsGlobal(numQuadPts*spaceDim);

  // Get cell information
  PetscDM dmMesh = fields->mesh().dmMesh();assert(dmMesh);
  assert(_materialIS);
  const PetscInt* cells = _materialIS->points();
  const PetscInt numCells = _materialIS->size();

  // Setup field visitors.
  scalar_array accCell(numBasis*spaceDim);
  topology::VecVisitorMesh accVisitor(fields->get("acceleration(t)"), "displacement");
  accVisitor.optimizeClosure();

  scalar_array velCell(numBasis*spaceDim);
  topology::VecVisitorMesh velVisitor(fields->get("velocity(t)"), "displacement");
  velVisitor.optimizeClosure();

  scalar_array dispCell(numBasis*spaceDim);
  scalar_array dispAdjCell(numBasis*spaceDim);
  topology::VecVisitorMesh dispVisitor(fields->get("disp(t)"), "displacement");
  dispVisitor.optimizeClosure();

  topology::VecVisitorMesh residualVisitor(residual, "displacement");
  residualVisitor.optimizeClosure();

  scalar_array coordsCell(numBasis*spaceDim); // :KULDGE: Update numBasis to numCorners after implementing higher order
  topology::CoordsVisitor coordsVisitor(dmMesh);

  _material->createPropsAndVarsVisitors();

  assert(_normalizer);
  const PylithScalar lengthScale = _normalizer->lengthScale();
  const PylithScalar gravityScale = _normalizer->pressureScale() / (_normalizer->lengthScale() * _normalizer->densityScale());

  const PylithScalar dt = _dt;assert(dt > 0);
  const PylithScalar viscosity = dt*_normViscosity;assert(_normViscosity >= 0.0);

  // Get parameters used in integration.
  scalar_array valuesIJ(numBasis);

  _logger->eventEnd(setupEvent);
#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventBegin(computeEvent);
#endif

  // Loop over cells
  for(PetscInt c = 0; c < numCells; ++c) {
    const PetscInt cell = cells[c];
    // Compute geometry information for current cell
#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(geometryEvent);
#endif
    coordsVisitor.getClosure(&coordsCell, cell);
    _quadrature->computeGeometry(&coordsCell[0], coordsCell.size(), cell);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(geometryEvent);
    _logger->eventBegin(stateVarsEvent);
#endif

    // Get state variables for cell.
    _material->retrievePropsAndVars(cell);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(stateVarsEvent);
    _logger->eventBegin(restrictEvent);
#endif

    // Reset element vector to zero
    _resetCellVector();

    // Restrict input fields to cell
    accVisitor.getClosure(&accCell, cell);
    velVisitor.getClosure(&velCell, cell);
    dispVisitor.getClosure(&dispCell, cell);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(restrictEvent);
    _logger->eventBegin(computeEvent);
#endif

    // Get cell geometry information that depends on cell
    const scalar_array& basis = _quadrature->basis();
    const scalar_array& basisDeriv = _quadrature->basisDeriv();
    const scalar_array& jacobianDet = _quadrature->jacobianDet();
    const scalar_array& quadPtsNondim = _quadrature->quadPts();

    // Compute body force vector if gravity is being used.
    if (_gravityField) {
      const spatialdata::geocoords::CoordSys* cs = fields->mesh().coordsys();assert(cs);

      // Get density at quadrature points for this cell
      const scalar_array& density = _material->calcDensity();

      quadPtsGlobal = quadPtsNondim;
      _normalizer->dimensionalize(&quadPtsGlobal[0], quadPtsGlobal.size(), lengthScale);

      // Compute action for element body forces
      spatialdata::spatialdb::SpatialDB* db = _gravityField;
      for (int iQuad = 0; iQuad < numQuadPts; ++iQuad) {
        const int err = db->query(&gravVec[0], gravVec.size(), &quadPtsGlobal[0], spaceDim, cs);
        if (err) {
          throw std::runtime_error("Unable to get gravity vector for point.");
        } // if
        _normalizer->nondimensionalize(&gravVec[0], gravVec.size(), gravityScale);
        const PylithScalar wt = quadWts[iQuad] * jacobianDet[iQuad] * density[iQuad];
        for (int iBasis = 0, iQ = iQuad * numBasis; iBasis < numBasis; ++iBasis) {
          const PylithScalar valI = wt * basis[iQ + iBasis];
          for (int iDim = 0; iDim < spaceDim; ++iDim) {
            _cellVector[iBasis*spaceDim+iDim] += valI * gravVec[iDim];
          } // for
        } // for
      } // for
      PetscLogFlops(numQuadPts * (2 + numBasis * (1 + 2 * spaceDim)));
    } // if

    // Compute action for inertial terms
    const scalar_array& density = _material->calcDensity();
    valuesIJ = 0.0;
    for (int iQuad = 0; iQuad < numQuadPts; ++iQuad) {
      const PylithScalar wt = quadWts[iQuad] * jacobianDet[iQuad] * density[iQuad];
      const int iQ = iQuad * numBasis;
      PylithScalar valJ = 0.0;
      for (int jBasis = 0; jBasis < numBasis; ++jBasis) {
	valJ += basis[iQ + jBasis];
      } // for
      valJ *= wt;
      for (int iBasis = 0; iBasis < numBasis; ++iBasis) {
	valuesIJ[iBasis] += basis[iQ + iBasis] * valJ;
      } // for
    } // for
    for (int iBasis = 0; iBasis < numBasis; ++iBasis) {
      for (int iDim = 0; iDim < spaceDim; ++iDim) {
	_cellVector[iBasis*spaceDim+iDim] -= valuesIJ[iBasis] * accCell[iBasis*spaceDim+iDim];
      } // for
    } // for

    // Numerical damping. Compute displacements adjusted by velocity
    // times normalized viscosity.
    for(PetscInt i = 0, dispSize = dispCell.size(); i < dispSize; ++i) {
      dispAdjCell[i] = dispCell[i] + viscosity * velCell[i];
    } // for

#if defined(DETAILED_EVENT_LOGGING)
    PetscLogFlops(numQuadPts*(4+numBasis*3));
    _logger->eventEnd(computeEvent);
    _logger->eventBegin(stressEvent);
#endif

    // Compute B(transpose) * sigma, first computing strains
    calcTotalStrainFn(&strainCell, basisDeriv, &dispAdjCell[0], numBasis, spaceDim, numQuadPts);
    const scalar_array& stressCell = _material->calcStress(strainCell, false);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(stressEvent);
    _logger->eventBegin(computeEvent);
#endif

    CALL_MEMBER_FN(*this, elasticityResidualFn)(stressCell);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(computeEvent);
    _logger->eventBegin(updateEvent);
#endif

    // Assemble cell contribution into field
    residualVisitor.setClosure(&_cellVector[0], _cellVector.size(), cell, ADD_VALUES);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(updateEvent);
#endif
  } // for
  _material->destroyPropsAndVarsVisitors();

#if !defined(DETAILED_EVENT_LOGGING)
  PetscLogFlops(numCells*numQuadPts*(4+numBasis*3));
  _logger->eventEnd(computeEvent);
#endif

  PYLITH_METHOD_END;
} // integrateResidualLumped

// ----------------------------------------------------------------------
// Compute matrix associated with operator.
void
pylith::feassemble::ElasticityExplicit::integrateJacobian(topology::Jacobian* jacobian,
							  const PylithScalar t,
							  topology::SolutionFields* fields)
{ // integrateJacobian
  PYLITH_METHOD_BEGIN;
  
  throw std::logic_error("ElasticityExplicit::integrateJacobian() not implemented. Use integrateJacobian(lumped) instead.");

  PYLITH_METHOD_END;
} // integrateJacobian

// ----------------------------------------------------------------------
// Compute matrix associated with operator.
void
pylith::feassemble::ElasticityExplicit::integrateJacobian(topology::Field* jacobian,
							  const PylithScalar t,
							  topology::SolutionFields* fields)
{ // integrateJacobian
  PYLITH_METHOD_BEGIN;
  
  assert(_quadrature);
  assert(_material);
  assert(jacobian);
  assert(fields);

  const int setupEvent = _logger->eventId("ElIJ setup");
  const int computeEvent = _logger->eventId("ElIJ compute");
#if defined(DETAILED_EVENT_LOGGING)
  const int geometryEvent = _logger->eventId("ElIJ geometry");
  const int restrictEvent = _logger->eventId("ElIJ restrict");
  const int stateVarsEvent = _logger->eventId("ElIJ stateVars");
  const int updateEvent = _logger->eventId("ElIJ update");
#endif

  _logger->eventBegin(setupEvent);

  // Get cell geometry information that doesn't depend on cell
  const int numQuadPts = _quadrature->numQuadPts();
  const scalar_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == size_t(numQuadPts));
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int cellDim = _quadrature->cellDim();
  if (cellDim != spaceDim)
    throw std::logic_error("Don't know how to integrate elasticity " \
			   "contribution to Jacobian matrix for cells with " \
			   "different dimensions than the spatial dimension.");

  // Get cell information
  PetscDM dmMesh = fields->mesh().dmMesh();assert(dmMesh);
  assert(_materialIS);
  const PetscInt* cells = _materialIS->points();
  const PetscInt numCells = _materialIS->size();

  // Get parameters used in integration.
  const PylithScalar dt = _dt;
  const PylithScalar dt2 = dt*dt;
  assert(dt > 0);

  // Setup field visitors.
  scalar_array valuesIJ(numBasis);

  topology::VecVisitorMesh jacobianVisitor(*jacobian, "displacement");
  // Don't optimize closure since we compute the Jacobian only once.

  _material->createPropsAndVarsVisitors();

  scalar_array coordsCell(numBasis*spaceDim); // :KLUDGE: numBasis to numCorners after switching to higher order
  topology::CoordsVisitor coordsVisitor(dmMesh);

  _logger->eventEnd(setupEvent);
#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventBegin(computeEvent);
#endif
  // Loop over cells
  for(PetscInt c = 0; c < numCells; ++c) {
    const PetscInt cell = cells[c];

    // Compute geometry information for current cell
#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(geometryEvent);
#endif
    coordsVisitor.getClosure(&coordsCell, cell);
    _quadrature->computeGeometry(&coordsCell[0], coordsCell.size(), cell);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(geometryEvent);
    _logger->eventBegin(stateVarsEvent);
#endif

    // Get state variables for cell.
    _material->retrievePropsAndVars(cell);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(stateVarsEvent);
    _logger->eventBegin(computeEvent);
#endif

    // Reset element matrix to zero
    _resetCellVector();

    // Get cell geometry information that depends on cell
    const scalar_array& basis = _quadrature->basis();
    const scalar_array& jacobianDet = _quadrature->jacobianDet();

    // Get material physical properties at quadrature points for this cell
    const scalar_array& density = _material->calcDensity();

    // Compute Jacobian for inertial terms
    valuesIJ = 0.0;
    for (int iQuad = 0; iQuad < numQuadPts; ++iQuad) {
      const PylithScalar wt = quadWts[iQuad] * jacobianDet[iQuad] * density[iQuad] / dt2;
      const int iQ = iQuad * numBasis;
      PylithScalar valJ = 0.0;
      for (int jBasis = 0; jBasis < numBasis; ++jBasis) {
        valJ += basis[iQ + jBasis];
      } // for
      valJ *= wt;
      for (int iBasis = 0; iBasis < numBasis; ++iBasis) {
        valuesIJ[iBasis] += basis[iQ + iBasis] * valJ;
      } // for
    } // for
    for (int iBasis = 0; iBasis < numBasis; ++iBasis) {
      for (int iDim = 0; iDim < spaceDim; ++iDim) {
        _cellVector[iBasis*spaceDim+iDim] += valuesIJ[iBasis];
      } // for
    } // for
    
#if defined(DETAILED_EVENT_LOGGING)
    PetscLogFlops(numQuadPts*(4 + numBasis*3) + numBasis*spaceDim);
    _logger->eventEnd(computeEvent);
    _logger->eventBegin(updateEvent);
#endif
    
    // Assemble cell contribution into lumped matrix.
    jacobianVisitor.setClosure(&_cellVector[0], _cellVector.size(), cell, ADD_VALUES);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(updateEvent);
#endif
  } // for
  _material->destroyPropsAndVarsVisitors();

#if !defined(DETAILED_EVENT_LOGGING)
  PetscLogFlops(numCells*(numQuadPts*(4 + numBasis*3) + numBasis*spaceDim));
  _logger->eventEnd(computeEvent);
#endif

  _needNewJacobian = false;
  _material->resetNeedNewJacobian();

  PYLITH_METHOD_END;
} // integrateJacobian


// End of file 
