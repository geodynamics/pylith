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

#include "ElasticityExplicitTri3.hh" // implementation of class methods

#include "Quadrature.hh" // USES Quadrature

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
const int pylith::feassemble::ElasticityExplicitTri3::_spaceDim = 2;
const int pylith::feassemble::ElasticityExplicitTri3::_cellDim = 2;
const int pylith::feassemble::ElasticityExplicitTri3::_tensorSize = 3;
const int pylith::feassemble::ElasticityExplicitTri3::_numBasis = 3;
const int pylith::feassemble::ElasticityExplicitTri3::_numCorners = 3;
const int pylith::feassemble::ElasticityExplicitTri3::_numQuadPts = 1;

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::ElasticityExplicitTri3::ElasticityExplicitTri3(void) :
  _dtm1(-1.0),
  _normViscosity(0.1)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::ElasticityExplicitTri3::~ElasticityExplicitTri3(void)
{ // destructor
  deallocate();
} // destructor
  
// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::feassemble::ElasticityExplicitTri3::deallocate(void)
{ // deallocate
  PYLITH_METHOD_BEGIN;

  IntegratorElasticity::deallocate();

  PYLITH_METHOD_END;
} // deallocate
  
// ----------------------------------------------------------------------
// Set time step for advancing from time t to time t+dt.
void
pylith::feassemble::ElasticityExplicitTri3::timeStep(const PylithScalar dt)
{ // timeStep
  PYLITH_METHOD_BEGIN;

  if (_dt != -1.0)
    _dtm1 = _dt;
  else
    _dtm1 = dt;
  _dt = dt;
  assert(_dt == _dtm1); // For now, don't allow variable time step
  if (_material)
    _material->timeStep(_dt);

  PYLITH_METHOD_END;
} // timeStep

// ----------------------------------------------------------------------
// Get stable time step for advancing from time t to time t+dt.
PylithScalar
pylith::feassemble::ElasticityExplicitTri3::stableTimeStep(const topology::Mesh& mesh) const
{ // stableTimeStep
  PYLITH_METHOD_BEGIN;

  assert(_material);
  PYLITH_METHOD_RETURN(_material->stableTimeStepExplicit(mesh, _quadrature));
} // stableTimeStep

// ----------------------------------------------------------------------
// Set normalized viscosity for numerical damping.
void
pylith::feassemble::ElasticityExplicitTri3::normViscosity(const PylithScalar viscosity)
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
pylith::feassemble::ElasticityExplicitTri3::integrateResidual(const topology::Field& residual,
							      const PylithScalar t,
							      topology::SolutionFields* const fields)
{ // integrateResidual
  PYLITH_METHOD_BEGIN;

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
  assert(_quadrature->numQuadPts() == _numQuadPts);
  assert(_quadrature->numBasis() == _numBasis);
  assert(_quadrature->spaceDim() == _spaceDim);
  assert(_quadrature->cellDim() == _cellDim);
  assert(_material->tensorSize() == _tensorSize);
  const int spaceDim = _spaceDim;
  const int cellDim = _cellDim;
  const int tensorSize = _tensorSize;
  const int numBasis = _numBasis;
  const int numCorners = _numCorners;
  const int numQuadPts = _numQuadPts;
  const int cellVectorSize = numBasis*spaceDim;
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

  scalar_array coordsCell(numCorners*spaceDim);
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
#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(restrictEvent);
#endif

    // Restrict input fields to cell
    accVisitor.getClosure(&accCell, cell);
    velVisitor.getClosure(&velCell, cell);
    dispVisitor.getClosure(&dispCell, cell);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(restrictEvent);
    _logger->eventBegin(geometryEvent);
#endif

    // Compute geometry information for current cell
    coordsVisitor.getClosure(&coordsCell, cell);
    const PylithScalar area = _area(coordsCell);assert(area > 0.0);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(geometryEvent);
    _logger->eventBegin(stateVarsEvent);
#endif

    // Get state variables for cell.
    _material->retrievePropsAndVars(cell);

    // Get density at quadrature points for this cell
    const scalar_array& density = _material->calcDensity();

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(stateVarsEvent);
    _logger->eventBegin(computeEvent);
#endif

    // Reset element vector to zero
    _resetCellVector();

    // Compute body force vector if gravity is being used.
    if (_gravityField) {
      const spatialdata::geocoords::CoordSys* cs = fields->mesh().coordsys();assert(cs);

      quadPtsGlobal = 0.0;
      for (int iBasis=0; iBasis < numBasis; ++iBasis) {
        for (int iDim=0; iDim < spaceDim; ++iDim) {
          quadPtsGlobal[iDim] += coordsCell[iBasis*spaceDim+iDim] / numBasis;
	} // for
      } // for
      _normalizer->dimensionalize(&quadPtsGlobal[0], quadPtsGlobal.size(), lengthScale);

      // Compute action for element body forces
      spatialdata::spatialdb::SpatialDB* db = _gravityField;
      const int err = db->query(&gravVec[0], gravVec.size(), &quadPtsGlobal[0], spaceDim, cs);
      if (err) {
        throw std::runtime_error("Unable to get gravity vector for point.");
      } // if
      _normalizer->nondimensionalize(&gravVec[0], gravVec.size(), gravityScale);
      const PylithScalar wtVertex = density[0] * area / 3.0;
      for (int iBasis=0; iBasis < numBasis; ++iBasis) {
        for (int iDim=0; iDim < spaceDim; ++iDim) {
            _cellVector[iBasis * spaceDim + iDim] += wtVertex * gravVec[iDim];
	} // for
      } // for
      PetscLogFlops(numBasis*spaceDim*2 + numBasis*spaceDim*2);
    } // if

    // Compute action for inertial terms
    const PylithScalar wtVertex = density[0] * area / 3.0;
    assert(cellVectorSize == dispCell.size());
    for(PetscInt i = 0; i < cellVectorSize; ++i) {
      _cellVector[i] -= wtVertex * accCell[i];
    } // for

#if defined(DETAILED_EVENT_LOGGING)
    PetscLogFlops(2 + numBasis*spaceDim*2);
    _logger->eventEnd(computeEvent);
    _logger->eventBegin(stressEvent);
#endif

    // Numerical damping. Compute displacements adjusted by velocity
    // times normalized viscosity.
    for(PetscInt i = 0; i < cellVectorSize; ++i) {
      dispAdjCell[i] = dispCell[i] + viscosity * velCell[i];
    } // for

    // Compute B(transpose) * sigma, first computing strains
    const PylithScalar x0 = coordsCell[0];
    const PylithScalar y0 = coordsCell[1];

    const PylithScalar x1 = coordsCell[2];
    const PylithScalar y1 = coordsCell[3];

    const PylithScalar x2 = coordsCell[4];
    const PylithScalar y2 = coordsCell[5];

    const PylithScalar scaleB = 2.0 * area;
    const PylithScalar b0 = (y1 - y2) / scaleB;
    const PylithScalar c0 = (x2 - x1) / scaleB;

    const PylithScalar b1 = (y2 - y0) / scaleB;
    const PylithScalar c1 = (x0 - x2) / scaleB;

    const PylithScalar b2 = (y0 - y1) / scaleB;
    const PylithScalar c2 = (x1 - x0) / scaleB;

    assert(strainCell.size() == 3);
    strainCell[0] = b2*dispAdjCell[4] + b1*dispAdjCell[2] + b0*dispAdjCell[0];
    
    strainCell[1] = c2*dispAdjCell[5] + c1*dispAdjCell[3] + c0*dispAdjCell[1];

    strainCell[2] = (b2*dispAdjCell[5] + c2*dispAdjCell[4] + b1*dispAdjCell[3] + 
		     c1*dispAdjCell[2] + b0*dispAdjCell[1] + c0*dispAdjCell[0]) / 2.0;

    const scalar_array& stressCell = _material->calcStress(strainCell, false);

#if defined(DETAILED_EVENT_LOGGING)
    PetscLogFlops(34);
    _logger->eventEnd(stressEvent);
    _logger->eventBegin(computeEvent);
#endif

    assert(_cellVector.size() == 6);
    assert(stressCell.size() == 3);
    _cellVector[0] -= (c0*stressCell[2] + b0*stressCell[0]) * area;
    _cellVector[1] -= (b0*stressCell[2] + c0*stressCell[1]) * area;
    _cellVector[2] -= (c1*stressCell[2] + b1*stressCell[0]) * area;
    _cellVector[3] -= (b1*stressCell[2] + c1*stressCell[1]) * area;
    _cellVector[4] -= (c2*stressCell[2] + b2*stressCell[0]) * area;
    _cellVector[5] -= (b2*stressCell[2] + c2*stressCell[1]) * area;

#if defined(DETAILED_EVENT_LOGGING)
    PetscLogFlops(30);
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
  PetscLogFlops(numCells*(2 + numBasis*spaceDim*2 + 34+30));
  _logger->eventEnd(computeEvent);
#endif

  PYLITH_METHOD_END;
} // integrateResidual

// ----------------------------------------------------------------------
// Compute matrix associated with operator.
void
pylith::feassemble::ElasticityExplicitTri3::integrateJacobian(topology::Jacobian* jacobian,
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
pylith::feassemble::ElasticityExplicitTri3::integrateJacobian(topology::Field* jacobian,
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
  assert(_quadrature->numBasis() == _numBasis);
  assert(_quadrature->spaceDim() == _spaceDim);
  assert(_quadrature->cellDim() == _cellDim);
  assert(_material->tensorSize() == _tensorSize);
  const int spaceDim = _spaceDim;
  const int cellDim = _cellDim;
  const int numCorners = _numCorners;
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

  // Setup visitors.
  topology::VecVisitorMesh jacobianVisitor(*jacobian, "displacement");
  // Don't optimize closure since we compute the Jacobian only once.

  _material->createPropsAndVarsVisitors();

  scalar_array coordsCell(numCorners*spaceDim);
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
    const PylithScalar area = _area(coordsCell);assert(area > 0.0);

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

    // Compute Jacobian for inertial terms
    const scalar_array& density = _material->calcDensity();
    _cellVector = density[0] * area / (3.0 * dt2);
    
#if defined(DETAILED_EVENT_LOGGING)
    PetscLogFlops(3);
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
  PetscLogFlops(numCells*3);
  _logger->eventEnd(computeEvent);
#endif

  _needNewJacobian = false;
  _material->resetNeedNewJacobian();

  PYLITH_METHOD_END;
} // integrateJacobian

// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::feassemble::ElasticityExplicitTri3::verifyConfiguration(const topology::Mesh& mesh) const
{ // verifyConfiguration
  PYLITH_METHOD_BEGIN;

  IntegratorElasticity::verifyConfiguration(mesh);

  assert(_quadrature);
  assert(_material);
  if (_spaceDim != _quadrature->spaceDim() || _cellDim != _quadrature->cellDim() || _numBasis != _quadrature->numBasis() ||  _numQuadPts != _quadrature->numQuadPts()) {
    std::ostringstream msg;
    msg << "User specified quadrature settings material '" << _material->label() << "' do not match ElasticityExplicitTri3 hardwired quadrature settings.\n"
	<< "  Space dim: " << _spaceDim << " (code), " << _quadrature->spaceDim() << " (user)\n"
	<< "  Cell dim: " << _cellDim << " (code), " << _quadrature->cellDim() << " (user)\n"
	<< "  # basis fns: " << _numBasis << " (code), " << _quadrature->numBasis() << " (user)\n"
	<< "  # quad points: " << _numQuadPts << " (code), " << _quadrature->numQuadPts() << " (user)";
    throw std::runtime_error(msg.str());
  } // if
} // verifyConfiguration

// ----------------------------------------------------------------------
// Compute area of triangular cell.
PylithScalar
pylith::feassemble::ElasticityExplicitTri3::_area(const scalar_array& coordinatesCell) const
{ // __area
  assert(6 == coordinatesCell.size());

  const PylithScalar x0 = coordinatesCell[0];
  const PylithScalar y0 = coordinatesCell[1];

  const PylithScalar x1 = coordinatesCell[2];
  const PylithScalar y1 = coordinatesCell[3];

  const PylithScalar x2 = coordinatesCell[4];
  const PylithScalar y2 = coordinatesCell[5];

  const PylithScalar area = 0.5*((x1-x0)*(y2-y0) - (x2-x0)*(y1-y0));
  PetscLogFlops(8);

  return area;  
} // _area


// End of file 
