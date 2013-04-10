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
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "ElasticityExplicitTet4.hh" // implementation of class methods

#include "Quadrature.hh" // USES Quadrature

#include "pylith/materials/ElasticMaterial.hh" // USES ElasticMaterial
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields
#include "pylith/topology/Jacobian.hh" // USES Jacobian
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor

#include "pylith/utils/array.hh" // USES scalar_array
#include "pylith/utils/macrodefs.h" // USES CALL_MEMBER_FN

#include "petscmat.h" // USES PetscMat
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimendional

#include "pylith/utils/petscerror.h" // USES CHECK_PETSC_ERROR
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error

//#define DETAILED_EVENT_LOGGING

// ----------------------------------------------------------------------
const int pylith::feassemble::ElasticityExplicitTet4::_spaceDim = 3;
const int pylith::feassemble::ElasticityExplicitTet4::_cellDim = 3;
const int pylith::feassemble::ElasticityExplicitTet4::_tensorSize = 6;
const int pylith::feassemble::ElasticityExplicitTet4::_numBasis = 4;
const int pylith::feassemble::ElasticityExplicitTet4::_numQuadPts = 1;

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::ElasticityExplicitTet4::ElasticityExplicitTet4(void) :
  _dtm1(-1.0),
  _normViscosity(0.1)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::ElasticityExplicitTet4::~ElasticityExplicitTet4(void)
{ // destructor
  deallocate();
} // destructor
  
// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::feassemble::ElasticityExplicitTet4::deallocate(void)
{ // deallocate
  PYLITH_METHOD_BEGIN;

  IntegratorElasticity::deallocate();

  PYLITH_METHOD_END;
} // deallocate
  
// ----------------------------------------------------------------------
// Set time step for advancing from time t to time t+dt.
void
pylith::feassemble::ElasticityExplicitTet4::timeStep(const PylithScalar dt)
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
pylith::feassemble::ElasticityExplicitTet4::stableTimeStep(const topology::Mesh& mesh) const
{ // stableTimeStep
  PYLITH_METHOD_BEGIN;
  
  assert(_material);
  PYLITH_METHOD_RETURN(_material->stableTimeStepExplicit(mesh, _quadrature));
} // stableTimeStep

// ----------------------------------------------------------------------
// Set normalized viscosity for numerical damping.
void
pylith::feassemble::ElasticityExplicitTet4::normViscosity(const PylithScalar viscosity)
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
pylith::feassemble::ElasticityExplicitTet4::integrateResidual(
        const topology::Field<topology::Mesh>& residual,
        const PylithScalar t,
        topology::SolutionFields* const fields)
{ // integrateResidual
  PYLITH_METHOD_BEGIN;
  
  /// Member prototype for _elasticityResidualXD()
  typedef void (pylith::feassemble::ElasticityExplicitTet4::*elasticityResidual_fn_type)
    (const scalar_array&);

  assert(_quadrature);
  assert(_material);
  assert(_logger);
  assert(fields);

  const int setupEvent = _logger->eventId("ElIR setup");
  const int geometryEvent = _logger->eventId("ElIR geometry");
  const int computeEvent = _logger->eventId("ElIR compute");
  const int restrictEvent = _logger->eventId("ElIR restrict");
  const int stateVarsEvent = _logger->eventId("ElIR stateVars");
  const int stressEvent = _logger->eventId("ElIR stress");
  const int updateEvent = _logger->eventId("ElIR update");

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
  topology::StratumIS materialIS(dmMesh, "material-id", _material->id());
  const PetscInt* cells = materialIS.points();
  const PetscInt numCells = materialIS.size();

  // Setup field visitors.
  topology::VecVisitorMesh accVisitor(fields->get("acceleration(t)"));
  PetscScalar* accCell = NULL;
  PetscInt accSize = 0;

  topology::VecVisitorMesh velVisitor(fields->get("velocity(t)"));
  PetscScalar* velCell = NULL;
  PetscInt velSize = 0;

  scalar_array dispAdjCell(numBasis*spaceDim);
  topology::VecVisitorMesh dispVisitor(fields->get("disp(t)"));
  PetscScalar* dispCell = NULL;
  PetscInt dispSize = 0;
  
  topology::VecVisitorMesh residualVisitor(residual);

  topology::CoordsVisitor coordsVisitor(dmMesh);
  PetscScalar *coordsCell = NULL;
  PetscInt coordsSize = 0;

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
    accVisitor.getClosure(&accCell, &accSize, cell);
    velVisitor.getClosure(&velCell, &velSize, cell);
    dispVisitor.getClosure(&dispCell, &dispSize, cell);
    assert(velSize == accSize);
    assert(dispSize == accSize);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(restrictEvent);
    _logger->eventBegin(geometryEvent);
#endif

    // Compute geometry information for current cell
    coordsVisitor.getClosure(&coordsCell, &coordsSize, cell);
    const PylithScalar volume = _volume(coordsCell, coordsSize);assert(volume > 0.0);

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
      const PylithScalar wtVertex = density[0] * volume / 4.0;
      for (int iBasis=0; iBasis < numBasis; ++iBasis) {
        for (int iDim=0; iDim < spaceDim; ++iDim) {
            _cellVector[iBasis * spaceDim + iDim] += wtVertex * gravVec[iDim];
	} // for
      } // for
      PetscLogFlops(numBasis*spaceDim*2 + numBasis*spaceDim*2);
    } // if

    // Compute action for inertial terms
    const PylithScalar wtVertex = density[0] * volume / 4.0;
    assert(cellVectorSize == dispSize);
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
    accVisitor.restoreClosure(&accCell, &accSize, cell);
    velVisitor.restoreClosure(&velCell, &velSize, cell);
    dispVisitor.restoreClosure(&dispCell, &dispSize, cell);

    // Compute B(transpose) * sigma, first computing strains
    const PylithScalar x0 = coordsCell[0];
    const PylithScalar y0 = coordsCell[1];
    const PylithScalar z0 = coordsCell[2];

    const PylithScalar x1 = coordsCell[3];
    const PylithScalar y1 = coordsCell[4];
    const PylithScalar z1 = coordsCell[5];

    const PylithScalar x2 = coordsCell[6];
    const PylithScalar y2 = coordsCell[7];
    const PylithScalar z2 = coordsCell[8];

    const PylithScalar x3 = coordsCell[9];
    const PylithScalar y3 = coordsCell[10];
    const PylithScalar z3 = coordsCell[11];

    const PylithScalar scaleB = 6.0 * volume;
    const PylithScalar b1 = (y1*(z3-z2)-y2*z3+y3*z2-(y3-y2)*z1) / scaleB;
    const PylithScalar c1 = (-x1*(z3-z2)+x2*z3-x3*z2-(x2-x3)*z1) / scaleB;
    const PylithScalar d1 = (-x2*y3-x1*(y2-y3)+x3*y2+(x2-x3)*y1) / scaleB;

    const PylithScalar b2 = (-y0*z3-y2*(z0-z3)+(y0-y3)*z2+y3*z0) / scaleB;
    const PylithScalar c2 = (x0*z3+x2*(z0-z3)+(x3-x0)*z2-x3*z0) / scaleB;
    const PylithScalar d2 = (x2*(y3-y0)-x0*y3-(x3-x0)*y2+x3*y0) / scaleB;

    const PylithScalar b3 = (-(y1-y0)*z3+y3*(z1-z0)-y0*z1+y1*z0) / scaleB;
    const PylithScalar c3 = (-(x0-x1)*z3-x3*(z1-z0)+x0*z1-x1*z0) / scaleB;
    const PylithScalar d3 = ((x0-x1)*y3-x0*y1-x3*(y0-y1)+x1*y0) / scaleB;

    const PylithScalar b4 = (-y0*(z2-z1)+y1*z2-y2*z1+(y2-y1)*z0) / scaleB;
    const PylithScalar c4 = (x0*(z2-z1)-x1*z2+x2*z1+(x1-x2)*z0) / scaleB;
    const PylithScalar d4 = (x1*y2+x0*(y1-y2)-x2*y1-(x1-x2)*y0) / scaleB;

    assert(strainCell.size() == 6);
    strainCell[0] = 
      b1 * dispAdjCell[0] + b2 * dispAdjCell[3] + 
      b3 * dispAdjCell[6] + b4 * dispAdjCell[9];
    strainCell[1] = 
      c3 * dispAdjCell[7] + c2 * dispAdjCell[4] + 
      c4 * dispAdjCell[10] + c1 * dispAdjCell[1];
    strainCell[2] = 
      d3 * dispAdjCell[8] + d2 * dispAdjCell[5] + 
      d1 * dispAdjCell[2] + d4 * dispAdjCell[11];
    strainCell[3] = 
      (c4 * dispAdjCell[9] + b3 * dispAdjCell[7] + 
       c3 * dispAdjCell[6] + b2 * dispAdjCell[4] + 
       c2 * dispAdjCell[3] + b4 * dispAdjCell[10] + 
       b1 * dispAdjCell[1] + c1 * dispAdjCell[0]) / 2.0;
    strainCell[4] = 
      (c3 * dispAdjCell[8] + d3 * dispAdjCell[7] + 
       c2 * dispAdjCell[5] + d2 * dispAdjCell[4] +
       c1 * dispAdjCell[2] + c4 * dispAdjCell[11] +
       d4 * dispAdjCell[10] + d1 * dispAdjCell[1]) / 2.0;
    strainCell[5] = 
      (d4 * dispAdjCell[9] + b3 * dispAdjCell[8] + 
       d3 * dispAdjCell[6] + b2 * dispAdjCell[5] +
       d2 * dispAdjCell[3] + b1 * dispAdjCell[2] + 
       b4 * dispAdjCell[11] + d1 * dispAdjCell[0]) / 2.0;

    const scalar_array& stressCell = _material->calcStress(strainCell, false);

#if defined(DETAILED_EVENT_LOGGING)
    PetscLogFlops(196);
    _logger->eventEnd(stressEvent);
    _logger->eventBegin(computeEvent);
#endif

    assert(_cellVector.size() == 12);
    assert(stressCell.size() == 6);
    _cellVector[ 0] -= (d1*stressCell[5]+c1*stressCell[3]+b1*stressCell[0]) * volume;
    _cellVector[ 1] -= (d1*stressCell[4]+b1*stressCell[3]+c1*stressCell[1]) * volume;
    _cellVector[ 2] -= (b1*stressCell[5]+c1*stressCell[4]+d1*stressCell[2]) * volume;
    _cellVector[ 3] -= (d2*stressCell[5]+c2*stressCell[3]+b2*stressCell[0]) * volume;
    _cellVector[ 4] -= (d2*stressCell[4]+b2*stressCell[3]+c2*stressCell[1]) * volume;
    _cellVector[ 5] -= (b2*stressCell[5]+c2*stressCell[4]+d2*stressCell[2]) * volume;
    _cellVector[ 6] -= (d3*stressCell[5]+c3*stressCell[3]+b3*stressCell[0]) * volume;
    _cellVector[ 7] -= (d3*stressCell[4]+b3*stressCell[3]+c3*stressCell[1]) * volume;
    _cellVector[ 8] -= (b3*stressCell[5]+c3*stressCell[4]+d3*stressCell[2]) * volume;
    _cellVector[ 9] -= (d4*stressCell[5]+c4*stressCell[3]+b4*stressCell[0]) * volume;
    _cellVector[10] -= (d4*stressCell[4]+b4*stressCell[3]+c4*stressCell[1]) * volume;
    _cellVector[11] -= (b4*stressCell[5]+c4*stressCell[4]+d4*stressCell[2]) * volume;

#if defined(DETAILED_EVENT_LOGGING)
    PetscLogFlops(84);
    _logger->eventEnd(computeEvent);
    _logger->eventBegin(updateEvent);
#endif

    coordsVisitor.restoreClosure(&coordsCell, &coordsSize, cell);

    // Assemble cell contribution into field
    residualVisitor.setClosure(&_cellVector[0], _cellVector.size(), cell, ADD_VALUES);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(updateEvent);
#endif
  } // for

#if !defined(DETAILED_EVENT_LOGGING)
  PetscLogFlops(numCells*(2 + numBasis*spaceDim*2 + 196+84));
  _logger->eventEnd(computeEvent);
#endif

  PYLITH_METHOD_END;
} // integrateResidual

// ----------------------------------------------------------------------
// Compute matrix associated with operator.
void
pylith::feassemble::ElasticityExplicitTet4::integrateJacobian(topology::Jacobian* jacobian,
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
pylith::feassemble::ElasticityExplicitTet4::integrateJacobian(topology::Field<topology::Mesh>* jacobian,
							      const PylithScalar t,
							      topology::SolutionFields* fields)
{ // integrateJacobian
  PYLITH_METHOD_BEGIN;
  
  assert(_quadrature);
  assert(_material);
  assert(jacobian);
  assert(fields);

  const int setupEvent = _logger->eventId("ElIJ setup");
  const int geometryEvent = _logger->eventId("ElIJ geometry");
  const int computeEvent = _logger->eventId("ElIJ compute");
  const int restrictEvent = _logger->eventId("ElIJ restrict");
  const int stateVarsEvent = _logger->eventId("ElIJ stateVars");
  const int updateEvent = _logger->eventId("ElIJ update");

  _logger->eventBegin(setupEvent);

  // Get cell geometry information that doesn't depend on cell
  assert(_quadrature->numBasis() == _numBasis);
  assert(_quadrature->spaceDim() == _spaceDim);
  assert(_quadrature->cellDim() == _cellDim);
  assert(_material->tensorSize() == _tensorSize);
  const int spaceDim = _spaceDim;
  const int cellDim = _cellDim;
  const int numBasis = _numBasis;
  if (cellDim != spaceDim)
    throw std::logic_error("Don't know how to integrate elasticity " \
			   "contribution to Jacobian matrix for cells with " \
			   "different dimensions than the spatial dimension.");

  // Get cell information
  PetscDM dmMesh = fields->mesh().dmMesh();assert(dmMesh);
  topology::StratumIS materialIS(dmMesh, "material-id", _material->id());
  const PetscInt* cells = materialIS.points();
  const PetscInt numCells = materialIS.size();

  // Get parameters used in integration.
  const PylithScalar dt = _dt;
  const PylithScalar dt2 = dt*dt;
  assert(dt > 0);

  // Setup visitors.
  topology::VecVisitorMesh jacobianVisitor(*jacobian);
  PetscScalar* jacobianCell = NULL;
  PetscInt jacobianSize = 0;

  topology::CoordsVisitor coordsVisitor(dmMesh);
  PetscScalar* coordsCell = NULL;
  PetscInt coordsSize = 0;

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
    coordsVisitor.getClosure(&coordsCell, &coordsSize, cell);
    const PylithScalar volume = _volume(coordsCell, coordsSize);assert(volume > 0.0);
    coordsVisitor.restoreClosure(&coordsCell, &coordsSize, cell);

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
    _cellVector = density[0] * volume / (4.0 * dt2);
    
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
pylith::feassemble::ElasticityExplicitTet4::verifyConfiguration(const topology::Mesh& mesh) const
{ // verifyConfiguration
  PYLITH_METHOD_BEGIN;
  
  IntegratorElasticity::verifyConfiguration(mesh);

  assert(_quadrature);
  assert(_material);
  if (_spaceDim != _quadrature->spaceDim() || _cellDim != _quadrature->cellDim() || _numBasis != _quadrature->numBasis() ||  _numQuadPts != _quadrature->numQuadPts()) {
    std::ostringstream msg;
    msg << "User specified quadrature settings material '" << _material->label() << "' do not match ElasticityExplicitTet4 hardwired quadrature settings.\n"
	<< "  Space dim: " << _spaceDim << " (code), " << _quadrature->spaceDim() << " (user)\n"
	<< "  Cell dim: " << _cellDim << " (code), " << _quadrature->cellDim() << " (user)\n"
	<< "  # basis fns: " << _numBasis << " (code), " << _quadrature->numBasis() << " (user)\n"
	<< "  # quad points: " << _numQuadPts << " (code), " << _quadrature->numQuadPts() << " (user)";
    throw std::runtime_error(msg.str());
  } // if

  PYLITH_METHOD_END;
} // verifyConfiguration

// ----------------------------------------------------------------------
// Compute volume of tetrahedral cell.
PylithScalar
pylith::feassemble::ElasticityExplicitTet4::_volume(const PylithScalar coordinatesCell[],
						    const int coordinatesSize) const
{ // __volume
  assert(12 == coordinatesSize);

  const PylithScalar x0 = coordinatesCell[0];
  const PylithScalar y0 = coordinatesCell[1];
  const PylithScalar z0 = coordinatesCell[2];

  const PylithScalar x1 = coordinatesCell[3];
  const PylithScalar y1 = coordinatesCell[4];
  const PylithScalar z1 = coordinatesCell[5];

  const PylithScalar x2 = coordinatesCell[6];
  const PylithScalar y2 = coordinatesCell[7];
  const PylithScalar z2 = coordinatesCell[8];

  const PylithScalar x3 = coordinatesCell[9];
  const PylithScalar y3 = coordinatesCell[10];
  const PylithScalar z3 = coordinatesCell[11];

  const PylithScalar det = 
    x1*(y2*z3-y3*z2)-y1*(x2*z3-x3*z2)+(x2*y3-x3*y2)*z1 - 
    x0*((y2*z3-y3*z2)-y1*(z3-z2)+(y3-y2)*z1) +
    y0*((x2*z3-x3*z2)-x1*(z3-z2)+(x3-x2)*z1) -
    z0*((x2*y3-x3*y2)-x1*(y3-y2)+(x3-x2)*y1);
    
  const PylithScalar volume = det / 6.0;
  PetscLogFlops(48);

  return volume;
} // _volume


// End of file 
