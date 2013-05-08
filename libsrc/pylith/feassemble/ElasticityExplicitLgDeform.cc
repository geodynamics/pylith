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
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "ElasticityExplicitLgDeform.hh" // implementation of class methods

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
#include "pylith/utils/macrodefs.h" // USES CALL_MEMBER_FN

#include "petscmat.h" // USES PetscMat
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimendional

#include "pylith/utils/error.h" // USES PYLITH_CHECK_ERROR
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::ElasticityExplicitLgDeform::ElasticityExplicitLgDeform(void) :
  _dtm1(-1.0),
  _normViscosity(0.1)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::ElasticityExplicitLgDeform::~ElasticityExplicitLgDeform(void)
{ // destructor
  deallocate();
} // destructor
  
// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::feassemble::ElasticityExplicitLgDeform::deallocate(void)
{ // deallocate
  PYLITH_METHOD_BEGIN;

  IntegratorElasticityLgDeform::deallocate();

  PYLITH_METHOD_END;
} // deallocate
  
// ----------------------------------------------------------------------
// Set time step for advancing from time t to time t+dt.
void
pylith::feassemble::ElasticityExplicitLgDeform::timeStep(const PylithScalar dt)
{ // timeStep
  PYLITH_METHOD_BEGIN;
  
  if (_dt != -1.0)
    _dtm1 = _dt;
  else
    _dtm1 = dt;
  _dt = dt;
  assert(_dt == _dtm1); // For now, don't allow variable time step
  if (0 != _material)
    _material->timeStep(_dt);

  PYLITH_METHOD_END;
} // timeStep

// ----------------------------------------------------------------------
// Get stable time step for advancing from time t to time t+dt.
PylithScalar
pylith::feassemble::ElasticityExplicitLgDeform::stableTimeStep(const topology::Mesh& mesh) const
{ // stableTimeStep
  PYLITH_METHOD_BEGIN;
  
  assert(_material);
  PYLITH_METHOD_RETURN(_material->stableTimeStepExplicit(mesh, _quadrature));
} // stableTimeStep

// ----------------------------------------------------------------------
// Set normalized viscosity for numerical damping.
void
pylith::feassemble::ElasticityExplicitLgDeform::normViscosity(const PylithScalar viscosity)
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
pylith::feassemble::ElasticityExplicitLgDeform::integrateResidual(const topology::Field<topology::Mesh>& residual,
								  const PylithScalar t,
								  topology::SolutionFields* const fields)
{ // integrateResidualLumped
  PYLITH_METHOD_BEGIN;
  
  /// Member prototype for _elasticityResidualXD()
  typedef void (pylith::feassemble::ElasticityExplicitLgDeform::*elasticityResidual_fn_type)
    (const scalar_array&, const scalar_array&);

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
  const int numQuadPts = _quadrature->numQuadPts();
  const scalar_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == numQuadPts);
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
  if (1 == cellDim) {
    elasticityResidualFn = 
      &pylith::feassemble::ElasticityExplicitLgDeform::_elasticityResidual1D;
    calcTotalStrainFn = 
      &pylith::feassemble::IntegratorElasticityLgDeform::_calcTotalStrain1D;
  } else if (2 == cellDim) {
    elasticityResidualFn = 
      &pylith::feassemble::ElasticityExplicitLgDeform::_elasticityResidual2D;
    calcTotalStrainFn = 
      &pylith::feassemble::IntegratorElasticityLgDeform::_calcTotalStrain2D;
  } else if (3 == cellDim) {
    elasticityResidualFn = 
      &pylith::feassemble::ElasticityExplicitLgDeform::_elasticityResidual3D;
    calcTotalStrainFn = 
      &pylith::feassemble::IntegratorElasticityLgDeform::_calcTotalStrain3D;
  } else
    assert(0);

  // Allocate vectors for cell values.
  scalar_array deformCell(numQuadPts*spaceDim*spaceDim);
  scalar_array strainCell(numQuadPts*tensorSize);
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

  scalar_array coordinatesCell(numBasis*spaceDim);
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
  _logger->eventBegin(computeEvent);

  // Loop over cells
  for(PetscInt c = 0; c < numCells; ++c) {
    const PetscInt cell = cells[c];
    // Compute geometry information for current cell
    coordsVisitor.getClosure(&coordsCell, &coordsSize, cell);
    _quadrature->computeGeometry(coordsCell, coordsSize, cell);

    // Get state variables for cell.
    _material->retrievePropsAndVars(cell);

    // Reset element vector to zero
    _resetCellVector();

    // Restrict input fields to cell
    accVisitor.getClosure(&accCell, &accSize, cell);
    velVisitor.getClosure(&velCell, &velSize, cell);
    dispVisitor.getClosure(&dispCell, &dispSize, cell);
    assert(velSize == accSize);
    assert(dispSize == accSize);

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
      for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
	const int err = db->query(&gravVec[0], gravVec.size(), &quadPtsGlobal[0], spaceDim, cs);
	if (err) {
	  throw std::runtime_error("Unable to get gravity vector for point.");
	} // if
	_normalizer->nondimensionalize(&gravVec[0], gravVec.size(), gravityScale);
	const PylithScalar wt = quadWts[iQuad] * jacobianDet[iQuad] * density[iQuad];
	for (int iBasis=0, iQ=iQuad*numBasis; iBasis < numBasis; ++iBasis) {
	  const PylithScalar valI = wt*basis[iQ+iBasis];
	  for (int iDim=0; iDim < spaceDim; ++iDim) {
	    _cellVector[iBasis*spaceDim+iDim] += valI*gravVec[iDim];
	  } // for
	} // for
      } // for
      PetscLogFlops(numQuadPts*(2+numBasis*(1+2*spaceDim)));
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
    PetscLogFlops(numQuadPts*(4+numBasis*3));

    // Numerical damping. Compute displacements adjusted by velocity
    // times normalized viscosity.
    for (PetscInt i = 0; i < dispSize; ++i) {
      dispAdjCell[i] = dispCell[i] + viscosity * velCell[i];
    } // for
    accVisitor.restoreClosure(&accCell, &accSize, cell);
    velVisitor.restoreClosure(&velCell, &velSize, cell);
    dispVisitor.restoreClosure(&dispCell, &dispSize, cell);

    // Compute B(transpose) * sigma, first computing strains
    _calcDeformation(&deformCell, basisDeriv, coordsCell, &dispAdjCell[0], numBasis, numQuadPts, spaceDim);
    calcTotalStrainFn(&strainCell, deformCell, numQuadPts);
    const scalar_array& stressCell = _material->calcStress(strainCell, true);

    CALL_MEMBER_FN(*this, elasticityResidualFn)(stressCell, dispAdjCell);
    
    // Assemble cell contribution into field
    residualVisitor.setClosure(&_cellVector[0], _cellVector.size(), cell, ADD_VALUES);

    coordsVisitor.restoreClosure(&coordsCell, &coordsSize, cell);
  } // for

  _logger->eventEnd(computeEvent);

  PYLITH_METHOD_END;
} // integrateResidual

// ----------------------------------------------------------------------
// Compute matrix associated with operator.
void
pylith::feassemble::ElasticityExplicitLgDeform::integrateJacobian(topology::Jacobian* jacobian,
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
pylith::feassemble::ElasticityExplicitLgDeform::integrateJacobian(topology::Field<topology::Mesh>* jacobian,
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
  const int numQuadPts = _quadrature->numQuadPts();
  const scalar_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == numQuadPts);
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int cellDim = _quadrature->cellDim();
  const int tensorSize = _material->tensorSize();
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

  // Setup field visitors.
  scalar_array valuesIJ(numBasis);
  topology::VecVisitorMesh jacobianVisitor(*jacobian);

  topology::CoordsVisitor coordsVisitor(dmMesh);
  PetscScalar* coordsCell = NULL;
  PetscInt coordsSize = 0;

  _logger->eventEnd(setupEvent);
  _logger->eventBegin(computeEvent);

  // Loop over cells
  for(PetscInt c = 0; c < numCells; ++c) {
    const PetscInt cell = cells[c];
    // Compute geometry information for current cell
    coordsVisitor.getClosure(&coordsCell, &coordsSize, cell);
    _quadrature->computeGeometry(coordsCell, coordsSize, cell);
    coordsVisitor.restoreClosure(&coordsCell, &coordsSize, cell);

    // Get state variables for cell.
    _material->retrievePropsAndVars(cell);

    // Reset element matrix to zero
    _resetCellMatrix();

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
  } // for

  _needNewJacobian = false;
  _material->resetNeedNewJacobian();

  _logger->eventEnd(computeEvent);

  PYLITH_METHOD_END;
} // integrateJacobian


// End of file 
