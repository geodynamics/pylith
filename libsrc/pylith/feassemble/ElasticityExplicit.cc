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

#include "ElasticityExplicit.hh" // implementation of class methods

#include "Quadrature.hh" // USES Quadrature
#include "CellGeometry.hh" // USES CellGeometry

#include "pylith/materials/ElasticMaterial.hh" // USES ElasticMaterial
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields
#include "pylith/topology/Jacobian.hh" // USES Jacobian

#include "pylith/utils/array.hh" // USES scalar_array
#include "pylith/utils/macrodefs.h" // USES CALL_MEMBER_FN
#include "pylith/utils/lapack.h" // USES LAPACKdgesvd

#include "petscmat.h" // USES PetscMat
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimendional

#include "pylith/utils/petscerror.h" // USES CHECK_PETSC_ERROR
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

//#define PRECOMPUTE_GEOMETRY
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
  IntegratorElasticity::deallocate();
} // deallocate
  
// ----------------------------------------------------------------------
// Set time step for advancing from time t to time t+dt.
void
pylith::feassemble::ElasticityExplicit::timeStep(const PylithScalar dt)
{ // timeStep
  if (_dt != -1.0)
    _dtm1 = _dt;
  else
    _dtm1 = dt;
  _dt = dt;
  assert(_dt == _dtm1); // For now, don't allow variable time step
  if (0 != _material)
    _material->timeStep(_dt);
} // timeStep

// ----------------------------------------------------------------------
// Get stable time step for advancing from time t to time t+dt.
PylithScalar
pylith::feassemble::ElasticityExplicit::stableTimeStep(const topology::Mesh& mesh) const
{ // stableTimeStep
  assert(_material);
  return _material->stableTimeStepExplicit(mesh, _quadrature);
} // stableTimeStep

// ----------------------------------------------------------------------
// Set normalized viscosity for numerical damping.
void
pylith::feassemble::ElasticityExplicit::normViscosity(const PylithScalar viscosity)
{ // normViscosity
  if (viscosity < 0.0) {
    std::ostringstream msg;
    msg << "Normalized viscosity (" << viscosity << ") must be nonnegative.";
    throw std::runtime_error(msg.str());
  } // if

  _normViscosity = viscosity;
} // normViscosity

// ----------------------------------------------------------------------
// Integrate constributions to residual term (r) for operator.
void
pylith::feassemble::ElasticityExplicit::integrateResidual(const topology::Field<topology::Mesh>& residual,
							  const PylithScalar t,
							  topology::SolutionFields* const fields)
{ // integrateResidual
  /// Member prototype for _elasticityResidualXD()
  typedef void (pylith::feassemble::ElasticityExplicit::*elasticityResidual_fn_type)
    (const scalar_array&);

  assert(0 != _quadrature);
  assert(0 != _material);
  assert(0 != _logger);
  assert(0 != fields);

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
      &pylith::feassemble::ElasticityExplicit::_elasticityResidual1D;
    calcTotalStrainFn =
      &pylith::feassemble::IntegratorElasticity::_calcTotalStrain1D;
  } else if (2 == cellDim) {
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
  DM              dmMesh = fields->mesh().dmMesh();
  IS              cellIS;
  const PetscInt *cells;
  PetscInt        numCells;
  PetscErrorCode  err;

  assert(dmMesh);
  err = DMPlexGetStratumIS(dmMesh, "material-id", _material->id(), &cellIS);CHECK_PETSC_ERROR(err);
  err = ISGetSize(cellIS, &numCells);CHECK_PETSC_ERROR(err);
  err = ISGetIndices(cellIS, &cells);CHECK_PETSC_ERROR(err);

  // Get sections
  topology::Field<topology::Mesh>& acc = fields->get("acceleration(t)");
  PetscSection accSection = acc.petscSection();
  Vec          accVec     = acc.localVector();
  assert(accSection);assert(accVec);

  topology::Field<topology::Mesh>& vel = fields->get("velocity(t)");
  PetscSection velSection = vel.petscSection();
  Vec          velVec     = vel.localVector();
  assert(velSection);assert(velVec);
  
  scalar_array dispAdjCell(numBasis*spaceDim);
  topology::Field<topology::Mesh>& dispT = fields->get("disp(t)");
  PetscSection dispTSection = dispT.petscSection();
  Vec          dispTVec     = dispT.localVector();
  assert(dispTSection);assert(dispTVec);
  
  PetscSection residualSection = residual.petscSection();
  Vec          residualVec     = residual.localVector();

  scalar_array coordinatesCell(numBasis*spaceDim);
  PetscSection coordSection;
  Vec          coordVec;
  err = DMPlexGetCoordinateSection(dmMesh, &coordSection);CHECK_PETSC_ERROR(err);
  err = DMGetCoordinatesLocal(dmMesh, &coordVec);CHECK_PETSC_ERROR(err);
  assert(coordSection);assert(coordVec);

  assert(0 != _normalizer);
  const PylithScalar lengthScale = _normalizer->lengthScale();
  const PylithScalar gravityScale =
    _normalizer->pressureScale() / (_normalizer->lengthScale() *
            _normalizer->densityScale());

  const PylithScalar dt = _dt;
  assert(_normViscosity >= 0.0);
  assert(dt > 0);
  const PylithScalar viscosity = dt*_normViscosity;

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
#if defined(PRECOMPUTE_GEOMETRY)
    _quadrature->retrieveGeometry(cell);
#else
    PetscScalar *coords;
    PetscInt     coordsSize;
    err = DMPlexVecGetClosure(dmMesh, coordSection, coordVec, cell, &coordsSize, &coords);CHECK_PETSC_ERROR(err);
    for(PetscInt i = 0; i < coordsSize; ++i) {coordinatesCell[i] = coords[i];}
    _quadrature->computeGeometry(coordinatesCell, cell);
    err = DMPlexVecRestoreClosure(dmMesh, coordSection, coordVec, cell, &coordsSize, &coords);CHECK_PETSC_ERROR(err);
#endif
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
    PetscScalar *accArray, *velArray, *dispTArray;
    PetscInt     accSize,   velSize,   dispTSize;
    err = DMPlexVecGetClosure(dmMesh, accSection,   accVec,   cell, &accSize,   &accArray);CHECK_PETSC_ERROR(err);
    err = DMPlexVecGetClosure(dmMesh, velSection,   velVec,   cell, &velSize,   &velArray);CHECK_PETSC_ERROR(err);
    err = DMPlexVecGetClosure(dmMesh, dispTSection, dispTVec, cell, &dispTSize, &dispTArray);CHECK_PETSC_ERROR(err);
    assert(velSize   == accSize);
    assert(dispTSize == accSize);

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
    if (0 != _gravityField) {
      const spatialdata::geocoords::CoordSys* cs = fields->mesh().coordsys();
      assert(0 != cs);

      // Get density at quadrature points for this cell
      const scalar_array& density = _material->calcDensity();

      quadPtsGlobal = quadPtsNondim;
      _normalizer->dimensionalize(&quadPtsGlobal[0], quadPtsGlobal.size(),
          lengthScale);

      // Compute action for element body forces
      spatialdata::spatialdb::SpatialDB* db = _gravityField;
      for (int iQuad = 0; iQuad < numQuadPts; ++iQuad) {
        const int err = db->query(&gravVec[0], gravVec.size(),
            &quadPtsGlobal[0], spaceDim, cs);
        if (err)
          throw std::runtime_error("Unable to get gravity vector for point.");
        _normalizer->nondimensionalize(&gravVec[0], gravVec.size(),
            gravityScale);
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
      for (int jBasis = 0; jBasis < numBasis; ++jBasis)
	valJ += basis[iQ + jBasis];
      valJ *= wt;
      for (int iBasis = 0; iBasis < numBasis; ++iBasis)
	valuesIJ[iBasis] += basis[iQ + iBasis] * valJ;
    } // for
    for (int iBasis = 0; iBasis < numBasis; ++iBasis)
      for (int iDim = 0; iDim < spaceDim; ++iDim)
	_cellVector[iBasis*spaceDim+iDim] -= valuesIJ[iBasis] *
	  accArray[iBasis*spaceDim+iDim];

#if defined(DETAILED_EVENT_LOGGING)
    PetscLogFlops(numQuadPts*(4+numBasis*3));
    _logger->eventEnd(computeEvent);
    _logger->eventBegin(stressEvent);
#endif

    // Numerical damping. Compute displacements adjusted by velocity
    // times normalized viscosity.
    for(PetscInt i = 0; i < dispTSize; ++i) {dispAdjCell[i] = dispTArray[i] + viscosity * velArray[i];}
    err = DMPlexVecRestoreClosure(dmMesh, accSection,   accVec,   cell, &accSize,   &accArray);CHECK_PETSC_ERROR(err);
    err = DMPlexVecRestoreClosure(dmMesh, velSection,   velVec,   cell, &velSize,   &velArray);CHECK_PETSC_ERROR(err);
    err = DMPlexVecRestoreClosure(dmMesh, dispTSection, dispTVec, cell, &dispTSize, &dispTArray);CHECK_PETSC_ERROR(err);

    // Compute B(transpose) * sigma, first computing strains
    calcTotalStrainFn(&strainCell, basisDeriv, dispAdjCell, numBasis, numQuadPts);

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
    err = DMPlexVecSetClosure(dmMesh, residualSection, residualVec, cell, &_cellVector[0], ADD_VALUES);CHECK_PETSC_ERROR(err);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(updateEvent);
#endif
  } // for
  err = ISRestoreIndices(cellIS, &cells);CHECK_PETSC_ERROR(err);
  err = ISDestroy(&cellIS);CHECK_PETSC_ERROR(err);

#if !defined(DETAILED_EVENT_LOGGING)
  PetscLogFlops(numCells*numQuadPts*(4+numBasis*3));
  _logger->eventEnd(computeEvent);
#endif
} // integrateResidualLumped

// ----------------------------------------------------------------------
// Compute matrix associated with operator.
void
pylith::feassemble::ElasticityExplicit::integrateJacobian(
					topology::Jacobian* jacobian,
					const PylithScalar t,
					topology::SolutionFields* fields)
{ // integrateJacobian
  throw std::logic_error("ElasticityExplicit::integrateJacobian() not implemented. Use integrateJacobian(lumped) instead.");
} // integrateJacobian

// ----------------------------------------------------------------------
// Compute matrix associated with operator.
void
pylith::feassemble::ElasticityExplicit::integrateJacobian(
			    topology::Field<topology::Mesh>* jacobian,
			    const PylithScalar t,
			    topology::SolutionFields* fields)
{ // integrateJacobian
  assert(0 != _quadrature);
  assert(0 != _material);
  assert(0 != jacobian);
  assert(0 != fields);

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
  if (cellDim != spaceDim)
    throw std::logic_error("Don't know how to integrate elasticity " \
			   "contribution to Jacobian matrix for cells with " \
			   "different dimensions than the spatial dimension.");

  // Get cell information
  DM              dmMesh = fields->mesh().dmMesh();
  IS              cellIS;
  const PetscInt *cells;
  PetscInt        numCells;
  PetscErrorCode  err;

  assert(dmMesh);
  err = DMPlexGetStratumIS(dmMesh, "material-id", _material->id(), &cellIS);CHECK_PETSC_ERROR(err);
  err = ISGetSize(cellIS, &numCells);CHECK_PETSC_ERROR(err);
  err = ISGetIndices(cellIS, &cells);CHECK_PETSC_ERROR(err);

  // Get parameters used in integration.
  const PylithScalar dt = _dt;
  const PylithScalar dt2 = dt*dt;
  assert(dt > 0);
  scalar_array valuesIJ(numBasis);

  // Get sections
  PetscSection jacSection = jacobian->petscSection();
  Vec          jacVec     = jacobian->localVector();

  scalar_array coordinatesCell(numBasis*spaceDim);
  PetscSection coordSection;
  Vec          coordVec;
  err = DMPlexGetCoordinateSection(dmMesh, &coordSection);CHECK_PETSC_ERROR(err);
  err = DMGetCoordinatesLocal(dmMesh, &coordVec);CHECK_PETSC_ERROR(err);
  assert(coordSection);assert(coordVec);

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
#if defined(PRECOMPUTE_GEOMETRY)
    _quadrature->retrieveGeometry(cell);
#else
    PetscScalar *coords;
    PetscInt     coordsSize;
    err = DMPlexVecGetClosure(dmMesh, coordSection, coordVec, cell, &coordsSize, &coords);CHECK_PETSC_ERROR(err);
    for(PetscInt i = 0; i < coordsSize; ++i) {coordinatesCell[i] = coords[i];}
    _quadrature->computeGeometry(coordinatesCell, cell);
    err = DMPlexVecRestoreClosure(dmMesh, coordSection, coordVec, cell, &coordsSize, &coords);CHECK_PETSC_ERROR(err);
#endif

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
      const PylithScalar wt = quadWts[iQuad] * jacobianDet[iQuad] * density[iQuad]
	/ dt2;
      const int iQ = iQuad * numBasis;
      PylithScalar valJ = 0.0;
      for (int jBasis = 0; jBasis < numBasis; ++jBasis)
        valJ += basis[iQ + jBasis];
      valJ *= wt;
      for (int iBasis = 0; iBasis < numBasis; ++iBasis) {
        valuesIJ[iBasis] += basis[iQ + iBasis] * valJ;
      } // for
    } // for
    for (int iBasis = 0; iBasis < numBasis; ++iBasis)
      for (int iDim = 0; iDim < spaceDim; ++iDim)
        _cellVector[iBasis*spaceDim+iDim] += valuesIJ[iBasis];
    
#if defined(DETAILED_EVENT_LOGGING)
    PetscLogFlops(numQuadPts*(4 + numBasis*3) + numBasis*spaceDim);
    _logger->eventEnd(computeEvent);
    _logger->eventBegin(updateEvent);
#endif
    
    // Assemble cell contribution into lumped matrix.
    err = DMPlexVecSetClosure(dmMesh, jacSection, jacVec, cell, &_cellVector[0], ADD_VALUES);CHECK_PETSC_ERROR(err);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(updateEvent);
#endif
  } // for
  err = ISRestoreIndices(cellIS, &cells);CHECK_PETSC_ERROR(err);
  err = ISDestroy(&cellIS);CHECK_PETSC_ERROR(err);

#if !defined(DETAILED_EVENT_LOGGING)
  PetscLogFlops(numCells*(numQuadPts*(4 + numBasis*3) + numBasis*spaceDim));
  _logger->eventEnd(computeEvent);
#endif

  _needNewJacobian = false;
  _material->resetNeedNewJacobian();
} // integrateJacobian


// End of file 
