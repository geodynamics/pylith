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

#include "pylith/utils/array.hh" // USES scalar_array
#include "pylith/utils/macrodefs.h" // USES CALL_MEMBER_FN
#include "pylith/utils/lapack.h" // USES LAPACKdgesvd

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
  _basisDerivArray.resize(_numQuadPts*_numBasis*_spaceDim);
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
  IntegratorElasticity::deallocate();
} // deallocate
  
// ----------------------------------------------------------------------
// Set time step for advancing from time t to time t+dt.
void
pylith::feassemble::ElasticityExplicitTet4::timeStep(const PylithScalar dt)
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
// Set normalized viscosity for numerical damping.
void
pylith::feassemble::ElasticityExplicitTet4::normViscosity(const PylithScalar viscosity)
{ // normViscosity
  if (viscosity < 0.0) {
    std::ostringstream msg;
    msg << "Normalized viscosity (" << viscosity << ") must be nonnegative.";
    throw std::runtime_error(msg.str());
  } // if

  _normViscosity = viscosity;
} // normViscosity

// ----------------------------------------------------------------------
// Set flag for setting constraints for total field solution or
// incremental field solution.
void
pylith::feassemble::ElasticityExplicitTet4::useSolnIncr(const bool flag)
{ // useSolnIncr
  _material->useElasticBehavior(false);
  if (!flag)
    throw std::logic_error("Non-incremental solution not supported for "
			   "explicit time integration of elasticity "
			   "equation.");
} // useSolnIncr

// ----------------------------------------------------------------------
// Integrate constributions to residual term (r) for operator.
void
pylith::feassemble::ElasticityExplicitTet4::integrateResidual(
			  const topology::Field<topology::Mesh>& residual,
			  const PylithScalar t,
			  topology::SolutionFields* const fields)
{ // integrateResidual
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
  assert(_quadrature->numQuadPts() == _numQuadPts);
  const scalar_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == _numQuadPts);
  assert(_quadrature->numBasis() == _numBasis);
  assert(_quadrature->spaceDim() == _spaceDim);
  assert(_quadrature->cellDim() == _cellDim);
  assert(_material->tensorSize() == _tensorSize);
  const int spaceDim = _spaceDim;
  const int cellDim = _cellDim;
  const int tensorSize = _tensorSize;
  const int numBasis = _numBasis;
  const int numQuadPts = _numQuadPts;
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
  DM              dmMesh = fields->mesh().dmMesh();
  IS              cellIS;
  const PetscInt *cells;
  PetscInt        numCells;
  PetscErrorCode  err;

  assert(dmMesh);
  err = DMComplexGetStratumIS(dmMesh, "material-id", _material->id(), &cellIS);CHECK_PETSC_ERROR(err);
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
  err = DMComplexGetCoordinateSection(dmMesh, &coordSection);CHECK_PETSC_ERROR(err);
  err = DMGetCoordinatesLocal(dmMesh, &coordVec);CHECK_PETSC_ERROR(err);
  assert(coordSection);assert(coordVec);

  assert(0 != _normalizer);
  const PylithScalar lengthScale = _normalizer->lengthScale();
  const PylithScalar gravityScale = 
    _normalizer->pressureScale() / (_normalizer->lengthScale() *
				    _normalizer->densityScale());

  const PylithScalar dt = _dt;
  assert(_normViscosity > 0.0);
  assert(dt > 0);
  const PylithScalar viscosity = dt*_normViscosity;

  _logger->eventEnd(setupEvent);
#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventBegin(computeEvent);
#endif

  // Loop over cells
  for(PetscInt c = 0; c < numCells; ++c) {
    const PetscInt cell = cells[c];
#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(geometryEvent);
#endif

    // Compute geometry information for current cell
    const PetscScalar *coords;
    PetscInt           coordsSize;
    err = DMComplexVecGetClosure(dmMesh, coordSection, coordVec, cell, &coordsSize, &coords);CHECK_PETSC_ERROR(err);
    for(PetscInt i = 0; i < coordsSize; ++i) {coordinatesCell[i] = coords[i];}
    const PylithScalar volume = _volume(coordinatesCell);
    assert(volume > 0.0);
    err = DMComplexVecRestoreClosure(dmMesh, coordSection, coordVec, cell, &coordsSize, &coords);CHECK_PETSC_ERROR(err);

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
    const PetscScalar *accArray, *velArray, *dispTArray;
    PetscInt           accSize,   velSize,   dispTSize;
    err = DMComplexVecGetClosure(dmMesh, accSection,   accVec,   cell, &accSize,   &accArray);CHECK_PETSC_ERROR(err);
    err = DMComplexVecGetClosure(dmMesh, velSection,   velVec,   cell, &velSize,   &velArray);CHECK_PETSC_ERROR(err);
    err = DMComplexVecGetClosure(dmMesh, dispTSection, dispTVec, cell, &dispTSize, &dispTArray);CHECK_PETSC_ERROR(err);
    assert(velSize   == accSize);
    assert(dispTSize == accSize);
    
#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(restrictEvent);
    _logger->eventBegin(computeEvent);
#endif

    const scalar_array& density = _material->calcDensity();
    assert(density.size() == 1);

    // Compute body force vector if gravity is being used.
    if (0 != _gravityField) {
      const spatialdata::geocoords::CoordSys* cs = fields->mesh().coordsys();
      assert(0 != cs);

      quadPtsGlobal = 0.0;
      for (int iBasis=0; iBasis < numBasis; ++iBasis)
        for (int iDim=0; iDim < spaceDim; ++iDim)
          quadPtsGlobal[iDim] += 
	    coordinatesCell[iBasis*spaceDim+iDim] / numBasis;
      _normalizer->dimensionalize(&quadPtsGlobal[0], quadPtsGlobal.size(),
          lengthScale);

      // Compute action for element body forces
      spatialdata::spatialdb::SpatialDB* db = _gravityField;
      const int err = db->query(&gravVec[0], gravVec.size(),
        &quadPtsGlobal[0], spaceDim, cs);
      if (err)
        throw std::runtime_error("Unable to get gravity vector for point.");
      _normalizer->nondimensionalize(&gravVec[0], gravVec.size(),
          gravityScale);
      const PylithScalar wtVertex = density[0] * volume / 4.0;
      for (int iBasis=0; iBasis < numBasis; ++iBasis)
        for (int iDim=0; iDim < spaceDim; ++iDim)
            _cellVector[iBasis * spaceDim + iDim] += wtVertex * gravVec[iDim];
      PetscLogFlops(numBasis*spaceDim*2 + numBasis*spaceDim*2);
    } // if

    // Compute action for inertial terms
    const PylithScalar wtVertex = density[0] * volume / 16.0;
    for (int iBasis = 0; iBasis < numBasis; ++iBasis)
      for (int jBasis = 0; jBasis < numBasis; ++jBasis)
        for (int iDim = 0; iDim < spaceDim; ++iDim)
            _cellVector[iBasis*spaceDim+iDim] -= 
	      wtVertex * accArray[jBasis*spaceDim+iDim];

#if defined(DETAILED_EVENT_LOGGING)
    PetscLogFlops(3 + numBasis*numBasis*spaceDim*2);
    _logger->eventEnd(computeEvent);
    _logger->eventBegin(stressEvent);
#endif

    // Numerical damping. Compute displacements adjusted by velocity
    // times normalized viscosity.
    for(PetscInt i = 0; i < dispTSize; ++i) {dispAdjCell[i] = dispTArray[i] + viscosity * velArray[i];}
    err = DMComplexVecRestoreClosure(dmMesh, accSection,   accVec,   cell, &accSize,   &accArray);CHECK_PETSC_ERROR(err);
    err = DMComplexVecRestoreClosure(dmMesh, velSection,   velVec,   cell, &velSize,   &velArray);CHECK_PETSC_ERROR(err);
    err = DMComplexVecRestoreClosure(dmMesh, dispTSection, dispTVec, cell, &dispTSize, &dispTArray);CHECK_PETSC_ERROR(err);

    // Compute B(transpose) * sigma, first computing strains
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
    _cellVector[0] -= (d1*stressCell[5]+c1*stressCell[3]+b1*stressCell[0]) * volume;
    _cellVector[1] -= (d1*stressCell[4]+b1*stressCell[3]+c1*stressCell[1]) * volume;
    _cellVector[2] -= (b1*stressCell[5]+c1*stressCell[4]+d1*stressCell[2]) * volume;
    _cellVector[3] -= (d2*stressCell[5]+c2*stressCell[3]+b2*stressCell[0]) * volume;
    _cellVector[4] -= (d2*stressCell[4]+b2*stressCell[3]+c2*stressCell[1]) * volume;
    _cellVector[5] -= (b2*stressCell[5]+c2*stressCell[4]+d2*stressCell[2]) * volume;
    _cellVector[6] -= (d3*stressCell[5]+c3*stressCell[3]+b3*stressCell[0]) * volume;
    _cellVector[7] -= (d3*stressCell[4]+b3*stressCell[3]+c3*stressCell[1]) * volume;
    _cellVector[8] -= (b3*stressCell[5]+c3*stressCell[4]+d3*stressCell[2]) * volume;
    _cellVector[9] -= (d4*stressCell[5]+c4*stressCell[3]+b4*stressCell[0]) * volume;
    _cellVector[10] -= (d4*stressCell[4]+b4*stressCell[3]+c4*stressCell[1]) * volume;
    _cellVector[11] -= (b4*stressCell[5]+c4*stressCell[4]+d4*stressCell[2]) * volume;

#if defined(DETAILED_EVENT_LOGGING)
    PetscLogFlops(84);
    _logger->eventEnd(computeEvent);
    _logger->eventBegin(updateEvent);
#endif

    // Assemble cell contribution into field
    err = DMComplexVecSetClosure(dmMesh, residualSection, residualVec, cell, &_cellVector[0], ADD_VALUES);CHECK_PETSC_ERROR(err);
#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(updateEvent);
#endif
  } // for

#if !defined(DETAILED_EVENT_LOGGING)
  PetscLogFlops(numCells*(3 + numBasis*numBasis*spaceDim*2 + 196+84));
  _logger->eventEnd(computeEvent);
#endif
} // integrateResidual

// ----------------------------------------------------------------------
// Integrate constributions to residual term (r) for operator.
void
pylith::feassemble::ElasticityExplicitTet4::integrateResidualLumped(
        const topology::Field<topology::Mesh>& residual,
        const PylithScalar t,
        topology::SolutionFields* const fields)
{ // integrateResidualLumped
  /// Member prototype for _elasticityResidualXD()
  typedef void (pylith::feassemble::ElasticityExplicitTet4::*elasticityResidual_fn_type)
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
  DM              dmMesh = fields->mesh().dmMesh();
  IS              cellIS;
  const PetscInt *cells;
  PetscInt        numCells;
  PetscErrorCode  err;

  assert(dmMesh);
  err = DMComplexGetStratumIS(dmMesh, "material-id", _material->id(), &cellIS);CHECK_PETSC_ERROR(err);
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
  err = DMComplexGetCoordinateSection(dmMesh, &coordSection);CHECK_PETSC_ERROR(err);
  err = DMGetCoordinatesLocal(dmMesh, &coordVec);CHECK_PETSC_ERROR(err);
  assert(coordSection);assert(coordVec);

  assert(0 != _normalizer);
  const PylithScalar lengthScale = _normalizer->lengthScale();
  const PylithScalar gravityScale =
    _normalizer->pressureScale() / (_normalizer->lengthScale() *
            _normalizer->densityScale());

  const PylithScalar dt = _dt;
  assert(_normViscosity > 0.0);
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
#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(restrictEvent);
#endif

    // Restrict input fields to cell
    const PetscScalar *accArray, *velArray, *dispTArray;
    PetscInt           accSize,   velSize,   dispTSize;
    err = DMComplexVecGetClosure(dmMesh, accSection,   accVec,   cell, &accSize,   &accArray);CHECK_PETSC_ERROR(err);
    err = DMComplexVecGetClosure(dmMesh, velSection,   velVec,   cell, &velSize,   &velArray);CHECK_PETSC_ERROR(err);
    err = DMComplexVecGetClosure(dmMesh, dispTSection, dispTVec, cell, &dispTSize, &dispTArray);CHECK_PETSC_ERROR(err);
    assert(velSize   == accSize);
    assert(dispTSize == accSize);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(restrictEvent);
    _logger->eventBegin(geometryEvent);
#endif

    // Compute geometry information for current cell
    const PetscScalar *coords;
    PetscInt           coordsSize;
    err = DMComplexVecGetClosure(dmMesh, coordSection, coordVec, cell, &coordsSize, &coords);CHECK_PETSC_ERROR(err);
    for(PetscInt i = 0; i < coordsSize; ++i) {coordinatesCell[i] = coords[i];}
    const PylithScalar volume = _volume(coordinatesCell);
    assert(volume > 0.0);
    err = DMComplexVecRestoreClosure(dmMesh, coordSection, coordVec, cell, &coordsSize, &coords);CHECK_PETSC_ERROR(err);

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
    if (0 != _gravityField) {
      const spatialdata::geocoords::CoordSys* cs = fields->mesh().coordsys();
      assert(0 != cs);

      quadPtsGlobal = 0.0;
      for (int iBasis=0; iBasis < numBasis; ++iBasis)
        for (int iDim=0; iDim < spaceDim; ++iDim)
          quadPtsGlobal[iDim] += 
	    coordinatesCell[iBasis*spaceDim+iDim] / numBasis;
      _normalizer->dimensionalize(&quadPtsGlobal[0], quadPtsGlobal.size(),
          lengthScale);

      // Compute action for element body forces
      spatialdata::spatialdb::SpatialDB* db = _gravityField;
      const int err = db->query(&gravVec[0], gravVec.size(),
        &quadPtsGlobal[0], spaceDim, cs);
      if (err)
        throw std::runtime_error("Unable to get gravity vector for point.");
      _normalizer->nondimensionalize(&gravVec[0], gravVec.size(),
          gravityScale);
      const PylithScalar wtVertex = density[0] * volume / 4.0;
      for (int iBasis=0; iBasis < numBasis; ++iBasis)
        for (int iDim=0; iDim < spaceDim; ++iDim)
            _cellVector[iBasis * spaceDim + iDim] += wtVertex * gravVec[iDim];
      PetscLogFlops(numBasis*spaceDim*2 + numBasis*spaceDim*2);
    } // if

    // Compute action for inertial terms
    const PylithScalar wtVertex = density[0] * volume / 4.0;
    for(PetscInt i = 0; i < dispTSize; ++i) {_cellVector[i] -= wtVertex * accArray[i];}

#if defined(DETAILED_EVENT_LOGGING)
    PetscLogFlops(2 + numBasis*spaceDim*2);
    _logger->eventEnd(computeEvent);
    _logger->eventBegin(stressEvent);
#endif

    // Numerical damping. Compute displacements adjusted by velocity
    // times normalized viscosity.
    for(PetscInt i = 0; i < dispTSize; ++i) {dispAdjCell[i] = dispTArray[i] + viscosity * velArray[i];}
    err = DMComplexVecRestoreClosure(dmMesh, accSection,   accVec,   cell, &accSize,   &accArray);CHECK_PETSC_ERROR(err);
    err = DMComplexVecRestoreClosure(dmMesh, velSection,   velVec,   cell, &velSize,   &velArray);CHECK_PETSC_ERROR(err);
    err = DMComplexVecRestoreClosure(dmMesh, dispTSection, dispTVec, cell, &dispTSize, &dispTArray);CHECK_PETSC_ERROR(err);

    // Compute B(transpose) * sigma, first computing strains
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
    _cellVector[0] -= 
      (d1*stressCell[5]+c1*stressCell[3]+b1*stressCell[0]) * volume;
    _cellVector[1] -= 
      (d1*stressCell[4]+b1*stressCell[3]+c1*stressCell[1]) * volume;
    _cellVector[2] -= 
      (b1*stressCell[5]+c1*stressCell[4]+d1*stressCell[2]) * volume;
    _cellVector[3] -= 
      (d2*stressCell[5]+c2*stressCell[3]+b2*stressCell[0]) * volume;
    _cellVector[4] -= 
      (d2*stressCell[4]+b2*stressCell[3]+c2*stressCell[1]) * volume;
    _cellVector[5] -= 
      (b2*stressCell[5]+c2*stressCell[4]+d2*stressCell[2]) * volume;
    _cellVector[6] -= 
      (d3*stressCell[5]+c3*stressCell[3]+b3*stressCell[0]) * volume;
    _cellVector[7] -= 
      (d3*stressCell[4]+b3*stressCell[3]+c3*stressCell[1]) * volume;
    _cellVector[8] -= 
      (b3*stressCell[5]+c3*stressCell[4]+d3*stressCell[2]) * volume;
    _cellVector[9] -= 
      (d4*stressCell[5]+c4*stressCell[3]+b4*stressCell[0]) * volume;
    _cellVector[10] -= 
      (d4*stressCell[4]+b4*stressCell[3]+c4*stressCell[1]) * volume;
    _cellVector[11] -= 
      (b4*stressCell[5]+c4*stressCell[4]+d4*stressCell[2]) * volume;

#if defined(DETAILED_EVENT_LOGGING)
    PetscLogFlops(84);
    _logger->eventEnd(computeEvent);
    _logger->eventBegin(updateEvent);
#endif

    // Assemble cell contribution into field
    err = DMComplexVecSetClosure(dmMesh, residualSection, residualVec, cell, &_cellVector[0], ADD_VALUES);CHECK_PETSC_ERROR(err);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(updateEvent);
#endif
  } // for
  err = ISRestoreIndices(cellIS, &cells);CHECK_PETSC_ERROR(err);
  err = ISDestroy(&cellIS);CHECK_PETSC_ERROR(err);

#if !defined(DETAILED_EVENT_LOGGING)
  PetscLogFlops(numCells*(2 + numBasis*spaceDim*2 + 196+84));
  _logger->eventEnd(computeEvent);
#endif
} // integrateResidualLumped

// ----------------------------------------------------------------------
// Compute matrix associated with operator.
void
pylith::feassemble::ElasticityExplicitTet4::integrateJacobian(
					topology::Jacobian* jacobian,
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
  const int tensorSize = _material->tensorSize();
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
  err = DMComplexGetStratumIS(dmMesh, "material-id", _material->id(), &cellIS);CHECK_PETSC_ERROR(err);
  err = ISGetSize(cellIS, &numCells);CHECK_PETSC_ERROR(err);
  err = ISGetIndices(cellIS, &cells);CHECK_PETSC_ERROR(err);

  // Get sections
  topology::Field<topology::Mesh>& dispT = fields->get("disp(t)");
  PetscSection dispTSection = dispT.petscSection();
  Vec          dispTVec     = dispT.localVector();
  assert(dispTSection);assert(dispTVec);

  // Get sparse matrix
  const PetscMat jacobianMat = jacobian->matrix();
  assert(0 != jacobianMat);

  // Get parameters used in integration.
  const PylithScalar dt = _dt;
  const PylithScalar dt2 = dt*dt;
  assert(dt > 0);

  scalar_array coordinatesCell(numBasis*spaceDim);
  PetscSection coordSection;
  Vec          coordVec;
  err = DMComplexGetCoordinateSection(dmMesh, &coordSection);CHECK_PETSC_ERROR(err);
  err = DMGetCoordinatesLocal(dmMesh, &coordVec);CHECK_PETSC_ERROR(err);
  assert(coordSection);assert(coordVec);

  _logger->eventEnd(setupEvent);
#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventBegin(computeEvent);
#endif

  // Loop over cells
  for(PetscInt c = 0; c < numCells; ++c) {
    const PetscInt cell = cells[c];
#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(geometryEvent);
#endif

    // Compute geometry information for current cell
    const PetscScalar *coords;
    PetscInt           coordsSize;
    err = DMComplexVecGetClosure(dmMesh, coordSection, coordVec, cell, &coordsSize, &coords);CHECK_PETSC_ERROR(err);
    for(PetscInt i = 0; i < coordsSize; ++i) {coordinatesCell[i] = coords[i];}
    const PylithScalar volume = _volume(coordinatesCell);
    assert(volume > 0.0);
    err = DMComplexVecRestoreClosure(dmMesh, coordSection, coordVec, cell, &coordsSize, &coords);CHECK_PETSC_ERROR(err);

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
    _resetCellMatrix();

    // Get material physical properties at quadrature points for this cell
    const scalar_array& density = _material->calcDensity();
    assert(density.size() == 1);

    // Compute Jacobian for inertial terms
    const PylithScalar wtVertex = density[0] * volume / (16.0 * dt2);
    for (int iBasis = 0; iBasis < numBasis; ++iBasis)
      for (int jBasis = 0; jBasis < numBasis; ++jBasis)
	for (int iDim=0; iDim < spaceDim; ++iDim) {
	  const int iBlock = (iBasis*spaceDim + iDim) * (numBasis*spaceDim);
	  const int jBlock = (jBasis*spaceDim + iDim);
	  _cellMatrix[iBlock+jBlock] += wtVertex;
	} // for
    
#if defined(DETAILED_EVENT_LOGGING)
    PetscLogFlops(numQuadPts*(3+numBasis*numBasis*spaceDim*1));
    _logger->eventEnd(computeEvent);
    _logger->eventBegin(updateEvent);
#endif
    
    // Assemble cell contribution into PETSc matrix.
    err = DMComplexMatSetClosure(dmMesh, dispTSection, PETSC_NULL, jacobianMat, cell, &_cellMatrix[0], ADD_VALUES);
    CHECK_PETSC_ERROR_MSG(err, "Update to PETSc Mat failed.");

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(updateEvent);
#endif
  } // for
  err = ISRestoreIndices(cellIS, &cells);CHECK_PETSC_ERROR(err);
  err = ISDestroy(&cellIS);CHECK_PETSC_ERROR(err);

#if !defined(DETAILED_EVENT_LOGGING)
  PetscLogFlops(numCells*(3+numBasis*numBasis*spaceDim*1));
  _logger->eventEnd(computeEvent);
#endif

  _needNewJacobian = false;
  _material->resetNeedNewJacobian();
} // integrateJacobian

// ----------------------------------------------------------------------
// Compute matrix associated with operator.
void
pylith::feassemble::ElasticityExplicitTet4::integrateJacobian(
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
  DM              dmMesh = fields->mesh().dmMesh();
  IS              cellIS;
  const PetscInt *cells;
  PetscInt        numCells;
  PetscErrorCode  err;

  assert(dmMesh);
  err = DMComplexGetStratumIS(dmMesh, "material-id", _material->id(), &cellIS);CHECK_PETSC_ERROR(err);
  err = ISGetSize(cellIS, &numCells);CHECK_PETSC_ERROR(err);
  err = ISGetIndices(cellIS, &cells);CHECK_PETSC_ERROR(err);

  // Get parameters used in integration.
  const PylithScalar dt = _dt;
  const PylithScalar dt2 = dt*dt;
  assert(dt > 0);

  // Get sections
  PetscSection jacSection = jacobian->petscSection();
  Vec          jacVec     = jacobian->localVector();

  scalar_array coordinatesCell(numBasis*spaceDim);
  PetscSection coordSection;
  Vec          coordVec;
  err = DMComplexGetCoordinateSection(dmMesh, &coordSection);CHECK_PETSC_ERROR(err);
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
    const PetscScalar *coords;
    PetscInt           coordsSize;
    err = DMComplexVecGetClosure(dmMesh, coordSection, coordVec, cell, &coordsSize, &coords);CHECK_PETSC_ERROR(err);
    for(PetscInt i = 0; i < coordsSize; ++i) {coordinatesCell[i] = coords[i];}
    const PylithScalar volume = _volume(coordinatesCell);
    assert(volume > 0.0);
    err = DMComplexVecRestoreClosure(dmMesh, coordSection, coordVec, cell, &coordsSize, &coords);CHECK_PETSC_ERROR(err);

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
    assert(volume > 0.0);
    _cellVector = density[0] * volume / (4.0 * dt2);
    
#if defined(DETAILED_EVENT_LOGGING)
    PetscLogFlops(3);
    _logger->eventEnd(computeEvent);
    _logger->eventBegin(updateEvent);
#endif
    
    // Assemble cell contribution into lumped matrix.
    err = DMComplexVecSetClosure(dmMesh, jacSection, jacVec, cell, &_cellVector[0], ADD_VALUES);CHECK_PETSC_ERROR(err);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(updateEvent);
#endif
  } // for
  err = ISRestoreIndices(cellIS, &cells);CHECK_PETSC_ERROR(err);
  err = ISDestroy(&cellIS);CHECK_PETSC_ERROR(err);

#if !defined(DETAILED_EVENT_LOGGING)
  PetscLogFlops(numCells*3);
  _logger->eventEnd(computeEvent);
#endif

  _needNewJacobian = false;
  _material->resetNeedNewJacobian();
} // integrateJacobian

// ----------------------------------------------------------------------
// Compute volume of tetrahedral cell.
PylithScalar
pylith::feassemble::ElasticityExplicitTet4::_volume(
			     const scalar_array& coordinatesCell) const
{ // __volume
  assert(12 == coordinatesCell.size());

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
