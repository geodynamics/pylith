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
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "AbsorbingDampers.hh" // implementation of object methods

#include "pylith/topology/Fields.hh" // HOLDSA Fields
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields
#include "pylith/topology/Jacobian.hh" // USES Jacobian
#include "pylith/feassemble/CellGeometry.hh" // USES CellGeometry

#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/petscerror.h" // USES CHECK_PETSC_ERROR
#include <cstring> // USES memcpy()
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

//#define DETAILED_EVENT_LOGGING

// ----------------------------------------------------------------------
// Default constructor.
pylith::bc::AbsorbingDampers::AbsorbingDampers(void) :
  _db(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::bc::AbsorbingDampers::~AbsorbingDampers(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void 
pylith::bc::AbsorbingDampers::deallocate(void)
{ // deallocate
  BCIntegratorSubMesh::deallocate();
  _db = 0; // :TODO: Use shared pointer
} // deallocate
  
// ----------------------------------------------------------------------
// Initialize boundary condition. Determine orienation and compute traction
// vector at integration points.
void
pylith::bc::AbsorbingDampers::initialize(const topology::Mesh& mesh,
					 const PylithScalar upDir[3])
{ // initialize
  assert(_boundaryMesh);
  assert(_quadrature);
  assert(_db);

  _initializeLogger();

  scalar_array up(3);
  for (int i=0; i < 3; ++i) {
    up[i] = upDir[i];
  } // for

  const int numCorners = _quadrature->numBasis();

  // Get 'surface' cells (1 dimension lower than top-level cells)
  PetscDM subMesh = _boundaryMesh->dmMesh();
  PetscInt cStart, cEnd;
  PetscErrorCode err;

  assert(subMesh);
  err = DMPlexGetHeightStratum(subMesh, 1, &cStart, &cEnd);CHECK_PETSC_ERROR(err);

  // Get damping constants at each quadrature point and rotate to
  // global coordinate frame using orientation information
  const feassemble::CellGeometry& cellGeometry = _quadrature->refGeometry();
  const int cellDim = _quadrature->cellDim() > 0 ? _quadrature->cellDim() : 1;
  const int numBasis = _quadrature->numBasis();
  const int numQuadPts = _quadrature->numQuadPts();
  const int spaceDim = cellGeometry.spaceDim();
  const int fiberDim = numQuadPts * spaceDim;
  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("BoundaryConditions");

  delete _parameters;
  _parameters = new topology::Fields<topology::Field<topology::SubMesh> >(*_boundaryMesh);
  assert(_parameters);
  _parameters->add("damping constants", "damping_constants", topology::FieldBase::FACES_FIELD, fiberDim);
  _parameters->get("damping constants").allocate();

  logger.stagePop();

  // Containers for orientation information
  const int orientationSize = spaceDim * spaceDim;
  const int jacobianSize = spaceDim * cellDim;
  scalar_array jacobian(jacobianSize);
  PylithScalar jacobianDet = 0;
  scalar_array orientation(orientationSize);

  // open database with material property information
  _db->open();
  int numValues = 0;
  if (cellDim > 0) {
    numValues = 3;
    const char* valueNames[] = { "density", "vp", "vs" };
    _db->queryVals(valueNames, numValues);
  } else {
    numValues = 2;
    const char* valueNames[] = { "density", "vp" };
    _db->queryVals(valueNames, numValues);
  } // else

  // Container for data returned in query of database
  scalar_array queryData(numValues);
  scalar_array quadPtRef(cellDim);
  scalar_array quadPtsGlobal(numQuadPts*spaceDim);

  // Container for damping constants for current cell
  scalar_array dampingConstsLocal(spaceDim);

  scalar_array coordinatesCell(numBasis*spaceDim);
  PetscSection coordSection;
  PetscVec coordVec;
  err = DMPlexGetCoordinateSection(subMesh, &coordSection);CHECK_PETSC_ERROR(err);
  err = DMGetCoordinatesLocal(subMesh, &coordVec);CHECK_PETSC_ERROR(err);
  assert(coordSection);assert(coordVec);

  assert(_normalizer);
  const PylithScalar lengthScale = _normalizer->lengthScale();
  const PylithScalar densityScale = _normalizer->densityScale();
  assert(_normalizer->timeScale() > 0);
  const PylithScalar velocityScale = _normalizer->lengthScale() / _normalizer->timeScale();

  PetscSection valueSection = _parameters->get("damping constants").petscSection();
  PetscVec valueVec = _parameters->get("damping constants").localVector();
  PetscScalar *dampingConstsArray;
  assert(valueSection);assert(valueVec);

  const spatialdata::geocoords::CoordSys* cs = _boundaryMesh->coordsys();

  // Compute quadrature information
  _quadrature->initializeGeometry();
#if defined(PRECOMPUTE_GEOMETRY)
  _quadrature->computeGeometry(*_boundaryMesh, cells);
#endif

  err = VecGetArray(valueVec, &dampingConstsArray);CHECK_PETSC_ERROR(err);
  for(PetscInt c = cStart; c < cEnd; ++c) {
    // Compute geometry information for current cell
#if defined(PRECOMPUTE_GEOMETRY)
    _quadrature->retrieveGeometry(*c_iter);
#else
    const PetscScalar *coords;
    PetscInt coordsSize;
    err = DMPlexVecGetClosure(subMesh, coordSection, coordVec, c, &coordsSize, &coords);CHECK_PETSC_ERROR(err);
    for (PetscInt i = 0; i < coordsSize; ++i) {
      coordinatesCell[i] = coords[i];
    } // for
    _quadrature->computeGeometry(coordinatesCell, c);
    err = DMPlexVecRestoreClosure(subMesh, coordSection, coordVec, c, &coordsSize, &coords);CHECK_PETSC_ERROR(err);
#endif
    PetscInt ddof, doff;

    err = PetscSectionGetDof(valueSection, c, &ddof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(valueSection, c, &doff);CHECK_PETSC_ERROR(err);
    assert(ddof == fiberDim);

    const scalar_array& quadPtsNondim = _quadrature->quadPts();
    const scalar_array& quadPtsRef = _quadrature->quadPtsRef();
    quadPtsGlobal = quadPtsNondim;
    _normalizer->dimensionalize(&quadPtsGlobal[0], quadPtsGlobal.size(), lengthScale);

    for(int iQuad = 0; iQuad < numQuadPts; ++iQuad) {
      // Compute damping constants in normal/tangential coordinates
      const int err = _db->query(&queryData[0], numValues, 
                                 &quadPtsGlobal[iQuad*spaceDim], spaceDim, cs);
      if (err) {
        std::ostringstream msg;
        msg << "Could not find parameters for physical properties at \n"
            << "(";
        for (int i=0; i < spaceDim; ++i)
          msg << "  " << quadPtsGlobal[iQuad*spaceDim+i];
        msg << ") for absorbing boundary condition " << _label << "\n"
            << "using spatial database " << _db->label() << ".";
        throw std::runtime_error(msg.str());
      } // if
      // Nondimensionalize damping constants
      const PylithScalar densityN = _normalizer->nondimensionalize(queryData[0], densityScale);
      const PylithScalar vpN = _normalizer->nondimensionalize(queryData[1], velocityScale);
      const PylithScalar vsN = (3 == numValues) ? _normalizer->nondimensionalize(queryData[2], velocityScale) : 0.0;
      
      const PylithScalar constTangential = densityN * vsN;
      const PylithScalar constNormal = densityN * vpN;
      const int numTangential = spaceDim-1;
      for (int iDim=0; iDim < numTangential; ++iDim)
        dampingConstsLocal[iDim] = constTangential;
      dampingConstsLocal[spaceDim-1] = constNormal;

      // Compute normal/tangential orientation
      memcpy(&quadPtRef[0], &quadPtsRef[iQuad*cellDim], cellDim*sizeof(PylithScalar));
      cellGeometry.jacobian(&jacobian, &jacobianDet, coordinatesCell, quadPtRef);
      cellGeometry.orientation(&orientation, jacobian, jacobianDet, up);
      assert(jacobianDet > 0.0);
      orientation /= jacobianDet;

      for (int iDim=0; iDim < spaceDim; ++iDim) {
        dampingConstsArray[doff+iQuad*spaceDim+iDim] = 0.0;
        for (int jDim=0; jDim < spaceDim; ++jDim)
          dampingConstsArray[doff+iQuad*spaceDim+iDim] += dampingConstsLocal[jDim]*orientation[jDim*spaceDim+iDim];
        // Ensure damping constants are positive
        dampingConstsArray[doff+iQuad*spaceDim+iDim] = fabs(dampingConstsArray[doff+iQuad*spaceDim+iDim]);
      } // for
    } // for
  } // for
  err = VecRestoreArray(valueVec, &dampingConstsArray);CHECK_PETSC_ERROR(err);

  _db->close();
} // initialize

// ----------------------------------------------------------------------
// Integrate contributions to residual term (r) for operator.
void
pylith::bc::AbsorbingDampers::integrateResidual(
			     const topology::Field<topology::Mesh>& residual,
			     const PylithScalar t,
			     topology::SolutionFields* const fields)
{ // integrateResidual
  assert(_quadrature);
  assert(_boundaryMesh);
  assert(_parameters);
  assert(fields);
  assert(_logger);

  const int setupEvent = _logger->eventId("AdIR setup");
  const int geometryEvent = _logger->eventId("AdIR geometry");
  const int computeEvent = _logger->eventId("AdIR compute");
  const int restrictEvent = _logger->eventId("AdIR restrict");
  const int updateEvent = _logger->eventId("AdIR update");

  _logger->eventBegin(setupEvent);

  // Get cell geometry information that doesn't depend on cell
  const int numQuadPts = _quadrature->numQuadPts();
  const scalar_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == numQuadPts);
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();

  // Allocate vectors for cell values.
  _initCellVector();

  // Get cell information
  PetscDM       subMesh = _boundaryMesh->dmMesh();
  PetscIS       subpointIS;
  PetscInt cStart, cEnd;
  PetscErrorCode err;

  assert(subMesh);
  err = DMPlexGetHeightStratum(subMesh, 1, &cStart, &cEnd);CHECK_PETSC_ERROR(err);
  err = DMPlexCreateSubpointIS(subMesh, &subpointIS);CHECK_PETSC_ERROR(err);

  // Get sections
  PetscSection dampingConstsSection = _parameters->get("damping constants").petscSection();
  PetscVec dampingConstsVec = _parameters->get("damping constants").localVector();
  PetscScalar *dampingConstsArray;
  assert(dampingConstsSection);assert(dampingConstsVec);

  // Use _cellVector for cell residual.
  PetscSection residualSection = residual.petscSection(), residualSubsection;
  PetscVec residualVec = residual.localVector();
  assert(residualSection);assert(residualVec);
  err = PetscSectionCreateSubmeshSection(residualSection, subpointIS, &residualSubsection);
  
  PetscSection velSection = fields->get("velocity(t)").petscSection(), velSubsection;
  PetscVec velVec = fields->get("velocity(t)").localVector();
  assert(velSection);assert(velVec);
  err = PetscSectionCreateSubmeshSection(velSection, subpointIS, &velSubsection);
  err = ISDestroy(&subpointIS);CHECK_PETSC_ERROR(err);
  
#if !defined(PRECOMPUTE_GEOMETRY)
  scalar_array coordinatesCell(numBasis*spaceDim);
  PetscSection coordSection;
  PetscVec coordVec;
  err = DMPlexGetCoordinateSection(subMesh, &coordSection);CHECK_PETSC_ERROR(err);
  err = DMGetCoordinatesLocal(subMesh, &coordVec);CHECK_PETSC_ERROR(err);
  assert(coordSection);assert(coordVec);
#endif

  _logger->eventEnd(setupEvent);
#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventBegin(computeEvent);
#endif

  err = VecGetArray(dampingConstsVec, &dampingConstsArray);CHECK_PETSC_ERROR(err);
  for(PetscInt c = cStart; c < cEnd; ++c) {
    // Get geometry information for current cell
#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(geometryEvent);
#endif
#if defined(PRECOMPUTE_GEOMETRY)
    _quadrature->retrieveGeometry(*c_iter);
#else
    const PetscScalar *coordsArray;
    PetscInt           coordsSize;
    err = DMPlexVecGetClosure(subMesh, coordSection, coordVec, c, &coordsSize, &coordsArray);CHECK_PETSC_ERROR(err);
    for(PetscInt i = 0; i < coordsSize; ++i) {coordinatesCell[i] = coordsArray[i];}
    _quadrature->computeGeometry(coordinatesCell, c);
    err = DMPlexVecRestoreClosure(subMesh, coordSection, coordVec, c, &coordsSize, &coordsArray);CHECK_PETSC_ERROR(err);
#endif

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(geometryEvent);
    _logger->eventBegin(restrictEvent);
#endif

    // Reset element vector to zero
    _resetCellVector();

    // Restrict input fields to cell
    const PetscScalar *velArray = NULL;
    PetscInt velSize;
    PetscInt ddof, doff;

    err = PetscSectionGetDof(dampingConstsSection, c, &ddof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dampingConstsSection, c, &doff);CHECK_PETSC_ERROR(err);
    assert(ddof == numQuadPts*spaceDim);
    err = DMPlexVecGetClosure(subMesh, velSubsection, velVec, c, &velSize, &velArray);CHECK_PETSC_ERROR(err);
    assert(velSize == numBasis*spaceDim);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(restrictEvent);
    _logger->eventBegin(computeEvent);
#endif

    // Get cell geometry information that depends on cell
    const scalar_array& basis = _quadrature->basis();
    const scalar_array& jacobianDet = _quadrature->jacobianDet();

    // Compute action for absorbing bc terms
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      const PylithScalar wt = quadWts[iQuad] * jacobianDet[iQuad];

      for (int iBasis=0; iBasis < numBasis; ++iBasis) {
        const PylithScalar valI = wt*basis[iQuad*numBasis+iBasis];
        for (int jBasis=0; jBasis < numBasis; ++jBasis) {
          const PylithScalar valIJ = valI * basis[iQuad*numBasis+jBasis];
          for (int iDim=0; iDim < spaceDim; ++iDim)
            _cellVector[iBasis*spaceDim+iDim] -= 
              dampingConstsArray[doff+iQuad*spaceDim+iDim] *
              valIJ * velArray[jBasis*spaceDim+iDim];
        } // for
      } // for
    } // for
    err = DMPlexVecRestoreClosure(subMesh, velSubsection, velVec, c, &velSize, &velArray);CHECK_PETSC_ERROR(err);
    err = DMPlexVecSetClosure(subMesh, residualSubsection, residualVec, c, &_cellVector[0], ADD_VALUES);CHECK_PETSC_ERROR(err);

#if defined(DETAILED_EVENT_LOGGING)
    PetscLogFlops(numQuadPts*(1+numBasis*(1+numBasis*(3*spaceDim))));
    _logger->eventEnd(computeEvent);
    _logger->eventBegin(updateEvent);
#endif

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(updateEvent);
#endif
  } // for
  err = VecRestoreArray(dampingConstsVec, &dampingConstsArray);CHECK_PETSC_ERROR(err);
  err = PetscSectionDestroy(&residualSubsection);CHECK_PETSC_ERROR(err);
  err = PetscSectionDestroy(&velSubsection);CHECK_PETSC_ERROR(err);

#if !defined(DETAILED_EVENT_LOGGING)
  PetscLogFlops((cEnd-cStart)*numQuadPts*(1+numBasis*(1+numBasis*(3*spaceDim))));
  _logger->eventEnd(computeEvent);
#endif
} // integrateResidual

// ----------------------------------------------------------------------
// Integrate contributions to residual term (r) for operator.
void
pylith::bc::AbsorbingDampers::integrateResidualLumped(
           const topology::Field<topology::Mesh>& residual,
           const PylithScalar t,
           topology::SolutionFields* const fields)
{ // integrateResidualLumped
  assert(_quadrature);
  assert(_boundaryMesh);
  assert(_parameters);
  assert(fields);
  assert(_logger);

  const int setupEvent = _logger->eventId("AdIR setup");
  const int geometryEvent = _logger->eventId("AdIR geometry");
  const int computeEvent = _logger->eventId("AdIR compute");
  const int restrictEvent = _logger->eventId("AdIR restrict");
  const int updateEvent = _logger->eventId("AdIR update");

  _logger->eventBegin(setupEvent);

  // Get cell geometry information that doesn't depend on cell
  const int numQuadPts = _quadrature->numQuadPts();
  const scalar_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == numQuadPts);
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();

  // Allocate vectors for cell values.
  _initCellVector();

  // Get cell information
  PetscDM       subMesh = _boundaryMesh->dmMesh();
  PetscIS       subpointIS;
  PetscInt cStart, cEnd;
  PetscErrorCode err;

  assert(subMesh);
  err = DMPlexGetHeightStratum(subMesh, 1, &cStart, &cEnd);CHECK_PETSC_ERROR(err);
  err = DMPlexCreateSubpointIS(subMesh, &subpointIS);CHECK_PETSC_ERROR(err);

  // Get sections
  PetscSection dampingConstsSection = _parameters->get("damping constants").petscSection();
  PetscVec dampingConstsVec = _parameters->get("damping constants").localVector();
  PetscScalar *dampingConstsArray;
  assert(dampingConstsSection);assert(dampingConstsVec);

  // Use _cellVector for cell values.
  PetscSection residualSection = residual.petscSection(), residualSubsection;
  PetscVec residualVec = residual.localVector();
  assert(residualSection);assert(residualVec);
  err = PetscSectionCreateSubmeshSection(residualSection, subpointIS, &residualSubsection);

  PetscSection velSection = fields->get("velocity(t)").petscSection(), velSubsection;
  PetscVec velVec = fields->get("velocity(t)").localVector();
  assert(velSection);assert(velVec);
  err = PetscSectionCreateSubmeshSection(velSection, subpointIS, &velSubsection);
  err = ISDestroy(&subpointIS);CHECK_PETSC_ERROR(err);

#if !defined(PRECOMPUTE_GEOMETRY)
  scalar_array coordinatesCell(numBasis*spaceDim);
  PetscSection coordSection;
  PetscVec coordVec;
  err = DMPlexGetCoordinateSection(subMesh, &coordSection);CHECK_PETSC_ERROR(err);
  err = DMGetCoordinatesLocal(subMesh, &coordVec);CHECK_PETSC_ERROR(err);
  assert(coordSection);assert(coordVec);
#endif

  _logger->eventEnd(setupEvent);
#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventBegin(computeEvent);
#endif

  err = VecGetArray(dampingConstsVec, &dampingConstsArray);CHECK_PETSC_ERROR(err);
  for(PetscInt c = cStart; c < cEnd; ++c) {
    // Get geometry information for current cell
#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(geometryEvent);
#endif
#if defined(PRECOMPUTE_GEOMETRY)
    _quadrature->retrieveGeometry(*c_iter);
#else
    const PetscScalar *coordsArray;
    PetscInt           coordsSize;
    err = DMPlexVecGetClosure(subMesh, coordSection, coordVec, c, &coordsSize, &coordsArray);CHECK_PETSC_ERROR(err);
    for(PetscInt i = 0; i < coordsSize; ++i) {coordinatesCell[i] = coordsArray[i];}
    _quadrature->computeGeometry(coordinatesCell, c);
    err = DMPlexVecRestoreClosure(subMesh, coordSection, coordVec, c, &coordsSize, &coordsArray);CHECK_PETSC_ERROR(err);
#endif

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(geometryEvent);
    _logger->eventBegin(restrictEvent);
#endif

    // Reset element vector to zero
    _resetCellVector();

    // Restrict input fields to cell
    const PetscScalar *velArray = NULL;
    PetscInt velSize;
    PetscInt ddof, doff;

    err = PetscSectionGetDof(dampingConstsSection, c, &ddof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dampingConstsSection, c, &doff);CHECK_PETSC_ERROR(err);
    assert(ddof == numQuadPts*spaceDim);
    err = DMPlexVecGetClosure(subMesh, velSubsection, velVec, c, &velSize, &velArray);CHECK_PETSC_ERROR(err);
    assert(velSize == numBasis*spaceDim);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(restrictEvent);
    _logger->eventBegin(computeEvent);
#endif

    // Get cell geometry information that depends on cell
    const scalar_array& basis = _quadrature->basis();
    const scalar_array& jacobianDet = _quadrature->jacobianDet();

    // Compute action for absorbing bc terms
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      const PylithScalar wt = quadWts[iQuad] * jacobianDet[iQuad];
      const int iQ = iQuad * numBasis;
      PylithScalar valJ = 0.0;
      for (int jBasis = 0; jBasis < numBasis; ++jBasis)
        valJ += basis[iQ + jBasis];
      valJ *= wt;
      for (int iBasis = 0; iBasis < numBasis; ++iBasis) {
        const PylithScalar valIJ = basis[iQ + iBasis] * valJ;
        for (int iDim = 0; iDim < spaceDim; ++iDim)
          _cellVector[iBasis*spaceDim+iDim] -= 
            dampingConstsArray[doff+iQuad*spaceDim+iDim] *
            valIJ * velArray[iBasis*spaceDim+iDim];
      } // for
    } // for
    err = DMPlexVecRestoreClosure(subMesh, velSubsection, velVec, c, &velSize, &velArray);CHECK_PETSC_ERROR(err);
    err = DMPlexVecSetClosure(subMesh, residualSubsection, residualVec, c, &_cellVector[0], ADD_VALUES);CHECK_PETSC_ERROR(err);
#if defined(DETAILED_EVENT_LOGGING)
    PetscLogFlops(numQuadPts*(1+numBasis+numBasis*(1+spaceDim*3)));
    _logger->eventEnd(computeEvent);
    _logger->eventBegin(updateEvent);
#endif

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(updateEvent);
#endif
  } // for
  err = VecRestoreArray(dampingConstsVec, &dampingConstsArray);CHECK_PETSC_ERROR(err);
  err = PetscSectionDestroy(&residualSubsection);CHECK_PETSC_ERROR(err);
  err = PetscSectionDestroy(&velSubsection);CHECK_PETSC_ERROR(err);

#if !defined(DETAILED_EVENT_LOGGING)
  PetscLogFlops((cEnd-cStart)*numQuadPts*(1+numBasis+numBasis*(1+spaceDim*3)));
  _logger->eventEnd(computeEvent);
#endif
} // integrateResidualLumped

// ----------------------------------------------------------------------
// Integrate contributions to Jacobian matrix (A) associated with
void
pylith::bc::AbsorbingDampers::integrateJacobian(
				      topology::Jacobian* jacobian,
				      const PylithScalar t,
				      topology::SolutionFields* const fields)
{ // integrateJacobian
  assert(_quadrature);
  assert(_boundaryMesh);
  assert(_logger);
  assert(jacobian);
  assert(fields);

  const int setupEvent = _logger->eventId("AdIJ setup");
  const int geometryEvent = _logger->eventId("AdIJ geometry");
  const int computeEvent = _logger->eventId("AdIJ compute");
  const int restrictEvent = _logger->eventId("AdIJ restrict");
  const int updateEvent = _logger->eventId("AdIJ update");

  _logger->eventBegin(setupEvent);

  // Get cell geometry information that doesn't depend on cell
  const int numQuadPts = _quadrature->numQuadPts();
  const scalar_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == numQuadPts);
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();

  // Get cell information
  PetscDM subMesh = _boundaryMesh->dmMesh();
  PetscIS subpointIS;
  PetscInt cStart, cEnd;
  PetscErrorCode err = 0;

  assert(subMesh);
  err = DMPlexGetHeightStratum(subMesh, 1, &cStart, &cEnd);CHECK_PETSC_ERROR(err);
  err = DMPlexCreateSubpointIS(subMesh, &subpointIS);CHECK_PETSC_ERROR(err);

  // Get sections
  PetscSection dampingConstsSection = _parameters->get("damping constants").petscSection();
  PetscVec dampingConstsVec = _parameters->get("damping constants").localVector();
  PetscScalar *dampingConstsArray;
  assert(dampingConstsSection);assert(dampingConstsVec);

  const topology::Field<topology::Mesh>& solution = fields->solution();
  PetscSection solutionSection = solution.petscSection(), solutionGlobalSection, solutionSubsection, solutionGlobalSubsection;
  PetscVec solutionVec = solution.localVector();
  PetscSF      sf;
  assert(solutionSection);assert(solutionVec);
  err = PetscSectionCreateSubmeshSection(solutionSection, subpointIS, &solutionSubsection);CHECK_PETSC_ERROR(err);
  err = DMGetPointSF(solution.dmMesh(), &sf);CHECK_PETSC_ERROR(err);
  err = PetscSectionCreateGlobalSection(solutionSection, sf, PETSC_FALSE, &solutionGlobalSection);CHECK_PETSC_ERROR(err);
  err = PetscSectionCreateSubmeshSection(solutionGlobalSection, subpointIS, &solutionGlobalSubsection);CHECK_PETSC_ERROR(err);
  err = ISDestroy(&subpointIS);CHECK_PETSC_ERROR(err);

  // Get sparse matrix
  const PetscMat jacobianMat = jacobian->matrix();
  assert(jacobianMat);

  // Get parameters used in integration.
  const PylithScalar dt = _dt;
  assert(dt > 0);

  // Allocate matrix for cell values.
  _initCellMatrix();

#if !defined(PRECOMPUTE_GEOMETRY)
  scalar_array coordinatesCell(numBasis*spaceDim);
  PetscSection coordSection;
  PetscVec coordVec;
  err = DMPlexGetCoordinateSection(subMesh, &coordSection);CHECK_PETSC_ERROR(err);
  err = DMGetCoordinatesLocal(subMesh, &coordVec);CHECK_PETSC_ERROR(err);
  assert(coordSection);assert(coordVec);
#endif

  _logger->eventEnd(setupEvent);
#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventBegin(computeEvent);
#endif

  err = VecGetArray(dampingConstsVec, &dampingConstsArray);CHECK_PETSC_ERROR(err);
  for(PetscInt c = cStart; c < cEnd; ++c) {
    // Compute geometry information for current cell
#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(geometryEvent);
#endif
#if defined(PRECOMPUTE_GEOMETRY)
    _quadrature->retrieveGeometry(*c_iter);
#else
    const PetscScalar *coordsArray;
    PetscInt           coordsSize;
    err = DMPlexVecGetClosure(subMesh, coordSection, coordVec, c, &coordsSize, &coordsArray);CHECK_PETSC_ERROR(err);
    for(PetscInt i = 0; i < coordsSize; ++i) {coordinatesCell[i] = coordsArray[i];}
    _quadrature->computeGeometry(coordinatesCell, c);
    err = DMPlexVecRestoreClosure(subMesh, coordSection, coordVec, c, &coordsSize, &coordsArray);CHECK_PETSC_ERROR(err);
#endif
#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(geometryEvent);
    _logger->eventBegin(restrictEvent);
#endif

    // Get damping constants
    PetscInt ddof, doff;

    err = PetscSectionGetDof(dampingConstsSection, c, &ddof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dampingConstsSection, c, &doff);CHECK_PETSC_ERROR(err);
    assert(ddof == numQuadPts*spaceDim);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(restrictEvent);
    _logger->eventBegin(computeEvent);
#endif

    // Reset element vector to zero
    _resetCellMatrix();

    // Get cell geometry information that depends on cell
    const scalar_array& basis = _quadrature->basis();
    const scalar_array& jacobianDet = _quadrature->jacobianDet();

    // Compute Jacobian for absorbing bc terms
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      const PylithScalar wt = 
        quadWts[iQuad] * jacobianDet[iQuad] / (2.0 * dt);
      for (int iBasis=0, iQ=iQuad*numBasis; iBasis < numBasis; ++iBasis) {
        const PylithScalar valI = wt*basis[iQ+iBasis];
        for (int jBasis=0; jBasis < numBasis; ++jBasis) {
          const PylithScalar valIJ = valI * basis[iQ+jBasis];
          for (int iDim=0; iDim < spaceDim; ++iDim) {
            const int iBlock = (iBasis*spaceDim + iDim) * (numBasis*spaceDim);
            const int jBlock = (jBasis*spaceDim + iDim);
            _cellMatrix[iBlock+jBlock] += valIJ * dampingConstsArray[doff+iQuad*spaceDim+iDim];
          } // for
        } // for
      } // for
    } // for
#if defined(DETAILED_EVENT_LOGGING)
    PetscLogFlops(numQuadPts*(3+numBasis*(1+numBasis*(1+2*spaceDim))));
    _logger->eventEnd(computeEvent);
    _logger->eventBegin(updateEvent);
#endif
    
    // Assemble cell contribution into PETSc Matrix
    err = DMPlexMatSetClosure(subMesh, solutionSubsection, solutionGlobalSubsection, jacobianMat, c, &_cellMatrix[0], ADD_VALUES);
    CHECK_PETSC_ERROR_MSG(err, "Update to PETSc Mat failed.");

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(updateEvent);
#endif
  } // for
  err = VecRestoreArray(dampingConstsVec, &dampingConstsArray);CHECK_PETSC_ERROR(err);
  err = PetscSectionDestroy(&solutionSubsection);CHECK_PETSC_ERROR(err);
  err = PetscSectionDestroy(&solutionGlobalSection);CHECK_PETSC_ERROR(err);
  err = PetscSectionDestroy(&solutionGlobalSubsection);CHECK_PETSC_ERROR(err);

#if !defined(DETAILED_EVENT_LOGGING)
  PetscLogFlops((cEnd-cStart)*numQuadPts*(3+numBasis*(1+numBasis*(1+2*spaceDim))));
  _logger->eventEnd(computeEvent);
#endif

  _needNewJacobian = false;
} // integrateJacobian

// ----------------------------------------------------------------------
// Integrate contributions to Jacobian matrix (A) associated with
void
pylith::bc::AbsorbingDampers::integrateJacobian(
			      topology::Field<topology::Mesh>* jacobian,
			      const PylithScalar t,
			      topology::SolutionFields* const fields)
{ // integrateJacobian
  assert(_quadrature);
  assert(_boundaryMesh);
  assert(_logger);
  assert(jacobian);
  assert(fields);

  const int setupEvent = _logger->eventId("AdIJ setup");
  const int geometryEvent = _logger->eventId("AdIJ geometry");
  const int computeEvent = _logger->eventId("AdIJ compute");
  const int restrictEvent = _logger->eventId("AdIJ restrict");
  const int updateEvent = _logger->eventId("AdIJ update");

  _logger->eventBegin(setupEvent);

  // Get cell geometry information that doesn't depend on cell
  const int numQuadPts = _quadrature->numQuadPts();
  const scalar_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == numQuadPts);
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();

  // Get cell information
  PetscDM subMesh = _boundaryMesh->dmMesh();
  PetscIS subpointIS;
  PetscInt cStart, cEnd;
  PetscErrorCode err;

  assert(subMesh);
  err = DMPlexGetHeightStratum(subMesh, 1, &cStart, &cEnd);CHECK_PETSC_ERROR(err);
  err = DMPlexCreateSubpointIS(subMesh, &subpointIS);CHECK_PETSC_ERROR(err);

  // Get parameters used in integration.
  const PylithScalar dt = _dt;
  assert(dt > 0);

  // Allocate matrix for cell values.
  _initCellMatrix();
  _initCellVector();

  // Get sections
  PetscSection dampingConstsSection = _parameters->get("damping constants").petscSection();
  PetscVec dampingConstsVec     = _parameters->get("damping constants").localVector();
  PetscScalar *dampingConstsArray;
  assert(dampingConstsSection);assert(dampingConstsVec);

  PetscSection jacobianSection = jacobian->petscSection(), jacobianSubsection;
  PetscVec jacobianVec     = jacobian->localVector();
  assert(jacobianSection);assert(jacobianVec);
  err = PetscSectionCreateSubmeshSection(jacobianSection, subpointIS, &jacobianSubsection);
  err = ISDestroy(&subpointIS);CHECK_PETSC_ERROR(err);

#if !defined(PRECOMPUTE_GEOMETRY)
  scalar_array coordinatesCell(numBasis*spaceDim);
  PetscSection coordSection;
  PetscVec coordVec;
  err = DMPlexGetCoordinateSection(subMesh, &coordSection);CHECK_PETSC_ERROR(err);
  err = DMGetCoordinatesLocal(subMesh, &coordVec);CHECK_PETSC_ERROR(err);
  assert(coordSection);assert(coordVec);
#endif

  _logger->eventEnd(setupEvent);
#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventBegin(computeEvent);
#endif

  err = VecGetArray(dampingConstsVec, &dampingConstsArray);CHECK_PETSC_ERROR(err);
  for(PetscInt c = cStart; c < cEnd; ++c) {
    // Compute geometry information for current cell
#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(geometryEvent);
#endif
#if defined(PRECOMPUTE_GEOMETRY)
    _quadrature->retrieveGeometry(*c_iter);
#else
    const PetscScalar *coordsArray;
    PetscInt           coordsSize;
    err = DMPlexVecGetClosure(subMesh, coordSection, coordVec, c, &coordsSize, &coordsArray);CHECK_PETSC_ERROR(err);
    for(PetscInt i = 0; i < coordsSize; ++i) {coordinatesCell[i] = coordsArray[i];}
    _quadrature->computeGeometry(coordinatesCell, c);
    err = DMPlexVecRestoreClosure(subMesh, coordSection, coordVec, c, &coordsSize, &coordsArray);CHECK_PETSC_ERROR(err);
#endif
#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(geometryEvent);
    _logger->eventBegin(restrictEvent);
#endif

    // Get damping constants
    PetscInt ddof, doff;

    err = PetscSectionGetDof(dampingConstsSection, c, &ddof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dampingConstsSection, c, &doff);CHECK_PETSC_ERROR(err);
    assert(ddof == numQuadPts*spaceDim);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(restrictEvent);
    _logger->eventBegin(computeEvent);
#endif

    // Reset element vector to zero
    _resetCellVector();

    // Get cell geometry information that depends on cell
    const scalar_array& basis = _quadrature->basis();
    const scalar_array& jacobianDet = _quadrature->jacobianDet();

    // Compute Jacobian for absorbing bc terms
    for (int iQuad = 0; iQuad < numQuadPts; ++iQuad) {
      const PylithScalar wt = quadWts[iQuad] * jacobianDet[iQuad] / (2.0 * dt);
      const int iQ = iQuad * numBasis;
      PylithScalar valJ = 0.0;
      for (int jBasis = 0; jBasis < numBasis; ++jBasis)
        valJ += basis[iQ + jBasis];
      valJ *= wt;
      for (int iBasis = 0; iBasis < numBasis; ++iBasis) {
        const PylithScalar valIJ = basis[iQ + iBasis] * valJ;
        for (int iDim = 0; iDim < spaceDim; ++iDim)
          _cellVector[iBasis * spaceDim + iDim] += valIJ
              * dampingConstsArray[doff+iQuad * spaceDim + iDim];
      } // for
    } // for
    err = DMPlexVecSetClosure(subMesh, jacobianSubsection, jacobianVec, c, &_cellVector[0], ADD_VALUES);CHECK_PETSC_ERROR(err);
#if defined(DETAILED_EVENT_LOGGING)
    PetscLogFlops(numQuadPts*(4+numBasis+numBasis*(1+spaceDim*2)));
    _logger->eventEnd(computeEvent);
    _logger->eventBegin(updateEvent);
#endif
#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(updateEvent);
#endif
  } // for
  err = VecRestoreArray(dampingConstsVec, &dampingConstsArray);CHECK_PETSC_ERROR(err);
  err = PetscSectionDestroy(&jacobianSubsection);CHECK_PETSC_ERROR(err);

#if !defined(DETAILED_EVENT_LOGGING)
  PetscLogFlops((cEnd-cStart)*numQuadPts*(4+numBasis+numBasis*(1+spaceDim*2)));
  _logger->eventEnd(computeEvent);
#endif

  _needNewJacobian = false;
} // integrateJacobian

// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::bc::AbsorbingDampers::verifyConfiguration(const topology::Mesh& mesh) const
{ // verifyConfiguration
  BCIntegratorSubMesh::verifyConfiguration(mesh);
} // verifyConfiguration

// ----------------------------------------------------------------------
// Initialize logger.
void
pylith::bc::AbsorbingDampers::_initializeLogger(void)
{ // initializeLogger
  delete _logger; _logger = new utils::EventLogger;
  assert(_logger);
  _logger->className("AbsorbingDampers");
  _logger->initialize();

  _logger->registerEvent("AdIR setup");
  _logger->registerEvent("AdIR geometry");
  _logger->registerEvent("AdIR compute");
  _logger->registerEvent("AdIR restrict");
  _logger->registerEvent("AdIR update");

  _logger->registerEvent("AdIJ setup");
  _logger->registerEvent("AdIJ geometry");
  _logger->registerEvent("AdIJ compute");
  _logger->registerEvent("AdIJ restrict");
  _logger->registerEvent("AdIJ update");
} // initializeLogger


// End of file 
