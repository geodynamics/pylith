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
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/topology/VisitorSubMesh.hh" // USES VecVisitorSubMesh, MatVisitorSubMesh
#include "pylith/topology/Stratum.hh" // USES Stratum

#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <cstring> // USES memcpy()
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

//#define PRECOMPUTE_GEOMETRY
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
  const PetscDM dmSubMesh = _boundaryMesh->dmMesh();assert(dmSubMesh);
  topology::Stratum heightStratum(dmSubMesh, topology::Stratum::HEIGHT, 1);
  const PetscInt cStart = heightStratum.begin();
  const PetscInt cEnd = heightStratum.end();

  // Get damping constants at each quadrature point and rotate to
  // global coordinate frame using orientation information
  const feassemble::CellGeometry& cellGeometry = _quadrature->refGeometry();
  const int cellDim = _quadrature->cellDim() > 0 ? _quadrature->cellDim() : 1;
  const int numBasis = _quadrature->numBasis();
  const int numQuadPts = _quadrature->numQuadPts();
  const int spaceDim = cellGeometry.spaceDim();
  const int fiberDim = numQuadPts * spaceDim;

  delete _parameters;
  _parameters = new topology::Fields<topology::Field<topology::SubMesh> >(*_boundaryMesh);
  assert(_parameters);
  _parameters->add("damping constants", "damping_constants", topology::FieldBase::FACES_FIELD, fiberDim);
  _parameters->get("damping constants").allocate();

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
  topology::Field<topology::SubMesh>& dampingConsts = _parameters->get("damping constants");
  topology::VecVisitorMesh dampingConstsVisitor(dampingConsts);

  scalar_array coordinatesCell(numBasis*spaceDim);
  PetscScalar *coordsCell = NULL;
  PetscInt coordsSize = 0;
  topology::CoordsVisitor coordsVisitor(dmSubMesh);

  assert(_normalizer);
  const PylithScalar lengthScale = _normalizer->lengthScale();
  const PylithScalar densityScale = _normalizer->densityScale();
  assert(_normalizer->timeScale() > 0);
  const PylithScalar velocityScale = _normalizer->lengthScale() / _normalizer->timeScale();

  const spatialdata::geocoords::CoordSys* cs = _boundaryMesh->coordsys();

  // Compute quadrature information
  _quadrature->initializeGeometry();

  PetscScalar* dampingConstsArray = dampingConstsVisitor.localArray();

  for(PetscInt c = cStart; c < cEnd; ++c) {
    // Compute geometry information for current cell
    coordsVisitor.getClosure(&coordsCell, &coordsSize, c);
    for (int i=0; i < coordsSize; ++i) { // :TODO: Remove copy
      coordinatesCell[i] = coordsCell[i];
    } // for
    _quadrature->computeGeometry(coordsCell, coordsSize, c);
    coordsVisitor.restoreClosure(&coordsCell, &coordsSize, c);

    const PetscInt doff = dampingConstsVisitor.sectionOffset(c);
    assert(fiberDim == dampingConstsVisitor.sectionDof(c));

    const scalar_array& quadPtsNondim = _quadrature->quadPts();
    const scalar_array& quadPtsRef = _quadrature->quadPtsRef();
    quadPtsGlobal = quadPtsNondim;
    _normalizer->dimensionalize(&quadPtsGlobal[0], quadPtsGlobal.size(), lengthScale);

    for(int iQuad = 0; iQuad < numQuadPts; ++iQuad) {
      // Compute damping constants in normal/tangential coordinates
      const int err = _db->query(&queryData[0], numValues, &quadPtsGlobal[iQuad*spaceDim], spaceDim, cs);
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
      for (int iDim=0; iDim < numTangential; ++iDim) {
        dampingConstsLocal[iDim] = constTangential;
      } // for
      dampingConstsLocal[spaceDim-1] = constNormal;

      // Compute normal/tangential orientation
      memcpy(&quadPtRef[0], &quadPtsRef[iQuad*cellDim], cellDim*sizeof(PylithScalar));
      cellGeometry.jacobian(&jacobian, &jacobianDet, coordinatesCell, quadPtRef);
      cellGeometry.orientation(&orientation, jacobian, jacobianDet, up);
      assert(jacobianDet > 0.0);
      orientation /= jacobianDet;

      for (int iDim=0; iDim < spaceDim; ++iDim) {
        dampingConstsArray[doff+iQuad*spaceDim+iDim] = 0.0;
        for (int jDim=0; jDim < spaceDim; ++jDim) {
          dampingConstsArray[doff+iQuad*spaceDim+iDim] += dampingConstsLocal[jDim]*orientation[jDim*spaceDim+iDim];
	} // for
        // Ensure damping constants are positive
        dampingConstsArray[doff+iQuad*spaceDim+iDim] = fabs(dampingConstsArray[doff+iQuad*spaceDim+iDim]);
      } // for
    } // for
  } // for

  _db->close();
} // initialize

// ----------------------------------------------------------------------
// Integrate contributions to residual term (r) for operator.
void
pylith::bc::AbsorbingDampers::integrateResidual(const topology::Field<topology::Mesh>& residual,
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

  // Get sections
  topology::Field<topology::SubMesh>& dampingConsts = _parameters->get("damping constants");
  topology::VecVisitorMesh dampingConstsVisitor(dampingConsts);
  PetscScalar* dampingConstsArray = dampingConstsVisitor.localArray();

  // Get subsections
  const PetscDM dmSubMesh = _boundaryMesh->dmMesh();assert(dmSubMesh);
  topology::SubMeshIS submeshIS(*_boundaryMesh);

  // Use _cellVector for cell residual.
  topology::VecVisitorSubMesh residualVisitor(residual, submeshIS);

  PetscScalar *velocityArray = NULL;
  PetscInt velocitySize = 0;
  topology::VecVisitorSubMesh velocityVisitor(fields->get("velocity(t)"), submeshIS);

  submeshIS.deallocate();
  
#if !defined(PRECOMPUTE_GEOMETRY)
  PetscScalar* coordsCell = NULL;
  PetscInt coordsSize = 0;
topology::CoordsVisitor coordsVisitor(dmSubMesh);
#endif

  // Get 'surface' cells (1 dimension lower than top-level cells)
  topology::Stratum heightStratum(dmSubMesh, topology::Stratum::HEIGHT, 1);
  const PetscInt cStart = heightStratum.begin();
  const PetscInt cEnd = heightStratum.end();

  _logger->eventEnd(setupEvent);
#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventBegin(computeEvent);
#endif

  for (PetscInt c = cStart; c < cEnd; ++c) {
    // Get geometry information for current cell
#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(geometryEvent);
#endif
#if defined(PRECOMPUTE_GEOMETRY)
#error("Code for PRECOMPUTE_GEOMETRY not implemented.")
#else
    coordsVisitor.getClosure(&coordsCell, &coordsSize, c);
    _quadrature->computeGeometry(coordsCell, coordsSize, c);
    coordsVisitor.restoreClosure(&coordsCell, &coordsSize, c);
#endif

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(geometryEvent);
    _logger->eventBegin(restrictEvent);
#endif

    // Reset element vector to zero
    _resetCellVector();

    // Restrict input fields to cell
    velocityVisitor.getClosure(&velocityArray, &velocitySize, c);
    assert(velocitySize == numBasis*spaceDim);

    const PetscInt doff = dampingConstsVisitor.sectionOffset(c);
    assert(numQuadPts*spaceDim == dampingConstsVisitor.sectionDof(c));

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
              valIJ * velocityArray[jBasis*spaceDim+iDim];
        } // for
      } // for
    } // for
    velocityVisitor.restoreClosure(&velocityArray, &velocitySize, c);

    residualVisitor.setClosure(&_cellVector[0], _cellVector.size(), c, ADD_VALUES);

#if defined(DETAILED_EVENT_LOGGING)
    PetscLogFlops(numQuadPts*(1+numBasis*(1+numBasis*(3*spaceDim))));
    _logger->eventEnd(computeEvent);
    _logger->eventBegin(updateEvent);
#endif

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(updateEvent);
#endif
  } // for

#if !defined(DETAILED_EVENT_LOGGING)
  PetscLogFlops((cEnd-cStart)*numQuadPts*(1+numBasis*(1+numBasis*(3*spaceDim))));
  _logger->eventEnd(computeEvent);
#endif
} // integrateResidual

// ----------------------------------------------------------------------
// Integrate contributions to residual term (r) for operator.
void
pylith::bc::AbsorbingDampers::integrateResidualLumped(const topology::Field<topology::Mesh>& residual,
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

  // Get 'surface' cells (1 dimension lower than top-level cells)
  const PetscDM dmSubMesh = _boundaryMesh->dmMesh();assert(dmSubMesh);
  topology::Stratum heightStratum(dmSubMesh, topology::Stratum::HEIGHT, 1);
  const PetscInt cStart = heightStratum.begin();
  const PetscInt cEnd = heightStratum.end();

  // Get sections
  topology::Field<topology::SubMesh>& dampingConsts = _parameters->get("damping constants");
  topology::VecVisitorMesh dampingConstsVisitor(dampingConsts);
  PetscScalar* dampingConstsArray = dampingConstsVisitor.localArray();

  // Get subsections
  topology::SubMeshIS submeshIS(*_boundaryMesh);
  // Use _cellVector for cell values.
  topology::VecVisitorSubMesh residualVisitor(residual, submeshIS);

  PetscScalar *velocityArray = NULL;
  PetscInt velocitySize;
  topology::VecVisitorSubMesh velocityVisitor(fields->get("velocity(t)"), submeshIS);

  submeshIS.deallocate();

#if !defined(PRECOMPUTE_GEOMETRY)
  PetscScalar *coordsCell = NULL;
  PetscInt coordsSize = 0;
  topology::CoordsVisitor coordsVisitor(dmSubMesh);
#endif

  _logger->eventEnd(setupEvent);
#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventBegin(computeEvent);
#endif

  for (PetscInt c=cStart; c < cEnd; ++c) {
    // Get geometry information for current cell
#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(geometryEvent);
#endif
#if defined(PRECOMPUTE_GEOMETRY)
#error("Code for PRECOMPUTE_GEOMETRY not implemented.")
#else
    coordsVisitor.getClosure(&coordsCell, &coordsSize, c);
    _quadrature->computeGeometry(coordsCell, coordsSize, c);
    coordsVisitor.restoreClosure(&coordsCell, &coordsSize, c);
#endif

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(geometryEvent);
    _logger->eventBegin(restrictEvent);
#endif

    // Reset element vector to zero
    _resetCellVector();

    // Restrict input fields to cell
    velocityVisitor.getClosure(&velocityArray, &velocitySize, c);
    assert(velocitySize == numBasis*spaceDim);

    const PetscInt doff = dampingConstsVisitor.sectionOffset(c);
    assert(numQuadPts*spaceDim == dampingConstsVisitor.sectionDof(c));

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
            valIJ * velocityArray[iBasis*spaceDim+iDim];
      } // for
    } // for
    velocityVisitor.restoreClosure(&velocityArray, &velocitySize, c);

    residualVisitor.setClosure(&_cellVector[0], _cellVector.size(), c, ADD_VALUES);

#if defined(DETAILED_EVENT_LOGGING)
    PetscLogFlops(numQuadPts*(1+numBasis+numBasis*(1+spaceDim*3)));
    _logger->eventEnd(computeEvent);
    _logger->eventBegin(updateEvent);
#endif

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(updateEvent);
#endif
  } // for

#if !defined(DETAILED_EVENT_LOGGING)
  PetscLogFlops((cEnd-cStart)*numQuadPts*(1+numBasis+numBasis*(1+spaceDim*3)));
  _logger->eventEnd(computeEvent);
#endif
} // integrateResidualLumped

// ----------------------------------------------------------------------
// Integrate contributions to Jacobian matrix (A) associated with
void
pylith::bc::AbsorbingDampers::integrateJacobian(topology::Jacobian* jacobian,
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

  // Get 'surface' cells (1 dimension lower than top-level cells)
  const PetscDM dmSubMesh = _boundaryMesh->dmMesh();assert(dmSubMesh);
  topology::Stratum heightStratum(dmSubMesh, topology::Stratum::HEIGHT, 1);
  const PetscInt cStart = heightStratum.begin();
  const PetscInt cEnd = heightStratum.end();

  // Get sections
  topology::Field<topology::SubMesh>& dampingConsts = _parameters->get("damping constants");
  topology::VecVisitorMesh dampingConstsVisitor(dampingConsts);
  PetscScalar* dampingConstsArray = dampingConstsVisitor.localArray();

  // Get sparse matrix
  const topology::Field<topology::Mesh>& solution = fields->solution();
  topology::SubMeshIS submeshIS(*_boundaryMesh);
  const PetscMat jacobianMat = jacobian->matrix();assert(jacobianMat);
  topology::MatVisitorSubMesh jacobianVisitor(jacobianMat, solution, submeshIS);
  submeshIS.deallocate();

  // Get parameters used in integration.
  const PylithScalar dt = _dt;
  assert(dt > 0);

  // Allocate matrix for cell values.
  _initCellMatrix();

#if !defined(PRECOMPUTE_GEOMETRY)
  PetscScalar *coordsCell = NULL;
  PetscInt coordsSize = 0;
  topology::CoordsVisitor coordsVisitor(dmSubMesh);
#endif

  _logger->eventEnd(setupEvent);
#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventBegin(computeEvent);
#endif

  for(PetscInt c = cStart; c < cEnd; ++c) {
    // Compute geometry information for current cell
#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(geometryEvent);
#endif
#if defined(PRECOMPUTE_GEOMETRY)
#error("Code for PRECOMPUTE_GEOMETRY not implemented")
#else
    coordsVisitor.getClosure(&coordsCell, &coordsSize, c);
    _quadrature->computeGeometry(coordsCell, coordsSize, c);
    coordsVisitor.restoreClosure(&coordsCell, &coordsSize, c);
#endif
#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(geometryEvent);
    _logger->eventBegin(restrictEvent);
#endif

    // Get damping constants
    const PetscInt doff = dampingConstsVisitor.sectionOffset(c);
    assert(numQuadPts*spaceDim == dampingConstsVisitor.sectionDof(c));

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
    jacobianVisitor.setClosure(&_cellMatrix[0], _cellMatrix.size(), c, ADD_VALUES);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(updateEvent);
#endif
  } // for

#if !defined(DETAILED_EVENT_LOGGING)
  PetscLogFlops((cEnd-cStart)*numQuadPts*(3+numBasis*(1+numBasis*(1+2*spaceDim))));
  _logger->eventEnd(computeEvent);
#endif

  _needNewJacobian = false;
} // integrateJacobian

// ----------------------------------------------------------------------
// Integrate contributions to Jacobian matrix (A) associated with
void
pylith::bc::AbsorbingDampers::integrateJacobian(topology::Field<topology::Mesh>* jacobian,
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

  // Get 'surface' cells (1 dimension lower than top-level cells)
  const PetscDM dmSubMesh = _boundaryMesh->dmMesh();assert(dmSubMesh);
  topology::Stratum heightStratum(dmSubMesh, topology::Stratum::HEIGHT, 1);
  const PetscInt cStart = heightStratum.begin();
  const PetscInt cEnd = heightStratum.end();

  // Get parameters used in integration.
  const PylithScalar dt = _dt;
  assert(dt > 0);

  // Allocate matrix for cell values.
  _initCellMatrix();
  _initCellVector();

  // Get sections
  topology::Field<topology::SubMesh>& dampingConsts = _parameters->get("damping constants");
  topology::VecVisitorMesh dampingConstsVisitor(dampingConsts);
  PetscScalar* dampingConstsArray = dampingConstsVisitor.localArray();

  topology::SubMeshIS submeshIS(*_boundaryMesh);
  topology::VecVisitorSubMesh jacobianVisitor(*jacobian, submeshIS);

  submeshIS.deallocate();
  
#if !defined(PRECOMPUTE_GEOMETRY)
  PetscScalar* coordsCell = NULL;
  PetscInt coordsSize = 0;
  topology::CoordsVisitor coordsVisitor(dmSubMesh);
#endif

  _logger->eventEnd(setupEvent);
#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventBegin(computeEvent);
#endif

  for(PetscInt c = cStart; c < cEnd; ++c) {
    // Compute geometry information for current cell
#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(geometryEvent);
#endif
#if defined(PRECOMPUTE_GEOMETRY)
#error("Code for PRECOMPUTE_GEOMETRY not implemented");
#else
    coordsVisitor.getClosure(&coordsCell, &coordsSize, c);
    _quadrature->computeGeometry(coordsCell, coordsSize, c);
    coordsVisitor.restoreClosure(&coordsCell, &coordsSize, c);
#endif
#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(geometryEvent);
    _logger->eventBegin(restrictEvent);
#endif

    // Get damping constants
    const PetscInt doff = dampingConstsVisitor.sectionOffset(c);
    assert(numQuadPts*spaceDim == dampingConstsVisitor.sectionDof(c));

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

    jacobianVisitor.setClosure(&_cellVector[0], _cellVector.size(), c, ADD_VALUES);

#if defined(DETAILED_EVENT_LOGGING)
    PetscLogFlops(numQuadPts*(4+numBasis+numBasis*(1+spaceDim*2)));
    _logger->eventEnd(computeEvent);
    _logger->eventBegin(updateEvent);
#endif
#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(updateEvent);
#endif
  } // for

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
