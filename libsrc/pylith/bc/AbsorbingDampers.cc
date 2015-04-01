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

#include "AbsorbingDampers.hh" // implementation of object methods

#include "pylith/topology/Fields.hh" // HOLDSA Fields
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields
#include "pylith/topology/Jacobian.hh" // USES Jacobian
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/topology/VisitorSubMesh.hh" // USES VecVisitorSubMesh, MatVisitorSubMesh
#include "pylith/topology/Stratum.hh" // USES Stratum

#include "pylith/feassemble/CellGeometry.hh" // USES CellGeometry
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature

#include "pylith/utils/EventLogger.hh" // USES EventLogger

#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

//#define DETAILED_EVENT_LOGGING

// ----------------------------------------------------------------------
// Default constructor.
pylith::bc::AbsorbingDampers::AbsorbingDampers(void) :
  _velocityVisitor(0),
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
  PYLITH_METHOD_BEGIN;

  delete _velocityVisitor; _velocityVisitor = 0;
  BCIntegratorSubMesh::deallocate();
  _db = 0; // :TODO: Use shared pointer

  PYLITH_METHOD_END;
} // deallocate
  
// ----------------------------------------------------------------------
// Initialize boundary condition. Determine orienation and compute traction
// vector at integration points.
void
pylith::bc::AbsorbingDampers::initialize(const topology::Mesh& mesh,
					 const PylithScalar upDir[3])
{ // initialize
  PYLITH_METHOD_BEGIN;

  assert(_boundaryMesh);
  assert(_quadrature);
  assert(_db);

  _initializeLogger();

  scalar_array up(3);
  for (int i=0; i < 3; ++i) {
    up[i] = upDir[i];
  } // for

  // Get 'surface' cells (1 dimension lower than top-level cells)
  const PetscDM dmSubMesh = _boundaryMesh->dmMesh();assert(dmSubMesh);
  topology::Stratum cellsStratum(dmSubMesh, topology::Stratum::HEIGHT, 1);
  const PetscInt cStart = cellsStratum.begin();
  const PetscInt cEnd = cellsStratum.end();

  // Get damping constants at each quadrature point and rotate to
  // global coordinate frame using orientation information
  const feassemble::CellGeometry& cellGeometry = _quadrature->refGeometry();
  const int cellDim = _quadrature->cellDim() > 0 ? _quadrature->cellDim() : 1;
  const int numBasis = _quadrature->numBasis();
  const int numQuadPts = _quadrature->numQuadPts();
  const int spaceDim = cellGeometry.spaceDim();
  const int fiberDim = numQuadPts * spaceDim;

  delete _parameters;
  _parameters = new topology::Fields(*_boundaryMesh);
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
  scalar_array quadPtsGlobal(numQuadPts*spaceDim);

  // Container for damping constants for current cell
  scalar_array dampingConstsLocal(spaceDim);
  topology::Field& dampingConsts = _parameters->get("damping constants");
  topology::VecVisitorMesh dampingConstsVisitor(dampingConsts);

  scalar_array coordsCell(numBasis*spaceDim); // :KULDGE: Update numBasis to numCorners after implementing higher order
  topology::CoordsVisitor coordsVisitor(dmSubMesh);

  assert(_normalizer);
  const PylithScalar lengthScale = _normalizer->lengthScale();
  const PylithScalar densityScale = _normalizer->densityScale();
  assert(_normalizer->timeScale() > 0);
  const PylithScalar velocityScale = _normalizer->lengthScale() / _normalizer->timeScale();

  const spatialdata::geocoords::CoordSys* cs = _boundaryMesh->coordsys();assert(cs);

  // Compute quadrature information
  _quadrature->initializeGeometry();

  // Optimize coordinate retrieval in closure
  topology::CoordsVisitor::optimizeClosure(dmSubMesh);

  PetscScalar* dampingConstsArray = dampingConstsVisitor.localArray();

  for(PetscInt c = cStart; c < cEnd; ++c) {
    // Compute geometry information for current cell
    coordsVisitor.getClosure(&coordsCell, c);
    _quadrature->computeGeometry(&coordsCell[0], coordsCell.size(), c);

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
        msg << ") for absorbing boundary condition '" << _label
            << "' using spatial database '" << _db->label() << "'.";
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
      cellGeometry.jacobian(&jacobian, &jacobianDet, &coordsCell[0], numBasis, spaceDim, &quadPtsRef[iQuad*cellDim], cellDim);
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

  PYLITH_METHOD_END;
} // initialize

// ----------------------------------------------------------------------
// Integrate contributions to residual term (r) for operator.
void
pylith::bc::AbsorbingDampers::integrateResidual(const topology::Field& residual,
						const PylithScalar t,
						topology::SolutionFields* const fields)
{ // integrateResidual
  PYLITH_METHOD_BEGIN;

  assert(_quadrature);
  assert(_boundaryMesh);
  assert(_parameters);
  assert(fields);
  assert(_logger);

  const int setupEvent = _logger->eventId("AdIR setup");
  const int computeEvent = _logger->eventId("AdIR compute");
#if defined(DETAILED_EVENT_LOGGING)
  const int geometryEvent = _logger->eventId("AdIR geometry");
  const int restrictEvent = _logger->eventId("AdIR restrict");
  const int updateEvent = _logger->eventId("AdIR update");
#endif

  _logger->eventBegin(setupEvent);

  // Get cell geometry information that doesn't depend on cell
  const int numQuadPts = _quadrature->numQuadPts();
  const scalar_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == size_t(numQuadPts));
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();

  // Allocate vectors for cell values.
  _initCellVector();

  // Get sections
  topology::Field& dampingConsts = _parameters->get("damping constants");
  topology::VecVisitorMesh dampingConstsVisitor(dampingConsts);
  PetscScalar* dampingConstsArray = dampingConstsVisitor.localArray();

  // Get subsections
  const PetscDM dmSubMesh = _boundaryMesh->dmMesh();assert(dmSubMesh);

  // Use _cellVector for cell residual.
  assert(_submeshIS);
  if (!_residualVisitor) {
    _residualVisitor = new topology::VecVisitorSubMesh(residual, *_submeshIS);assert(_residualVisitor);
  } // if
  if (!_velocityVisitor) {
    _velocityVisitor = new topology::VecVisitorSubMesh(fields->get("velocity(t)"), *_submeshIS);assert(_velocityVisitor);
  } // if
  scalar_array velocityCell(numBasis*spaceDim);
  
  scalar_array coordsCell(numBasis*spaceDim); // :KULDGE: Update numBasis to numCorners after implementing higher order
  topology::CoordsVisitor coordsVisitor(dmSubMesh);

  // Get 'surface' cells (1 dimension lower than top-level cells)
  topology::Stratum cellsStratum(dmSubMesh, topology::Stratum::HEIGHT, 1);
  const PetscInt cStart = cellsStratum.begin();
  const PetscInt cEnd = cellsStratum.end();

  _logger->eventEnd(setupEvent);
#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventBegin(computeEvent);
#endif

  for (PetscInt c = cStart; c < cEnd; ++c) {
    // Get geometry information for current cell
#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(geometryEvent);
#endif
    coordsVisitor.getClosure(&coordsCell, c);
    _quadrature->computeGeometry(&coordsCell[0], coordsCell.size(), c);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(geometryEvent);
    _logger->eventBegin(restrictEvent);
#endif

    // Reset element vector to zero
    _resetCellVector();

    // Restrict input fields to cell
    _velocityVisitor->getClosure(&velocityCell, c);

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
              valIJ * velocityCell[jBasis*spaceDim+iDim];
        } // for
      } // for
    } // for

    _residualVisitor->setClosure(&_cellVector[0], _cellVector.size(), c, ADD_VALUES);

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

  PYLITH_METHOD_END;
} // integrateResidual

// ----------------------------------------------------------------------
// Integrate contributions to residual term (r) for operator.
void
pylith::bc::AbsorbingDampers::integrateResidualLumped(const topology::Field& residual,
						      const PylithScalar t,
						      topology::SolutionFields* const fields)
{ // integrateResidualLumped
  PYLITH_METHOD_BEGIN;

  assert(_quadrature);
  assert(_boundaryMesh);
  assert(_parameters);
  assert(fields);
  assert(_logger);

  const int setupEvent = _logger->eventId("AdIR setup");
  const int computeEvent = _logger->eventId("AdIR compute");
#if defined(DETAILED_EVENT_LOGGING)
  const int geometryEvent = _logger->eventId("AdIR geometry");
  const int restrictEvent = _logger->eventId("AdIR restrict");
  const int updateEvent = _logger->eventId("AdIR update");
#endif

  _logger->eventBegin(setupEvent);

  // Get cell geometry information that doesn't depend on cell
  const int numQuadPts = _quadrature->numQuadPts();
  const scalar_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == size_t(numQuadPts));
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();

  // Allocate vectors for cell values.
  _initCellVector();

  // Get 'surface' cells (1 dimension lower than top-level cells)
  const PetscDM dmSubMesh = _boundaryMesh->dmMesh();assert(dmSubMesh);
  topology::Stratum cellsStratum(dmSubMesh, topology::Stratum::HEIGHT, 1);
  const PetscInt cStart = cellsStratum.begin();
  const PetscInt cEnd = cellsStratum.end();

  // Get sections
  topology::Field& dampingConsts = _parameters->get("damping constants");
  topology::VecVisitorMesh dampingConstsVisitor(dampingConsts);
  PetscScalar* dampingConstsArray = dampingConstsVisitor.localArray();

  // Get subsections
  // Use _cellVector for cell values.
  assert(_submeshIS);
  if (!_residualVisitor) {
    _residualVisitor = new topology::VecVisitorSubMesh(residual, *_submeshIS);assert(_residualVisitor);
  } // if
  if (!_velocityVisitor) {
    _velocityVisitor = new topology::VecVisitorSubMesh(fields->get("velocity(t)"), *_submeshIS);assert(_velocityVisitor);
  } // if
  scalar_array velocityCell(numBasis*spaceDim);

  scalar_array coordsCell(numBasis*spaceDim); // :KULDGE: Update numBasis to numCorners after implementing higher order
  topology::CoordsVisitor coordsVisitor(dmSubMesh);

  _logger->eventEnd(setupEvent);
#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventBegin(computeEvent);
#endif

  for (PetscInt c=cStart; c < cEnd; ++c) {
    // Get geometry information for current cell
#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(geometryEvent);
#endif
    coordsVisitor.getClosure(&coordsCell, c);
    _quadrature->computeGeometry(&coordsCell[0], coordsCell.size(), c);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(geometryEvent);
    _logger->eventBegin(restrictEvent);
#endif

    // Reset element vector to zero
    _resetCellVector();

    // Restrict input fields to cell
    _velocityVisitor->getClosure(&velocityCell, c);

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
            valIJ * velocityCell[iBasis*spaceDim+iDim];
      } // for
    } // for

    _residualVisitor->setClosure(&_cellVector[0], _cellVector.size(), c, ADD_VALUES);

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

  PYLITH_METHOD_END;
} // integrateResidualLumped

// ----------------------------------------------------------------------
// Integrate contributions to Jacobian matrix (A) associated with
void
pylith::bc::AbsorbingDampers::integrateJacobian(topology::Jacobian* jacobian,
						const PylithScalar t,
						topology::SolutionFields* const fields)
{ // integrateJacobian
  PYLITH_METHOD_BEGIN;

  assert(_quadrature);
  assert(_boundaryMesh);
  assert(_logger);
  assert(jacobian);
  assert(fields);

  const int setupEvent = _logger->eventId("AdIJ setup");
  const int computeEvent = _logger->eventId("AdIJ compute");
#if defined(DETAILED_EVENT_LOGGING)
  const int geometryEvent = _logger->eventId("AdIJ geometry");
  const int restrictEvent = _logger->eventId("AdIJ restrict");
  const int updateEvent = _logger->eventId("AdIJ update");
#endif

  _logger->eventBegin(setupEvent);

  // Get cell geometry information that doesn't depend on cell
  const int numQuadPts = _quadrature->numQuadPts();
  const scalar_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == size_t(numQuadPts));
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();

  // Get 'surface' cells (1 dimension lower than top-level cells)
  const PetscDM dmSubMesh = _boundaryMesh->dmMesh();assert(dmSubMesh);
  topology::Stratum cellsStratum(dmSubMesh, topology::Stratum::HEIGHT, 1);
  const PetscInt cStart = cellsStratum.begin();
  const PetscInt cEnd = cellsStratum.end();

  // Get sections
  topology::Field& dampingConsts = _parameters->get("damping constants");
  topology::VecVisitorMesh dampingConstsVisitor(dampingConsts);
  PetscScalar* dampingConstsArray = dampingConstsVisitor.localArray();

  // Get sparse matrix
  const topology::Field& solution = fields->solution();
  const PetscMat jacobianMat = jacobian->matrix();assert(jacobianMat);
  if (!_jacobianMatVisitor) {
    assert(_submeshIS);
    _jacobianMatVisitor = new topology::MatVisitorSubMesh(jacobianMat, solution, *_submeshIS);assert(_jacobianMatVisitor);
  } // if

  // Get parameters used in integration.
  const PylithScalar dt = _dt;
  assert(dt > 0);

  // Allocate matrix for cell values.
  _initCellMatrix();

  scalar_array coordsCell(numBasis*spaceDim); // :KULDGE: Update numBasis to numCorners after implementing higher order
  topology::CoordsVisitor coordsVisitor(dmSubMesh);

  _logger->eventEnd(setupEvent);
#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventBegin(computeEvent);
#endif

  for(PetscInt c = cStart; c < cEnd; ++c) {
    // Compute geometry information for current cell
#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(geometryEvent);
#endif
    coordsVisitor.getClosure(&coordsCell, c);
    _quadrature->computeGeometry(&coordsCell[0], coordsCell.size(), c);

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
    _jacobianMatVisitor->setClosure(&_cellMatrix[0], _cellMatrix.size(), c, ADD_VALUES);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(updateEvent);
#endif
  } // for

#if !defined(DETAILED_EVENT_LOGGING)
  PetscLogFlops((cEnd-cStart)*numQuadPts*(3+numBasis*(1+numBasis*(1+2*spaceDim))));
  _logger->eventEnd(computeEvent);
#endif

  _needNewJacobian = false;

  PYLITH_METHOD_END;
} // integrateJacobian

// ----------------------------------------------------------------------
// Integrate contributions to Jacobian matrix (A) associated with
void
pylith::bc::AbsorbingDampers::integrateJacobian(topology::Field* jacobian,
						const PylithScalar t,
						topology::SolutionFields* const fields)
{ // integrateJacobian
  PYLITH_METHOD_BEGIN;

  assert(_quadrature);
  assert(_boundaryMesh);
  assert(_logger);
  assert(jacobian);
  assert(fields);

  const int setupEvent = _logger->eventId("AdIJ setup");
  const int computeEvent = _logger->eventId("AdIJ compute");
#if defined(DETAILED_EVENT_LOGGING)
  const int geometryEvent = _logger->eventId("AdIJ geometry");
  const int restrictEvent = _logger->eventId("AdIJ restrict");
  const int updateEvent = _logger->eventId("AdIJ update");
#endif

  _logger->eventBegin(setupEvent);

  // Get cell geometry information that doesn't depend on cell
  const int numQuadPts = _quadrature->numQuadPts();
  const scalar_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == size_t(numQuadPts));
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();

  // Get 'surface' cells (1 dimension lower than top-level cells)
  const PetscDM dmSubMesh = _boundaryMesh->dmMesh();assert(dmSubMesh);
  topology::Stratum cellsStratum(dmSubMesh, topology::Stratum::HEIGHT, 1);
  const PetscInt cStart = cellsStratum.begin();
  const PetscInt cEnd = cellsStratum.end();

  // Get parameters used in integration.
  const PylithScalar dt = _dt;
  assert(dt > 0);

  // Allocate matrix for cell values.
  _initCellMatrix();
  _initCellVector();

  // Get sections
  topology::Field& dampingConsts = _parameters->get("damping constants");
  topology::VecVisitorMesh dampingConstsVisitor(dampingConsts);
  PetscScalar* dampingConstsArray = dampingConstsVisitor.localArray();

  if (!_jacobianMatVisitor) {
    assert(_submeshIS);
    _jacobianVecVisitor = new topology::VecVisitorSubMesh(*jacobian, *_submeshIS);assert(_jacobianVecVisitor);
  } // if
  
  scalar_array coordsCell(numBasis*spaceDim); // :KULDGE: Update numBasis to numCorners after implementing higher order
  topology::CoordsVisitor coordsVisitor(dmSubMesh);

  _logger->eventEnd(setupEvent);
#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventBegin(computeEvent);
#endif

  for(PetscInt c = cStart; c < cEnd; ++c) {
    // Compute geometry information for current cell
#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(geometryEvent);
#endif
    coordsVisitor.getClosure(&coordsCell, c);
    _quadrature->computeGeometry(&coordsCell[0], coordsCell.size(), c);

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

    _jacobianVecVisitor->setClosure(&_cellVector[0], _cellVector.size(), c, ADD_VALUES);

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

  PYLITH_METHOD_END;
} // integrateJacobian

// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::bc::AbsorbingDampers::verifyConfiguration(const topology::Mesh& mesh) const
{ // verifyConfiguration
  PYLITH_METHOD_BEGIN;

  BCIntegratorSubMesh::verifyConfiguration(mesh);

  PYLITH_METHOD_END;
} // verifyConfiguration

// ----------------------------------------------------------------------
// Initialize logger.
void
pylith::bc::AbsorbingDampers::_initializeLogger(void)
{ // initializeLogger
   PYLITH_METHOD_BEGIN;

 delete _logger; _logger = new utils::EventLogger;assert(_logger);
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

  PYLITH_METHOD_END;
} // initializeLogger


// End of file 
