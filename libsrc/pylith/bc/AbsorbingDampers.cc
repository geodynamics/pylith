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

#include "pylith/topology/FieldsNew.hh" // HOLDSA FieldsNew
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

//#define PRECOMPUTE_GEOMETRY
//#define DETAILED_EVENT_LOGGING

// ----------------------------------------------------------------------
typedef pylith::topology::SubMesh::SieveMesh SieveSubMesh;
typedef pylith::topology::SubMesh::RealUniformSection SubRealUniformSection;
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::RealSection RealSection;
typedef pylith::topology::Field<pylith::topology::SubMesh>::RestrictVisitor RestrictVisitor;
typedef pylith::topology::Field<pylith::topology::SubMesh>::UpdateAddVisitor UpdateAddVisitor;
typedef ALE::ISieveVisitor::IndicesVisitor<RealSection,SieveMesh::order_type,PylithInt> IndicesVisitor;

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
  _db = 0; // :TODO: Use shared pointer
} // deallocate
  
// ----------------------------------------------------------------------
// Initialize boundary condition. Determine orienation and compute traction
// vector at integration points.
void
pylith::bc::AbsorbingDampers::initialize(const topology::Mesh& mesh,
					 const PylithScalar upDir[3])
{ // initialize
  assert(0 != _boundaryMesh);
  assert(0 != _quadrature);
  assert(0 != _db);

  _initializeLogger();

  scalar_array up(3);
  for (int i=0; i < 3; ++i)
    up[i] = upDir[i];

  const int numCorners = _quadrature->numBasis();

  // Get 'surface' cells (1 dimension lower than top-level cells)
  const ALE::Obj<SieveSubMesh>& sieveSubMesh = _boundaryMesh->sieveMesh();
  assert(!sieveSubMesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& cells = 
    sieveSubMesh->heightStratum(1);
  assert(!cells.isNull());
  const SieveSubMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveSubMesh::label_sequence::iterator cellsEnd = cells->end();

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
  _parameters = 
    new topology::FieldsNew<topology::SubMesh>(*_boundaryMesh);
  assert(0 != _parameters);
  _parameters->add("damping constants", "damping_constants", fiberDim);
  _parameters->allocate(cells);

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
  scalar_array dampingConstsLocal(fiberDim);
  scalar_array dampingConstsGlobal(fiberDim);

  scalar_array coordinatesCell(numCorners*spaceDim);
  const ALE::Obj<RealSection>& coordinates = 
    sieveSubMesh->getRealSection("coordinates");
  assert(!coordinates.isNull());
  RestrictVisitor coordsVisitor(*coordinates, 
				coordinatesCell.size(), &coordinatesCell[0]);

  assert(0 != _normalizer);
  const PylithScalar lengthScale = _normalizer->lengthScale();
  const PylithScalar densityScale = _normalizer->densityScale();
  assert(_normalizer->timeScale() > 0);
  const PylithScalar velocityScale = 
    _normalizer->lengthScale() / _normalizer->timeScale();

  const ALE::Obj<SubRealUniformSection>& parametersSection = 
    _parameters->section();
  assert(!parametersSection.isNull());

  const spatialdata::geocoords::CoordSys* cs = _boundaryMesh->coordsys();

  // Compute quadrature information
  _quadrature->initializeGeometry();
#if defined(PRECOMPUTE_GEOMETRY)
  _quadrature->computeGeometry(*_boundaryMesh, cells);
#endif

  for(SieveSubMesh::label_sequence::iterator c_iter = cellsBegin;
      c_iter != cellsEnd;
      ++c_iter) {
    // Compute geometry information for current cell
#if defined(PRECOMPUTE_GEOMETRY)
    _quadrature->retrieveGeometry(*c_iter);
#else
    coordsVisitor.clear();
    sieveSubMesh->restrictClosure(*c_iter, coordsVisitor);
    _quadrature->computeGeometry(coordinatesCell, *c_iter);
#endif

    const scalar_array& quadPtsNondim = _quadrature->quadPts();
    const scalar_array& quadPtsRef = _quadrature->quadPtsRef();
    quadPtsGlobal = quadPtsNondim;
    _normalizer->dimensionalize(&quadPtsGlobal[0], quadPtsGlobal.size(), 
				lengthScale);

    dampingConstsGlobal = 0.0;
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
      const PylithScalar densityN = 
	_normalizer->nondimensionalize(queryData[0], densityScale);
      const PylithScalar vpN = 
	_normalizer->nondimensionalize(queryData[1], velocityScale);
      const PylithScalar vsN = (3 == numValues) ?
	_normalizer->nondimensionalize(queryData[2], velocityScale) :
	0.0;
      
      const PylithScalar constTangential = densityN * vsN;
      const PylithScalar constNormal = densityN * vpN;
      const int numTangential = spaceDim-1;
      for (int iDim=0; iDim < numTangential; ++iDim)
	dampingConstsLocal[iDim] = constTangential;
      dampingConstsLocal[spaceDim-1] = constNormal;

      // Compute normal/tangential orientation
      memcpy(&quadPtRef[0], &quadPtsRef[iQuad*cellDim], 
	     cellDim*sizeof(PylithScalar));
#if defined(PRECOMPUTE_GEOMETRY)
      coordsVisitor.clear();
      sieveSubMesh->restrictClosure(*c_iter, coordsVisitor);
#endif
      cellGeometry.jacobian(&jacobian, &jacobianDet, coordinatesCell, quadPtRef);
      cellGeometry.orientation(&orientation, jacobian, jacobianDet, up);
      assert(jacobianDet > 0.0);
      orientation /= jacobianDet;

      for (int iDim=0; iDim < spaceDim; ++iDim) {
	for (int jDim=0; jDim < spaceDim; ++jDim)
	  dampingConstsGlobal[iQuad*spaceDim+iDim] += 
	    dampingConstsLocal[jDim]*orientation[jDim*spaceDim+iDim];
	// Ensure damping constants are positive
	dampingConstsGlobal[iQuad*spaceDim+iDim] = 
	  fabs(dampingConstsGlobal[iQuad*spaceDim+iDim]);
      } // for
    } // for
    parametersSection->updatePoint(*c_iter, &dampingConstsGlobal[0]);
  } // for

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
  assert(0 != _quadrature);
  assert(0 != _boundaryMesh);
  assert(0 != _parameters);
  assert(0 != fields);
  assert(0 != _logger);

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
  const ALE::Obj<SieveSubMesh>& sieveSubMesh = _boundaryMesh->sieveMesh();
  assert(!sieveSubMesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& cells = 
    sieveSubMesh->heightStratum(1);
  assert(!cells.isNull());
  const SieveSubMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveSubMesh::label_sequence::iterator cellsEnd = cells->end();

  // Get sections
  const ALE::Obj<SubRealUniformSection>& parametersSection = 
    _parameters->section();
  assert(!parametersSection.isNull());

  // Use _cellVector for cell residual.
  PetscSection residualSection = residual.petscSection();
  Vec          residualVec     = residual.localVector();
  assert(residualSection);assert(residualVec);
  UpdateAddVisitor residualVisitor(*residualSection, &_cellVector[0]);
  
  scalar_array velCell(numBasis*spaceDim);
  PetscSection velSection = fields->get("velocity(t)").petscSection();
  PetscSection velVec     = fields->get("velocity(t)").localVector();
  assert(velSection);assert(velVec);
  RestrictVisitor velVisitor(*velSection, velCell.size(), &velCell[0]);
  
#if !defined(PRECOMPUTE_GEOMETRY)
  scalar_array coordinatesCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& coordinates = 
    sieveSubMesh->getRealSection("coordinates");
  RestrictVisitor coordsVisitor(*coordinates,
				coordinatesCell.size(), &coordinatesCell[0]);
#endif

  _logger->eventEnd(setupEvent);
#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventBegin(computeEvent);
#endif

  for (SieveSubMesh::label_sequence::iterator c_iter=cellsBegin;
       c_iter != cellsEnd;
       ++c_iter) {
    // Get geometry information for current cell
#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(geometryEvent);
#endif
#if defined(PRECOMPUTE_GEOMETRY)
    _quadrature->retrieveGeometry(*c_iter);
#else
    coordsVisitor.clear();
    sieveSubMesh->restrictClosure(*c_iter, coordsVisitor);
    _quadrature->computeGeometry(coordinatesCell, *c_iter);
#endif

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(geometryEvent);
    _logger->eventBegin(restrictEvent);
#endif

    // Reset element vector to zero
    _resetCellVector();

    // Restrict input fields to cell
    velVisitor.clear();
    sieveSubMesh->restrictClosure(*c_iter, velVisitor);
    assert(numQuadPts*spaceDim == 
	   parametersSection->getFiberDimension(*c_iter));
    const PylithScalar* dampersCell = parametersSection->restrictPoint(*c_iter);

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
	      dampersCell[iQuad*spaceDim+iDim] *
	      valIJ * velCell[jBasis*spaceDim+iDim];
        } // for
      } // for
    } // for

#if defined(DETAILED_EVENT_LOGGING)
    PetscLogFlops(numQuadPts*(1+numBasis*(1+numBasis*(3*spaceDim))));
    _logger->eventEnd(computeEvent);
    _logger->eventBegin(updateEvent);
#endif

    // Assemble cell contribution into field
    residualVisitor.clear();
    sieveSubMesh->updateClosure(*c_iter, residualVisitor);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(updateEvent);
#endif
  } // for

#if !defined(DETAILED_EVENT_LOGGING)
  PetscLogFlops(cells->size()*numQuadPts*(1+numBasis*(1+numBasis*(3*spaceDim))));
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
  assert(0 != _quadrature);
  assert(0 != _boundaryMesh);
  assert(0 != _parameters);
  assert(0 != fields);
  assert(0 != _logger);

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
  const ALE::Obj<SieveSubMesh>& sieveSubMesh = _boundaryMesh->sieveMesh();
  assert(!sieveSubMesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& cells =
    sieveSubMesh->heightStratum(1);
  assert(!cells.isNull());
  const SieveSubMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveSubMesh::label_sequence::iterator cellsEnd = cells->end();

  // Get sections
  const ALE::Obj<SubRealUniformSection>& parametersSection = 
    _parameters->section();
  assert(!parametersSection.isNull());

  // Use _cellVector for cell values.
  const ALE::Obj<RealSection>& residualSection = residual.section();
  assert(!residualSection.isNull());
  UpdateAddVisitor residualVisitor(*residualSection, &_cellVector[0]);

  scalar_array velCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& velSection = 
    fields->get("velocity(t)").section();
  assert(!velSection.isNull());
  RestrictVisitor velVisitor(*velSection, velCell.size(), &velCell[0]);

#if !defined(PRECOMPUTE_GEOMETRY)
  scalar_array coordinatesCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& coordinates =
    sieveSubMesh->getRealSection("coordinates");
  RestrictVisitor coordsVisitor(*coordinates,
				coordinatesCell.size(), &coordinatesCell[0]);
#endif

  _logger->eventEnd(setupEvent);
#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventBegin(computeEvent);
#endif

  for (SieveSubMesh::label_sequence::iterator c_iter=cellsBegin;
       c_iter != cellsEnd;
       ++c_iter) {
    // Get geometry information for current cell
#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(geometryEvent);
#endif
#if defined(PRECOMPUTE_GEOMETRY)
    _quadrature->retrieveGeometry(*c_iter);
#else
    coordsVisitor.clear();
    sieveSubMesh->restrictClosure(*c_iter, coordsVisitor);
    _quadrature->computeGeometry(coordinatesCell, *c_iter);
#endif

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(geometryEvent);
    _logger->eventBegin(restrictEvent);
#endif

    // Reset element vector to zero
    _resetCellVector();

    // Restrict input fields to cell
    velVisitor.clear();
    sieveSubMesh->restrictClosure(*c_iter, velVisitor);
    assert(numQuadPts*spaceDim == 
	   parametersSection->getFiberDimension(*c_iter));
    const PylithScalar* dampersCell = parametersSection->restrictPoint(*c_iter);

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
          _cellVector[iBasis*spaceDim+iDim] -= valIJ * 
	    dampersCell[iQuad*spaceDim+iDim] * velCell[iBasis*spaceDim+iDim];
      } // for
    } // for
#if defined(DETAILED_EVENT_LOGGING)
    PetscLogFlops(numQuadPts*(1+numBasis+numBasis*(1+spaceDim*3)));
    _logger->eventEnd(computeEvent);
    _logger->eventBegin(updateEvent);
#endif

    // Assemble cell contribution into field
    residualVisitor.clear();
    sieveSubMesh->updateClosure(*c_iter, residualVisitor);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(updateEvent);
#endif
  } // for
#if !defined(DETAILED_EVENT_LOGGING)
  PetscLogFlops(cells->size()*numQuadPts*(1+numBasis+numBasis*(1+spaceDim*3)));
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
  assert(0 != _quadrature);
  assert(0 != _boundaryMesh);
  assert(0 != _logger);
  assert(0 != jacobian);
  assert(0 != fields);

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
  const ALE::Obj<SieveSubMesh>& sieveSubMesh = _boundaryMesh->sieveMesh();
  assert(!sieveSubMesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& cells = 
    sieveSubMesh->heightStratum(1);
  assert(!cells.isNull());
  const SieveSubMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveSubMesh::label_sequence::iterator cellsEnd = cells->end();

  // Get sections
  const ALE::Obj<SubRealUniformSection>& parametersSection = 
    _parameters->section();
  assert(!parametersSection.isNull());

  const topology::Field<topology::Mesh>& solution = fields->solution();
  const ALE::Obj<SieveMesh>& sieveMesh = solution.mesh().sieveMesh();
  const ALE::Obj<RealSection>& solutionSection = solution.section();
  assert(!solutionSection.isNull());
  const ALE::Obj<SieveMesh::order_type>& globalOrder = 
    sieveMesh->getFactory()->getGlobalOrder(sieveMesh, "default", 
					    solutionSection);
  assert(!globalOrder.isNull());

  const ALE::Obj<SieveMesh::sieve_type>& sieve = sieveMesh->getSieve();
  assert(!sieve.isNull());
  const int closureSize = 
    int(pow(sieve->getMaxConeSize(), sieveMesh->depth()));
  IndicesVisitor jacobianVisitor(*solutionSection, *globalOrder,
				 closureSize*spaceDim);

  // Get sparse matrix
  const PetscMat jacobianMat = jacobian->matrix();
  assert(0 != jacobianMat);

  // Get parameters used in integration.
  const PylithScalar dt = _dt;
  assert(dt > 0);

  // Allocate matrix for cell values.
  _initCellMatrix();

#if !defined(PRECOMPUTE_GEOMETRY)
  scalar_array coordinatesCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& coordinates =
    sieveSubMesh->getRealSection("coordinates");
  RestrictVisitor coordsVisitor(*coordinates,
				coordinatesCell.size(), &coordinatesCell[0]);
#endif

  _logger->eventEnd(setupEvent);
#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventBegin(computeEvent);
#endif

  for (SieveSubMesh::label_sequence::iterator c_iter=cellsBegin;
       c_iter != cellsEnd;
       ++c_iter) {
    // Compute geometry information for current cell
#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(geometryEvent);
#endif
#if defined(PRECOMPUTE_GEOMETRY)
    _quadrature->retrieveGeometry(*c_iter);
#else
    coordsVisitor.clear();
    sieveSubMesh->restrictClosure(*c_iter, coordsVisitor);
    _quadrature->computeGeometry(coordinatesCell, *c_iter);
#endif
#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(geometryEvent);
    _logger->eventBegin(restrictEvent);
#endif

    // Get damping constants
    assert(numQuadPts*spaceDim == parametersSection->getFiberDimension(*c_iter));
    const PylithScalar* dampingConstsCell = parametersSection->restrictPoint(*c_iter);

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
            _cellMatrix[iBlock+jBlock] += 
	      valIJ * dampingConstsCell[iQuad*spaceDim+iDim];
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
    jacobianVisitor.clear();
    PetscErrorCode err = updateOperator(jacobianMat, *sieveSubMesh->getSieve(), 
					jacobianVisitor, *c_iter,
					&_cellMatrix[0], ADD_VALUES);
    CHECK_PETSC_ERROR_MSG(err, "Update to PETSc Mat failed.");

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(updateEvent);
#endif
  } // for

#if !defined(DETAILED_EVENT_LOGGING)
  PetscLogFlops(cells->size()*numQuadPts*(3+numBasis*(1+numBasis*(1+2*spaceDim))));
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
  assert(0 != _quadrature);
  assert(0 != _boundaryMesh);
  assert(0 != _logger);
  assert(0 != jacobian);
  assert(0 != fields);

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
  const ALE::Obj<SieveSubMesh>& sieveSubMesh = _boundaryMesh->sieveMesh();
  assert(!sieveSubMesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& cells = 
    sieveSubMesh->heightStratum(1);
  assert(!cells.isNull());
  const SieveSubMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveSubMesh::label_sequence::iterator cellsEnd = cells->end();

  // Get parameters used in integration.
  const PylithScalar dt = _dt;
  assert(dt > 0);

  // Allocate matrix for cell values.
  _initCellMatrix();
  _initCellVector();

  // Get sections
  const ALE::Obj<SubRealUniformSection>& parametersSection = 
    _parameters->section();
  assert(!parametersSection.isNull());

  const topology::Field<topology::Mesh>& solution = fields->solution();
  const ALE::Obj<SieveMesh>& sieveMesh = solution.mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<RealSection>& solutionSection = solution.section();
  assert(!solutionSection.isNull());

  const ALE::Obj<RealSection>& jacobianSection = jacobian->section();
  assert(!jacobianSection.isNull());
  UpdateAddVisitor jacobianVisitor(*jacobianSection, &_cellVector[0]);

#if !defined(PRECOMPUTE_GEOMETRY)
  scalar_array coordinatesCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& coordinates = 
    sieveSubMesh->getRealSection("coordinates");
  RestrictVisitor coordsVisitor(*coordinates,
				coordinatesCell.size(), &coordinatesCell[0]);
#endif

  _logger->eventEnd(setupEvent);
#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventBegin(computeEvent);
#endif

  for (SieveSubMesh::label_sequence::iterator c_iter=cellsBegin;
       c_iter != cellsEnd;
       ++c_iter) {
    // Compute geometry information for current cell
#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(geometryEvent);
#endif
#if defined(PRECOMPUTE_GEOMETRY)
    _quadrature->retrieveGeometry(*c_iter);
#else
    coordsVisitor.clear();
    sieveSubMesh->restrictClosure(*c_iter, coordsVisitor);
    _quadrature->computeGeometry(coordinatesCell, *c_iter);
#endif
#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(geometryEvent);
    _logger->eventBegin(restrictEvent);
#endif

    // Get damping constants
    assert(numQuadPts*spaceDim == parametersSection->getFiberDimension(*c_iter));
    const PylithScalar* dampingConstsCell = parametersSection->restrictPoint(*c_iter);

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
              * dampingConstsCell[iQuad * spaceDim + iDim];
      } // for
    } // for
#if defined(DETAILED_EVENT_LOGGING)
    PetscLogFlops(numQuadPts*(4+numBasis+numBasis*(1+spaceDim*2)));
    _logger->eventEnd(computeEvent);
    _logger->eventBegin(updateEvent);
#endif

    // Assemble cell contribution into lumped matrix.
    jacobianVisitor.clear();
    sieveSubMesh->updateClosure(*c_iter, jacobianVisitor);
#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(updateEvent);
#endif
  } // for

#if !defined(DETAILED_EVENT_LOGGING)
  PetscLogFlops(cells->size()*numQuadPts*(4+numBasis+numBasis*(1+spaceDim*2)));
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
  assert(0 != _logger);
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
