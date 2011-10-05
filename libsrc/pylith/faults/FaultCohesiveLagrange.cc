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
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "FaultCohesiveLagrange.hh" // implementation of object methods

#include "CohesiveTopology.hh" // USES CohesiveTopology

#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/feassemble/CellGeometry.hh" // USES CellGeometry
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/SubMesh.hh" // USES SubMesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/Jacobian.hh" // USES Jacobian
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields
#include "pylith/utils/macrodefs.h" // USES CALL_MEMBER_FN

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB

#include <cmath> // USES pow(), sqrt()
#include <strings.h> // USES strcasecmp()
#include <cstring> // USES strlen()
#include <cstdlib> // USES atoi()
#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

#include <petscblaslapack.h> // USES svd and dgemm

//#define PRECOMPUTE_GEOMETRY
//#define DETAILED_EVENT_LOGGING

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::RealSection RealSection;
typedef pylith::topology::SubMesh::SieveMesh SieveSubMesh;
typedef pylith::topology::Field<pylith::topology::SubMesh>::RestrictVisitor RestrictVisitor;
typedef pylith::topology::Field<pylith::topology::SubMesh>::UpdateAddVisitor UpdateAddVisitor;
typedef ALE::ISieveVisitor::IndicesVisitor<RealSection,SieveSubMesh::order_type,PetscInt> IndicesVisitor;

// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::FaultCohesiveLagrange::FaultCohesiveLagrange(void)
{ // constructor
  _useLagrangeConstraints = true;
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::faults::FaultCohesiveLagrange::~FaultCohesiveLagrange(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::faults::FaultCohesiveLagrange::deallocate(void)
{ // deallocate
  FaultCohesive::deallocate();

} // deallocate

// ----------------------------------------------------------------------
// Initialize fault. Determine orientation and setup boundary
void
pylith::faults::FaultCohesiveLagrange::initialize(const topology::Mesh& mesh,
					     const double upDir[3])
{ // initialize
  assert(0 != upDir);
  assert(0 != _quadrature);
  assert(0 != _normalizer);

  _initializeLogger();

  const spatialdata::geocoords::CoordSys* cs = mesh.coordsys();
  assert(0 != cs);

  delete _faultMesh;
  _faultMesh = new topology::SubMesh();
  CohesiveTopology::createFaultParallel(_faultMesh, mesh, id(),
					_useLagrangeConstraints);
  _initializeCohesiveInfo(mesh);

  delete _fields;
  _fields = new topology::Fields<topology::Field<topology::SubMesh> >(
    *_faultMesh);

  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  //logger.stagePush("Fault");

  // Allocate dispRel field
  const ALE::Obj<SieveSubMesh>& faultSieveMesh = _faultMesh->sieveMesh();
  assert(!faultSieveMesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& vertices =
      faultSieveMesh->depthStratum(0);
  assert(!vertices.isNull());
  _fields->add("relative disp", "relative_disp");
  topology::Field<topology::SubMesh>& dispRel = _fields->get("relative disp");
  dispRel.newSection(vertices, cs->spaceDim());
  dispRel.allocate();
  dispRel.vectorFieldType(topology::FieldBase::VECTOR);
  dispRel.scale(_normalizer->lengthScale());

  const ALE::Obj<SieveSubMesh::label_sequence>& cells =
      faultSieveMesh->heightStratum(0);
  assert(!cells.isNull());
  const SieveSubMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveSubMesh::label_sequence::iterator cellsEnd = cells->end();
  _quadrature->initializeGeometry();
#if defined(PRECOMPUTE_GEOMETRY)
  _quadrature->computeGeometry(*_faultMesh, cells);
#endif

  _fields->add("distribution", "distribution", pylith::topology::FieldBase::CELLS_FIELD, 1);
  topology::Field<topology::SubMesh>& dist = _fields->get("distribution");
  dist.allocate();
  const ALE::Obj<RealSection>& distSection = dist.section();
  assert(!distSection.isNull());
  const double rank = (double) distSection->commRank();

  // Loop over cells in fault mesh, set distribution
  for (SieveSubMesh::label_sequence::iterator c_iter = cellsBegin; c_iter
      != cellsEnd; ++c_iter) {
    distSection->updatePoint(*c_iter, &rank);
  } // for

  // Compute orientation at vertices in fault mesh.
  _calcOrientation(upDir);

  //logger.stagePop();
} // initialize

// ----------------------------------------------------------------------
void
pylith::faults::FaultCohesiveLagrange::splitField(topology::Field<
    topology::Mesh>* field)
{ // splitField
  assert(0 != field);

  const ALE::Obj<RealSection>& section = field->section();
  assert(!section.isNull());
  if (0 == section->getNumSpaces())
    return; // Return if there are no fibrations

  const int spaceDim = field->mesh().dimension();
  const int fibrationLagrange = spaceDim;

  // Add space for Lagrange multipliers if it does not yet exist.
  if (spaceDim == section->getNumSpaces())
    section->addSpace();

  const int numVertices = _cohesiveVertices.size();
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_lagrange = _cohesiveVertices[iVertex].lagrange;
    assert(spaceDim == section->getFiberDimension(v_lagrange));
    // Reset displacement fibration fiber dimension to zero.
    for (int fibration=0; fibration < spaceDim; ++fibration)
      section->setFiberDimension(v_lagrange, 0, fibration);
    // Set Lagrange fibration fiber dimension.
    section->setFiberDimension(v_lagrange, spaceDim, fibrationLagrange);
  } // for
} // splitField

// ----------------------------------------------------------------------
// Integrate contribution of cohesive cells to residual term.
void
pylith::faults::FaultCohesiveLagrange::integrateResidual(
			 const topology::Field<topology::Mesh>& residual,
			 const double t,
			 topology::SolutionFields* const fields)
{ // integrateResidual
  assert(0 != fields);
  assert(0 != _fields);
  assert(0 != _logger);

  // Cohesive cells with conventional vertices N and P, and constraint
  // vertex L make contributions to the assembled residual:
  //
  // DOF P: \int_{S_f^+} \tensor{N}_m^T \cdot \tensor{N}_p \cdot \vec{l}_p dS
  // DOF N: -\int_{S_f^+} \tensor{N}_m^T \cdot \tensor{N}_p \cdot \vec{l}_p dS
  // DOF L: \int_S_f \tensor{N}_p^T ( \tensor{R} \cdot \vec{d} 
  //                 -\tensor{N}_{n^+} \cdot \vec{u}_{n^+}
  //                 +\tensor{N}_{n^-} \cdot \vec{u}_{n^-} dS

  // ASK MATT: Why does fault mesh not use Lagrange vertices?

  const int setupEvent = _logger->eventId("FaIR setup");
  const int geometryEvent = _logger->eventId("FaIR geometry");
  const int computeEvent = _logger->eventId("FaIR compute");
  const int restrictEvent = _logger->eventId("FaIR restrict");
  const int updateEvent = _logger->eventId("FaIR update");

  _logger->eventBegin(setupEvent);

  // Get cell geometry information that doesn't depend on cell
  const int numQuadPts = _quadrature->numQuadPts();
  const double_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == numQuadPts);
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int cellDim = _quadrature->cellDim();
  assert(cellDim == spaceDim-1);

  // Get cohesive cell information
  const ALE::Obj<SieveMesh>& sieveMesh = fields->mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& cells = 
    sieveMesh->getLabelStratum("material-id", id());
  assert(!cells.isNull());
  const SieveMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveMesh::label_sequence::iterator cellsEnd = cells->end();
  const int numCells = cells->size();

  // Get sections associated with cohesive cells
  double_array residualCell(3*numBasis*spaceDim);
  const ALE::Obj<RealSection>& residualSection = residual.section();
  assert(!residualSection.isNull());
  UpdateAddVisitor residualVisitor(*residualSection, &residualCell[0]);

  double_array dispTCell(3*numBasis*spaceDim);
  topology::Field<topology::Mesh>& dispT = fields->get("disp(t)");
  const ALE::Obj<RealSection>& dispTSection = dispT.section();
  assert(!dispTSection.isNull());
  RestrictVisitor dispTVisitor(*dispTSection, dispTCell.size(), &dispTCell[0]);

  double_array dispTIncrCell(3*numBasis*spaceDim);
  topology::Field<topology::Mesh>& dispTIncr = fields->get("dispIncr(t->t+dt)");
  const ALE::Obj<RealSection>& dispTIncrSection = dispTIncr.section();
  assert(!dispTIncrSection.isNull());
  RestrictVisitor dispTIncrVisitor(*dispTIncrSection, 
				   dispTIncrCell.size(), &dispTIncrCell[0]);

  double_array dispTpdtCell(3*numBasis*spaceDim);

  // Get fault cell information
  const ALE::Obj<SieveMesh>& faultSieveMesh = _faultMesh->sieveMesh();
  assert(!faultSieveMesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& faultCells =
    faultSieveMesh->heightStratum(0);
  assert(!faultCells.isNull());
  assert(faultCells->size() == cells->size());

  // Get sections associated with fault cells
  double_array coordinatesCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& coordinates = 
    faultSieveMesh->getRealSection("coordinates");
  assert(!coordinates.isNull());
  RestrictVisitor coordsVisitor(*coordinates, 
				coordinatesCell.size(), &coordinatesCell[0]);

  double_array dispRelCell(numBasis*spaceDim);
  topology::Field<topology::SubMesh>& dispRel = _fields->get("relative disp");
  const ALE::Obj<RealSection>& dispRelSection = dispRel.section();
  assert(!dispRelSection.isNull());
  RestrictVisitor dispRelVisitor(*dispRelSection, 
				 dispRelCell.size(), &dispRelCell[0]);

  _logger->eventEnd(setupEvent);
#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventBegin(computeEvent);
#endif

  // Loop over cells
  for (SieveMesh::label_sequence::iterator c_iter=cellsBegin;
	 c_iter != cellsEnd;
       ++c_iter) {
    topology::SubMesh::SieveMesh::point_type c_fault = 
      _cohesiveToFault[*c_iter];

    // Compute geometry information for current cell
#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(geometryEvent);
#endif
#if defined(PRECOMPUTE_GEOMETRY)
    _quadrature->retrieveGeometry(c_fault);
#else
    coordsVisitor.clear();
    faultSieveMesh->restrictClosure(c_fault, coordsVisitor);
    _quadrature->computeGeometry(coordinatesCell, c_fault);
#endif

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(geometryEvent);
    _logger->eventBegin(restrictEvent);
#endif

    // Restrict input fields to cell
    dispTVisitor.clear();
    sieveMesh->restrictClosure(*c_iter, dispTVisitor);
    dispTIncrVisitor.clear();
    sieveMesh->restrictClosure(*c_iter, dispTIncrVisitor);

    dispRelVisitor.clear();
    faultSieveMesh->restrictClosure(c_fault, dispRelVisitor);

    // Get cell geometry information that depends on cell
    const double_array& basis = _quadrature->basis();
    const double_array& basisDeriv = _quadrature->basisDeriv();
    const double_array& jacobianDet = _quadrature->jacobianDet();

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(restrictEvent);
    _logger->eventBegin(computeEvent);
#endif

    // Compute current estimate of displacement at time t+dt using
    // solution increment.
    dispTpdtCell = dispTCell + dispTIncrCell;

    residualCell = 0.0;

    // Compute action for positive side of fault and Lagrange constraints
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      const double wt = quadWts[iQuad] * jacobianDet[iQuad];
      const int iQ = iQuad * numBasis;
      for (int iBasis=0; iBasis < numBasis; ++iBasis) {
        const double valI = wt*basis[iQ+iBasis];

	// Add entries to residual at
	// iN = DOF at negative node
	// iP = DOF at positive node
	// iL = DOF at constraint node

	// Indices for negative vertex
	const int iBN = 0*numBasis*spaceDim + iBasis*spaceDim;
	
	// Indices for positive vertex
	const int iBP = 1*numBasis*spaceDim + iBasis*spaceDim;
	
	// Indices for Lagrange vertex
	const int iBL = 2*numBasis*spaceDim + iBasis*spaceDim;

        for (int jBasis=0; jBasis < numBasis; ++jBasis) {
          const double valIJ = valI * basis[iQ+jBasis];

	  // Indices for negative vertex
	  const int jBN = 0*numBasis*spaceDim + jBasis*spaceDim;

	  // Indices for positive vertex
	  const int jBP = 1*numBasis*spaceDim + jBasis*spaceDim;

	  // Indices for Lagrange vertex
	  const int jBL = 2*numBasis*spaceDim + jBasis*spaceDim;

          for (int iDim=0; iDim < spaceDim; ++iDim) {
	    // negative side of the fault
            residualCell[iBN + iDim] += valIJ * dispTpdtCell[jBL + iDim];
	    
	    // positive side of the fault
            residualCell[iBP + iDim] -= valIJ * dispTpdtCell[jBL + iDim];
	    
	    // Lagrange constraint
	    residualCell[iBL + iDim] -= valIJ * 
	      (dispTpdtCell[jBP + iDim] - dispTpdtCell[jBN + iDim] -
	       dispRelCell[jBasis*spaceDim+iDim]);

#if 0
	    std::cout << "iBasis: " << iBasis
		      << ", jBasis: " << jBasis
		      << ", iDim: " << iDim
		      << ", valIJ: " << valIJ
		      << ", jacobianDet: " << jacobianDet[iQuad]
		      << ", dispRel: " << dispRelCell[jBasis*spaceDim+iDim]
		      << ", dispP: " << dispTpdtCell[jBP + iDim]
		      << ", dispN: " << dispTpdtCell[jBN + iDim]
		      << ", dispL: " << dispTpdtCell[jBL + iDim]
		      << ", residualN: " << residualCell[iBN + iDim]
		      << ", residualP: " << residualCell[iBP + iDim]
		      << ", residualL: " << residualCell[iBL + iDim]
		      << std::endl;
#endif

	  } // for
        } // for
      } // for
    } // for

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(computeEvent);
    _logger->eventBegin(updateEvent);
#endif

    // Assemble cell contribution into field
    residualVisitor.clear();
    sieveMesh->updateClosure(*c_iter, residualVisitor);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(updateEvent);
#endif
  } // for
  PetscLogFlops(numCells*spaceDim*spaceDim*8);

#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventEnd(computeEvent);
#endif
} // integrateResidual

// ----------------------------------------------------------------------
// Compute Jacobian matrix (A) associated with operator.
void
pylith::faults::FaultCohesiveLagrange::integrateJacobian(
				   topology::Jacobian* jacobian,
				   const double t,
				   topology::SolutionFields* const fields)
{ // integrateJacobian
  assert(0 != jacobian);
  assert(0 != fields);
  assert(0 != _fields);
  assert(0 != _logger);

  const int setupEvent = _logger->eventId("FaIJ setup");
  const int geometryEvent = _logger->eventId("FaIJ geometry");
  const int computeEvent = _logger->eventId("FaIJ compute");
  const int restrictEvent = _logger->eventId("FaIJ restrict");
  const int updateEvent = _logger->eventId("FaIJ update");

  _logger->eventBegin(setupEvent);

  // Add constraint information to Jacobian matrix; Entries are
  // associated with vertices ik, jk, ki, and kj.

  // Get cell geometry information that doesn't depend on cell
  const int numQuadPts = _quadrature->numQuadPts();
  const double_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == numQuadPts);
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int cellDim = _quadrature->cellDim();
  assert (cellDim == spaceDim-1);

  // Get cohesive cell information
  const ALE::Obj<SieveMesh>& sieveMesh = fields->mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& cells = 
    sieveMesh->getLabelStratum("material-id", id());
  assert(!cells.isNull());
  const SieveMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveMesh::label_sequence::iterator cellsEnd = cells->end();

  // Get sections
  const ALE::Obj<RealSection>& solnSection = fields->solution().section();
  assert(!solnSection.isNull());

  // Get fault cell information
  const ALE::Obj<SieveMesh>& faultSieveMesh = _faultMesh->sieveMesh();
  assert(!faultSieveMesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& faultCells =
    faultSieveMesh->heightStratum(0);
  assert(!faultCells.isNull());
  assert(faultCells->size() == cells->size());

  // Get sections associated with fault cells
  double_array coordinatesCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& coordinates = 
    faultSieveMesh->getRealSection("coordinates");
  assert(!coordinates.isNull());
  RestrictVisitor coordsVisitor(*coordinates, 
				coordinatesCell.size(), &coordinatesCell[0]);

  // Get sparse matrix
  double_array jacobianCell(3*numBasis*spaceDim * 3*numBasis*spaceDim);
  const PetscMat jacobianMat = jacobian->matrix();
  assert(0 != jacobianMat);

  const ALE::Obj<SieveMesh::order_type>& globalOrder = 
    sieveMesh->getFactory()->getGlobalOrder(sieveMesh, "default", solnSection);
  assert(!globalOrder.isNull());
  // We would need to request unique points here if we had an interpolated mesh
  IndicesVisitor jacobianVisitor(*solnSection, *globalOrder,
				 (int) pow(sieveMesh->getSieve()->getMaxConeSize(),
					   sieveMesh->depth())*spaceDim);

  _logger->eventEnd(setupEvent);
#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventBegin(computeEvent);
#endif

  // Loop over cells
  for (SieveMesh::label_sequence::iterator c_iter=cellsBegin;
       c_iter != cellsEnd;
       ++c_iter) {
    topology::SubMesh::SieveMesh::point_type c_fault = 
      _cohesiveToFault[*c_iter];

    // Compute geometry information for current cell
#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(geometryEvent);
#endif
#if defined(PRECOMPUTE_GEOMETRY)
    _quadrature->retrieveGeometry(c_fault);
#else
    coordsVisitor.clear();
    faultSieveMesh->restrictClosure(c_fault, coordsVisitor);
    _quadrature->computeGeometry(coordinatesCell, c_fault);
#endif

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(geometryEvent);
    _logger->eventBegin(computeEvent);
#endif

    // Get cell geometry information that depends on cell
    const double_array& basis = _quadrature->basis();
    const double_array& jacobianDet = _quadrature->jacobianDet();

    jacobianCell = 0.0;

    const int rowSize = 3*numBasis*spaceDim;

    // Compute Jacobian for constrants
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      const double wt = quadWts[iQuad] * jacobianDet[iQuad];
      for (int iBasis=0, iQ=iQuad*numBasis; iBasis < numBasis; ++iBasis) {
        const double valI = wt*basis[iQ+iBasis];

	// First index for positive, negative, and Lagrange vertices
	const int iBN = 0*numBasis*spaceDim + iBasis*spaceDim;
	const int iBP = 1*numBasis*spaceDim + iBasis*spaceDim;
	const int iBL = 2*numBasis*spaceDim + iBasis*spaceDim;

        for (int jBasis=0; jBasis < numBasis; ++jBasis) {
          const double valIJ = valI * basis[iQ+jBasis];

	  // First index for positive, negative, and Lagrange vertices
	  const int jBN = 0*numBasis*spaceDim + jBasis*spaceDim;
	  const int jBP = 1*numBasis*spaceDim + jBasis*spaceDim;
	  const int jBL = 2*numBasis*spaceDim + jBasis*spaceDim;

	  for (int iDim=0; iDim < spaceDim; ++iDim) {
	    // Add entries to Jacobian at (i,j) where
	    // iN = DOF at negative node, jL = DOF at constraint node
	    // iL = DOF at constraint node, jN = DOF at negative node
	    // iP = DOF at positive node, jL = DOF at constraint node
	    // iL = DOF at constraint node, jP = DOF at positive node

	    // Indices for negative vertex
            const int iN = (iBN + iDim) * rowSize; // row
            const int jN = (jBN + iDim); // col

	    // Indices for positive vertex
            const int iP = (iBP + iDim) * rowSize; // row
            const int jP = (jBP + iDim); // col

	    // Indices for Lagrange vertex
            const int iL = (iBL + iDim) * rowSize; // row
            const int jL = (jBL + iDim); // col

            jacobianCell[iN + jL] -= valIJ;
            jacobianCell[iL + jN] -= valIJ;

            jacobianCell[iP + jL] += valIJ;
            jacobianCell[iL + jP] += valIJ;
          } // for
        } // for
      } // for
    } // for

#if defined(DETAILED_EVENT_LOGGING)
    PetscLogFlops(spaceDim*spaceDim*2);
    _logger->eventEnd(computeEvent);
    _logger->eventBegin(updateEvent);
#endif

    // Assemble cell contribution into PETSc matrix.
    jacobianVisitor.clear();
    PetscErrorCode err = updateOperator(jacobianMat, *sieveMesh->getSieve(),
					jacobianVisitor, *c_iter,
					&jacobianCell[0], ADD_VALUES);
    CHECK_PETSC_ERROR_MSG(err, "Update to PETSc Mat failed.");

  } // for

#if 0 // DEBUGGING
  sieveMesh->getSendOverlap()->view("Send domain overlap");
  sieveMesh->getRecvOverlap()->view("Receive domain overlap");
#endif

  _needNewJacobian = false;
} // integrateJacobian

// ----------------------------------------------------------------------
// Compute Jacobian matrix (A) associated with operator.
void
pylith::faults::FaultCohesiveLagrange::integrateJacobian(
	                  topology::Field<topology::Mesh>* jacobian,
			  const double t,
			  topology::SolutionFields* const fields)
{ // integrateJacobian
  assert(0 != jacobian);
  assert(0 != fields);
  assert(0 != _fields);
  assert(0 != _logger);

  const int setupEvent = _logger->eventId("FaIJ setup");
  const int geometryEvent = _logger->eventId("FaIJ geometry");
  const int computeEvent = _logger->eventId("FaIJ compute");
  const int restrictEvent = _logger->eventId("FaIJ restrict");
  const int updateEvent = _logger->eventId("FaIJ update");

  _logger->eventBegin(setupEvent);

  // Integrate constraint information.

  // Get cell geometry information that doesn't depend on cell
  const int numQuadPts = _quadrature->numQuadPts();
  const double_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == numQuadPts);
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int cellDim = _quadrature->cellDim();
  assert (cellDim == spaceDim-1);

  // Get cohesive cell information
  const ALE::Obj<SieveMesh>& sieveMesh = fields->mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& cells = 
    sieveMesh->getLabelStratum("material-id", id());
  assert(!cells.isNull());
  const SieveMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveMesh::label_sequence::iterator cellsEnd = cells->end();
  const int numCells = cells->size();

  // Get sections
  const ALE::Obj<RealSection>& solnSection = fields->solution().section();
  assert(!solnSection.isNull());

  double_array jacobianCell(3*numBasis*spaceDim);
  const ALE::Obj<RealSection>& jacobianSection = jacobian->section();
  assert(!jacobianSection.isNull());
  UpdateAddVisitor jacobianVisitor(*jacobianSection, &jacobianCell[0]);

  // Get fault cell information
  const ALE::Obj<SieveMesh>& faultSieveMesh = _faultMesh->sieveMesh();
  assert(!faultSieveMesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& faultCells =
    faultSieveMesh->heightStratum(0);
  assert(!faultCells.isNull());
  const SieveMesh::label_sequence::iterator faultCellsBegin = 
    faultCells->begin();
  const SieveMesh::label_sequence::iterator faultCellsEnd = 
    faultCells->end();

  // Get sections associated with fault cells
  double_array coordinatesCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& coordinates = 
    faultSieveMesh->getRealSection("coordinates");
  assert(!coordinates.isNull());
  RestrictVisitor coordsVisitor(*coordinates, 
				coordinatesCell.size(), &coordinatesCell[0]);

  _logger->eventEnd(setupEvent);
#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventBegin(computeEvent);
#endif

  // Loop over cells
  for (SieveMesh::label_sequence::iterator c_iter=cellsBegin;
       c_iter != cellsEnd;
       ++c_iter) {
    topology::SubMesh::SieveMesh::point_type c_fault = 
      _cohesiveToFault[*c_iter];

    // Compute geometry information for current cell
#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(geometryEvent);
#endif
#if defined(PRECOMPUTE_GEOMETRY)
    _quadrature->retrieveGeometry(c_fault);
#else
    coordsVisitor.clear();
    faultSieveMesh->restrictClosure(c_fault, coordsVisitor);
    _quadrature->computeGeometry(coordinatesCell, c_fault);
#endif

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(geometryEvent);
    _logger->eventBegin(restrictEvent);
#endif

    // Get cell geometry information that depends on cell
    const double_array& basis = _quadrature->basis();
    const double_array& basisDeriv = _quadrature->basisDeriv();
    const double_array& jacobianDet = _quadrature->jacobianDet();

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(restrictEvent);
    _logger->eventBegin(computeEvent);
#endif

    jacobianCell = 0.0;

    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      const double wt = quadWts[iQuad] * jacobianDet[iQuad];
      const int iQ = iQuad * numBasis;
      double valJ = 0.0;
      for (int jBasis = 0; jBasis < numBasis; ++jBasis)
	valJ += basis[iQ + jBasis];
      valJ *= wt;
      for (int iBasis=0; iBasis < numBasis; ++iBasis) {
	const double valIJ = valJ * basis[iQ+iBasis];
	for (int iDim=0; iDim < spaceDim; ++iDim) {
	  // Lagrange constraints
	  jacobianCell[2*numBasis*spaceDim+iBasis*spaceDim+iDim] += valIJ;
	} // for
      } // for
    } // for


#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(computeEvent);
    _logger->eventBegin(updateEvent);
#endif

    // Assemble cell contribution into field
    jacobianVisitor.clear();
    sieveMesh->updateClosure(*c_iter, jacobianVisitor);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(updateEvent);
#endif
  } // for
  PetscLogFlops(numCells*numQuadPts*(2+numBasis) +
		numCells*numQuadPts*numBasis*spaceDim*3);

#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventEnd(computeEvent);
#endif

  _needNewJacobian = false;
} // integrateJacobian

// ----------------------------------------------------------------------
// Integrate contributions to Jacobian matrix (A) associated with
// operator.
void
pylith::faults::FaultCohesiveLagrange::calcPreconditioner(
				   PetscMat* const precondMatrix,
				   topology::Jacobian* const jacobian,
				   topology::SolutionFields* const fields)
{ // calcPreconditioner
  assert(0 != precondMatrix);
  assert(0 != jacobian);
  assert(0 != fields);
  assert(0 != _fields);
  assert(0 != _logger);

  /** We have A = [K L^T]
   *              [L   0]
   *
   * Compute P^(-1) = -( [L] [K]^(-1) [L]^T ) for each cohesive cell
   * using the diagonal of K to approximate [K]^(-1).
   *
   * Decompose [K] into [Kn] and [Kp], where [Kn] contains the terms
   * for vertices on the negative side of the fault and [Kp] contains
   * the terms for vertices on the positive side of the fault.
   *
   * P^(-1) = L_{ik} (1.0/Kn_{kk} + 1.0/Kp{kk}) L_{jk}
   *
   * L_{ij} = sum_cells sum_{qi} N_i N_j w_qi |J|_{qi}
   *
   * :KLUDGE: There is a difference between P^(-1) using the full
   * system matrices and assemblying P^(-1) from the cell
   * matrices. The square of the sums of the basis functions (system
   * matrices) is greater than the sums of the squares of the basis
   * functions (cell matrices). Empirical experiments suggests that
   * multiplying the cell matrices by a factor of about 5.0 improves
   * the rate of convergence.
   */
  const double assemblyFactor = 5.0;

  //#define USE_DIAGONAL_PRECONDITIONER // Extract only the diagonal terms of P.

  const int setupEvent = _logger->eventId("FaPr setup");
  const int computeEvent = _logger->eventId("FaPr compute");
  const int restrictEvent = _logger->eventId("FaPr restrict");
  const int updateEvent = _logger->eventId("FaPr update");

  _logger->eventBegin(setupEvent);

  // Get cell information and setup storage for cell data
  const int numQuadPts = _quadrature->numQuadPts();
  const double_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == numQuadPts);
  const int spaceDim = _quadrature->spaceDim();
  const int numBasis = _quadrature->numBasis();
  const int nrowsF = numBasis*spaceDim; // number of rows/cols in fault matrix

  // Size of fault preconditioner matrix for cell
  const int matrixSizeF = nrowsF * nrowsF;

  // Allocate vectors for vertex values
  int_array indicesN(nrowsF);
  int_array indicesP(nrowsF);
  int_array indicesLagrange(nrowsF);
#if defined(USE_DIAGONAL_PRECONDITIONER)
  double_array preconditionerCell(nrowsF);
#else
  double_array preconditionerCell(matrixSizeF);
#endif
  double_array jacobianCellP(matrixSizeF);
  double_array jacobianCellN(matrixSizeF);
  double_array basisProducts(numBasis*numBasis);

  // Get sections
  const ALE::Obj<RealSection>& solutionSection = fields->solution().section();
  assert(!solutionSection.isNull());

  const int numConstraintVert = numBasis;
  const int numCorners = 3 * numConstraintVert; // cohesive cell

  // Get cohesive cells
  const ALE::Obj<SieveMesh>& sieveMesh = fields->mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::order_type>& lagrangeGlobalOrder =
      sieveMesh->getFactory()->getGlobalOrder(sieveMesh, "faultDefault",
        solutionSection, spaceDim);
  assert(!lagrangeGlobalOrder.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& cellsCohesive =
    sieveMesh->getLabelStratum("material-id", id());
  assert(!cellsCohesive.isNull());
  const SieveMesh::label_sequence::iterator cellsCohesiveBegin =
    cellsCohesive->begin();
  const SieveMesh::label_sequence::iterator cellsCohesiveEnd =
    cellsCohesive->end();

  const PetscMat jacobianMatrix = jacobian->matrix();
  assert(0 != jacobianMatrix);
  const ALE::Obj<SieveMesh::order_type>& globalOrder =
    sieveMesh->getFactory()->getGlobalOrder(sieveMesh, "default", 
					    solutionSection);
  assert(!globalOrder.isNull());
  // We would need to request unique points here if we had an interpolated mesh
  IndicesVisitor jacobianVisitor(*solutionSection, *globalOrder,
				 (int) pow(sieveMesh->getSieve()->getMaxConeSize(),
					   sieveMesh->depth())*spaceDim);

  const ALE::Obj<SieveMesh::sieve_type>& sieve = sieveMesh->getSieve();
  assert(!sieve.isNull());
  ALE::ISieveVisitor::NConeRetriever<SieveMesh::sieve_type> ncV(*sieve,
      (size_t) pow(sieve->getMaxConeSize(), std::max(0, sieveMesh->depth())));

  // Get fault cell information
  const ALE::Obj<SieveMesh>& faultSieveMesh = _faultMesh->sieveMesh();
  assert(!faultSieveMesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& faultCells =
    faultSieveMesh->heightStratum(0);
  assert(!faultCells.isNull());
  assert(faultCells->size() == cellsCohesive->size());

  // Get sections associated with fault cells
  double_array coordinatesCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& coordinates = 
    faultSieveMesh->getRealSection("coordinates");
  assert(!coordinates.isNull());
  RestrictVisitor coordsVisitor(*coordinates, 
				coordinatesCell.size(), &coordinatesCell[0]);

  _logger->eventEnd(setupEvent);
#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventBegin(computeEvent);
#endif

  for (SieveMesh::label_sequence::iterator c_iter=cellsCohesiveBegin;
       c_iter != cellsCohesiveEnd;
       ++c_iter) {
    topology::SubMesh::SieveMesh::point_type c_fault = 
      _cohesiveToFault[*c_iter];

    // Compute geometry information for current cell
#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(geometryEvent);
#endif
#if defined(PRECOMPUTE_GEOMETRY)
    _quadrature->retrieveGeometry(c_fault);
#else
    coordsVisitor.clear();
    faultSieveMesh->restrictClosure(c_fault, coordsVisitor);
    _quadrature->computeGeometry(coordinatesCell, c_fault);
#endif

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(geometryEvent);
    _logger->eventBegin(restrictEvent);
#endif

    // Get cone for cohesive cell
    ncV.clear();
    ALE::ISieveTraversal<SieveMesh::sieve_type>::orientedClosure(*sieve,
								 *c_iter, ncV);
    const int coneSize = ncV.getSize();
    assert(coneSize == numCorners);
    const SieveMesh::point_type* cohesiveCone = ncV.getPoints();
    assert(0 != cohesiveCone);

    // Get indices and Jacobian values
    jacobianCellP = 0.0;
    jacobianCellN = 0.0;
    for (int iBasis = 0; iBasis < numBasis; ++iBasis) {
      // constraint vertex k
      const int v_negative = cohesiveCone[0*numBasis+iBasis];
      const int v_positive = cohesiveCone[1*numBasis+iBasis];
      const int v_lagrange = cohesiveCone[2*numBasis+iBasis];
      
      for (int iDim=0, iB=iBasis*spaceDim; iDim < spaceDim; ++iDim) {
	if (globalOrder->isLocal(v_negative))
	  indicesN[iB+iDim] = globalOrder->getIndex(v_negative) + iDim;
	else
	  indicesN[iB+iDim] = -1;
	if (globalOrder->isLocal(v_positive))
	  indicesP[iB+iDim] = globalOrder->getIndex(v_positive) + iDim;
	else
	  indicesP[iB+iDim] = -1;
	if (globalOrder->isLocal(v_lagrange))
	  indicesLagrange[iB+iDim] = lagrangeGlobalOrder->getIndex(v_lagrange) + iDim;
	else
	  indicesLagrange[iB+iDim] = -1;

	// Set matrix diagonal entries to 1.0 (used when vertex is not local).
	jacobianCellN[iB+iDim] = 1.0;
	jacobianCellP[iB+iDim] = 1.0;
      } // for
    } // for
    
    // Get values from Jacobian matrix.
    PetscErrorCode err = 0;
    err = MatGetValues(jacobianMatrix, 
		       indicesN.size(), &indicesN[0],
		       indicesN.size(), &indicesN[0],
		       &jacobianCellN[0]);
    CHECK_PETSC_ERROR_MSG(err, "Restrict from PETSc Mat failed.");
    
    err = MatGetValues(jacobianMatrix, 
		       indicesP.size(), &indicesP[0],
		       indicesP.size(), &indicesP[0],
		       &jacobianCellP[0]);
    CHECK_PETSC_ERROR_MSG(err, "Restrict from PETSc Mat failed.");
   
#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(restrictEvent);
    _logger->eventBegin(computeEvent);
#endif

    // Integrate preconditioner terms
    
    // Get cell geometry information that depends on cell
    const double_array& basis = _quadrature->basis();
    const double_array& jacobianDet = _quadrature->jacobianDet();

    // Compute product of basis functions.
    // Want values summed over quadrature points
    basisProducts = 0.0;
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      const double wt = quadWts[iQuad] * jacobianDet[iQuad];

      for (int iBasis=0, iQ=iQuad*numBasis; iBasis < numBasis; ++iBasis) {
        const double valI = wt*basis[iQ+iBasis];

        for (int jBasis=0; jBasis < numBasis; ++jBasis) {

	  basisProducts[iBasis*numBasis+jBasis] += valI*basis[iQ+jBasis];
	} // for
      } // for
    } // for

    preconditionerCell = 0.0;

#if defined(USE_DIAGONAL_PRECONDITIONER)
    for (int iBasis=0; iBasis < numBasis; ++iBasis) {
      // First index Lagrange vertices
      const int iB = iBasis*spaceDim;

      const int jBasis = iBasis; // Compute only diagonal entries
	  
      // First index for Lagrange vertices
      const int jB = jBasis*spaceDim;

      for (int kBasis=0; kBasis < numBasis; ++kBasis) {
	const int kB = kBasis*spaceDim;

	const double valIK = basisProducts[iBasis*numBasis+kBasis];
	const double valJK = basisProducts[jBasis*numBasis+kBasis];

	for (int iDim=0; iDim < spaceDim; ++iDim) {
	  const int k = (kB+iDim)*nrowsF+(kB+iDim);
	  const double jacobianNPInv = 
	    1.0/jacobianCellN[k] + 1.0/jacobianCellP[k];

	  preconditionerCell[iB+iDim] -= valIK*valJK*jacobianNPInv;
        } // for
      } // for
    } // for
#else
    for (int iBasis=0; iBasis < numBasis; ++iBasis) {
      // First index Lagrange vertices
      const int iB = iBasis*spaceDim;
      
      for (int jBasis=0; jBasis < numBasis; ++jBasis) {
	
	// First index for Lagrange vertices
	const int jB = jBasis*spaceDim;
	
	for (int kBasis=0; kBasis < numBasis; ++kBasis) {
	  const int kB = kBasis*spaceDim;
	  
	  const double valIK = basisProducts[iBasis*numBasis+kBasis];
	  const double valJK = basisProducts[jBasis*numBasis+kBasis];
	  
	  for (int iDim=0; iDim < spaceDim; ++iDim) {
	    // Add entries to Jacobian at (i,j) where
	    // iL = DOF at constraint node, jL = DOF at constraint node
	    const int iL = (iB + iDim) * nrowsF; // row
	    const int jL = (jB + iDim); // col
	    
	    const int k = (kB+iDim)*nrowsF+(kB+iDim);
	    const double jacobianNPInv = 
	      1.0/jacobianCellN[k] + 1.0/jacobianCellP[k];
	    
	    preconditionerCell[iL + jL] -= valIK*valJK*jacobianNPInv;
	  } // for
        } // for
      } // for
    } // for
#endif
    preconditionerCell *= assemblyFactor;

#if 0 // DEBUGGING
#if defined(USE_DIAGONAL_PRECONDITIONER)
    std::cout << "1/P_cell " << *c_iter << std::endl;
    for(int i = 0; i < nrowsF; ++i) {
      std::cout << "  " << preconditionerCell[i] << std::endl;
    }
#else
    std::cout << "1/P_cell " << *c_iter << std::endl;
    for(int i = 0; i < nrowsF; ++i) {
      for(int j = 0; j < nrowsF; ++j) {
        std::cout << "  " << preconditionerCell[i*nrowsF+j];
      }
      std::cout << std::endl;
    }
#endif
#endif

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(computeEvent);
    _logger->eventBegin(updateEvent);
#endif

#if defined(USE_DIAGONAL_PRECONDITIONER)
    for (int i=0; i < nrowsF; ++i)
      MatSetValue(*precondMatrix,
		  indicesLagrange[i],
		  indicesLagrange[i],
		  preconditionerCell[i],
		  ADD_VALUES);
#else
    err = MatSetValues(*precondMatrix,
                 indicesLagrange.size(), &indicesLagrange[0],
                 indicesLagrange.size(), &indicesLagrange[0],
                 &preconditionerCell[0],
                 ADD_VALUES);
#endif
    CHECK_PETSC_ERROR_MSG(err, "Setting values in fault preconditioner failed.");

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(updateEvent);
#endif

  } // for

#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventEnd(computeEvent);
#endif

} // calcPreconditioner

// ----------------------------------------------------------------------
// Adjust solution from solver with lumped Jacobian to match Lagrange
// multiplier constraints.
void
pylith::faults::FaultCohesiveLagrange::adjustSolnLumped(topology::SolutionFields* const fields,
                                                        const topology::Field<
							topology::Mesh>& jacobian)
{ // adjustSolnLumped
  assert(0 != fields);
  assert(0 != _quadrature);

  // Cohesive cells with conventional vertices N and P, and constraint
  // vertex L require 2 adjustments to the solution:
  //
  //   * DOF L: Compute increment in Lagrange multipliers
  //
  //     d\vec{l}_p = \tensor{S}^{-1} \cdot (-\vec{r}_p^{*} + \tensor{L}_p 
  //                   \cdot (d\vec{u}_{n+}^* - d\vec{u}_{n-}^*))
  //     \tensor{S} = \tensor{L}_p \cdot 
  //                   (\tensor{K}_{n+n+}{-1} + \tensor{K}_{n-n-})
  //                   \cdot \tensor{L}_p^T
  //
  //   * DOF i and j: Adjust displacement increment (solution) to create slip
  //     consistent with Lagrange multiplier constraints
  //
  //     d\vec{u}_{n+} = d\vec{u}_{n+}^{*} - \tensor{K}_{n+n+}^{-1} \cdot 
  //                     \tensor{L}_p^T \cdot d\vec{l}_p
  //     d\vec{u}_{n-} = d\vec{u}_{n-}^{*} + \tensor{K}_{n-n-}^{-1} \cdot 
  //                     \tensor{L}_p^T \cdot d\vec{l}_p

  const int setupEvent = _logger->eventId("FaAS setup");
  const int geometryEvent = _logger->eventId("FaAS geometry");
  const int computeEvent = _logger->eventId("FaAS compute");
  const int restrictEvent = _logger->eventId("FaAS restrict");
  const int updateEvent = _logger->eventId("FaAS update");

  _logger->eventBegin(setupEvent);

  // Get cell information and setup storage for cell data
  const int spaceDim = _quadrature->spaceDim();

  // Get section information
  const ALE::Obj<RealSection>& jacobianSection = jacobian.section();
  assert(!jacobianSection.isNull());

  const ALE::Obj<RealSection>& residualSection =
      fields->get("residual").section();
  assert(!residualSection.isNull());

  double_array dispTIncrVertexN(spaceDim);
  double_array dispTIncrVertexP(spaceDim);
  double_array dispTIncrVertexL(spaceDim); // Lagrange multiplier
  const ALE::Obj<RealSection>& dispTIncrSection = fields->get(
    "dispIncr(t->t+dt)").section();
  assert(!dispTIncrSection.isNull());

  const ALE::Obj<RealSection>& dispTIncrAdjSection = fields->get(
    "dispIncr adjust").section();
  assert(!dispTIncrAdjSection.isNull());

  const ALE::Obj<SieveMesh>& sieveMesh = fields->mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::order_type>& globalOrder =
    sieveMesh->getFactory()->getGlobalOrder(sieveMesh, "default", 
					    jacobianSection);
  assert(!globalOrder.isNull());

  _logger->eventEnd(setupEvent);

#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventBegin(computeEvent);
#endif

  const int numVertices = _cohesiveVertices.size();
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_lagrange = _cohesiveVertices[iVertex].lagrange;
    const int v_fault = _cohesiveVertices[iVertex].fault;
    const int v_negative = _cohesiveVertices[iVertex].negative;
    const int v_positive = _cohesiveVertices[iVertex].positive;

    // Compute contribution only if Lagrange constraint is local.
    if (!globalOrder->isLocal(v_lagrange))
      continue;

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(restrictEvent);
#endif

    // Get residual at cohesive cell's vertices.
    assert(spaceDim == residualSection->getFiberDimension(v_lagrange));
    const double* residualVertexL = residualSection->restrictPoint(v_lagrange);
    assert(residualVertexL);

    // Get jacobian at cohesive cell's vertices.
    assert(spaceDim == jacobianSection->getFiberDimension(v_negative));
    const double* jacobianVertexN = jacobianSection->restrictPoint(v_negative);
    assert(jacobianVertexN);

    assert(spaceDim == jacobianSection->getFiberDimension(v_positive));
    const double* jacobianVertexP = jacobianSection->restrictPoint(v_positive);
    assert(jacobianVertexP);

    assert(spaceDim == jacobianSection->getFiberDimension(v_lagrange));
    const double* jacobianVertexL = jacobianSection->restrictPoint(v_lagrange);
    assert(jacobianVertexL);

    // Get dispIncr(t) at cohesive cell's vertices.
    assert(spaceDim == dispTIncrSection->getFiberDimension(v_negative));
    dispTIncrSection->restrictPoint(v_negative, &dispTIncrVertexN[0],
				    dispTIncrVertexN.size());

    assert(spaceDim == dispTIncrSection->getFiberDimension(v_positive));
    dispTIncrSection->restrictPoint(v_positive, &dispTIncrVertexP[0],
				    dispTIncrVertexP.size());

    assert(spaceDim == dispTIncrSection->getFiberDimension(v_lagrange));
    dispTIncrSection->restrictPoint(v_lagrange, &dispTIncrVertexL[0],
				    dispTIncrVertexL.size());

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(restrictEvent);
    _logger->eventBegin(computeEvent);
#endif

    for (int iDim=0; iDim < spaceDim; ++iDim) {
      assert(jacobianVertexL[iDim] > 0.0);
      const double S = (1.0/jacobianVertexP[iDim] + 1.0/jacobianVertexN[iDim]) *
	jacobianVertexL[iDim]*jacobianVertexL[iDim];
      dispTIncrVertexL[iDim] = 1.0/S * 
	(-residualVertexL[iDim] +
	 jacobianVertexL[iDim] * 
	 (dispTIncrVertexP[iDim] - dispTIncrVertexN[iDim]));

      assert(jacobianVertexN[iDim] > 0.0);
      dispTIncrVertexN[iDim] = 
	+jacobianVertexL[iDim]/jacobianVertexN[iDim]*dispTIncrVertexL[iDim];

      assert(jacobianVertexP[iDim] > 0.0);
      dispTIncrVertexP[iDim] = 
	-jacobianVertexL[iDim]/jacobianVertexP[iDim]*dispTIncrVertexL[iDim];

    } // for

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(computeEvent);
    _logger->eventBegin(updateEvent);
#endif

    // Adjust displacements to account for Lagrange multiplier values
    // (assumed to be zero in preliminary solve).
    assert(dispTIncrVertexN.size() == 
	   dispTIncrAdjSection->getFiberDimension(v_negative));
    dispTIncrAdjSection->updateAddPoint(v_negative, &dispTIncrVertexN[0]);

    assert(dispTIncrVertexP.size() == 
	   dispTIncrAdjSection->getFiberDimension(v_positive));
    dispTIncrAdjSection->updateAddPoint(v_positive, &dispTIncrVertexP[0]);

    // Set Lagrange multiplier value. Value from preliminary solve is
    // bogus due to artificial diagonal entry.
    assert(dispTIncrVertexL.size() == 
	   dispTIncrSection->getFiberDimension(v_lagrange));
    dispTIncrSection->updatePoint(v_lagrange, &dispTIncrVertexL[0]);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(updateEvent);
#endif
  } // for

#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventEnd(computeEvent);
#endif
} // adjustSolnLumped

// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::faults::FaultCohesiveLagrange::verifyConfiguration(const topology::Mesh& mesh) const
{ // verifyConfiguration
  assert(0 != _quadrature);

  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  assert(!sieveMesh.isNull());

  if (!sieveMesh->hasIntSection(label())) {
    std::ostringstream msg;
    msg << "Mesh missing group of vertices '" << label()
        << " for boundary condition.";
    throw std::runtime_error(msg.str());
  } // if  

  // check compatibility of mesh and quadrature scheme
  const int dimension = mesh.dimension() - 1;
  if (_quadrature->cellDim() != dimension) {
    std::ostringstream msg;
    msg << "Dimension of reference cell in quadrature scheme ("
        << _quadrature->cellDim()
        << ") does not match dimension of cells in mesh (" << dimension
        << ") for fault '" << label() << "'.";
    throw std::runtime_error(msg.str());
  } // if

  // Check quadrature against mesh
  const int numCorners = _quadrature->numBasis();
  const ALE::Obj<SieveMesh::label_sequence>& cells =
      sieveMesh->getLabelStratum("material-id", id());
  assert(!cells.isNull());
  const SieveMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveMesh::label_sequence::iterator cellsEnd = cells->end();
  for (SieveMesh::label_sequence::iterator c_iter = cellsBegin; c_iter
      != cellsEnd; ++c_iter) {
    const int cellNumCorners = sieveMesh->getNumCellCorners(*c_iter);
    if (3 * numCorners != cellNumCorners) {
      std::ostringstream msg;
      msg << "Number of vertices in reference cell (" << numCorners
          << ") is not compatible with number of vertices (" << cellNumCorners
          << ") in cohesive cell " << *c_iter << " for fault '" << label()
          << "'.";
      throw std::runtime_error(msg.str());
    } // if
  } // for
} // verifyConfiguration

// ----------------------------------------------------------------------
// Verify constraints are acceptable.
void
pylith::faults::FaultCohesiveLagrange::checkConstraints(const topology::Field<topology::Mesh>& solution) const
{ // checkConstraints Check to make sure no vertices connected to the
  // fault are constrained.

  const ALE::Obj<RealSection>& section = solution.section();
  assert(!section.isNull());

  const int spaceDim = solution.mesh().dimension();

  const int numVertices = _cohesiveVertices.size();
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_negative = _cohesiveVertices[iVertex].negative;    
    const int fiberDimN = section->getFiberDimension(v_negative);
    assert(spaceDim == fiberDimN);
    const int numConstraintsN = section->getConstraintDimension(v_negative);
    if (numConstraintsN > 0) {
      std::ostringstream msg;
      msg << "Vertex with label '" << v_negative << "' on negative side "
	  << "of fault '" << label() << "' is constrained.\n"
	  << "Fault vertices cannot be constrained.";
      throw std::runtime_error(msg.str());
    } // if
    
    const int v_positive = _cohesiveVertices[iVertex].positive;
    const int fiberDimP = section->getFiberDimension(v_positive);
    assert(spaceDim == fiberDimP);
    const int numConstraintsP = section->getConstraintDimension(v_positive);
    if (numConstraintsP > 0) {
      std::ostringstream msg;
      msg << "Vertex with label '" << v_positive << "' on positive side "
	  << "of fault '" << label() << "' is constrained.\n"
	  << "Fault vertices cannot be constrained.";
      throw std::runtime_error(msg.str());
    } // if
  } // for
} // checkConstraints

// ----------------------------------------------------------------------
// Initialize auxiliary cohesive cell information.
void pylith::faults::FaultCohesiveLagrange::_initializeCohesiveInfo(const topology::Mesh& mesh)
{ // _initializeCohesiveInfo
  assert(0 != _quadrature);

  // Get cohesive cells
  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& cells =
      sieveMesh->getLabelStratum("material-id", id());
  assert(!cells.isNull());
  const SieveMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveMesh::label_sequence::iterator cellsEnd = cells->end();

  const int numConstraintVert = _quadrature->numBasis();
  const int numCorners = 3 * numConstraintVert; // cohesive cell

  const ALE::Obj<SieveSubMesh>& faultSieveMesh = _faultMesh->sieveMesh();
  assert(!faultSieveMesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& faultCells =
    faultSieveMesh->heightStratum(0);
  assert(!faultCells.isNull());
  SieveSubMesh::label_sequence::iterator f_iter = faultCells->begin();

  SieveSubMesh::renumbering_type& renumbering = faultSieveMesh->getRenumbering();
  const SieveSubMesh::renumbering_type::const_iterator renumberingEnd =
    renumbering.end();
  const ALE::Obj<SieveSubMesh::label_sequence>& vertices =
      faultSieveMesh->depthStratum(0);

  _cohesiveToFault.clear();
  typedef std::map<int, int> indexmap_type;
  indexmap_type indexMap;
  _cohesiveVertices.resize(vertices->size());
  int index = 0;

  const ALE::Obj<SieveMesh::sieve_type>& sieve = mesh.sieveMesh()->getSieve();
  assert(!sieve.isNull());

  ALE::ISieveVisitor::NConeRetriever<SieveMesh::sieve_type> ncV(*sieve,
      (size_t) pow(sieve->getMaxConeSize(), std::max(0, sieveMesh->depth())));

  for (SieveMesh::label_sequence::iterator c_iter=cellsBegin;
      c_iter != cellsEnd;
      ++c_iter, ++f_iter) {
    _cohesiveToFault[*c_iter] = *f_iter;

    // Get oriented closure
    ncV.clear();
    ALE::ISieveTraversal<SieveMesh::sieve_type>::orientedClosure(*sieve, *c_iter, ncV);
    const int coneSize = ncV.getSize();
    assert(coneSize == numCorners);
    const SieveMesh::point_type *cone = ncV.getPoints();
    assert(0 != cone);

    for (int iConstraint = 0; iConstraint < numConstraintVert; ++iConstraint) {
      // normal cohesive vertices i and j and constraint vertex k
      const int indexI = iConstraint;
      const int indexJ = iConstraint + numConstraintVert;
      const int indexK = iConstraint + 2 * numConstraintVert;

      const int v_lagrange = cone[indexK];
      const int v_negative = cone[indexI];
      const int v_positive = cone[indexJ];
      const int v_fault = renumbering[v_lagrange];

      if (indexMap.end() == indexMap.find(v_lagrange)) {
        _cohesiveVertices[index].lagrange = v_lagrange;
        _cohesiveVertices[index].positive = v_positive;
        _cohesiveVertices[index].negative = v_negative;
        _cohesiveVertices[index].fault = v_fault;
#if 0
	std::cout << "cohesiveVertices[" << index << "]: "
		  << "l: " << v_lagrange
		  << ", p: " << v_positive
		  << ", n: " << v_negative
		  << ", f: " << v_fault
		  << std::endl;
#endif
        indexMap[v_lagrange] = index; // add index to map
        ++index;
      } // if
    } // for
  } // for
} // _initializeCohesiveInfo

// ----------------------------------------------------------------------
// Initialize logger.
void
pylith::faults::FaultCohesiveLagrange::_initializeLogger(void)
{ // initializeLogger
  delete _logger;
  _logger = new utils::EventLogger;
  assert(0 != _logger);
  _logger->className("FaultCohesiveLagrange");
  _logger->initialize();

  _logger->registerEvent("FaAS setup");
  _logger->registerEvent("FaAS geometry");
  _logger->registerEvent("FaAS compute");
  _logger->registerEvent("FaAS restrict");
  _logger->registerEvent("FaAS update");

  _logger->registerEvent("FaIR setup");
  _logger->registerEvent("FaIR geometry");
  _logger->registerEvent("FaIR compute");
  _logger->registerEvent("FaIR restrict");
  _logger->registerEvent("FaIR update");

  _logger->registerEvent("FaIJ setup");
  _logger->registerEvent("FaIJ geometry");
  _logger->registerEvent("FaIJ compute");
  _logger->registerEvent("FaIJ restrict");
  _logger->registerEvent("FaIJ update");

  _logger->registerEvent("FaPr setup");
  _logger->registerEvent("FaPr geometry");
  _logger->registerEvent("FaPr compute");
  _logger->registerEvent("FaPr restrict");
  _logger->registerEvent("FaPr update");
} // initializeLogger

// ----------------------------------------------------------------------
// Calculate orientation at fault vertices.
void
pylith::faults::FaultCohesiveLagrange::_calcOrientation(const double upDir[3])
{ // _calcOrientation
  assert(0 != upDir);
  assert(0 != _faultMesh);
  assert(0 != _fields);

  double_array upDirArray(upDir, 3);

  // Get vertices in fault mesh.
  const ALE::Obj<SieveSubMesh>& faultSieveMesh = _faultMesh->sieveMesh();
  assert(!faultSieveMesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& vertices =
      faultSieveMesh->depthStratum(0);
  const SieveSubMesh::label_sequence::iterator verticesBegin =
      vertices->begin();
  const SieveSubMesh::label_sequence::iterator verticesEnd = vertices->end();

  // Containers for orientation information.
  const int cohesiveDim = _faultMesh->dimension();
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int orientationSize = spaceDim * spaceDim;
  const feassemble::CellGeometry& cellGeometry = _quadrature->refGeometry();
  const double_array& verticesRef = cellGeometry.vertices();
  const int jacobianSize = (cohesiveDim > 0) ? spaceDim * cohesiveDim : 1;
  const double_array& quadWts = _quadrature->quadWts();
  double_array jacobian(jacobianSize);
  double jacobianDet = 0;
  double_array orientationVertex(orientationSize);
  double_array coordinatesCell(numBasis * spaceDim);
  double_array refCoordsVertex(cohesiveDim);

  // Allocate orientation field.
  _fields->add("orientation", "orientation");
  topology::Field<topology::SubMesh>& orientation = _fields->get("orientation");
  const topology::Field<topology::SubMesh>& dispRel = 
    _fields->get("relative disp");
  orientation.newSection(dispRel, orientationSize);
  const ALE::Obj<RealSection>& orientationSection = orientation.section();
  assert(!orientationSection.isNull());
  // Create subspaces for along-strike, up-dip, and normal directions
  for (int iDim = 0; iDim <= cohesiveDim; ++iDim)
    orientationSection->addSpace();
  for (int iDim = 0; iDim <= cohesiveDim; ++iDim)
    orientationSection->setFiberDimension(vertices, spaceDim, iDim);
  orientation.allocate();
  orientation.zero();

  // Get fault cells.
  const ALE::Obj<SieveSubMesh::label_sequence>& cells =
      faultSieveMesh->heightStratum(0);
  assert(!cells.isNull());
  const SieveSubMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveSubMesh::label_sequence::iterator cellsEnd = cells->end();

  // Compute orientation of fault at constraint vertices

  // Get section containing coordinates of vertices
  const ALE::Obj<RealSection>& coordinatesSection =
      faultSieveMesh->getRealSection("coordinates");
  assert(!coordinatesSection.isNull());
  RestrictVisitor coordinatesVisitor(*coordinatesSection,
				     coordinatesCell.size(),
				     &coordinatesCell[0]);

  // Set orientation function
  assert(cohesiveDim == _quadrature->cellDim());
  assert(spaceDim == _quadrature->spaceDim());

  // Loop over cohesive cells, computing orientation weighted by
  // jacobian at constraint vertices

  const ALE::Obj<SieveSubMesh::sieve_type>& sieve = faultSieveMesh->getSieve();
  assert(!sieve.isNull());

  ALE::ISieveVisitor::NConeRetriever<SieveMesh::sieve_type>
      ncV(*sieve, (size_t) pow(sieve->getMaxConeSize(), std::max(0,
        faultSieveMesh->depth())));

  for (SieveSubMesh::label_sequence::iterator c_iter = cellsBegin; c_iter
      != cellsEnd; ++c_iter) {
    // Get orientations at fault cell's vertices.
    coordinatesVisitor.clear();
    faultSieveMesh->restrictClosure(*c_iter, coordinatesVisitor);

    ncV.clear();
    ALE::ISieveTraversal<SieveSubMesh::sieve_type>::orientedClosure(*sieve,
      *c_iter, ncV);
    const int coneSize = ncV.getSize();
    const SieveSubMesh::point_type *cone = ncV.getPoints();

    for (int v = 0; v < coneSize; ++v) {
      // Compute Jacobian and determinant of Jacobian at vertex
      memcpy(&refCoordsVertex[0], &verticesRef[v * cohesiveDim], cohesiveDim
          * sizeof(double));
      cellGeometry.jacobian(&jacobian, &jacobianDet, coordinatesCell,
        refCoordsVertex);

      // Compute orientation
      cellGeometry.orientation(&orientationVertex, jacobian, jacobianDet,
        upDirArray);

      // Update orientation
      orientationSection->updateAddPoint(cone[v], &orientationVertex[0]);
    } // for
  } // for

  //orientation.view("ORIENTATION BEFORE COMPLETE");

  // Assemble orientation information
  orientation.complete();

  // Loop over vertices, make orientation information unit magnitude
  double_array vertexDir(orientationSize);
  int count = 0;
  for (SieveSubMesh::label_sequence::iterator v_iter = verticesBegin; v_iter
      != verticesEnd; ++v_iter, ++count) {
    orientationVertex = 0.0;
    orientationSection->restrictPoint(*v_iter, &orientationVertex[0],
      orientationVertex.size());
    for (int iDim = 0; iDim < spaceDim; ++iDim) {
      double mag = 0;
      for (int jDim = 0, index = iDim * spaceDim; jDim < spaceDim; ++jDim)
        mag += pow(orientationVertex[index + jDim], 2);
      mag = sqrt(mag);
      assert(mag > 0.0);
      for (int jDim = 0, index = iDim * spaceDim; jDim < spaceDim; ++jDim)
        orientationVertex[index + jDim] /= mag;
    } // for

    orientationSection->updatePoint(*v_iter, &orientationVertex[0]);
  } // for
  PetscLogFlops(count * orientationSize * 4);

  if (1 == cohesiveDim && vertices->size() > 0) {
    // Default sense of positive slip is left-lateral and
    // fault-opening.
    // 
    // If fault is dipping, then we use the up-dir to make sure the
    // sense of positive slip is reverse and fault-opening.
    //
    // Check orientation of first vertex, (1) if dot product of the
    // normal-dir with preferred up-dir is positive, then we want dot
    // product of shear-dir and preferred up-dir to be positive and
    // (2) if the dot product of the normal-dir with preferred up-dir
    // is negative, then we want the dot product of the shear-dir and
    // preferred up-dir to be negative.
    //
    // When we flip the shear direction, we create a left-handed
    // coordinate system, but it gives the correct sense of slip. In
    // reality the shear/normal directions that are used are the
    // opposite of what we would want, but we cannot flip the fault
    // normal direction because it is tied to how the cohesive cells
    // are created.
    assert(vertices->size() > 0);
    orientationSection->restrictPoint(*vertices->begin(),
      &orientationVertex[0], orientationVertex.size());

    assert(2 == spaceDim);
    const double* shearDirVertex = &orientationVertex[0];
    const double* normalDirVertex = &orientationVertex[2];
    const double shearDirDot = 
      upDir[0] * shearDirVertex[0] + upDir[1] * shearDirVertex[1];
    const double normalDirDot = 
      upDir[0] * normalDirVertex[0] + upDir[1] * normalDirVertex[1];

    const int ishear = 0;
    const int inormal = 2;
    if (normalDirDot * shearDirDot < 0.0) {
      // Flip shear direction
      for (SieveSubMesh::label_sequence::iterator v_iter=verticesBegin; 
	   v_iter != verticesEnd;
	   ++v_iter) {
        orientationSection->restrictPoint(*v_iter, &orientationVertex[0],
					  orientationVertex.size());
        assert(4 == orientationSection->getFiberDimension(*v_iter));
        for (int iDim = 0; iDim < 2; ++iDim) // flip shear
          orientationVertex[ishear + iDim] *= -1.0;
	
        // Update orientation
        orientationSection->updatePoint(*v_iter, &orientationVertex[0]);
      } // for
      PetscLogFlops(3 + count * 2);
    } // if

  } else if (2 == cohesiveDim && vertices->size() > 0) {
    // Check orientation of first vertex, if dot product of fault
    // normal with preferred normal is negative, flip up/down dip
    // direction.
    //
    // Check orientation of first vertex, (1) if dot product of the
    // normal-dir with preferred up-dir is positive, then we want dot
    // product of shear-dir and preferred up-dir to be positive and
    // (2) if the dot product of the normal-dir with preferred up-dir
    // is negative, then we want the dot product of the shear-dir and
    // preferred up-dir to be negative.
    //
    // When we flip the up/down dip direction, we create a left-handed
    // strike/dip/normal coordinate system, but it gives the correct
    // sense of slip. In reality the strike/dip/normal directions that
    // are used are the opposite of what we would want, but we cannot
    // flip the fault normal direction because it is tied to how the
    // cohesive cells are created.

    assert(vertices->size() > 0);
    orientationSection->restrictPoint(*vertices->begin(),
      &orientationVertex[0], orientationVertex.size());

    assert(3 == spaceDim);
    const double* dipDirVertex = &orientationVertex[3];
    const double* normalDirVertex = &orientationVertex[6];
    const double dipDirDot = 
      upDir[0]*dipDirVertex[0] + 
      upDir[1]*dipDirVertex[1] + 
      upDir[2]*dipDirVertex[2];
    const double normalDirDot = 
      upDir[0]*normalDirVertex[0] +
      upDir[1]*normalDirVertex[1] +
      upDir[2]*normalDirVertex[2];

    const int istrike = 0;
    const int idip = 3;
    const int inormal = 6;
    if (dipDirDot * normalDirDot < 0.0 ||
	fabs(normalDirVertex[2] + 1.0) < 0.001) {
      // if fault normal is (0,0,+-1) then up-dir dictates reverse
      // motion for case with normal (0,0,1), so we reverse the dip-dir
      // if we have (0,0,-1).

      // Flip dip direction
      for (SieveSubMesh::label_sequence::iterator v_iter = verticesBegin; v_iter
	     != verticesEnd; ++v_iter) {
        orientationSection->restrictPoint(*v_iter, &orientationVertex[0],
					  orientationVertex.size());
        assert(9 == orientationSection->getFiberDimension(*v_iter));
        for (int iDim = 0; iDim < 3; ++iDim) // flip dip
          orientationVertex[idip + iDim] *= -1.0;
	
        // Update direction
        orientationSection->updatePoint(*v_iter, &orientationVertex[0]);
      } // for
      PetscLogFlops(5 + count * 3);
    } // if
  } // if

  //orientation.view("ORIENTATION");
} // _calcOrientation

// ----------------------------------------------------------------------
// Compute change in tractions on fault surface using solution.
void
pylith::faults::FaultCohesiveLagrange::_calcTractionsChange(
    topology::Field<topology::SubMesh>* tractions,
    const topology::Field<topology::Mesh>& dispT)
{ // _calcTractionsChange
  assert(0 != tractions);
  assert(0 != _faultMesh);
  assert(0 != _fields);
  assert(0 != _normalizer);

  tractions->label("traction_change");
  tractions->scale(_normalizer->pressureScale());

  // Fiber dimension of tractions matches spatial dimension.
  const int spaceDim = _quadrature->spaceDim();
  double_array tractionsVertex(spaceDim);

  // Get sections.
  const ALE::Obj<RealSection>& dispTSection = dispT.section();
  assert(!dispTSection.isNull());

  const ALE::Obj<RealSection>& orientationSection =
    _fields->get("orientation").section();
  assert(!orientationSection.isNull());

  // Allocate buffer for tractions field (if necessary).
  const ALE::Obj<RealSection>& tractionsSection = tractions->section();
  if (tractionsSection.isNull()) {
    ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
    //logger.stagePush("Fault");

    const topology::Field<topology::SubMesh>& dispRel = 
      _fields->get("relative disp");
    tractions->cloneSection(dispRel);

    //logger.stagePop();
  } // if
  assert(!tractionsSection.isNull());
  tractions->zero();

  const double pressureScale = tractions->scale();

  const int numVertices = _cohesiveVertices.size();
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_lagrange = _cohesiveVertices[iVertex].lagrange;
    const int v_fault = _cohesiveVertices[iVertex].fault;

    assert(spaceDim == dispTSection->getFiberDimension(v_lagrange));
    const double* dispTVertex = dispTSection->restrictPoint(v_lagrange);
    assert(0 != dispTVertex);

    assert(spaceDim*spaceDim == orientationSection->getFiberDimension(v_fault));
    const double* orientationVertex = 
      orientationSection->restrictPoint(v_fault);
    assert(orientationVertex);

    // Rotate from global coordinate system to fault (orientation)
    tractionsVertex = 0.0;
    for (int iDim=0; iDim < spaceDim; ++iDim)
      for (int jDim=0; jDim < spaceDim; ++jDim)
	tractionsVertex[iDim] += 
	  orientationVertex[iDim*spaceDim+jDim] * dispTVertex[jDim];

    assert(tractionsVertex.size() == 
	   tractionsSection->getFiberDimension(v_fault));
    tractionsSection->updatePoint(v_fault, &tractionsVertex[0]);
  } // for

  PetscLogFlops(numVertices * (1 + spaceDim) );

#if 0 // DEBUGGING
  tractions->view("TRACTIONS");
#endif
} // _calcTractionsChange

// ----------------------------------------------------------------------
// Transform field from local (fault) coordinate system to
// global coordinate system.
void
pylith::faults::FaultCohesiveLagrange::_faultToGlobal(topology::Field<topology::SubMesh>* field)
{ // _faultToGlobal
  assert(field);
  assert(0 != _faultMesh);
  assert(0 != _fields);

  // Fiber dimension of vector field matches spatial dimension.
  const int spaceDim = _quadrature->spaceDim();
  double_array fieldVertexGlobal(spaceDim);

  // Get sections.
  const ALE::Obj<RealSection>& fieldSection = field->section();
  assert(!fieldSection.isNull());

  const ALE::Obj<RealSection>& orientationSection =
    _fields->get("orientation").section();
  assert(!orientationSection.isNull());

  const int numVertices = _cohesiveVertices.size();
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_fault = _cohesiveVertices[iVertex].fault;

    assert(spaceDim == fieldSection->getFiberDimension(v_fault));
    const double* fieldVertexFault = fieldSection->restrictPoint(v_fault);
    assert(fieldVertexFault);

    assert(spaceDim*spaceDim == orientationSection->getFiberDimension(v_fault));
    const double* orientationVertex = orientationSection->restrictPoint(v_fault);
    assert(orientationVertex);

    // Rotate from fault to global coordinate system (transpose orientation)
    fieldVertexGlobal = 0.0;
    for (int iDim=0; iDim < spaceDim; ++iDim)
      for (int jDim=0; jDim < spaceDim; ++jDim)
	fieldVertexGlobal[iDim] += 
	  orientationVertex[jDim*spaceDim+iDim] * fieldVertexFault[jDim];

    assert(fieldVertexGlobal.size() == 
	   fieldSection->getFiberDimension(v_fault));
    fieldSection->updatePoint(v_fault, &fieldVertexGlobal[0]);
  } // for
  
  PetscLogFlops(numVertices * (2*spaceDim*spaceDim) );

#if 0 // DEBUGGING
  field->view("FIELD (GLOBAL)");
#endif
} // _faultToGlobal

// ----------------------------------------------------------------------
// Transform field from global coordinate system to local (fault)
// coordinate system.
void
pylith::faults::FaultCohesiveLagrange::_globalToFault(topology::Field<topology::SubMesh>* field)
{ // _globalToFault
  assert(field);
  assert(0 != _faultMesh);
  assert(0 != _fields);

  // Fiber dimension of vector field matches spatial dimension.
  const int spaceDim = _quadrature->spaceDim();
  double_array fieldVertexFault(spaceDim);

  // Get sections.
  const ALE::Obj<RealSection>& fieldSection = field->section();
  assert(!fieldSection.isNull());

  const ALE::Obj<RealSection>& orientationSection =
    _fields->get("orientation").section();
  assert(!orientationSection.isNull());

  const int numVertices = _cohesiveVertices.size();
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_fault = _cohesiveVertices[iVertex].fault;

    assert(spaceDim == fieldSection->getFiberDimension(v_fault));
    const double* fieldVertexGlobal = fieldSection->restrictPoint(v_fault);
    assert(fieldVertexGlobal);

    assert(spaceDim*spaceDim == orientationSection->getFiberDimension(v_fault));
    const double* orientationVertex = orientationSection->restrictPoint(v_fault);
    assert(orientationVertex);

    // Rotate from global coordinate system to fault (orientation)
    fieldVertexFault = 0.0;
    for (int iDim=0; iDim < spaceDim; ++iDim)
      for (int jDim=0; jDim < spaceDim; ++jDim)
	fieldVertexFault[iDim] += 
	  orientationVertex[iDim*spaceDim+jDim] * fieldVertexGlobal[jDim];

    assert(fieldVertexFault.size() == 
	   fieldSection->getFiberDimension(v_fault));
    fieldSection->updatePoint(v_fault, &fieldVertexFault[0]);
  } // for
  
  PetscLogFlops(numVertices * (2*spaceDim*spaceDim) );

#if 0 // DEBUGGING
  field->view("FIELD (FAULT)");
#endif
} // _faultToGlobal

// ----------------------------------------------------------------------
// Allocate buffer for vector field.
void
pylith::faults::FaultCohesiveLagrange::_allocateBufferVectorField(void)
{ // _allocateBufferVectorField
  assert(0 != _fields);
  if (_fields->hasField("buffer (vector)"))
    return;

  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("Output");

  // Create vector field; use same shape/chart as relative
  // displacement field.
  assert(0 != _faultMesh);
  _fields->add("buffer (vector)", "buffer");
  topology::Field<topology::SubMesh>& buffer = _fields->get("buffer (vector)");
  const topology::Field<topology::SubMesh>& dispRel = 
    _fields->get("relative disp");
  buffer.cloneSection(dispRel);
  buffer.zero();
  assert(buffer.vectorFieldType() == topology::FieldBase::VECTOR);

  logger.stagePop();
} // _allocateBufferVectorField

// ----------------------------------------------------------------------
// Allocate buffer for scalar field.
void
pylith::faults::FaultCohesiveLagrange::_allocateBufferScalarField(void)
{ // _allocateBufferScalarField
  assert(0 != _fields);
  if (_fields->hasField("buffer (scalar)"))
    return;

  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("Output");

  // Create vector field; use same shape/chart as area field.
  assert(0 != _faultMesh);
  _fields->add("buffer (scalar)", "buffer");
  topology::Field<topology::SubMesh>& buffer = _fields->get("buffer (scalar)");
  const ALE::Obj<SieveSubMesh>& faultSieveMesh = _faultMesh->sieveMesh();
  assert(!faultSieveMesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& vertices =
      faultSieveMesh->depthStratum(0);
  assert(!vertices.isNull());
  buffer.newSection(vertices, 1);
  buffer.allocate();
  buffer.vectorFieldType(topology::FieldBase::SCALAR);
  buffer.scale(1.0);
  buffer.zero();
  assert(buffer.vectorFieldType() == topology::FieldBase::SCALAR);

  logger.stagePop();
} // _allocateBufferScalarField

// ----------------------------------------------------------------------
//  Get submatrix of Jacobian matrix associated with the negative and
//  positive sides of the fault.
void
pylith::faults::FaultCohesiveLagrange::_getJacobianSubmatrixNP(
				PetscMat* jacobianSub,
				std::map<int,int>* indicesMatToSubmat,
				const topology::Jacobian& jacobian,
				const topology::SolutionFields& fields)
{ // _getJacobianSubmatrixNP
  assert(jacobianSub);
  assert(indicesMatToSubmat);

  // Get global order
  const ALE::Obj<SieveMesh>& sieveMesh = fields.mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<RealSection>& solutionSection = fields.solution().section();
  assert(!solutionSection.isNull());
  const ALE::Obj<SieveMesh::order_type>& globalOrder =
    sieveMesh->getFactory()->getGlobalOrder(sieveMesh, "default",
					    solutionSection);
  assert(!globalOrder.isNull());

  // Get Jacobian matrix
  const PetscMat jacobianMatrix = jacobian.matrix();
  assert(0 != jacobianMatrix);

  const spatialdata::geocoords::CoordSys* cs = fields.mesh().coordsys();
  assert(0 != cs);
  const int spaceDim = cs->spaceDim();

  const int numVertices = _cohesiveVertices.size();
  int numIndicesNP = 0;
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_lagrange = _cohesiveVertices[iVertex].lagrange;
    if (globalOrder->isLocal(v_lagrange))
      numIndicesNP += 2;
  } // for
  int_array indicesNP(numIndicesNP*spaceDim);

  for (int iVertex=0, indexNP=0; iVertex < numVertices; ++iVertex) {
    const int v_lagrange = _cohesiveVertices[iVertex].lagrange;
    const int v_negative = _cohesiveVertices[iVertex].negative;
    const int v_positive = _cohesiveVertices[iVertex].positive;
    
    // Compute contribution only if Lagrange constraint is local.
    if (!globalOrder->isLocal(v_lagrange))
      continue;
    
    // Set global order indices
    for(int iDim=0; iDim < spaceDim; ++iDim)
      indicesNP[indexNP*spaceDim+iDim] = 
	globalOrder->getIndex(v_negative) + iDim;
    ++indexNP;
    for(int iDim=0; iDim < spaceDim; ++iDim)
      indicesNP[indexNP*spaceDim+iDim] = 
	globalOrder->getIndex(v_positive) + iDim;
    ++indexNP;
  } // for
  
  // MatGetSubMatrices requires sorted indices
  std::sort(&indicesNP[0], &indicesNP[indicesNP.size()]);  
  
  PetscMat* subMat[1];
  IS indicesIS[1];
  PetscErrorCode err = 0;
  err = ISCreateGeneral(PETSC_COMM_SELF, indicesNP.size(), &indicesNP[0],
			PETSC_USE_POINTER, &indicesIS[0]);
  CHECK_PETSC_ERROR(err);
  err = MatGetSubMatrices(jacobianMatrix, 1, indicesIS,
			  indicesIS, MAT_INITIAL_MATRIX, subMat);
  CHECK_PETSC_ERROR(err);
  err = ISDestroy(&indicesIS[0]); CHECK_PETSC_ERROR(err);

  *jacobianSub = *subMat[0];
  err = PetscObjectReference((PetscObject) *subMat[0]); 
  CHECK_PETSC_ERROR(err);
  err = MatDestroyMatrices(1, &subMat[0]); CHECK_PETSC_ERROR(err);

  // Create map from global indices to local indices (using only the
  // first index as to match the global order.
  indicesMatToSubmat->clear();
  const int indicesNPSize = indicesNP.size();
  for (int i=0; i < indicesNPSize; i+=spaceDim)
    (*indicesMatToSubmat)[indicesNP[i]] = i;

} // _getJacobianSubmatrixNP

// ----------------------------------------------------------------------
// Adjust solution in lumped formulation to match slip for 1-D.
void
pylith::faults::FaultCohesiveLagrange::_adjustSolnLumped1D(
					  double_array* lagrangeTIncr,
					  double_array* dispTIncrN,
					  double_array* dispTIncrP,
					  const double_array& slip,
					  const double_array& orientation,
					  const double_array& dispTN,
					  const double_array& dispTP,
					  const double_array& residualN,
					  const double_array& residualP,
					  const double_array& jacobianN,
					  const double_array& jacobianP)
{ // _adjustSoln1D
  assert(0 != lagrangeTIncr);
  assert(0 != dispTIncrN);
  assert(0 != dispTIncrP);

  assert(jacobianN[0] > 0.0);
  assert(jacobianP[0] > 0.0);
  
  const double Sinv = jacobianN[0] * jacobianP[0]
    / (jacobianN[0] + jacobianP[0]);
  
  // Aru = A_i^{-1} r_i - A_j^{-1} r_j + u_i - u_j
  const double Aru = residualN[0] / jacobianN[0]
    - residualP[0] / jacobianP[0] + dispTN[0]
    - dispTP[0];
  
  // dl_k = D^{-1}( - C_{ki} Aru - d_k)
  const double Aruslip = -Aru - slip[0];
  const double dlp = Sinv * Aruslip;
  
  // Update displacements at negative vertex
  (*dispTIncrN)[0] = +1.0 / jacobianN[0] * dlp;
  
  // Update displacements at positive vertex
  (*dispTIncrP)[0] = -1.0 / jacobianP[0] * dlp;
  
  // Update Lagrange multiplier.
  (*lagrangeTIncr)[0] = dlp;

  PetscLogFlops(17);
} // _adjustSoln1D

// ----------------------------------------------------------------------
// Adjust solution in lumped formulation to match slip for 2-D.
void
pylith::faults::FaultCohesiveLagrange::_adjustSolnLumped2D(
					  double_array* lagrangeTIncr,
					  double_array* dispTIncrN,
					  double_array* dispTIncrP,
					  const double_array& slip,
					  const double_array& orientation,
					  const double_array& dispTN,
					  const double_array& dispTP,
					  const double_array& residualN,
					  const double_array& residualP,
					  const double_array& jacobianN,
					  const double_array& jacobianP)
{ // _adjustSoln2D
  assert(0 != lagrangeTIncr);
  assert(0 != dispTIncrN);
  assert(0 != dispTIncrP);

  assert(jacobianN[0] > 0.0);
  assert(jacobianN[1] > 0.0);
  assert(jacobianP[0] > 0.0);
  assert(jacobianP[1] > 0.0);
  
  const double Cpx = orientation[0];
  const double Cpy = orientation[1];
  const double Cqx = orientation[2];
  const double Cqy = orientation[3];
  
  // Check to make sure Jacobian is same at all DOF for
  // vertices i and j (means S is diagonal with equal enties).
  assert(jacobianN[0] == jacobianN[1]);
  assert(jacobianP[0] == jacobianP[1]);
  
  const double Sinv = jacobianN[0] * jacobianP[0]
    / (jacobianN[0] + jacobianP[0]);
  
  // Aru = A_i^{-1} r_i - A_j^{-1} r_j + u_i - u_j
  const double Arux = residualN[0] / jacobianN[0]
    - residualP[0] / jacobianP[0] + dispTN[0]
    - dispTP[0];
  const double Aruy = residualN[1] / jacobianN[1]
    - residualP[1] / jacobianP[1] + dispTN[1]
    - dispTP[1];
  
  // dl_k = S^{-1}(-C_{ki} Aru - d_k)
  const double Arup = Cpx * Arux + Cpy * Aruy;
  const double Aruq = Cqx * Arux + Cqy * Aruy;
  const double Arupslip = -Arup - slip[0];
  const double Aruqslip = -Aruq - slip[1];
  const double dlp = Sinv * Arupslip;
  const double dlq = Sinv * Aruqslip;
  
  const double dlx = Cpx * dlp + Cqx * dlq;
  const double dly = Cpy * dlp + Cqy * dlq;
  
  // Update displacements at negative vertex.
  (*dispTIncrN)[0] = dlx / jacobianN[0];
  (*dispTIncrN)[1] = dly / jacobianN[1];
  
  // Update displacements at positive vertex.
  (*dispTIncrP)[0] = -dlx / jacobianP[0];
  (*dispTIncrP)[1] = -dly / jacobianP[0];
  
  // Update Lagrange multiplier.
  (*lagrangeTIncr)[0] = dlp;
  (*lagrangeTIncr)[1] = dlq;

    PetscLogFlops(41);
} // _adjustSoln2D

// ----------------------------------------------------------------------
// Adjust solution in lumped formulation to match slip for 3-D.
void
pylith::faults::FaultCohesiveLagrange::_adjustSolnLumped3D(
					  double_array* lagrangeTIncr,
					  double_array* dispTIncrN,
					  double_array* dispTIncrP,
					  const double_array& slip,
					  const double_array& orientation,
					  const double_array& dispTN,
					  const double_array& dispTP,
					  const double_array& residualN,
					  const double_array& residualP,
					  const double_array& jacobianN,
					  const double_array& jacobianP)
{ // _adjustSoln3D
  assert(0 != lagrangeTIncr);
  assert(0 != dispTIncrN);
  assert(0 != dispTIncrP);

  assert(jacobianN[0] > 0.0);
  assert(jacobianN[1] > 0.0);
  assert(jacobianN[2] > 0.0);
  assert(jacobianP[0] > 0.0);
  assert(jacobianP[1] > 0.0);
  assert(jacobianP[2] > 0.0);

  const double Cpx = orientation[0];
  const double Cpy = orientation[1];
  const double Cpz = orientation[2];
  const double Cqx = orientation[3];
  const double Cqy = orientation[4];
  const double Cqz = orientation[5];
  const double Crx = orientation[6];
  const double Cry = orientation[7];
  const double Crz = orientation[8];

  // Check to make sure Jacobian is same at all DOF for
  // vertices i and j (means S is diagonal with equal enties).
  assert(jacobianN[0] == jacobianN[1] && 
	 jacobianN[0] == jacobianN[2]);
  assert(jacobianP[0] == jacobianP[1] && 
	 jacobianP[0] == jacobianP[2]);

  const double Sinv = jacobianN[0] * jacobianP[0]
    / (jacobianN[0] + jacobianP[0]);

  // Aru = A_i^{-1} r_i - A_j^{-1} r_j + u_i - u_j
  const double Arux = residualN[0] / jacobianN[0]
    - residualP[0] / jacobianP[0] + dispTN[0]
    - dispTP[0];
  const double Aruy = residualN[1] / jacobianN[1]
    - residualP[1] / jacobianP[1] + dispTN[1]
    - dispTP[1];
  const double Aruz = residualN[2] / jacobianN[2]
    - residualP[2] / jacobianP[2] + dispTN[2]
    - dispTP[2];

  // dl_k = D^{-1}( -C_{ki} Aru - d_k)
  const double Arup = Cpx * Arux + Cpy * Aruy + Cpz * Aruz;
  const double Aruq = Cqx * Arux + Cqy * Aruy + Cqz * Aruz;
  const double Arur = Crx * Arux + Cry * Aruy + Crz * Aruz;
  const double Arupslip = -Arup - slip[0];
  const double Aruqslip = -Aruq - slip[1];
  const double Arurslip = -Arur - slip[2];
  const double dlp = Sinv * Arupslip;
  const double dlq = Sinv * Aruqslip;
  const double dlr = Sinv * Arurslip;

  const double dlx = Cpx * dlp + Cqx * dlq + Crx * dlr;
  const double dly = Cpy * dlp + Cqy * dlq + Cry * dlr;
  const double dlz = Cpz * dlp + Cqz * dlq + Crz * dlr;

  // Update displacements at negative vertex.
  (*dispTIncrN)[0] = dlx / jacobianN[0];
  (*dispTIncrN)[1] = dly / jacobianN[1];
  (*dispTIncrN)[2] = dlz / jacobianN[2];

  // Update displacements at positive vertex.
  (*dispTIncrP)[0] = -dlx / jacobianP[0];
  (*dispTIncrP)[1] = -dly / jacobianP[1];
  (*dispTIncrP)[2] = -dlz / jacobianP[2];

  // Update Lagrange multiplier.
  (*lagrangeTIncr)[0] = dlp;
  (*lagrangeTIncr)[1] = dlq;
  (*lagrangeTIncr)[2] = dlr;

  PetscLogFlops(72);
} // _adjustSoln3D

// ----------------------------------------------------------------------
// Get cell field associated with integrator.
const pylith::topology::Field<pylith::topology::SubMesh>&
pylith::faults::FaultCohesiveLagrange::cellField(const char* name,
                                                 const topology::SolutionFields* fields)
{ // cellField
  if (0 == strcasecmp("distribution", name)) {
    const topology::Field<topology::SubMesh>& dist = _fields->get("distribution");
    return dist;
  }
  // Should not reach this point if requested field was found
  std::ostringstream msg;
  msg << "Request for unknown cell field '" << name << "' for fault '"
      << label() << ".";
  throw std::runtime_error(msg.str());

  // Satisfy return values
  assert(0 != _fields);
  const topology::Field<topology::SubMesh>& buffer = _fields->get(
    "buffer (vector)");
  return buffer;
} // cellField


// End of file 
