// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
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

//#define PRECOMPUTE_GEOMETRY
//#define DETAILED_EVENT_LOGGING

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::RealSection RealSection;
typedef pylith::topology::SubMesh::SieveMesh SieveSubMesh;

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
					     const double upDir[3],
					     const double normalDir[3])
{ // initialize
  assert(0 != upDir);
  assert(0 != normalDir);
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

  // Allocate slip field
  const ALE::Obj<SieveSubMesh>& faultSieveMesh = _faultMesh->sieveMesh();
  assert(!faultSieveMesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& vertices =
      faultSieveMesh->depthStratum(0);
  assert(!vertices.isNull());
  _fields->add("slip", "slip");
  topology::Field<topology::SubMesh>& slip = _fields->get("slip");
  slip.newSection(vertices, cs->spaceDim());
  slip.allocate();
  slip.vectorFieldType(topology::FieldBase::VECTOR);
  slip.scale(_normalizer->lengthScale());

  const ALE::Obj<SieveSubMesh::label_sequence>& cells =
      faultSieveMesh->heightStratum(0);
  assert(!cells.isNull());
  const SieveSubMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveSubMesh::label_sequence::iterator cellsEnd = cells->end();
  _quadrature->initializeGeometry();
#if defined(PRECOMPUTE_GEOMETRY)
  _quadrature->computeGeometry(*_faultMesh, cells);
#endif

  // Compute orientation at vertices in fault mesh.
  _calcOrientation(upDir, normalDir);

  // Compute tributary area for each vertex in fault mesh.
  _calcArea();

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
// Integrate contribution of cohesive cells to residual term that do
// not require assembly across cells, vertices, or processors.
void
pylith::faults::FaultCohesiveLagrange::integrateResidualAssembled(const topology::Field<
                                                                      topology::Mesh>& residual,
                                                                  const double t,
                                                                  topology::SolutionFields* const fields) { // integrateResidualAssembled
  assert(0 != fields);
  assert(0 != _fields);
  assert(0 != _logger);

  // Cohesive cells with normal vertices i and j, and constraint
  // vertex k make contributions to the assembled residual:
  //
  //   * DOF i and j: internal forces in soln field associated with
  //                  slip  -[C]^T{L(t)+dL(t)}
  //   * DOF k: slip values  -[C]{u(t)+dt(t)}
  //   * DOF k: slip values {D(t+dt)}

  const int setupEvent = _logger->eventId("FaIR setup");
  const int geometryEvent = _logger->eventId("FaIR geometry");
  const int computeEvent = _logger->eventId("FaIR compute");
  const int restrictEvent = _logger->eventId("FaIR restrict");
  const int updateEvent = _logger->eventId("FaIR update");

  _logger->eventBegin(setupEvent);

  // Get cell information and setup storage for cell data
  const int spaceDim = _quadrature->spaceDim();
  const int orientationSize = spaceDim * spaceDim;

  // Allocate vectors for vertex values
  double_array slipVertex(spaceDim);
  double_array orientationVertex(orientationSize);
  double_array dispTVertexN(spaceDim);
  double_array dispTVertexP(spaceDim);
  double_array dispTVertexL(spaceDim);
  double_array dispTIncrVertexN(spaceDim);
  double_array dispTIncrVertexP(spaceDim);
  double_array dispTIncrVertexL(spaceDim);
  double_array dispTpdtVertexN(spaceDim);
  double_array dispTpdtVertexP(spaceDim);
  double_array dispTpdtVertexL(spaceDim);
  double_array residualVertexN(spaceDim);
  double_array residualVertexP(spaceDim);
  double_array residualVertexL(spaceDim);

  // Get sections
  topology::Field<topology::SubMesh>& slip = _fields->get("slip");
  const ALE::Obj<RealSection>& slipSection = slip.section();
  assert(!slipSection.isNull());

  const ALE::Obj<RealSection>& residualSection = residual.section();
  assert(!residualSection.isNull());

  const ALE::Obj<RealSection>& orientationSection =
      _fields->get("orientation").section();
  assert(!orientationSection.isNull());

  topology::Field<topology::Mesh>& dispT = fields->get("disp(t)");
  const ALE::Obj<RealSection>& dispTSection = dispT.section();
  assert(!dispTSection.isNull());

  topology::Field<topology::Mesh>& dispTIncr = fields->get("dispIncr(t->t+dt)");
  const ALE::Obj<RealSection>& dispTIncrSection = dispTIncr.section();
  assert(!dispTIncrSection.isNull());

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

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(restrictEvent);
#endif

    // Get slip at fault vertex.
    slipSection->restrictPoint(v_fault, &slipVertex[0], slipVertex.size());

    // Get orientations at fault vertex.
    orientationSection->restrictPoint(v_fault, &orientationVertex[0],
				      orientationVertex.size());

    // Get disp(t) at conventional vertices and Lagrange vertex.
    dispTSection->restrictPoint(v_negative, &dispTVertexN[0],
				dispTVertexN.size());
    dispTSection->restrictPoint(v_positive, &dispTVertexP[0],
				dispTVertexP.size());
    dispTSection->restrictPoint(v_lagrange, &dispTVertexL[0],
				dispTVertexL.size());

    // Get dispIncr(t->t+dt) at conventional vertices and Lagrange vertex.
    dispTIncrSection->restrictPoint(v_negative, &dispTIncrVertexN[0],
				    dispTIncrVertexN.size());
    dispTIncrSection->restrictPoint(v_positive, &dispTIncrVertexP[0],
				    dispTIncrVertexP.size());
    dispTIncrSection->restrictPoint(v_lagrange, &dispTIncrVertexL[0],
				    dispTIncrVertexL.size());

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(restrictEvent);
    _logger->eventBegin(computeEvent);
#endif

    // Compute current estimate of displacement at time t+dt using
    // solution increment.
    dispTpdtVertexN = dispTVertexN + dispTIncrVertexN;
    dispTpdtVertexP = dispTVertexP + dispTIncrVertexP;
    dispTpdtVertexL = dispTVertexL + dispTIncrVertexL;

    // Entries associated with constraint forces applied at negative vertex
    residualVertexN = 0.0;
    for (int iDim = 0; iDim < spaceDim; ++iDim)
      for (int kDim = 0; kDim < spaceDim; ++kDim)
        residualVertexN[iDim] -= 
	  dispTpdtVertexL[kDim] * -orientationVertex[kDim*spaceDim+iDim];

    // Entries associated with constraint forces applied at positive vertex
    residualVertexP = -residualVertexN;

    // Entries associated with relative displacements between positive
    // and negative vertices for Lagrange vertex.
    residualVertexL = slipVertex;
    for (int kDim = 0; kDim < spaceDim; ++kDim)
      for (int iDim = 0; iDim < spaceDim; ++iDim)
        residualVertexL[kDim] -= (dispTpdtVertexP[iDim] - dispTpdtVertexN[iDim])
            * orientationVertex[kDim*spaceDim+iDim];

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(computeEvent);
    _logger->eventBegin(updateEvent);
#endif

    assert(residualVertexN.size() == 
	   residualSection->getFiberDimension(v_negative));
    residualSection->updateAddPoint(v_negative, &residualVertexN[0]);

    assert(residualVertexP.size() == 
	   residualSection->getFiberDimension(v_positive));
    residualSection->updateAddPoint(v_positive, &residualVertexP[0]);

    assert(residualVertexL.size() == 
	   residualSection->getFiberDimension(v_lagrange));
    residualSection->updateAddPoint(v_lagrange, &residualVertexL[0]);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(updateEvent);
#endif
  } // for
  PetscLogFlops(numVertices*spaceDim*spaceDim*8);

#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventEnd(computeEvent);
#endif
} // integrateResidualAssembled

// ----------------------------------------------------------------------
// Compute Jacobian matrix (A) associated with operator that do not
// require assembly across cells, vertices, or processors.
void
pylith::faults::FaultCohesiveLagrange::integrateJacobianAssembled(topology::Jacobian* jacobian,
                                                                  const double t,
                                                                  topology::SolutionFields* const fields)
{ // integrateJacobianAssembled
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

  // Add constraint information to Jacobian matrix; these are the
  // direction cosines. Entries are associated with vertices ik, jk,
  // ki, and kj.

  // Get cell information and setup storage for cell data
  const int spaceDim = _quadrature->spaceDim();
  const int orientationSize = spaceDim * spaceDim;

  // Allocate vectors for vertex values
  double_array orientationVertex(orientationSize);
  double_array jacobianVertex(spaceDim*spaceDim);
  int_array indicesL(spaceDim);
  int_array indicesN(spaceDim);
  int_array indicesP(spaceDim);
  int_array indicesRel(spaceDim);
  for (int i=0; i < spaceDim; ++i)
    indicesRel[i] = i;

  // Get sections
  const ALE::Obj<RealSection>& solutionSection = fields->solution().section();
  assert(!solutionSection.isNull());

  const ALE::Obj<RealSection>& orientationSection =
      _fields->get("orientation").section();
  assert(!orientationSection.isNull());

  const ALE::Obj<SieveMesh>& sieveMesh = fields->mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::order_type>& globalOrder =
      sieveMesh->getFactory()->getGlobalOrder(sieveMesh, "default",
        solutionSection);
  assert(!globalOrder.isNull());

  const PetscMat jacobianMatrix = jacobian->matrix();
  assert(0 != jacobianMatrix);

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

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(restrictEvent);
#endif

    // Get orientations at fault cell's vertices.
    orientationSection->restrictPoint(v_fault, &orientationVertex[0], orientationVertex.size());

    // Set global order indices
    indicesL = indicesRel + globalOrder->getIndex(v_lagrange);
    indicesN = indicesRel + globalOrder->getIndex(v_negative);
    indicesP = indicesRel + globalOrder->getIndex(v_positive);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(restrictEvent);
    _logger->eventBegin(updateEvent);
#endif

    // Values associated with [C]
    // Values at positive vertex, entry L,P in Jacobian
    jacobianVertex = orientationVertex;
    MatSetValues(jacobianMatrix,
                 indicesL.size(), &indicesL[0],
                 indicesP.size(), &indicesP[0],
                 &jacobianVertex[0],
                 INSERT_VALUES);

    // Values at negative vertex, entry L,N in Jacobian
    jacobianVertex *= -1.0;
    MatSetValues(jacobianMatrix,
                 indicesL.size(), &indicesL[0],
                 indicesN.size(), &indicesN[0],
                 &jacobianVertex[0],
                 INSERT_VALUES);

    // Values associated with [C]^T
    // Transpose orientation matrix
    for (int iDim=0; iDim < spaceDim; ++iDim)
      for (int jDim=0; jDim < spaceDim; ++jDim)
        jacobianVertex[iDim*spaceDim+jDim] = orientationVertex[jDim*spaceDim+iDim];
    // Values at positive vertex, entry P,L in Jacobian
    MatSetValues(jacobianMatrix,
                 indicesP.size(), &indicesP[0],
                 indicesL.size(), &indicesL[0],
                 &jacobianVertex[0],
                 INSERT_VALUES);

    // Values at negative vertex, entry L,N in Jacobian
    jacobianVertex *= -1.0;
    MatSetValues(jacobianMatrix,
                 indicesN.size(), &indicesN[0],
                 indicesL.size(), &indicesL[0],
                 &jacobianVertex[0],
                 INSERT_VALUES);

    // Values at Lagrange vertex, entry L,L in Jacobian
    jacobianVertex = 0.0;
    MatSetValues(jacobianMatrix,
                 indicesL.size(), &indicesL[0],
                 indicesL.size(), &indicesL[0],
                 &jacobianVertex[0],
                 INSERT_VALUES);

#if defined(DETAILED_EVENT_LOGGING)
    PetscLogFlops(spaceDim*spaceDim*2);
    _logger->eventBegin(updateEvent);
#endif

  } // for

#if !defined(DETAILED_EVENT_LOGGING)
  PetscLogFlops(numVertices*(spaceDim*spaceDim*2));
  _logger->eventEnd(computeEvent);
#endif

  _needNewJacobian = false;
} // integrateJacobianAssembled

// ----------------------------------------------------------------------
// Compute Jacobian matrix (A) associated with operator that do not
// require assembly across cells, vertices, or processors.
void
pylith::faults::FaultCohesiveLagrange::integrateJacobianAssembled(
	                  topology::Field<topology::Mesh>* jacobian,
			  const double t,
			  topology::SolutionFields* const fields)
{ // integrateJacobianAssembled
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

  // Add ones to diagonal Jacobian matrix (as field) for
  // convenience. Instead of including the constraints in the Jacobian
  // matrix, we adjust the solution to account for the Lagrange
  // multipliers as part of the solve.

  const int spaceDim = _quadrature->spaceDim();
  double_array jacobianVertex(spaceDim);
  jacobianVertex = 1.0;
  const ALE::Obj<RealSection>& jacobianSection = jacobian->section();
  assert(!jacobianSection.isNull());

  _logger->eventEnd(setupEvent);
#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventBegin(computeEvent);
#endif

  const int numVertices = _cohesiveVertices.size();
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_lagrange = _cohesiveVertices[iVertex].lagrange;

    assert(jacobianSection->getFiberDimension(v_lagrange) == spaceDim);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(updateEvent);
#endif
    jacobianSection->updatePoint(v_lagrange, &jacobianVertex[0]);
#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(updateEvent);
#endif
  } // for
  PetscLogFlops(0);

#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventEnd(computeEvent);
#endif

  _needNewJacobian = false;
} // integrateJacobianAssembled

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

  /** We have J = [A C^T]
   *              [C   0]
   *
   * We want to approximate C A^(-1) C^T.
   *
   * Consider Lagrange vertex L that constrains the relative
   * displacement between vertex N on the negative side of the fault
   * and vertex P on the positive side of the fault.
   *
   * If we approximate A(-1) by 1/diag(A), then we can write 
   * C A^(-1) C^T for a 2-D case as
   *
   * [-R00 -R01  R00 R01][Ai_nx 0      0     0    ][-R00 -R10]
   * [-R10 -R11  R10 R11][      Ai_ny  0     0    ][-R01 -R11]
   *                     [      0      Ai_px 0    ][ R00  R10]
   *                     [                   Ai_py][ R01  R11]
   *
   * where
   *
   * Ai_nx is the inverse of the diag(A) for DOF x of vertex N
   * Ai_ny is the inverse of the diag(A) for DOF y of vertex N
   * Ai_px is the inverse of the diag(A) for DOF x of vertex P
   * Ai_py is the inverse of the diag(A) for DOF y of vertex P
   *
   * If Ai_nx == Ai_ny and Ai_px == Ai_py, then the result is
   * diagonal. Otherwise, the off-diagonal terms will be nonzero,
   * but we expect them to be small. Since we already approximate
   * the inverse of A by the inverse of the diagonal, we drop the
   * off-diagonal terms of C A^(-1) C^T:
   *
   * Term for DOF x of vertex L is: 
   * R00^2 (Ai_nx + Ai_px) + R01^2 (Ai_ny + Ai_py)
   *
   * Term for DOF y of vertex L is: 
   * R10^2 (Ai_nx + Ai_px) + R11^2 (Ai_ny + Ai_py)
   *
   * Translate DOF for global vertex L into DOF in split field for
   * preconditioner.
   */

  const int setupEvent = _logger->eventId("FaPr setup");
  const int computeEvent = _logger->eventId("FaPr compute");
  const int restrictEvent = _logger->eventId("FaPr restrict");
  const int updateEvent = _logger->eventId("FaPr update");

  _logger->eventBegin(setupEvent);

  // Get cell information and setup storage for cell data
  const int spaceDim = _quadrature->spaceDim();
  const int orientationSize = spaceDim * spaceDim;

  // Allocate vectors for vertex values
  double_array orientationVertex(orientationSize);
  double_array jacobianVertexP(spaceDim*spaceDim);
  double_array jacobianVertexN(spaceDim*spaceDim);
  double_array jacobianInvVertexP(spaceDim);
  double_array jacobianInvVertexN(spaceDim);
  double_array precondVertexL(spaceDim);
  int_array indicesN(spaceDim);
  int_array indicesP(spaceDim);
  int_array indicesRel(spaceDim);
  for (int i=0; i < spaceDim; ++i)
    indicesRel[i] = i;

  // Get sections
  const ALE::Obj<RealSection>& solutionSection = fields->solution().section();
  assert(!solutionSection.isNull());

  const ALE::Obj<RealSection>& orientationSection =
      _fields->get("orientation").section();
  assert(!orientationSection.isNull());

  const ALE::Obj<SieveMesh>& sieveMesh = fields->mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::order_type>& globalOrder =
      sieveMesh->getFactory()->getGlobalOrder(sieveMesh, "default",
        solutionSection);
  assert(!globalOrder.isNull());

  // Get Jacobian matrix
  const PetscMat jacobianMatrix = jacobian->matrix();
  assert(0 != jacobianMatrix);

  const ALE::Obj<SieveMesh::order_type>& lagrangeGlobalOrder =
      sieveMesh->getFactory()->getGlobalOrder(sieveMesh, "faultDefault",
        solutionSection, spaceDim);
  assert(!lagrangeGlobalOrder.isNull());

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

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(restrictEvent);
#endif

    // Get orientations at fault cell's vertices.
    orientationSection->restrictPoint(v_fault, &orientationVertex[0],
				      orientationVertex.size());

    // Set global order indices
    indicesN = indicesRel + globalOrder->getIndex(v_negative);
    indicesP = indicesRel + globalOrder->getIndex(v_positive);

    PetscErrorCode err = 0;

    err = MatGetValues(jacobianMatrix,
		       indicesN.size(), &indicesN[0],
		       indicesN.size(), &indicesN[0],
		       &jacobianVertexN[0]); CHECK_PETSC_ERROR(err);
    err = MatGetValues(jacobianMatrix,
		       indicesP.size(), &indicesP[0],
		       indicesP.size(), &indicesP[0],
		       &jacobianVertexP[0]); CHECK_PETSC_ERROR(err);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(restrictEvent);
    _logger->eventBegin(computeEvent);
#endif

    // Compute inverse of Jacobian diagonals
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      jacobianInvVertexN[iDim] = 1.0/jacobianVertexN[iDim*spaceDim+iDim];
      jacobianInvVertexP[iDim] = 1.0/jacobianVertexP[iDim*spaceDim+iDim];
    } // for

    // Compute [C] [Adiag]^(-1) [C]^T
    precondVertexL = 0.0;
    for (int kDim=0; kDim < spaceDim; ++kDim)
      for (int iDim=0; iDim < spaceDim; ++iDim)
	precondVertexL[kDim] += 
	  orientationVertex[kDim*spaceDim+iDim] * 
	  orientationVertex[kDim*spaceDim+iDim] * 
	  (jacobianInvVertexN[iDim] + jacobianInvVertexP[iDim]);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(computeEvent);
    _logger->eventBegin(updateEvent);
#endif

    // Set global preconditioner index associated with Lagrange constraint vertex.
    const int indexLprecond = lagrangeGlobalOrder->getIndex(v_lagrange);
    
    // Set diagonal entries in preconditioned matrix.
    for (int iDim=0; iDim < spaceDim; ++iDim)
      MatSetValue(*precondMatrix,
		  indexLprecond + iDim,
		  indexLprecond + iDim,
		  precondVertexL[iDim],
		  INSERT_VALUES);
    
#if defined(DETAILED_EVENT_LOGGING)
    PetscLogFlops(spaceDim*spaceDim*4);
    _logger->eventBegin(updateEvent);
#endif

  } // for

#if !defined(DETAILED_EVENT_LOGGING)
  PetscLogFlops(numVertices*(spaceDim*spaceDim*4));
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
  /// Member prototype for _adjustSolnLumpedXD()
  typedef void (pylith::faults::FaultCohesiveLagrange::*adjustSolnLumped_fn_type)
    (double_array*, double_array*, double_array*,
     const double_array&, const double_array&, 
     const double_array&, const double_array&, 
     const double_array&, const double_array&, 
     const double_array&, const double_array&);

  assert(0 != fields);
  assert(0 != _quadrature);

  // Cohesive cells with conventional vertices i and j, and constraint
  // vertex k require 2 adjustments to the solution:
  //
  //   * DOF k: Compute increment in Lagrange multipliers
  //            dl_k = S^{-1} (-C_ki (A_i^{-1} r_i - C_kj A_j^{-1} r_j + u_i - u_j) - d_k)
  //            S = C_ki (A_i^{-1} + A_j^{-1}) C_ki^T
  //
  //   * DOF i and j: Adjust displacement increment (solution) to create slip
  //     consistent with Lagrange multiplier constraints
  //            du_i = +A_i^-1 C_ki^T dlk
  //            du_j = -A_j^-1 C_kj^T dlk

  const int setupEvent = _logger->eventId("FaAS setup");
  const int geometryEvent = _logger->eventId("FaAS geometry");
  const int computeEvent = _logger->eventId("FaAS compute");
  const int restrictEvent = _logger->eventId("FaAS restrict");
  const int updateEvent = _logger->eventId("FaAS update");

  _logger->eventBegin(setupEvent);

  // Get cell information and setup storage for cell data
  const int spaceDim = _quadrature->spaceDim();
  const int orientationSize = spaceDim * spaceDim;

  // Get section information
  double_array orientationVertex(orientationSize);
  const ALE::Obj<RealSection>& orientationSection =
      _fields->get("orientation").section();
  assert(!orientationSection.isNull());

  double_array slipVertex(spaceDim);
  const ALE::Obj<RealSection>& slipSection = _fields->get("slip").section();
  assert(!slipSection.isNull());

  double_array jacobianVertexN(spaceDim);
  double_array jacobianVertexP(spaceDim);
  const ALE::Obj<RealSection>& jacobianSection = jacobian.section();
  assert(!jacobianSection.isNull());

  double_array residualVertexN(spaceDim);
  double_array residualVertexP(spaceDim);
  const ALE::Obj<RealSection>& residualSection =
      fields->get("residual").section();

  double_array dispTVertexN(spaceDim);
  double_array dispTVertexP(spaceDim);
  const ALE::Obj<RealSection>& dispTSection = fields->get("disp(t)").section();

  double_array dispTIncrVertexN(spaceDim);
  double_array dispTIncrVertexP(spaceDim);
  double_array lagrangeIncrVertex(spaceDim);
  const ALE::Obj<RealSection>& dispTIncrSection = fields->get(
    "dispIncr(t->t+dt)").section();

  adjustSolnLumped_fn_type adjustSolnLumpedFn;
  switch (spaceDim) { // switch
  case 1:
    adjustSolnLumpedFn = 
      &pylith::faults::FaultCohesiveLagrange::_adjustSolnLumped1D;
    break;
  case 2: 
    adjustSolnLumpedFn = 
      &pylith::faults::FaultCohesiveLagrange::_adjustSolnLumped2D;
    break;
  case 3:
    adjustSolnLumpedFn = 
      &pylith::faults::FaultCohesiveLagrange::_adjustSolnLumped3D;
    break;
  default :
    assert(0);
    throw std::logic_error("Unknown spatial dimension in "
			   "FaultCohesiveLagrange::adjustSolnLumped.");
  } // switch

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

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(restrictEvent);
#endif

    // Get orientations at fault cell's vertices.
    orientationSection->restrictPoint(v_fault, &orientationVertex[0],
				      orientationVertex.size());

    // Get slip at fault cell's vertices.
    slipSection->restrictPoint(v_fault, &slipVertex[0], slipVertex.size());

    // Get residual at cohesive cell's vertices.
    residualSection->restrictPoint(v_negative, &residualVertexN[0],
			   residualVertexN.size());
    residualSection->restrictPoint(v_positive, &residualVertexP[0], 
				   residualVertexP.size());
    
    // Get jacobian at cohesive cell's vertices.
    jacobianSection->restrictPoint(v_negative, &jacobianVertexN[0], 
				   jacobianVertexN.size());
    jacobianSection->restrictPoint(v_positive, &jacobianVertexP[0], 
				   jacobianVertexP.size());

    // Get disp(t) at cohesive cell's vertices.
    dispTSection->restrictPoint(v_negative, &dispTVertexN[0],
				dispTVertexN.size());
    dispTSection->restrictPoint(v_positive, &dispTVertexP[0],
				dispTVertexP.size());

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(restrictEvent);
    _logger->eventBegin(computeEvent);
#endif

    CALL_MEMBER_FN(*this, 
		   adjustSolnLumpedFn)(&lagrangeIncrVertex, 
				       &dispTIncrVertexN, &dispTIncrVertexP,
				       slipVertex, orientationVertex, 
				       dispTVertexN, dispTVertexP,
				       residualVertexN, residualVertexP,
				       jacobianVertexN, jacobianVertexP);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(computeEvent);
    _logger->eventBegin(updateEvent);
#endif

    // Adjust displacements to account for Lagrange multiplier values
    // (assumed to be zero in perliminary solve).
    assert(dispTIncrVertexN.size() == 
	   dispTIncrSection->getFiberDimension(v_negative));
    dispTIncrSection->updateAddPoint(v_negative, &dispTIncrVertexN[0]);

    assert(dispTIncrVertexP.size() == 
	   dispTIncrSection->getFiberDimension(v_positive));
    dispTIncrSection->updateAddPoint(v_positive, &dispTIncrVertexP[0]);

    // Set Lagrange multiplier value. Value from preliminary solve is
    // bogus due to artificial diagonal entry of 1.0.
    assert(lagrangeIncrVertex.size() == 
	   dispTIncrSection->getFiberDimension(v_lagrange));
    dispTIncrSection->updatePoint(v_lagrange, &lagrangeIncrVertex[0]);

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
  const int numCorners = _quadrature->refGeometry().numCorners();
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

  SubMesh::renumbering_type& renumbering = faultSieveMesh->getRenumbering();
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
  typedef ALE::SieveAlg<SieveMesh> SieveAlg;

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
    const Mesh::point_type *cone = ncV.getPoints();
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
pylith::faults::FaultCohesiveLagrange::_calcOrientation(const double upDir[3],
						   const double normalDir[3])
{ // _calcOrientation
  assert(0 != upDir);
  assert(0 != normalDir);
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
  const topology::Field<topology::SubMesh>& slip = _fields->get("slip");
  orientation.newSection(slip, orientationSize);
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
  topology::Mesh::RestrictVisitor coordinatesVisitor(*coordinatesSection,
    coordinatesCell.size(), &coordinatesCell[0]);

  // Set orientation function
  assert(cohesiveDim == _quadrature->cellDim());
  assert(spaceDim == _quadrature->spaceDim());

  // Loop over cohesive cells, computing orientation weighted by
  // jacobian at constraint vertices

  const ALE::Obj<SieveSubMesh::sieve_type>& sieve = faultSieveMesh->getSieve();
  assert(!sieve.isNull());
  typedef ALE::SieveAlg<SieveSubMesh> SieveAlg;

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
    const Mesh::point_type *cone = ncV.getPoints();

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

  if (2 == cohesiveDim && vertices->size() > 0) {
    // Check orientation of first vertex, if dot product of fault
    // normal with preferred normal is negative, flip up/down dip
    // direction.
    //
    // If the user gives the correct normal direction (points from
    // footwall to ahanging wall), we should end up with
    // left-lateral-slip, reverse-slip, and fault-opening for positive
    // slip values.
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
    double_array normalDirVertex(&orientationVertex[6], 3);
    const double normalDot = normalDir[0] * normalDirVertex[0] + normalDir[1]
        * normalDirVertex[1] + normalDir[2] * normalDirVertex[2];

    const int istrike = 0;
    const int idip = 3;
    const int inormal = 6;
    if (normalDot < 0.0) {
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
void
pylith::faults::FaultCohesiveLagrange::_calcArea(void)
{ // _calcArea
  assert(0 != _faultMesh);
  assert(0 != _fields);

  // Containers for area information
  const int cellDim = _quadrature->cellDim();
  const int numBasis = _quadrature->numBasis();
  const int numQuadPts = _quadrature->numQuadPts();
  const int spaceDim = _quadrature->spaceDim();
  const feassemble::CellGeometry& cellGeometry = _quadrature->refGeometry();
  const double_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == numQuadPts);
  double jacobianDet = 0;
  double_array areaCell(numBasis);

  // Get vertices in fault mesh.
  const ALE::Obj<SieveSubMesh>& faultSieveMesh = _faultMesh->sieveMesh();
  assert(!faultSieveMesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& vertices =
      faultSieveMesh->depthStratum(0);
  const SieveSubMesh::label_sequence::iterator verticesBegin =
      vertices->begin();
  const SieveSubMesh::label_sequence::iterator verticesEnd = vertices->end();

  // Allocate area field.
  _fields->add("area", "area");

  topology::Field<topology::SubMesh>& area = _fields->get("area");
  const topology::Field<topology::SubMesh>& slip = _fields->get("slip");
  area.newSection(slip, 1);
  area.allocate();
  area.zero();
  const ALE::Obj<RealSection>& areaSection = area.section();
  assert(!areaSection.isNull());
  topology::Mesh::UpdateAddVisitor areaVisitor(*areaSection, &areaCell[0]);

  double_array coordinatesCell(numBasis * spaceDim);
  const ALE::Obj<RealSection>& coordinates = faultSieveMesh->getRealSection(
    "coordinates");
  assert(!coordinates.isNull());
  topology::Mesh::RestrictVisitor coordsVisitor(*coordinates,
    coordinatesCell.size(), &coordinatesCell[0]);

  const ALE::Obj<SieveSubMesh::label_sequence>& cells =
      faultSieveMesh->heightStratum(0);
  assert(!cells.isNull());
  const SieveSubMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveSubMesh::label_sequence::iterator cellsEnd = cells->end();

  // Loop over cells in fault mesh, compute area
  for (SieveSubMesh::label_sequence::iterator c_iter = cellsBegin; c_iter
      != cellsEnd; ++c_iter) {
    areaCell = 0.0;

    // Compute geometry information for current cell
#if defined(PRECOMPUTE_GEOMETRY)
    _quadrature->retrieveGeometry(*c_iter);
#else
    coordsVisitor.clear();
    faultSieveMesh->restrictClosure(*c_iter, coordsVisitor);
    _quadrature->computeGeometry(coordinatesCell, *c_iter);
#endif

    // Get cell geometry information that depends on cell
    const double_array& basis = _quadrature->basis();
    const double_array& jacobianDet = _quadrature->jacobianDet();

    // Compute area
    for (int iQuad = 0; iQuad < numQuadPts; ++iQuad) {
      const double wt = quadWts[iQuad] * jacobianDet[iQuad];
      for (int iBasis = 0; iBasis < numBasis; ++iBasis) {
        const double dArea = wt * basis[iQuad * numBasis + iBasis];
        areaCell[iBasis] += dArea;
      } // for
    } // for
    areaVisitor.clear();
    faultSieveMesh->updateClosure(*c_iter, areaVisitor);

    PetscLogFlops( numQuadPts*(1+numBasis*2) );
  } // for

  // Assemble area information
  area.complete();

#if 0 // DEBUGGING
  area.view("AREA");
  //_faultMesh->getSendOverlap()->view("Send fault overlap");
  //_faultMesh->getRecvOverlap()->view("Receive fault overlap");
#endif
} // _calcArea

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
  const ALE::Obj<RealSection>& areaSection = _fields->get("area").section();
  assert(!areaSection.isNull());
  const ALE::Obj<RealSection>& dispTSection = dispT.section();
  assert(!dispTSection.isNull());

  // Allocate buffer for tractions field (if necessary).
  const ALE::Obj<RealSection>& tractionsSection = tractions->section();
  if (tractionsSection.isNull()) {
    ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
    //logger.stagePush("Fault");

    const topology::Field<topology::SubMesh>& slip = _fields->get("slip");
    tractions->cloneSection(slip);

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
    assert(spaceDim == tractionsSection->getFiberDimension(v_fault));
    assert(1 == areaSection->getFiberDimension(v_fault));

    const double* dispTVertex = dispTSection->restrictPoint(v_lagrange);
    assert(0 != dispTVertex);
    const double* areaVertex = areaSection->restrictPoint(v_fault);
    assert(0 != areaVertex);

    for (int i = 0; i < spaceDim; ++i)
      tractionsVertex[i] = dispTVertex[i] / areaVertex[0];

    assert(tractionsVertex.size() == tractionsSection->getFiberDimension(v_fault));
    tractionsSection->updatePoint(v_fault, &tractionsVertex[0]);
  } // for

  PetscLogFlops(numVertices * (1 + spaceDim) );

#if 0 // DEBUGGING
  tractions->view("TRACTIONS");
#endif
} // _calcTractionsChange

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

  // Create vector field; use same shape/chart as slip field.
  assert(0 != _faultMesh);
  _fields->add("buffer (vector)", "buffer");
  topology::Field<topology::SubMesh>& buffer = _fields->get("buffer (vector)");
  const topology::Field<topology::SubMesh>& slip = _fields->get("slip");
  buffer.cloneSection(slip);
  buffer.zero();

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
  const topology::Field<topology::SubMesh>& area = _fields->get("area");
  buffer.cloneSection(area);
  buffer.zero();

  logger.stagePop();
} // _allocateBufferScalarField

// ----------------------------------------------------------------------
// Adjust solution in lumped formulation to match slip for 1-D.
void
pylith::faults::FaultCohesiveLagrange::_adjustSolnLumped1D(
					  double_array* lagrangeIncr,
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
  assert(0 != lagrangeIncr);
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
  (*lagrangeIncr)[0] = dlp;

  PetscLogFlops(17);
} // _adjustSoln1D

// ----------------------------------------------------------------------
// Adjust solution in lumped formulation to match slip for 2-D.
void
pylith::faults::FaultCohesiveLagrange::_adjustSolnLumped2D(
					  double_array* lagrangeIncr,
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
  assert(0 != lagrangeIncr);
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
  (*lagrangeIncr)[0] = dlp;
  (*lagrangeIncr)[1] = dlq;

    PetscLogFlops(41);
} // _adjustSoln2D

// ----------------------------------------------------------------------
// Adjust solution in lumped formulation to match slip for 3-D.
void
pylith::faults::FaultCohesiveLagrange::_adjustSolnLumped3D(
					  double_array* lagrangeIncr,
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
  assert(0 != lagrangeIncr);
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
  (*lagrangeIncr)[0] = dlp;
  (*lagrangeIncr)[1] = dlq;
  (*lagrangeIncr)[2] = dlr;

  PetscLogFlops(72);
} // _adjustSoln3D


// End of file 
