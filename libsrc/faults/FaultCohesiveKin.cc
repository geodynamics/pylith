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

#include "FaultCohesiveKin.hh" // implementation of object methods
#include "EqKinSrc.hh" // USES EqKinSrc
#include "CohesiveTopology.hh" // USES CohesiveTopology
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/feassemble/CellGeometry.hh" // USES CellGeometry
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/SubMesh.hh" // USES SubMesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/Jacobian.hh" // USES Jacobian
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields
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

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::RealSection RealSection;
typedef pylith::topology::SubMesh::SieveMesh SieveSubMesh;

// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::FaultCohesiveKin::FaultCohesiveKin(void) { // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::faults::FaultCohesiveKin::~FaultCohesiveKin(void) { // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void pylith::faults::FaultCohesiveKin::deallocate(void) { // deallocate
  FaultCohesive::deallocate();

  // :TODO: Use shared pointers for earthquake sources
} // deallocate

// ----------------------------------------------------------------------
// Set kinematic earthquake source.
void pylith::faults::FaultCohesiveKin::eqsrcs(const char* const * names,
                                              const int numNames,
                                              EqKinSrc** sources,
                                              const int numSources) { // eqsrcs
  assert(numNames == numSources);

  // :TODO: Use shared pointers for earthquake sources
  _eqSrcs.clear();
  for (int i = 0; i < numSources; ++i) {
    if (0 == sources[i])
      throw std::runtime_error("Null earthquake source.");
    _eqSrcs[std::string(names[i])] = sources[i];
  } // for
} // eqsrcs

// ----------------------------------------------------------------------
// Initialize fault. Determine orientation and setup boundary
void pylith::faults::FaultCohesiveKin::initialize(const topology::Mesh& mesh,
                                                  const double upDir[3],
                                                  const double normalDir[3]) { // initialize
  assert(0 != upDir);
  assert(0 != normalDir);
  assert(0 != _quadrature);
  assert(0 != _normalizer);

  _initializeLogger();

  const spatialdata::geocoords::CoordSys* cs = mesh.coordsys();
  assert(0 != cs);

  delete _faultMesh;
  _faultMesh = new topology::SubMesh();
  CohesiveTopology::createFaultParallel(_faultMesh, &_cohesiveToFault, mesh, id(), useLagrangeConstraints());

  _initializeCohesiveInfo(mesh);

  delete _fields;
  _fields = new topology::Fields<topology::Field<topology::SubMesh> >(
    *_faultMesh);

  const srcs_type::const_iterator srcsEnd = _eqSrcs.end();
  for (srcs_type::iterator s_iter = _eqSrcs.begin(); s_iter != srcsEnd; ++s_iter) {
    EqKinSrc* src = s_iter->second;
    assert(0 != src);
    src->initialize(*_faultMesh, *_normalizer);
  } // for

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
void pylith::faults::FaultCohesiveKin::splitField(topology::Field<
    topology::Mesh>* field) { // splitFields
  assert(0 != field);

  const ALE::Obj<RealSection>& section = field->section();
  assert(!section.isNull());
  if (0 == section->getNumSpaces())
    return; // Return if there are no fibrations

  const int fibrationDisp = 0;
  const int fibrationLagrange = 1;

  // Get domain Sieve mesh
  const ALE::Obj<SieveMesh>& sieveMesh = field->mesh().sieveMesh();
  assert(!sieveMesh.isNull());

  // Get fault Sieve mesh
  const ALE::Obj<SieveSubMesh>& faultSieveMesh = _faultMesh->sieveMesh();
  assert(!faultSieveMesh.isNull());

  const ALE::Obj<SieveMesh::label_sequence>& vertices =
      sieveMesh->depthStratum(0);
  assert(!vertices.isNull());
  const SieveSubMesh::label_sequence::iterator verticesBegin =
      vertices->begin();
  const SieveSubMesh::label_sequence::iterator verticesEnd = vertices->end();
  SieveSubMesh::renumbering_type& renumbering =
      faultSieveMesh->getRenumbering();
  const SieveSubMesh::renumbering_type::const_iterator renumberingEnd =
      renumbering.end();
  for (SieveSubMesh::label_sequence::iterator v_iter = verticesBegin; v_iter
      != verticesEnd; ++v_iter)
    if (renumbering.find(*v_iter) != renumberingEnd) {
      const int vertexFault = renumbering[*v_iter];
      const int vertexMesh = *v_iter;
      const int fiberDim = section->getFiberDimension(vertexMesh);
      assert(fiberDim > 0);
      // Reset displacement fibration fiber dimension to zero.
      section->setFiberDimension(vertexMesh, 0, fibrationDisp);
      // Set Langrange fibration fiber dimension.
      section->setFiberDimension(vertexMesh, fiberDim, fibrationLagrange);
    } // if
} // splitFields

// ----------------------------------------------------------------------
// Integrate contribution of cohesive cells to residual term that
// require assembly across processors.
void pylith::faults::FaultCohesiveKin::integrateResidual(const topology::Field<
                                                             topology::Mesh>& residual,
                                                         const double t,
                                                         topology::SolutionFields* const fields) { // integrateResidual
  assert(0 != fields);
  assert(0 != _quadrature);
  assert(0 != _fields);
  assert(0 != _logger);

  // Cohesive cells with normal vertices i and j, and constraint
  // vertex k make 2 contributions to the residual:
  //
  //   * DOF i and j: internal forces in soln field associated with 
  //                  slip  -[C]^T{L(t)+dL(t)}
  //   * DOF k: slip values  -[C]{u(t)+dt(t)}

  const int setupEvent = _logger->eventId("FkIR setup");
  const int geometryEvent = _logger->eventId("FkIR geometry");
  const int computeEvent = _logger->eventId("FkIR compute");
  const int restrictEvent = _logger->eventId("FkIR restrict");
  const int updateEvent = _logger->eventId("FkIR update");

  _logger->eventBegin(setupEvent);

  // Get cell information and setup storage for cell data
  const int spaceDim = _quadrature->spaceDim();
  const int orientationSize = spaceDim * spaceDim;
  const int numBasis = _quadrature->numBasis();
  const int numConstraintVert = numBasis;
  const int numCorners = 3 * numConstraintVert; // cohesive cell
  const int numQuadPts = _quadrature->numQuadPts();
  const double_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == numQuadPts);

  // Allocate vectors for cell values
  double_array orientationCell(numConstraintVert * orientationSize);
  double_array dispTCell(numCorners * spaceDim);
  double_array dispTIncrCell(numCorners * spaceDim);
  double_array dispTpdtCell(numCorners * spaceDim);
  double_array residualCell(numCorners * spaceDim);

  // Tributary area for the current for each vertex.
  double_array areaVertexCell(numConstraintVert);

  // Total fault area associated with each vertex (assembled over all cells).
  double_array areaCell(numConstraintVert);

  // Get cohesive cells
  const ALE::Obj<SieveMesh>& sieveMesh = residual.mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& cellsCohesive =
      sieveMesh->getLabelStratum("material-id", id());
  assert(!cellsCohesive.isNull());
  const SieveMesh::label_sequence::iterator cellsCohesiveBegin =
      cellsCohesive->begin();
  const SieveMesh::label_sequence::iterator cellsCohesiveEnd =
      cellsCohesive->end();
  const int cellsCohesiveSize = cellsCohesive->size();

  // Get fault Sieve mesh
  const ALE::Obj<SieveSubMesh>& faultSieveMesh = _faultMesh->sieveMesh();
  assert(!faultSieveMesh.isNull());

  // Get section information
  const ALE::Obj<RealSection>& orientationSection =
      _fields->get("orientation").section();
  assert(!orientationSection.isNull());
  topology::Mesh::RestrictVisitor orientationVisitor(*orientationSection,
    orientationCell.size(), &orientationCell[0]);

  const ALE::Obj<RealSection>& areaSection = _fields->get("area").section();
  assert(!areaSection.isNull());
  topology::Mesh::RestrictVisitor areaVisitor(*areaSection, areaCell.size(),
    &areaCell[0]);

  topology::Field<topology::Mesh>& dispT = fields->get("disp(t)");
  const ALE::Obj<RealSection>& dispTSection = dispT.section();
  assert(!dispTSection.isNull());
  topology::Mesh::RestrictVisitor dispTVisitor(*dispTSection, dispTCell.size(),
    &dispTCell[0]);

  topology::Field<topology::Mesh>& dispTIncr = fields->get("dispIncr(t->t+dt)");
  const ALE::Obj<RealSection>& dispTIncrSection = dispTIncr.section();
  assert(!dispTIncrSection.isNull());
  topology::Mesh::RestrictVisitor dispTIncrVisitor(*dispTIncrSection,
    dispTIncrCell.size(), &dispTIncrCell[0]);

  const ALE::Obj<RealSection>& residualSection = residual.section();
  topology::Mesh::UpdateAddVisitor residualVisitor(*residualSection,
    &residualCell[0]);

  double_array coordinatesCell(numBasis * spaceDim);
  const ALE::Obj<RealSection>& coordinates = faultSieveMesh->getRealSection(
    "coordinates");
  assert(!coordinates.isNull());
  topology::Mesh::RestrictVisitor coordsVisitor(*coordinates,
    coordinatesCell.size(), &coordinatesCell[0]);

  _logger->eventEnd(setupEvent);

  for (SieveMesh::label_sequence::iterator c_iter = cellsCohesiveBegin; c_iter
      != cellsCohesiveEnd; ++c_iter) {
    const SieveMesh::point_type c_fault = _cohesiveToFault[*c_iter];
    areaVertexCell = 0.0;
    residualCell = 0.0;

    // Compute geometry information for current cell
    _logger->eventBegin(geometryEvent);
#if defined(PRECOMPUTE_GEOMETRY)
    _quadrature->retrieveGeometry(c_fault);
#else
    coordsVisitor.clear();
    faultSieveMesh->restrictClosure(c_fault, coordsVisitor);
    _quadrature->computeGeometry(coordinatesCell, c_fault);
#endif
    _logger->eventEnd(geometryEvent);

    // Get cell geometry information that depends on cell
    const double_array& basis = _quadrature->basis();
    const double_array& jacobianDet = _quadrature->jacobianDet();

    _logger->eventBegin(restrictEvent);

    // Get orientations at fault cell's vertices.
    orientationVisitor.clear();
    faultSieveMesh->restrictClosure(c_fault, orientationVisitor);

    // Get area at fault cell's vertices.
    areaVisitor.clear();
    faultSieveMesh->restrictClosure(c_fault, areaVisitor);

    // Get disp(t) at cohesive cell's vertices.
    dispTVisitor.clear();
    sieveMesh->restrictClosure(*c_iter, dispTVisitor);

    // Get dispIncr(t) at cohesive cell's vertices.
    dispTIncrVisitor.clear();
    sieveMesh->restrictClosure(*c_iter, dispTIncrVisitor);

    _logger->eventEnd(restrictEvent);
    _logger->eventBegin(computeEvent);

    // Compute current estimate of displacement at time t+dt using
    // solution increment.
    dispTpdtCell = dispTCell + dispTIncrCell;

    // Compute contributory area for cell (to weight contributions)
    for (int iQuad = 0; iQuad < numQuadPts; ++iQuad) {
      const double wt = quadWts[iQuad] * jacobianDet[iQuad];
      for (int iBasis = 0; iBasis < numBasis; ++iBasis) {
        const double dArea = wt * basis[iQuad * numBasis + iBasis];
        areaVertexCell[iBasis] += dArea;
      } // for
    } // for

    for (int iConstraint = 0; iConstraint < numConstraintVert; ++iConstraint) {
      // Blocks in cell matrix associated with normal cohesive
      // vertices i and j and constraint vertex k
      const int indexI = iConstraint;
      const int indexJ = iConstraint + numConstraintVert;
      const int indexK = iConstraint + 2 * numConstraintVert;

      assert(areaCell[iConstraint] > 0);
      const double wt = areaVertexCell[iConstraint] / areaCell[iConstraint];

      // Get orientation at constraint vertex
      const double* orientationVertex = &orientationCell[iConstraint
          * orientationSize];
      assert(0 != orientationVertex);

      // Entries associated with constraint forces applied at node i
      for (int iDim = 0; iDim < spaceDim; ++iDim) {
        for (int kDim = 0; kDim < spaceDim; ++kDim)
          residualCell[indexI * spaceDim + iDim] -= dispTpdtCell[indexK
              * spaceDim + kDim] * -orientationVertex[kDim * spaceDim + iDim]
              * wt;
      } // for

      // Entries associated with constraint forces applied at node j
      for (int jDim = 0; jDim < spaceDim; ++jDim) {
        for (int kDim = 0; kDim < spaceDim; ++kDim)
          residualCell[indexJ * spaceDim + jDim] -= dispTpdtCell[indexK
              * spaceDim + kDim] * orientationVertex[kDim * spaceDim + jDim]
              * wt;
      } // for

      // Entries associated with relative displacements between node i
      // and node j for constraint node k
      for (int kDim = 0; kDim < spaceDim; ++kDim) {
        for (int iDim = 0; iDim < spaceDim; ++iDim)
          residualCell[indexK * spaceDim + kDim] -= (dispTpdtCell[indexJ
              * spaceDim + iDim] - dispTpdtCell[indexI * spaceDim + iDim])
              * orientationVertex[kDim * spaceDim + iDim] * wt;
      } // for
    } // for
    _logger->eventEnd(computeEvent);

#if 0 // DEBUGGING
    std::cout << "Updating fault residual for cell " << *c_iter << std::endl;
    for(int i = 0; i < numCorners*spaceDim; ++i) {
      std::cout << "  dispTpdt["<<i<<"]: " << dispTpdtCell[i] << std::endl;
    }
    for(int i = 0; i < numCorners*spaceDim; ++i) {
      std::cout << "  dispT["<<i<<"]: " << dispTCell[i] << std::endl;
    }
    for(int i = 0; i < numCorners*spaceDim; ++i) {
      std::cout << "  dispIncr["<<i<<"]: " << dispTIncrCell[i] << std::endl;
    }
    for(int i = 0; i < numCorners*spaceDim; ++i) {
      std::cout << "  v["<<i<<"]: " << residualCell[i] << std::endl;
    }
#endif

    _logger->eventBegin(updateEvent);
    residualVisitor.clear();
    sieveMesh->updateClosure(*c_iter, residualVisitor);
    _logger->eventEnd(updateEvent);
  } // for

  // FIX THIS
  PetscLogFlops(cellsCohesiveSize*numConstraintVert*spaceDim*spaceDim*7);
} // integrateResidual

// ----------------------------------------------------------------------
// Integrate contribution of cohesive cells to residual term that do
// not require assembly across cells, vertices, or processors.
void pylith::faults::FaultCohesiveKin::integrateResidualAssembled(const topology::Field<
                                                                      topology::Mesh>& residual,
                                                                  const double t,
                                                                  topology::SolutionFields* const fields) { // integrateResidualAssembled
  assert(0 != fields);
  assert(0 != _fields);
  assert(0 != _logger);

  // Cohesive cells with normal vertices i and j, and constraint
  // vertex k make contributions to the assembled residual:
  //
  //   * DOF k: slip values {D(t+dt)}

  const int setupEvent = _logger->eventId("FkIR setup");
  const int geometryEvent = _logger->eventId("FkIR geometry");
  const int computeEvent = _logger->eventId("FkIR compute");
  const int restrictEvent = _logger->eventId("FkIR restrict");
  const int updateEvent = _logger->eventId("FkIR update");

  _logger->eventBegin(setupEvent);

  topology::Field<topology::SubMesh>& slip = _fields->get("slip");
  slip.zero();
  // Compute slip field at current time step
  const srcs_type::const_iterator srcsEnd = _eqSrcs.end();
  for (srcs_type::iterator s_iter = _eqSrcs.begin(); s_iter != srcsEnd; ++s_iter) {
    EqKinSrc* src = s_iter->second;
    assert(0 != src);
    if (t >= src->originTime())
      src->slip(&slip, t);
  } // for

  const int spaceDim = _quadrature->spaceDim();

  // Get sections
  const ALE::Obj<SieveMesh>& sieveMesh = residual.mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveSubMesh>& faultSieveMesh = _faultMesh->sieveMesh();
  assert(!faultSieveMesh.isNull());
  const ALE::Obj<RealSection>& slipSection = slip.section();
  assert(!slipSection.isNull());
  const ALE::Obj<RealSection>& residualSection = residual.section();
  assert(!residualSection.isNull());

  const ALE::Obj<SieveMesh::label_sequence>& vertices =
      sieveMesh->depthStratum(0);
  assert(!vertices.isNull());
  const SieveSubMesh::label_sequence::iterator verticesBegin =
      vertices->begin();
  const SieveSubMesh::label_sequence::iterator verticesEnd = vertices->end();
  SieveSubMesh::renumbering_type& renumbering =
      faultSieveMesh->getRenumbering();
  const SieveSubMesh::renumbering_type::const_iterator renumberingEnd =
      renumbering.end();

  _logger->eventEnd(setupEvent);

  for (SieveSubMesh::label_sequence::iterator v_iter = verticesBegin; v_iter
      != verticesEnd; ++v_iter)
    if (renumbering.find(*v_iter) != renumberingEnd) {
      const int vertexFault = renumbering[*v_iter];
      const int vertexMesh = *v_iter;

      _logger->eventBegin(restrictEvent);
      const double* slipVertex = slipSection->restrictPoint(vertexFault);
      _logger->eventEnd(restrictEvent);

      assert(spaceDim == slipSection->getFiberDimension(vertexFault));
      assert(spaceDim == residualSection->getFiberDimension(vertexMesh));
      assert(0 != slipVertex);

      _logger->eventBegin(updateEvent);
      residualSection->updateAddPoint(vertexMesh, slipVertex);
      _logger->eventEnd(updateEvent);
    } // if
} // integrateResidualAssembled

// ----------------------------------------------------------------------
// Compute Jacobian matrix (A) associated with operator that do not
// require assembly across cells, vertices, or processors.
void pylith::faults::FaultCohesiveKin::integrateJacobianAssembled(topology::Jacobian* jacobian,
                                                                  const double t,
                                                                  topology::SolutionFields* const fields) { // integrateJacobianAssembled
  assert(0 != jacobian);
  assert(0 != fields);
  assert(0 != _fields);
  assert(0 != _logger);

  const int setupEvent = _logger->eventId("FkIJ setup");
  const int geometryEvent = _logger->eventId("FkIJ geometry");
  const int computeEvent = _logger->eventId("FkIJ compute");
  const int restrictEvent = _logger->eventId("FkIJ restrict");
  const int updateEvent = _logger->eventId("FkIJ update");

  _logger->eventBegin(setupEvent);

  // Add constraint information to Jacobian matrix; these are the
  // direction cosines. Entries are associated with vertices ik, jk,
  // ki, and kj.

  // Get cohesive cells
  const ALE::Obj<SieveMesh>& sieveMesh = fields->mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& cellsCohesive =
      sieveMesh->getLabelStratum("material-id", id());
  assert(!cellsCohesive.isNull());
  const SieveMesh::label_sequence::iterator cellsCohesiveBegin =
      cellsCohesive->begin();
  const SieveMesh::label_sequence::iterator cellsCohesiveEnd =
      cellsCohesive->end();
  const int cellsCohesiveSize = cellsCohesive->size();

  const int spaceDim = _quadrature->spaceDim();
  const int orientationSize = spaceDim * spaceDim;

  const int numConstraintVert = _quadrature->numBasis();
  const int numCorners = 3 * numConstraintVert; // cohesive cell
  double_array matrixCell(numCorners * spaceDim * numCorners * spaceDim);
  double_array orientationCell(numConstraintVert * orientationSize);

  // Get section information
  const ALE::Obj<SieveSubMesh>& faultSieveMesh = _faultMesh->sieveMesh();
  assert(!faultSieveMesh.isNull());
  const ALE::Obj<RealSection>& solutionSection = fields->solution().section();
  assert(!solutionSection.isNull());
  const ALE::Obj<RealSection>& orientationSection =
      _fields->get("orientation").section();
  assert(!orientationSection.isNull());
  topology::Mesh::RestrictVisitor orientationVisitor(*orientationSection,
    orientationCell.size(), &orientationCell[0]);

#if 0 // DEBUGGING
  // Check that fault cells match cohesive cells
  ALE::ISieveVisitor::PointRetriever<sieve_type> cV(std::max(1, mesh->getSieve()->getMaxConeSize()));
  ALE::ISieveVisitor::PointRetriever<sieve_type> cV2(std::max(1, _faultMesh->getSieve()->getMaxConeSize()));
  Mesh::renumbering_type& fRenumbering = _faultMesh->getRenumbering();
  const int rank = mesh->commRank();

  for (Mesh::label_sequence::iterator c_iter = cellsCohesiveBegin;
      c_iter != cellsCohesiveEnd;
      ++c_iter) {
    mesh->getSieve()->cone(*c_iter, cV);
    const int coneSize = cV.getSize();
    const Mesh::point_type *cone = cV.getPoints();
    const int faceSize = coneSize / 3;
    const Mesh::point_type face = _cohesiveToFault[*c_iter];
    _faultMesh->getSieve()->cone(face, cV2);
    const int fConeSize = cV2.getSize();
    const Mesh::point_type *fCone = cV2.getPoints();

    assert(0 == coneSize % faceSize);
    assert(faceSize == fConeSize);
    // Use last vertices (contraints) for fault mesh
    for(int i = 2*faceSize, j = 0; i < 3*faceSize; ++i, ++j) {
      assert(fRenumbering[cone[i]] == fCone[j]);
    }
    cV.clear();
    cV2.clear();
  }
#endif

  const PetscMat jacobianMatrix = jacobian->matrix();
  assert(0 != jacobianMatrix);
  const ALE::Obj<SieveMesh::order_type>& globalOrder =
      sieveMesh->getFactory()->getGlobalOrder(sieveMesh, "default",
        solutionSection);
  assert(!globalOrder.isNull());
  // We would need to request unique points here if we had an interpolated mesh
  topology::Mesh::IndicesVisitor jacobianVisitor(*solutionSection,
    *globalOrder, (int) pow(sieveMesh->getSieve()->getMaxConeSize(),
      sieveMesh->depth()) * spaceDim);

  _logger->eventEnd(setupEvent);

  for (SieveMesh::label_sequence::iterator c_iter = cellsCohesiveBegin; c_iter
      != cellsCohesiveEnd; ++c_iter) {
    const SieveMesh::point_type c_fault = _cohesiveToFault[*c_iter];

    _logger->eventBegin(restrictEvent);

    matrixCell = 0.0;
    // Get orientations at fault cell's vertices.
    orientationVisitor.clear();
    faultSieveMesh->restrictClosure(c_fault, orientationVisitor);

    _logger->eventEnd(restrictEvent);
    _logger->eventBegin(computeEvent);

    for (int iConstraint = 0; iConstraint < numConstraintVert; ++iConstraint) {
      // Blocks in cell matrix associated with normal cohesive
      // vertices i and j and constraint vertex k
      const int indexI = iConstraint;
      const int indexJ = iConstraint + numConstraintVert;
      const int indexK = iConstraint + 2 * numConstraintVert;

      // Get orientation at constraint vertex
      const double* orientationVertex = &orientationCell[iConstraint
          * orientationSize];
      assert(0 != orientationVertex);

      // Entries associated with constraint forces applied at node i
      for (int iDim = 0; iDim < spaceDim; ++iDim)
        for (int kDim = 0; kDim < spaceDim; ++kDim) {
          const int row = indexI * spaceDim + iDim;
          const int col = indexK * spaceDim + kDim;
          matrixCell[row * numCorners * spaceDim + col]
              = -orientationVertex[kDim * spaceDim + iDim];
          matrixCell[col * numCorners * spaceDim + row]
              = -orientationVertex[kDim * spaceDim + iDim];
        } // for

      // Entries associated with constraint forces applied at node j
      for (int jDim = 0; jDim < spaceDim; ++jDim)
        for (int kDim = 0; kDim < spaceDim; ++kDim) {
          const int row = indexJ * spaceDim + jDim;
          const int col = indexK * spaceDim + kDim;
          matrixCell[row * numCorners * spaceDim + col]
              = orientationVertex[kDim * spaceDim + jDim];
          matrixCell[col * numCorners * spaceDim + row]
              = orientationVertex[kDim * spaceDim + jDim];
        } // for
    } // for

    _logger->eventEnd(computeEvent);

    // Insert cell contribution into PETSc Matrix
    _logger->eventBegin(updateEvent);
    jacobianVisitor.clear();
    PetscErrorCode err = updateOperator(jacobianMatrix, *sieveMesh->getSieve(),
      jacobianVisitor, *c_iter, &matrixCell[0], INSERT_VALUES);
    CHECK_PETSC_ERROR_MSG(err, "Update to PETSc Mat failed.");
    _logger->eventEnd(updateEvent);

  } // for
  PetscLogFlops(cellsCohesiveSize*numConstraintVert*spaceDim*spaceDim*4);
  _needNewJacobian = false;
} // integrateJacobianAssembled

// ----------------------------------------------------------------------
// Compute Jacobian matrix (A) associated with operator that do not
// require assembly across cells, vertices, or processors.
void pylith::faults::FaultCohesiveKin::integrateJacobianAssembled(topology::Field<
                                                                      topology::Mesh>* jacobian,
                                                                  const double t,
                                                                  topology::SolutionFields* const fields) { // integrateJacobianAssembled
  assert(0 != jacobian);
  assert(0 != fields);
  assert(0 != _fields);
  assert(0 != _logger);

  const int setupEvent = _logger->eventId("FkIJ setup");
  const int geometryEvent = _logger->eventId("FkIJ geometry");
  const int computeEvent = _logger->eventId("FkIJ compute");
  const int restrictEvent = _logger->eventId("FkIJ restrict");
  const int updateEvent = _logger->eventId("FkIJ update");

  _logger->eventBegin(setupEvent);

  // Add ones to diagonal Jacobian matrix (as field) for
  // convenience. Instead of including the constraints in the Jacobian
  // matrix, we adjust the solution to account for the Lagrange
  // multipliers as part of the solve.

  // Get domain Sieve mesh
  const ALE::Obj<SieveMesh>& sieveMesh = fields->mesh().sieveMesh();
  assert(!sieveMesh.isNull());

  // Get fault Sieve mesh
  const ALE::Obj<SieveSubMesh>& faultSieveMesh = _faultMesh->sieveMesh();
  assert(!faultSieveMesh.isNull());

  const ALE::Obj<SieveMesh::label_sequence>& vertices =
      sieveMesh->depthStratum(0);
  assert(!vertices.isNull());
  const SieveSubMesh::label_sequence::iterator verticesBegin =
      vertices->begin();
  const SieveSubMesh::label_sequence::iterator verticesEnd = vertices->end();
  SieveSubMesh::renumbering_type& renumbering =
      faultSieveMesh->getRenumbering();
  const SieveSubMesh::renumbering_type::const_iterator renumberingEnd =
      renumbering.end();

  // Get section information
  const int spaceDim = _quadrature->spaceDim();
  double_array jacobianVertex(spaceDim);
  jacobianVertex = 1.0;
  const ALE::Obj<RealSection>& jacobianSection = jacobian->section();
  assert(!jacobianSection.isNull());

  _logger->eventEnd(setupEvent);

  for (SieveSubMesh::label_sequence::iterator v_iter = verticesBegin; v_iter
      != verticesEnd; ++v_iter)
    if (renumbering.find(*v_iter) != renumberingEnd) {
      assert(jacobianSection->getFiberDimension(*v_iter) == spaceDim);

      _logger->eventBegin(updateEvent);
      jacobianSection->updatePoint(*v_iter, &jacobianVertex[0]);
      _logger->eventEnd(updateEvent);
    } // if

  PetscLogFlops(0);

  _needNewJacobian = false;
} // integrateJacobianAssembled

// ----------------------------------------------------------------------
// Update state variables as needed.
void pylith::faults::FaultCohesiveKin::updateStateVars(const double t,
                                                       topology::SolutionFields* const fields) { // updateStateVars
  assert(0 != fields);
  assert(0 != _fields);

} // updateStateVars

// ----------------------------------------------------------------------
// Adjust solution from solver with lumped Jacobian to match Lagrange
// multiplier constraints.
void pylith::faults::FaultCohesiveKin::adjustSolnLumped(topology::SolutionFields* const fields,
                                                        const topology::Field<
                                                            topology::Mesh>& jacobian) { // adjustSolnLumped
  assert(0 != fields);
  assert(0 != _quadrature);

  // Cohesive cells with conventional vertices i and j, and constraint
  // vertex k require 2 adjustments to the solution:
  //
  //   * DOF k: Compute increment in Lagrange multipliers
  //            dl_k = S^{-1} (-C_ki (A_i^{-1} r_i - C_kj A_j^{-1} r_j + u_i - u_j) - d_k)
  //            S = C_ki (A_i^{-1} + A_j^{-1}) C_ki^T
  //   * DOF i and j: Adjust displacement increment (solution) to account
  //            for Lagrange multiplier constraints
  //            du_i = +A_i^-1 C_ki^T dlk
  //            du_j = -A_j^-1 C_kj^T dlk

  const int setupEvent = _logger->eventId("FkAS setup");
  const int geometryEvent = _logger->eventId("FkAS geometry");
  const int computeEvent = _logger->eventId("FkAS compute");
  const int restrictEvent = _logger->eventId("FkAS restrict");
  const int updateEvent = _logger->eventId("FkAS update");

  _logger->eventBegin(setupEvent);

  // Get cell information and setup storage for cell data
  const int spaceDim = _quadrature->spaceDim();
  const int orientationSize = spaceDim * spaceDim;
  const int numBasis = _quadrature->numBasis();
  const int numConstraintVert = numBasis;
  const int numCorners = 3 * numConstraintVert; // cohesive cell
  const int numQuadPts = _quadrature->numQuadPts();
  const double_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == numQuadPts);

  // Get cohesive cells
  const ALE::Obj<SieveMesh>& sieveMesh =
      fields->get("dispIncr(t->t+dt)").mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& cellsCohesive =
      sieveMesh->getLabelStratum("material-id", id());
  assert(!cellsCohesive.isNull());
  const SieveMesh::label_sequence::iterator cellsCohesiveBegin =
      cellsCohesive->begin();
  const SieveMesh::label_sequence::iterator cellsCohesiveEnd =
      cellsCohesive->end();
  const int cellsCohesiveSize = cellsCohesive->size();

  // Get fault Sieve mesh
  const ALE::Obj<SieveSubMesh>& faultSieveMesh = _faultMesh->sieveMesh();
  assert(!faultSieveMesh.isNull());

  // Get section information
  double_array orientationCell(numConstraintVert * orientationSize);
  const ALE::Obj<RealSection>& orientationSection =
      _fields->get("orientation").section();
  assert(!orientationSection.isNull());
  topology::Mesh::RestrictVisitor orientationVisitor(*orientationSection,
    orientationCell.size(), &orientationCell[0]);

  double_array slipCell(numConstraintVert * spaceDim);
  const ALE::Obj<RealSection>& slipSection = _fields->get("slip").section();
  assert(!slipSection.isNull());
  topology::Mesh::RestrictVisitor slipVisitor(*slipSection, slipCell.size(),
    &slipCell[0]);

  // Tributary area for the current for each vertex.
  double_array areaVertexCell(numConstraintVert);
  // Total fault area associated with each vertex (assembled over all cells).
  double_array areaCell(numConstraintVert);
  const ALE::Obj<RealSection>& areaSection = _fields->get("area").section();
  assert(!areaSection.isNull());
  topology::Mesh::RestrictVisitor areaVisitor(*areaSection, areaCell.size(),
    &areaCell[0]);

  double_array jacobianCell(numCorners * spaceDim);
  const ALE::Obj<RealSection>& jacobianSection = jacobian.section();
  assert(!jacobianSection.isNull());
  topology::Mesh::RestrictVisitor jacobianVisitor(*jacobianSection,
    jacobianCell.size(), &jacobianCell[0]);

  double_array residualCell(numCorners * spaceDim);
  const ALE::Obj<RealSection>& residualSection =
      fields->get("residual").section();
  topology::Mesh::RestrictVisitor residualVisitor(*residualSection,
    residualCell.size(), &residualCell[0]);

  double_array dispTCell(numCorners * spaceDim);
  const ALE::Obj<RealSection>& dispTSection = fields->get("disp(t)").section();
  topology::Mesh::RestrictVisitor dispTVisitor(*dispTSection, dispTCell.size(),
    &dispTCell[0]);

  double_array solutionCell(numCorners * spaceDim);
  const ALE::Obj<RealSection>& solutionSection = fields->get(
    "dispIncr(t->t+dt)").section();
  topology::Mesh::UpdateAddVisitor solutionVisitor(*solutionSection,
    &solutionCell[0]);

  double_array coordinatesCell(numBasis * spaceDim);
  const ALE::Obj<RealSection>& coordinates = faultSieveMesh->getRealSection(
    "coordinates");
  assert(!coordinates.isNull());
  topology::Mesh::RestrictVisitor coordsVisitor(*coordinates,
    coordinatesCell.size(), &coordinatesCell[0]);

  _logger->eventEnd(setupEvent);

  for (SieveMesh::label_sequence::iterator c_iter = cellsCohesiveBegin; c_iter
      != cellsCohesiveEnd; ++c_iter) {
    const SieveMesh::point_type c_fault = _cohesiveToFault[*c_iter];
    areaVertexCell = 0.0;
    solutionCell = 0.0;

    // Compute geometry information for current cell
    _logger->eventBegin(geometryEvent);
#if defined(PRECOMPUTE_GEOMETRY)
    _quadrature->retrieveGeometry(c_fault);
#else
    coordsVisitor.clear();
    faultSieveMesh->restrictClosure(c_fault, coordsVisitor);
    _quadrature->computeGeometry(coordinatesCell, c_fault);
#endif
    _logger->eventEnd(geometryEvent);

    // Get cell geometry information that depends on cell
    const double_array& basis = _quadrature->basis();
    const double_array& jacobianDet = _quadrature->jacobianDet();

    _logger->eventBegin(restrictEvent);

    // Get orientations at fault cell's vertices.
    orientationVisitor.clear();
    faultSieveMesh->restrictClosure(c_fault, orientationVisitor);

    // Get area at fault cell's vertices.
    areaVisitor.clear();
    faultSieveMesh->restrictClosure(c_fault, areaVisitor);

    // Get slip at fault cell's vertices.
    slipVisitor.clear();
    faultSieveMesh->restrictClosure(c_fault, slipVisitor);

    // Get residual at cohesive cell's vertices.
    residualVisitor.clear();
    sieveMesh->restrictClosure(*c_iter, residualVisitor);

    // Get jacobian at cohesive cell's vertices.
    jacobianVisitor.clear();
    sieveMesh->restrictClosure(*c_iter, jacobianVisitor);

    // Get disp(t) at cohesive cell's vertices.
    dispTVisitor.clear();
    sieveMesh->restrictClosure(*c_iter, dispTVisitor);

    _logger->eventEnd(restrictEvent);
    _logger->eventBegin(computeEvent);

    // Compute contributory area for cell (to weight contributions)
    for (int iQuad = 0; iQuad < numQuadPts; ++iQuad) {
      const double wt = quadWts[iQuad] * jacobianDet[iQuad];
      for (int iBasis = 0; iBasis < numBasis; ++iBasis) {
        const double dArea = wt * basis[iQuad * numBasis + iBasis];
        areaVertexCell[iBasis] += dArea;
      } // for
    } // for

    for (int iConstraint = 0; iConstraint < numConstraintVert; ++iConstraint) {
      // Blocks in cell matrix associated with normal cohesive
      // vertices i and j and constraint vertex k
      const int indexI = iConstraint;
      const int indexJ = iConstraint + numConstraintVert;
      const int indexK = iConstraint + 2 * numConstraintVert;

      assert(areaCell[iConstraint] > 0);
      const double wt = areaVertexCell[iConstraint] / areaCell[iConstraint];

      // Get orientation at constraint vertex
      const double* orientationVertex = &orientationCell[iConstraint
          * orientationSize];
      assert(0 != orientationVertex);

      // Get slip at constraint vertex
      const double* slipVertex = &slipCell[iConstraint * spaceDim];
      assert(0 != slipVertex);

      // Get Jacobian at conventional vertices i and j
      const double* Ai = &jacobianCell[indexI * spaceDim];
      assert(0 != Ai);
      const double* Aj = &jacobianCell[indexJ * spaceDim];
      assert(0 != Aj);

      // Get residual at conventional vertices i and j
      const double* ri = &residualCell[indexI * spaceDim];
      assert(0 != ri);
      const double* rj = &residualCell[indexJ * spaceDim];
      assert(0 != rj);

      // Get disp(t) at conventional vertices i and j
      const double* ui = &dispTCell[indexI * spaceDim];
      assert(0 != ui);
      const double* uj = &dispTCell[indexJ * spaceDim];
      assert(0 != uj);

      switch (spaceDim) { // switch
      case 1: {
        assert(Ai[0] > 0.0);
        assert(Aj[0] > 0.0);

        const double Sinv = Ai[0] * Aj[0] / (Ai[0] + Aj[0]);

        // Aru = A_i^{-1} r_i - A_j^{-1} r_j + u_i - u_j
        const double Aru = ri[0] / Ai[0] - rj[0] / Aj[0] + ui[0] - uj[0];

        // dl_k = D^{-1}( - C_{ki} Aru - d_k)
        const double Aruslip = -Aru - slipVertex[0];
        const double dlp = wt * Sinv * Aruslip;

        // Update displacements at node I
        solutionCell[indexI * spaceDim + 0] = +1.0 / Ai[0] * dlp;

        // Update displacements at node J
        solutionCell[indexJ * spaceDim + 0] = -1.0 / Aj[0] * dlp;

        // Update Lagrange multipliers
        solutionCell[indexK * spaceDim + 0] = dlp;

        PetscLogFlops(17);
        break;
      } // case 1
      case 2: {
        assert(Ai[0] > 0.0);
        assert(Ai[1] > 0.0);
        assert(Aj[0] > 0.0);
        assert(Aj[1] > 0.0);

        const double Cpx = orientationVertex[0];
        const double Cpy = orientationVertex[1];
        const double Cqx = orientationVertex[2];
        const double Cqy = orientationVertex[3];

        // Check to make sure Jacobian is same at all DOF for
        // vertices i and j (means S is diagonal with equal enties).
        assert(Ai[0] == Ai[1]);
        assert(Aj[0] == Aj[1]);

        const double Sinv = Ai[0] * Aj[0] / (Ai[0] + Aj[0]);

        // Aru = A_i^{-1} r_i - A_j^{-1} r_j + u_i - u_j
        const double Arux = ri[0] / Ai[0] - rj[0] / Aj[0] + ui[0] - uj[0];
        const double Aruy = ri[1] / Ai[1] - rj[1] / Aj[1] + ui[1] - uj[1];

        // dl_k = S^{-1}(-C_{ki} Aru - d_k)
        const double Arup = Cpx * Arux + Cpy * Aruy;
        const double Aruq = Cqx * Arux + Cqy * Aruy;
        const double Arupslip = -Arup - slipVertex[0];
        const double Aruqslip = -Aruq - slipVertex[1];
        const double dlp = Sinv * Arupslip;
        const double dlq = Sinv * Aruqslip;

        const double dlx = wt * (Cpx * dlp + Cqx * dlq);
        const double dly = wt * (Cpy * dlp + Cqy * dlq);

        // Update displacements at node I
        solutionCell[indexI * spaceDim + 0] = dlx / Ai[0];
        solutionCell[indexI * spaceDim + 1] = dly / Ai[1];

        // Update displacements at node J
        solutionCell[indexJ * spaceDim + 0] = -dlx / Aj[0];
        solutionCell[indexJ * spaceDim + 1] = -dly / Aj[0];

        // Update Lagrange multipliers
        solutionCell[indexK * spaceDim + 0] = wt * dlp;
        solutionCell[indexK * spaceDim + 1] = wt * dlq;

        PetscLogFlops(41);
        break;
      } // case 2
      case 3: {
        assert(Ai[0] > 0.0);
        assert(Ai[1] > 0.0);
        assert(Ai[2] > 0.0);
        assert(Aj[0] > 0.0);
        assert(Aj[1] > 0.0);
        assert(Aj[2] > 0.0);

        const double Cpx = orientationVertex[0];
        const double Cpy = orientationVertex[1];
        const double Cpz = orientationVertex[2];
        const double Cqx = orientationVertex[3];
        const double Cqy = orientationVertex[4];
        const double Cqz = orientationVertex[5];
        const double Crx = orientationVertex[6];
        const double Cry = orientationVertex[7];
        const double Crz = orientationVertex[8];

        // Check to make sure Jacobian is same at all DOF for
        // vertices i and j (means S is diagonal with equal enties).
        assert(Ai[0] == Ai[1] && Ai[0] == Ai[2]);
        assert(Aj[0] == Aj[1] && Aj[0] == Aj[2]);

        const double Sinv = Ai[0] * Aj[0] / (Ai[0] + Aj[0]);

        // Aru = A_i^{-1} r_i - A_j^{-1} r_j + u_i - u_j
        const double Arux = ri[0] / Ai[0] - rj[0] / Aj[0] + ui[0] - uj[0];
        const double Aruy = ri[1] / Ai[1] - rj[1] / Aj[1] + ui[1] - uj[1];
        const double Aruz = ri[2] / Ai[2] - rj[2] / Aj[2] + ui[2] - uj[2];

        // dl_k = D^{-1}( -C_{ki} Aru - d_k)
        const double Arup = Cpx * Arux + Cpy * Aruy + Cpz * Aruz;
        const double Aruq = Cqx * Arux + Cqy * Aruy + Cqz * Aruz;
        const double Arur = Crx * Arux + Cry * Aruy + Crz * Aruz;
        const double Arupslip = -Arup - slipVertex[0];
        const double Aruqslip = -Aruq - slipVertex[1];
        const double Arurslip = -Arur - slipVertex[2];
        const double dlp = Sinv * Arupslip;
        const double dlq = Sinv * Aruqslip;
        const double dlr = Sinv * Arurslip;

        const double dlx = wt * (Cpx * dlp + Cqx * dlq + Crx * dlr);
        const double dly = wt * (Cpy * dlp + Cqy * dlq + Cry * dlr);
        const double dlz = wt * (Cpz * dlp + Cqz * dlq + Crz * dlr);

        // Update displacements at node I
        solutionCell[indexI * spaceDim + 0] = dlx / Ai[0];
        solutionCell[indexI * spaceDim + 1] = dly / Ai[1];
        solutionCell[indexI * spaceDim + 2] = dlz / Ai[2];

        // Update displacements at node J
        solutionCell[indexJ * spaceDim + 0] = -dlx / Aj[0];
        solutionCell[indexJ * spaceDim + 1] = -dly / Aj[1];
        solutionCell[indexJ * spaceDim + 2] = -dlz / Aj[2];

        // Update Lagrange multipliers
        solutionCell[indexK * spaceDim + 0] = wt * dlp;
        solutionCell[indexK * spaceDim + 1] = wt * dlq;
        solutionCell[indexK * spaceDim + 2] = wt * dlr;

        PetscLogFlops(72);
        break;
      } // case 3
      default:
        assert(0);
        throw std::logic_error("Unknown spatial dimension.");
      } // switch
    } // for

    _logger->eventEnd(computeEvent);

#if 0 // DEBUGGING
    std::cout << "Adjusting lumped solution for cell " << *c_iter << std::endl;
    for(int i = 0; i < numCorners*spaceDim; ++i) {
      std::cout << "  dispTpdt["<<i<<"]: " << dispTpdtCell[i] << std::endl;
    }
    for(int i = 0; i < numCorners*spaceDim; ++i) {
      std::cout << "  dispT["<<i<<"]: " << dispTCell[i] << std::endl;
    }
    for(int i = 0; i < numCorners*spaceDim; ++i) {
      std::cout << "  dispIncr["<<i<<"]: " << dispTIncrCell[i] << std::endl;
    }
    for(int i = 0; i < numCorners*spaceDim; ++i) {
      std::cout << "  v["<<i<<"]: " << residualCell[i] << std::endl;
    }
#endif

    _logger->eventBegin(updateEvent);
    solutionVisitor.clear();
    sieveMesh->updateClosure(*c_iter, solutionVisitor);
    _logger->eventEnd(updateEvent);
  } // for

  // FIX THIS
  PetscLogFlops(0);
} // adjustSolnLumped

// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void pylith::faults::FaultCohesiveKin::verifyConfiguration(const topology::Mesh& mesh) const { // verifyConfiguration
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
// Get vertex field associated with integrator.
const pylith::topology::Field<pylith::topology::SubMesh>&
pylith::faults::FaultCohesiveKin::vertexField(const char* name,
                                              const topology::SolutionFields* fields) { // vertexField
  assert(0 != _faultMesh);
  assert(0 != _quadrature);
  assert(0 != _normalizer);
  assert(0 != _fields);

  const int cohesiveDim = _faultMesh->dimension();
  const int spaceDim = _quadrature->spaceDim();

  const int slipStrLen = strlen("final_slip");
  const int timeStrLen = strlen("slip_time");

  double scale = 0.0;
  int fiberDim = 0;
  if (0 == strcasecmp("slip", name)) {
    const topology::Field<topology::SubMesh>& slip = _fields->get("slip");
    return slip;

  } else if (cohesiveDim > 0 && 0 == strcasecmp("strike_dir", name)) {
    const ALE::Obj<RealSection>& orientationSection = _fields->get(
      "orientation").section();
    assert(!orientationSection.isNull());
    const ALE::Obj<RealSection>& dirSection = orientationSection->getFibration(
      0);
    assert(!dirSection.isNull());
    _allocateBufferVectorField();
    topology::Field<topology::SubMesh>& buffer =
        _fields->get("buffer (vector)");
    buffer.copy(dirSection);
    buffer.label("strike_dir");
    buffer.scale(1.0);
    return buffer;

  } else if (2 == cohesiveDim && 0 == strcasecmp("dip_dir", name)) {
    const ALE::Obj<RealSection>& orientationSection = _fields->get(
      "orientation").section();
    assert(!orientationSection.isNull());
    const ALE::Obj<RealSection>& dirSection = orientationSection->getFibration(
      1);
    _allocateBufferVectorField();
    topology::Field<topology::SubMesh>& buffer =
        _fields->get("buffer (vector)");
    buffer.copy(dirSection);
    buffer.label("dip_dir");
    buffer.scale(1.0);
    return buffer;

  } else if (0 == strcasecmp("normal_dir", name)) {
    const ALE::Obj<RealSection>& orientationSection = _fields->get(
      "orientation").section();
    assert(!orientationSection.isNull());
    const int space = (0 == cohesiveDim) ? 0 : (1 == cohesiveDim) ? 1 : 2;
    const ALE::Obj<RealSection>& dirSection = orientationSection->getFibration(
      space);
    assert(!dirSection.isNull());
    _allocateBufferVectorField();
    topology::Field<topology::SubMesh>& buffer =
        _fields->get("buffer (vector)");
    buffer.copy(dirSection);
    buffer.label("normal_dir");
    buffer.scale(1.0);
    return buffer;

  } else if (0 == strncasecmp("final_slip_X", name, slipStrLen)) {
    const std::string value = std::string(name).substr(slipStrLen + 1);

    const srcs_type::const_iterator s_iter = _eqSrcs.find(value);
    assert(s_iter != _eqSrcs.end());
    return s_iter->second->finalSlip();

  } else if (0 == strncasecmp("slip_time_X", name, timeStrLen)) {
    const std::string value = std::string(name).substr(timeStrLen + 1);
    const srcs_type::const_iterator s_iter = _eqSrcs.find(value);
    assert(s_iter != _eqSrcs.end());
    return s_iter->second->slipTime();

  } else if (0 == strcasecmp("traction_change", name)) {
    assert(0 != fields);
    const topology::Field<topology::Mesh>& dispT = fields->get("disp(t)");
    _allocateBufferVectorField();
    topology::Field<topology::SubMesh>& buffer =
        _fields->get("buffer (vector)");
    _calcTractionsChange(&buffer, dispT);
    return buffer;

  } else {
    std::ostringstream msg;
    msg << "Request for unknown vertex field '" << name << "' for fault '"
        << label() << "'.";
    throw std::runtime_error(msg.str());
  } // else


  // Satisfy return values
  assert(0 != _fields);
  const topology::Field<topology::SubMesh>& buffer = _fields->get(
    "buffer (vector)");
  return buffer;
} // vertexField

// ----------------------------------------------------------------------
// Get cell field associated with integrator.
const pylith::topology::Field<pylith::topology::SubMesh>&
pylith::faults::FaultCohesiveKin::cellField(const char* name,
                                            const topology::SolutionFields* fields) { // cellField
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

// ----------------------------------------------------------------------
// Initialize auxiliary cohesive cell information.
void pylith::faults::FaultCohesiveKin::_initializeCohesiveInfo(const topology::Mesh& mesh) { // _initializeCohesiveInfo
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
  SubMesh::renumbering_type& renumbering = faultSieveMesh->getRenumbering();
  const SieveSubMesh::renumbering_type::const_iterator renumberingEnd =
    renumbering.end();

  typedef std::map<int, int> indexmap_type;
  indexmap_type indexMap;
  _cohesiveVertices.resize(renumbering.size());
  int index = 0;

  const ALE::Obj<SieveMesh::sieve_type>& sieve = mesh.sieveMesh()->getSieve();
  assert(!sieve.isNull());
  typedef ALE::SieveAlg<SieveMesh> SieveAlg;

  ALE::ISieveVisitor::NConeRetriever<SieveMesh::sieve_type> ncV(*sieve,
      (size_t) pow(sieve->getMaxConeSize(), std::max(0, sieveMesh->depth())));

  for (SieveMesh::label_sequence::iterator c_iter=cellsBegin;
      c_iter != cellsEnd;
      ++c_iter) {
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
	std::cout << "cohesiveVertices[" << index << "]: "
		  << "l: " << v_lagrange
		  << ", p: " << v_positive
		  << ", n: " << v_negative
		  << ", f: " << v_fault
		  << std::endl;
        indexMap[v_lagrange] = index; // add index to map
        ++index;
      } // if
    } // for
  } // for
} // _initializeCohesiveInfo

// ----------------------------------------------------------------------
// Initialize logger.
void pylith::faults::FaultCohesiveKin::_initializeLogger(void) { // initializeLogger
  delete _logger;
  _logger = new utils::EventLogger;
  assert(0 != _logger);
  _logger->className("FaultCohesiveKin");
  _logger->initialize();

  _logger->registerEvent("FkAS setup");
  _logger->registerEvent("FkAS geometry");
  _logger->registerEvent("FkAS compute");
  _logger->registerEvent("FkAS restrict");
  _logger->registerEvent("FkAS update");

  _logger->registerEvent("FkIR setup");
  _logger->registerEvent("FkIR geometry");
  _logger->registerEvent("FkIR compute");
  _logger->registerEvent("FkIR restrict");
  _logger->registerEvent("FkIR update");

  _logger->registerEvent("FkIJ setup");
  _logger->registerEvent("FkIJ geometry");
  _logger->registerEvent("FkIJ compute");
  _logger->registerEvent("FkIJ restrict");
  _logger->registerEvent("FkIJ update");
} // initializeLogger

// ----------------------------------------------------------------------
// Calculate orientation at fault vertices.
void pylith::faults::FaultCohesiveKin::_calcOrientation(const double upDir[3],
                                                        const double normalDir[3]) { // _calcOrientation
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
void pylith::faults::FaultCohesiveKin::_calcArea(void) { // _calcArea
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
// NOTE: We must convert vertex labels to fault vertex labels
void pylith::faults::FaultCohesiveKin::_calcTractionsChange(topology::Field<
    topology::SubMesh>* tractions, const topology::Field<topology::Mesh>& dispT) { // _calcTractionsChange
  assert(0 != tractions);
  assert(0 != _faultMesh);
  assert(0 != _fields);
  assert(0 != _normalizer);

  tractions->label("traction_change");
  tractions->scale(_normalizer->pressureScale());

  // Get vertices from mesh of domain.
  const ALE::Obj<SieveMesh>& sieveMesh = dispT.mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& vertices =
      sieveMesh->depthStratum(0);
  assert(!vertices.isNull());
  const SieveMesh::label_sequence::iterator verticesBegin = vertices->begin();
  const SieveMesh::label_sequence::iterator verticesEnd = vertices->end();

  // Get fault vertices
  const ALE::Obj<SieveMesh>& faultSieveMesh = _faultMesh->sieveMesh();
  assert(!faultSieveMesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& fvertices =
      faultSieveMesh->depthStratum(0);
  const SieveSubMesh::label_sequence::iterator fverticesBegin =
      fvertices->begin();
  const SieveSubMesh::label_sequence::iterator fverticesEnd = fvertices->end();

  // Get sections.
  const ALE::Obj<RealSection>& areaSection = _fields->get("area").section();
  assert(!areaSection.isNull());
  const ALE::Obj<RealSection>& dispTSection = dispT.section();
  assert(!dispTSection.isNull());

  const int numFaultVertices = fvertices->size();
  Mesh::renumbering_type& renumbering = faultSieveMesh->getRenumbering();
  const SieveSubMesh::renumbering_type::const_iterator renumberingEnd =
      renumbering.end();

#if 0 // DEBUGGING, MOVE TO SEPARATE CHECK METHOD
  // Check fault mesh and volume mesh coordinates
  const ALE::Obj<RealSection>& coordinates = mesh->getRealSection("coordinates");
  const ALE::Obj<RealSection>& fCoordinates = _faultMesh->getRealSection("coordinates");

  for (Mesh::label_sequence::iterator v_iter = verticesBegin;
      v_iter != verticesEnd;
      ++v_iter) {
    if (renumbering.find(*v_iter) != renumberingEnd) {
      const int v = *v_iter;
      const int dim = coordinates->getFiberDimension(*v_iter);
      const double *a = coordinates->restrictPoint(*v_iter);
      const int fv = renumbering[*v_iter];
      const int fDim = fCoordinates->getFiberDimension(fv);
      const double *fa = fCoordinates->restrictPoint(fv);

      if (dim != fDim) throw ALE::Exception("Coordinate fiber dimensions do not match");
      for(int d = 0; d < dim; ++d) {
        if (a[d] != fa[d]) throw ALE::Exception("Coordinate values do not match");
      }
    }
  }
#endif

  // Fiber dimension of tractions matches spatial dimension.
  const int fiberDim = _quadrature->spaceDim();
  double_array tractionsVertex(fiberDim);

  // Allocate buffer for tractions field (if nec.).
  const ALE::Obj<RealSection>& tractionsSection = tractions->section();
  if (tractionsSection.isNull()) {
    ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
    //logger.stagePush("Fault");

    const topology::Field<topology::SubMesh>& slip = _fields->get("slip");
    tractions->newSection(slip, fiberDim);
    tractions->allocate();

    //logger.stagePop();
  } // if
  assert(!tractionsSection.isNull());
  tractions->zero();

  const double pressureScale = tractions->scale();
  for (SieveMesh::label_sequence::iterator v_iter = verticesBegin; v_iter
      != verticesEnd; ++v_iter)
    if (renumbering.find(*v_iter) != renumberingEnd) {
      const int vertexMesh = *v_iter;
      const int vertexFault = renumbering[*v_iter];
      assert(fiberDim == dispTSection->getFiberDimension(vertexMesh));
      assert(fiberDim == tractionsSection->getFiberDimension(vertexFault));
      assert(1 == areaSection->getFiberDimension(vertexFault));

      const double* dispTVertex = dispTSection->restrictPoint(vertexMesh);
      assert(0 != dispTVertex);
      const double* areaVertex = areaSection->restrictPoint(vertexFault);
      assert(0 != areaVertex);

      for (int i = 0; i < fiberDim; ++i)
        tractionsVertex[i] = dispTVertex[i] / areaVertex[0];

      tractionsSection->updatePoint(vertexFault, &tractionsVertex[0]);
    } // if

  PetscLogFlops(numFaultVertices * (1 + fiberDim) );

#if 0 // DEBUGGING
  tractions->view("TRACTIONS");
#endif
} // _calcTractionsChange

// ----------------------------------------------------------------------
// Allocate buffer for vector field.
void pylith::faults::FaultCohesiveKin::_allocateBufferVectorField(void) { // _allocateBufferVectorField
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
void pylith::faults::FaultCohesiveKin::_allocateBufferScalarField(void) { // _allocateBufferScalarField
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


// End of file 
