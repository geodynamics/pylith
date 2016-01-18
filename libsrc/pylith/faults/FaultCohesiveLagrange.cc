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

#include "FaultCohesiveLagrange.hh" // implementation of object methods

#include "CohesiveTopology.hh" // USES CohesiveTopology

#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/feassemble/CellGeometry.hh" // USES CellGeometry
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/Jacobian.hh" // USES Jacobian
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields
#include "pylith/topology/VisitorMesh.hh" // USES VisitorMesh
#include "pylith/topology/VisitorSubMesh.hh" // USES SubMeshIS
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor

#include "pylith/utils/EventLogger.hh" // USES EventLogger
#include "pylith/utils/macrodefs.h" // USES CALL_MEMBER_FN

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <cmath> // USES pow(), sqrt()
#include <strings.h> // USES strcasecmp()
#include <cstring> // USES strlen()
#include <cstdlib> // USES atoi()
#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

#include <iostream>

//#define DETAILED_EVENT_LOGGING

// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::FaultCohesiveLagrange::FaultCohesiveLagrange(void) :
  _cohesiveIS(0)
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
  PYLITH_METHOD_BEGIN;
  
  FaultCohesive::deallocate();
  delete _cohesiveIS; _cohesiveIS = 0;

  PYLITH_METHOD_END;
} // deallocate

// ----------------------------------------------------------------------
// Initialize fault. Determine orientation and setup boundary
void
pylith::faults::FaultCohesiveLagrange::initialize(const topology::Mesh& mesh,
						  const PylithScalar upDir[3])
{ // initialize
  PYLITH_METHOD_BEGIN;

  assert(upDir);
  assert(_quadrature);
  assert(_normalizer);

  _initializeLogger();

  const spatialdata::geocoords::CoordSys* cs = mesh.coordsys();assert(cs);

  const bool isSubMesh = true;
  delete _faultMesh; _faultMesh = new topology::Mesh(isSubMesh);assert(_faultMesh);
  CohesiveTopology::createFaultParallel(_faultMesh, mesh, id(), label(), _useLagrangeConstraints); // :TODO: Obsolete?
  
  topology::MeshOps::checkTopology(*_faultMesh);

  // Optimize coordinate retrieval in closure
  PetscDM faultDMMesh = _faultMesh->dmMesh();assert(faultDMMesh);
  topology::CoordsVisitor::optimizeClosure(faultDMMesh);

  _initializeCohesiveInfo(mesh);

  delete _fields; _fields = new topology::Fields(*_faultMesh);assert(_fields);

  // Allocate dispRel field
  const int spaceDim = cs->spaceDim();
  _fields->add("relative disp", "relative_disp");
  topology::Field& dispRel = _fields->get("relative disp");
  dispRel.newSection(topology::FieldBase::VERTICES_FIELD, spaceDim); // :TODO: Update?

  topology::SubMeshIS faultMeshIS(*_faultMesh);

  topology::Stratum verticesStratum(faultDMMesh, topology::Stratum::DEPTH, 0);

  const int numVertices = _cohesiveVertices.size();
  PetscErrorCode err;

  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int e_lagrange = _cohesiveVertices[iVertex].lagrange;
    const int v_fault    = _cohesiveVertices[iVertex].fault;

    if (e_lagrange < 0) {err = PetscSectionSetConstraintDof(dispRel.localSection(), -v_fault, spaceDim);PYLITH_CHECK_ERROR(err);}
  }
  dispRel.allocate();
  {
    PetscInt *ind = (spaceDim > 0) ? new PetscInt[spaceDim] : 0;

    for (int i = 0; i < spaceDim; ++i) ind[i] = i;
    for (int iVertex = 0; iVertex < numVertices; ++iVertex) {
      const int e_lagrange = _cohesiveVertices[iVertex].lagrange;
      const int v_fault    = _cohesiveVertices[iVertex].fault;

      if (e_lagrange < 0) {err = PetscSectionSetConstraintIndices(dispRel.localSection(), -v_fault, ind);PYLITH_CHECK_ERROR(err);}
    }
    delete[] ind; ind = 0;
  }
  dispRel.vectorFieldType(topology::FieldBase::VECTOR);
  dispRel.scale(_normalizer->lengthScale());

  _quadrature->initializeGeometry();

  // Compute orientation at vertices in fault mesh.
  _calcOrientation(upDir);

  // Compute tributary area for each vertex in fault mesh.
  _calcArea();

  PYLITH_METHOD_END;
} // initialize

// ----------------------------------------------------------------------
void
pylith::faults::FaultCohesiveLagrange::setupSolnDof(topology::Field* field)
{ // setupSolnDof
  PYLITH_METHOD_BEGIN;

  assert(field);

  const int indexDisp = field->subfieldInfo("displacement").index;
  const int indexLagrange = field->subfieldInfo("lagrange_multiplier").index;

  PetscDM dmMesh = field->dmMesh();assert(dmMesh);
  PetscSection fieldSection  = field->localSection();assert(fieldSection);
  PetscErrorCode err;

  assert(_quadrature);
  const int spaceDim = _quadrature->spaceDim();

  const int numVertices = _cohesiveVertices.size();
  for(PetscInt iVertex = 0; iVertex < numVertices; ++iVertex) {
    const int v_positive = _cohesiveVertices[iVertex].positive;
    const int v_negative = _cohesiveVertices[iVertex].negative;
    const int e_lagrange = _cohesiveVertices[iVertex].lagrange;

    if (e_lagrange < 0) continue;
    // Set DOF in section (all subfields)
    err = PetscSectionSetDof(fieldSection, e_lagrange, spaceDim);PYLITH_CHECK_ERROR(err);
    err = PetscSectionSetDof(fieldSection, v_positive, spaceDim);PYLITH_CHECK_ERROR(err);
    err = PetscSectionSetDof(fieldSection, v_negative, spaceDim);PYLITH_CHECK_ERROR(err);

    // Set DOF in displacement subfield
    err = PetscSectionSetFieldDof(fieldSection, v_positive, indexDisp, spaceDim);PYLITH_CHECK_ERROR(err);
    err = PetscSectionSetFieldDof(fieldSection, v_negative, indexDisp, spaceDim);PYLITH_CHECK_ERROR(err);

    // Set DOF in Lagrange multiplier subfield
    err = PetscSectionSetFieldDof(fieldSection, e_lagrange, indexLagrange, spaceDim);PYLITH_CHECK_ERROR(err);
  } // for

  PYLITH_METHOD_END;
} // setupSolnDof



// ----------------------------------------------------------------------
// Integrate contribution of cohesive cells to residual term.
void
pylith::faults::FaultCohesiveLagrange::integrateResidual(const topology::Field& residual,
							 const PylithScalar t,
							 topology::SolutionFields* const fields)
{ // integrateResidual
  PYLITH_METHOD_BEGIN;

  assert(fields);
  assert(_fields);
  assert(_logger);

  // Cohesive cells with conventional vertices N and P, and constraint
  // vertex L make contributions to the assembled residual:
  //
  // DOF P: \int_{S_f^+} \tensor{N}_m^T \cdot \tensor{N}_p \cdot \vec{l}_p dS
  // DOF N: -\int_{S_f^+} \tensor{N}_m^T \cdot \tensor{N}_p \cdot \vec{l}_p dS
  // DOF L: \int_S_f \tensor{N}_p^T ( \tensor{R} \cdot \vec{d} 
  //                 -\tensor{N}_{n^+} \cdot \vec{u}_{n^+}
  //                 +\tensor{N}_{n^-} \cdot \vec{u}_{n^-} dS

  const int setupEvent = _logger->eventId("FaIR setup");
  const int computeEvent = _logger->eventId("FaIR compute");
#if defined(DETAILED_EVENT_LOGGING)
  const int geometryEvent = _logger->eventId("FaIR geometry");
  const int restrictEvent = _logger->eventId("FaIR restrict");
  const int updateEvent = _logger->eventId("FaIR update");
#endif

  _logger->eventBegin(setupEvent);

  // Get cell geometry information that doesn't depend on cell
  const int spaceDim = _quadrature->spaceDim();

  // Get sections associated with cohesive cells
  PetscSection residualSection = residual.localSection();assert(residualSection);
  PetscSection residualGlobalSection = residual.globalSection();assert(residualGlobalSection);

  topology::VecVisitorMesh residualVisitor(residual);
  PetscScalar* residualArray = residualVisitor.localArray();

  topology::Field& dispT = fields->get("disp(t)");
  topology::VecVisitorMesh dispTVisitor(dispT);
  PetscScalar* dispTArray = dispTVisitor.localArray();

  topology::Field& dispTIncr = fields->get("dispIncr(t->t+dt)");
  topology::VecVisitorMesh dispTIncrVisitor(dispTIncr);
  PetscScalar* dispTIncrArray = dispTIncrVisitor.localArray();

  topology::Field& dispRel = _fields->get("relative disp");
  topology::VecVisitorMesh dispRelVisitor(dispRel);
  PetscScalar* dispRelArray = dispRelVisitor.localArray();

  topology::Field& area = _fields->get("area");
  topology::VecVisitorMesh areaVisitor(area);
  PetscScalar* areaArray = areaVisitor.localArray();

  // Get fault information
  PetscDM dmMesh = fields->mesh().dmMesh();assert(dmMesh);

  _logger->eventEnd(setupEvent);
#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventBegin(computeEvent);
#endif

  // Loop over fault vertices
  const int numVertices = _cohesiveVertices.size();
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int e_lagrange = _cohesiveVertices[iVertex].lagrange;
    const int v_fault = _cohesiveVertices[iVertex].fault;
    const int v_negative = _cohesiveVertices[iVertex].negative;
    const int v_positive = _cohesiveVertices[iVertex].positive;

    if (e_lagrange < 0) { // Skip clamped edges.
      continue;
    } // if

    // Compute contribution only if Lagrange constraint is local.
    PetscInt goff  = 0;
    PetscErrorCode err = PetscSectionGetOffset(residualGlobalSection, e_lagrange, &goff);PYLITH_CHECK_ERROR(err);
    if (goff < 0)
      continue;

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(restrictEvent);
#endif

    // Get relative dislplacement at fault vertex.
    const PetscInt droff = dispRelVisitor.sectionOffset(v_fault);
    assert(spaceDim == dispRelVisitor.sectionDof(v_fault));

    // Get area associated with fault vertex.
    const PetscInt aoff = areaVisitor.sectionOffset(v_fault);
    assert(1 == areaVisitor.sectionDof(v_fault));
    const PylithScalar areaValue = areaArray[aoff];

    // Get disp(t) at conventional vertices and Lagrange vertex.
    const PetscInt dtnoff = dispTVisitor.sectionOffset(v_negative);
    assert(spaceDim == dispTVisitor.sectionDof(v_negative));

    const PetscInt dtpoff = dispTVisitor.sectionOffset(v_positive);
    assert(spaceDim == dispTVisitor.sectionDof(v_positive));

    const PetscInt dtloff = dispTVisitor.sectionOffset(e_lagrange);
    assert(spaceDim == dispTVisitor.sectionDof(e_lagrange));

    // Get dispIncr(t->t+dt) at conventional vertices and Lagrange vertex.
    const PetscInt dinoff = dispTIncrVisitor.sectionOffset(v_negative);
    assert(spaceDim == dispTIncrVisitor.sectionDof(v_negative));

    const PetscInt dipoff = dispTIncrVisitor.sectionOffset(v_positive);
    assert(spaceDim == dispTIncrVisitor.sectionDof(v_positive));

    const PetscInt diloff = dispTIncrVisitor.sectionOffset(e_lagrange);
    assert(spaceDim == dispTIncrVisitor.sectionDof(e_lagrange));

    // Assemble contributions into field
    const PetscInt rnoff = residualVisitor.sectionOffset(v_negative);
    assert(spaceDim == residualVisitor.sectionDof(v_negative));

    const PetscInt rpoff = residualVisitor.sectionOffset(v_positive);
    assert(spaceDim == residualVisitor.sectionDof(v_positive));

    const PetscInt rloff = residualVisitor.sectionOffset(e_lagrange);
    assert(spaceDim == residualVisitor.sectionDof(e_lagrange));

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(restrictEvent);
    _logger->eventBegin(updateEvent);
#endif

    for(PetscInt d = 0; d < spaceDim; ++d) {
      const PylithScalar residualN = areaValue * (dispTArray[dtloff+d] + dispTIncrArray[diloff+d]);
      residualArray[rnoff+d] += +residualN;
      residualArray[rpoff+d] += -residualN;
      residualArray[rloff+d] += -areaValue * (dispTArray[dtpoff+d] + dispTIncrArray[dipoff+d] - dispTArray[dtnoff+d] - dispTIncrArray[dinoff+d] - dispRelArray[droff+d]);
    } // for

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(updateEvent);
#endif
  } // for
  PetscLogFlops(numVertices*spaceDim*10);

#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventEnd(computeEvent);
#endif

  PYLITH_METHOD_END;
} // integrateResidual

// ----------------------------------------------------------------------
// Compute Jacobian matrix (A) associated with operator.
void
pylith::faults::FaultCohesiveLagrange::integrateJacobian(topology::Jacobian* jacobian,
							 const PylithScalar t,
							 topology::SolutionFields* const fields)
{ // integrateJacobian
  PYLITH_METHOD_BEGIN;

  assert(jacobian);
  assert(fields);
  assert(_fields);
  assert(_logger);

  const int setupEvent = _logger->eventId("FaIJ setup");
  const int computeEvent = _logger->eventId("FaIJ compute");
#if defined(DETAILED_EVENT_LOGGING)
  const int geometryEvent = _logger->eventId("FaIJ geometry");
  const int restrictEvent = _logger->eventId("FaIJ restrict");
  const int updateEvent = _logger->eventId("FaIJ update");
#endif

  _logger->eventBegin(setupEvent);

  // Add constraint information to Jacobian matrix; Entries are
  // associated with vertices ik, jk, ki, and kj.

  // Get cell geometry information that doesn't depend on cell
  const int spaceDim = _quadrature->spaceDim();

  // Get fields.
  topology::Field& area = _fields->get("area");
  topology::VecVisitorMesh areaVisitor(area);
  const PetscScalar* areaArray = areaVisitor.localArray();
  
  PetscSection solnSection = fields->solution().localSection();assert(solnSection);
  PetscSection solnGlobalSection = fields->solution().globalSection();assert(solnGlobalSection);

  // Get fault information
  PetscDM dmMesh = fields->mesh().dmMesh();assert(dmMesh);

  // Allocate vectors for vertex values
  scalar_array jacobianVertex(spaceDim*spaceDim);
  int_array indicesL(spaceDim);
  int_array indicesN(spaceDim);
  int_array indicesP(spaceDim);
  int_array indicesRel(spaceDim);
  for (int i=0; i < spaceDim; ++i) {
    indicesRel[i] = i;
  } // for

  // Get sparse matrix
  const PetscMat jacobianMatrix = jacobian->matrix();assert(jacobianMatrix);

  _logger->eventEnd(setupEvent);
#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventBegin(computeEvent);
#endif

  PetscErrorCode err = 0;
  const int numVertices = _cohesiveVertices.size();
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int e_lagrange = _cohesiveVertices[iVertex].lagrange;
    const int v_fault = _cohesiveVertices[iVertex].fault;
    const int v_negative = _cohesiveVertices[iVertex].negative;
    const int v_positive = _cohesiveVertices[iVertex].positive;

    if (e_lagrange < 0) { // Skip clamped edges.
      continue;
    } // if

    // Compute contribution only if Lagrange constraint is local.
    PetscInt gloff = 0;
    err = PetscSectionGetOffset(solnGlobalSection, e_lagrange, &gloff);PYLITH_CHECK_ERROR(err);
    if (gloff < 0)
      continue;

    PetscInt gnoff = 0;
    err = PetscSectionGetOffset(solnGlobalSection, v_negative, &gnoff);PYLITH_CHECK_ERROR(err);
    gnoff = gnoff < 0 ? -(gnoff+1) : gnoff;

    PetscInt gpoff = 0;
    err = PetscSectionGetOffset(solnGlobalSection, v_positive, &gpoff);PYLITH_CHECK_ERROR(err);
    gpoff = gpoff < 0 ? -(gpoff+1) : gpoff;

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(restrictEvent);
#endif

    // Get area associated with fault vertex.
    const PetscInt aoff = areaVisitor.sectionOffset(v_fault);
    assert(1 == areaVisitor.sectionDof(v_fault));

    // Set global order indices
    indicesL = indicesRel + gloff;
    indicesN = indicesRel + gnoff;
    indicesP = indicesRel + gpoff;
    PetscInt cdof;
    err = PetscSectionGetConstraintDof(solnSection, v_negative, &cdof);PYLITH_CHECK_ERROR(err);assert(0 == cdof);
    err = PetscSectionGetConstraintDof(solnSection, v_positive, &cdof);PYLITH_CHECK_ERROR(err);assert(0 == cdof);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(restrictEvent);
    _logger->eventBegin(updateEvent);
#endif

    // Set diagonal entries of Jacobian at positive vertex to area
    // associated with vertex.
    for (int iDim=0; iDim < spaceDim; ++iDim)
      jacobianVertex[iDim*spaceDim+iDim] = areaArray[aoff];

    // Values at positive vertex, entry L,P in Jacobian
    err = MatSetValues(jacobianMatrix, 
		       indicesL.size(), &indicesL[0], 
		       indicesP.size(), &indicesP[0], 
		       &jacobianVertex[0], ADD_VALUES);PYLITH_CHECK_ERROR(err);

    // Values at positive vertex, entry P,L in Jacobian
    err = MatSetValues(jacobianMatrix, 
		       indicesP.size(), &indicesP[0], 
		       indicesL.size(), &indicesL[0], 
		       &jacobianVertex[0], ADD_VALUES);PYLITH_CHECK_ERROR(err);
    
    // Values at negative vertex, entry L,N in Jacobian
    jacobianVertex *= -1.0;
    err = MatSetValues(jacobianMatrix,
		       indicesL.size(), &indicesL[0],
		       indicesN.size(), &indicesN[0],
		       &jacobianVertex[0], ADD_VALUES);PYLITH_CHECK_ERROR(err);

    // Values at negative vertex, entry N,L in Jacobian
    err = MatSetValues(jacobianMatrix,
		       indicesN.size(), &indicesN[0],
		       indicesL.size(), &indicesL[0],
		       &jacobianVertex[0], ADD_VALUES);PYLITH_CHECK_ERROR(err);

    // Values at Lagrange vertex, entry L,L in Jacobian
    // We must have entries on the diagonal.
    jacobianVertex = 0.0;
    err = MatSetValues(jacobianMatrix,
		       indicesL.size(), &indicesL[0],
		       indicesL.size(), &indicesL[0],
		       &jacobianVertex[0], ADD_VALUES);PYLITH_CHECK_ERROR(err);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(updateEvent);
#endif

  } // for
  PetscLogFlops(numVertices*spaceDim*2);

#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventEnd(computeEvent);
#endif

  _needNewJacobian = false;

  PYLITH_METHOD_END;
} // integrateJacobian

// ----------------------------------------------------------------------
// Compute Jacobian matrix (A) associated with operator.
void
pylith::faults::FaultCohesiveLagrange::integrateJacobian(topology::Field* jacobian,
							 const PylithScalar t,
							 topology::SolutionFields* const fields)
{ // integrateJacobian
  PYLITH_METHOD_BEGIN;

  assert(jacobian);
  assert(fields);
  assert(_fields);
  assert(_logger);

  const int setupEvent = _logger->eventId("FaIJ setup");
  const int computeEvent = _logger->eventId("FaIJ compute");
#if defined(DETAILED_EVENT_LOGGING)
  const int geometryEvent = _logger->eventId("FaIJ geometry");
  const int restrictEvent = _logger->eventId("FaIJ restrict");
  const int updateEvent = _logger->eventId("FaIJ update");
#endif

  _logger->eventBegin(setupEvent);

  // Add ones to diagonal Jacobian matrix (as field) for
  // convenience. Instead of including the constraints in the Jacobian
  // matrix, we adjust the solution to account for the Lagrange
  // multipliers as part of the solve.

  const int spaceDim  = _quadrature->spaceDim();

  PetscSection jacobianGlobalSection = jacobian->globalSection();assert(jacobianGlobalSection);

  topology::VecVisitorMesh jacobianVisitor(*jacobian);
  PetscScalar* jacobianArray = jacobianVisitor.localArray();

  _logger->eventEnd(setupEvent);
#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventBegin(computeEvent);
#endif

  PetscErrorCode err = 0;
  const int numVertices = _cohesiveVertices.size();
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int e_lagrange = _cohesiveVertices[iVertex].lagrange;

    if (e_lagrange < 0) { // Skip clamped edges
      continue;
    } // if

    // Compute contribution only if Lagrange constraint is local.
    PetscInt goff = 0;
    err = PetscSectionGetOffset(jacobianGlobalSection, e_lagrange, &goff);PYLITH_CHECK_ERROR(err);
    if (goff < 0) 
      continue;

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(updateEvent);
#endif
    const PetscInt off = jacobianVisitor.sectionOffset(e_lagrange);
    assert(spaceDim == jacobianVisitor.sectionDof(e_lagrange));

    for(PetscInt d = 0; d < spaceDim; ++d) {
      jacobianArray[off+d] = 1.0;
    } // for

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(updateEvent);
#endif
  } // for
  PetscLogFlops(0);

#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventEnd(computeEvent);
#endif

  _needNewJacobian = false;

  PYLITH_METHOD_END;
} // integrateJacobian

// ----------------------------------------------------------------------
// Integrate contributions to Jacobian matrix (A) associated with
// operator.
void
pylith::faults::FaultCohesiveLagrange::calcPreconditioner(PetscMat* const precondMatrix,
							  topology::Jacobian* const jacobian,
							  topology::SolutionFields* const fields)
{ // calcPreconditioner
  PYLITH_METHOD_BEGIN;

  assert(precondMatrix);
  assert(jacobian);
  assert(fields);
  assert(_fields);
  assert(_logger);

  /** We have A = [K L^T]
   *              [L   0]
   *
   * Compute Pmat = -( [L] [K]^(-1) [L]^T ) using the diagonal of K to
   * approximate [K]^(-1).
   *
   * Decompose [K] into [Kn] and [Kp], where [Kn] contains the terms
   * for vertices on the negative side of the fault and [Kp] contains
   * the terms for vertices on the positive side of the fault.
   *
   * Pmat = L_{ik} (1.0/Kn_{kk} + 1.0/Kp{kk}) L_{ik}
   *
   * Because we use quadrature points located at the vertices,
   * L_{ii} = area, L_{ij} = 0 if i != j
   */

  const int setupEvent = _logger->eventId("FaPr setup");
  const int computeEvent = _logger->eventId("FaPr compute");
#if defined(DETAILED_EVENT_LOGGING)
  const int restrictEvent = _logger->eventId("FaPr restrict");
  const int updateEvent = _logger->eventId("FaPr update");
#endif

  _logger->eventBegin(setupEvent);

  // Get cell information and setup storage for cell data
  const int spaceDim = _quadrature->spaceDim();

  // Allocate vectors for vertex values
  scalar_array jacobianVertexP(spaceDim*spaceDim);
  scalar_array jacobianVertexN(spaceDim*spaceDim);
  scalar_array jacobianInvVertexP(spaceDim);
  scalar_array jacobianInvVertexN(spaceDim);
  scalar_array precondVertexL(spaceDim);
  int_array indicesN(spaceDim);
  int_array indicesP(spaceDim);
  int_array indicesRel(spaceDim);
  for (int i=0; i < spaceDim; ++i) {
    indicesRel[i] = i;
  } // for

  // Get fields
  topology::Field& area = _fields->get("area");
  topology::VecVisitorMesh areaVisitor(area);
  const PetscScalar* areaArray = areaVisitor.localArray();

  PetscSection solnGlobalSection = fields->solution().globalSection();assert(solnGlobalSection);

  PetscDM lagrangeDM = fields->solution().subfieldInfo("lagrange_multiplier").dm;assert(lagrangeDM);
  PetscSection lagrangeGlobalSection = NULL;
  PetscErrorCode err = DMGetDefaultGlobalSection(lagrangeDM, &lagrangeGlobalSection);PYLITH_CHECK_ERROR(err);

  _logger->eventEnd(setupEvent);
#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventBegin(computeEvent);
#endif

  PetscMat jacobianNP;
  std::map<int, int> indicesMatToSubmat;
  _getJacobianSubmatrixNP(&jacobianNP, &indicesMatToSubmat, *jacobian, *fields);
  
  const int numVertices = _cohesiveVertices.size();
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int e_lagrange = _cohesiveVertices[iVertex].lagrange;
    const int v_fault = _cohesiveVertices[iVertex].fault;
    const int v_negative = _cohesiveVertices[iVertex].negative;
    const int v_positive = _cohesiveVertices[iVertex].positive;

    if (e_lagrange < 0) { // Skip clamped edges.
      continue;
    } // if

    // Compute contribution only if Lagrange constraint is local.
    PetscInt gloff = 0;
    err = PetscSectionGetOffset(solnGlobalSection, e_lagrange, &gloff);PYLITH_CHECK_ERROR(err);
    if (gloff < 0) {
      continue;
    } // if

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(restrictEvent);
#endif

    // Get area associated with fault vertex.
    const PetscInt aoff = areaVisitor.sectionOffset(v_fault);
    assert(1 == areaVisitor.sectionDof(v_fault));

    PetscInt gnoff = 0;
    err = PetscSectionGetOffset(solnGlobalSection, v_negative, &gnoff);PYLITH_CHECK_ERROR(err);
    indicesN = indicesRel + indicesMatToSubmat[gnoff];
    err = MatGetValues(jacobianNP,
                       indicesN.size(), &indicesN[0], indicesN.size(), &indicesN[0],
                       &jacobianVertexN[0]);PYLITH_CHECK_ERROR(err);

    PetscInt gpoff = 0;
    err = PetscSectionGetOffset(solnGlobalSection, v_positive, &gpoff);PYLITH_CHECK_ERROR(err);
    indicesP = indicesRel + indicesMatToSubmat[gpoff];
    err = MatGetValues(jacobianNP,
                       indicesP.size(), &indicesP[0], indicesP.size(), &indicesP[0],
                       &jacobianVertexP[0]);PYLITH_CHECK_ERROR(err);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(restrictEvent);
    _logger->eventBegin(computeEvent);
#endif

    // Compute inverse of Jacobian diagonals
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      jacobianInvVertexN[iDim] = 1.0/jacobianVertexN[iDim*spaceDim+iDim];
      jacobianInvVertexP[iDim] = 1.0/jacobianVertexP[iDim*spaceDim+iDim];
    } // for

    // Compute -[L] [Adiag]^(-1) [L]^T
    //   L_{ii} = L^T{ii} = areaVertex
    //   Adiag^{-1}_{ii} = jacobianInvVertexN[i] + jacobianInvVertexP[i]
    precondVertexL = 0.0;
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      precondVertexL[iDim] -= areaArray[aoff] * areaArray[aoff] * 
        (jacobianInvVertexN[iDim] + jacobianInvVertexP[iDim]);
    } // for
    

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(computeEvent);
    _logger->eventBegin(updateEvent);
#endif

    // Set diagonal entries in preconditioned matrix.
    PetscInt poff = 0;
    err = PetscSectionGetOffset(lagrangeGlobalSection, e_lagrange, &poff);PYLITH_CHECK_ERROR(err);

    for (int iDim=0; iDim < spaceDim; ++iDim) {
      err = MatSetValue(*precondMatrix, poff+iDim, poff+iDim, precondVertexL[iDim], INSERT_VALUES);PYLITH_CHECK_ERROR(err);
    } // for

#if 0 // DEBUGGING
    std::cout << "1/P_vertex " << e_lagrange << ", poff: " << poff << std::endl;
    for(int iDim = 0; iDim < spaceDim; ++iDim) {
      std::cout << "  " << precondVertexL[iDim] << std::endl;
    } // for
#endif

    
#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(updateEvent);
#endif
  } // for
  err = MatDestroy(&jacobianNP);PYLITH_CHECK_ERROR(err);
  PetscLogFlops(numVertices*spaceDim*6);

#if !defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(computeEvent);
#endif

  PYLITH_METHOD_END;
} // calcPreconditioner

// ----------------------------------------------------------------------
// Adjust solution from solver with lumped Jacobian to match Lagrange
// multiplier constraints.
void
pylith::faults::FaultCohesiveLagrange::adjustSolnLumped(topology::SolutionFields* const fields,
							const PylithScalar t,
                                                        const topology::Field& jacobian)
{ // adjustSolnLumped
  PYLITH_METHOD_BEGIN;

  assert(fields);
  assert(_quadrature);

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
  const int computeEvent = _logger->eventId("FaAS compute");
#if defined(DETAILED_EVENT_LOGGING)
  const int geometryEvent = _logger->eventId("FaAS geometry");
  const int restrictEvent = _logger->eventId("FaAS restrict");
#endif

  _logger->eventBegin(setupEvent);

  // Get cell information and setup storage for cell data
  const int spaceDim = _quadrature->spaceDim();

  // Get fields
  topology::Field& area = _fields->get("area");
  topology::VecVisitorMesh areaVisitor(area);
  const PetscScalar* areaArray = areaVisitor.localArray();

  topology::VecVisitorMesh jacobianVisitor(jacobian);
  const PetscScalar* jacobianArray = jacobianVisitor.localArray();

  topology::Field& residual = fields->get("residual");
  topology::VecVisitorMesh residualVisitor(residual);
  const PetscScalar* residualArray = residualVisitor.localArray();

  topology::Field& dispTIncr = fields->get("dispIncr(t->t+dt)");
  topology::VecVisitorMesh dispTIncrVisitor(dispTIncr);
  PetscScalar* dispTIncrArray = dispTIncrVisitor.localArray();

  topology::Field& dispTIncrAdj = fields->get("dispIncr adjust");
  topology::VecVisitorMesh dispTIncrAdjVisitor(dispTIncrAdj);
  PetscScalar* dispTIncrAdjArray = dispTIncrAdjVisitor.localArray();

  PetscSection jacobianGlobalSection = jacobian.globalSection();assert(jacobianGlobalSection);

  _logger->eventEnd(setupEvent);

#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventBegin(computeEvent);
#endif

  PetscErrorCode err = 0;
  const int numVertices = _cohesiveVertices.size();
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int e_lagrange = _cohesiveVertices[iVertex].lagrange;
    const int v_fault = _cohesiveVertices[iVertex].fault;
    const int v_negative = _cohesiveVertices[iVertex].negative;
    const int v_positive = _cohesiveVertices[iVertex].positive;

    if (e_lagrange < 0) { // Skip clamped edges
      continue;
    } // if

    // Set Lagrange multiplier value. Value from preliminary solve is
    // bogus due to artificial diagonal entry.
    const PetscInt dtloff = dispTIncrVisitor.sectionOffset(e_lagrange);
    assert(spaceDim == dispTIncrVisitor.sectionDof(e_lagrange));
    for(PetscInt d = 0; d < spaceDim; ++d) {
      dispTIncrArray[dtloff+d] = 0.0;
    } // for

    // Compute contribution only if Lagrange constraint is local.
    PetscInt goff;
    err = PetscSectionGetOffset(jacobianGlobalSection, e_lagrange, &goff);PYLITH_CHECK_ERROR(err);
    if (goff < 0) {
      continue;
    } // if

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(restrictEvent);
#endif

    // Residual at Lagrange vertex.
    const PetscInt rloff = residualVisitor.sectionOffset(e_lagrange);
    assert(spaceDim == residualVisitor.sectionDof(e_lagrange));

    // Get jacobian at cohesive cell's vertices.
    const PetscInt jnoff = jacobianVisitor.sectionOffset(v_negative);
    assert(spaceDim == jacobianVisitor.sectionDof(v_negative));

    const PetscInt jpoff = jacobianVisitor.sectionOffset(v_positive);
    assert(spaceDim == jacobianVisitor.sectionDof(v_positive));

    // Area at fault vertex.
    const PetscInt aoff = areaVisitor.sectionOffset(v_fault);
    assert(1 == areaVisitor.sectionDof(v_fault));assert(areaArray[aoff] > 0.0);

    // Get dispIncr(t) at cohesive cell vertices.
    const PetscInt dtnoff = dispTIncrVisitor.sectionOffset(v_negative);
    assert(spaceDim == dispTIncrVisitor.sectionDof(v_negative));

    const PetscInt dtpoff = dispTIncrVisitor.sectionOffset(v_positive);
    assert(spaceDim == dispTIncrVisitor.sectionDof(v_positive));

    // Get dispIncrAdj at cohesive cell vertices.
    const PetscInt danoff = dispTIncrAdjVisitor.sectionOffset(v_negative);
    assert(spaceDim == dispTIncrAdjVisitor.sectionDof(v_negative));

    const PetscInt dapoff = dispTIncrAdjVisitor.sectionOffset(v_positive);
    assert(spaceDim == dispTIncrAdjVisitor.sectionDof(v_positive));

    const PetscInt daloff = dispTIncrAdjVisitor.sectionOffset(e_lagrange);
    assert(spaceDim == dispTIncrAdjVisitor.sectionDof(e_lagrange));
    
#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(restrictEvent);
    _logger->eventBegin(computeEvent);
#endif

    const PetscScalar areaVertex = areaArray[aoff];
    for(PetscInt d = 0; d < spaceDim; ++d) {
      const PylithScalar S = (1.0/jacobianArray[jpoff+d] + 1.0/jacobianArray[jnoff+d]) * areaVertex * areaVertex;
      // Set Lagrange multiplier value (value from preliminary solve is bogus due to artificial diagonal entry)
      dispTIncrAdjArray[daloff+d] = 1.0/S * (-residualArray[rloff+d] + areaArray[aoff] * (dispTIncrArray[dtpoff+d] - dispTIncrArray[dtnoff+d]));

      // Adjust displacements to account for Lagrange multiplier values (assumed to be zero in preliminary solve).
      assert(jacobianArray[jnoff+d] > 0.0);
      dispTIncrAdjArray[danoff+d] +=  +areaVertex / jacobianArray[jnoff+d] * dispTIncrAdjArray[dtloff+d];

      assert(jacobianArray[jpoff+d] > 0.0);
      dispTIncrAdjArray[dapoff+d] += -areaVertex / jacobianArray[jpoff+d] * dispTIncrAdjArray[dtloff+d];
    } // for

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(computeEvent);
#endif

  } // for

#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventEnd(computeEvent);
#endif

  PYLITH_METHOD_END;
} // adjustSolnLumped

// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::faults::FaultCohesiveLagrange::verifyConfiguration(const topology::Mesh& mesh) const
{ // verifyConfiguration
  PYLITH_METHOD_BEGIN;

  assert(_quadrature);
  PetscErrorCode err = 0;

  PetscDM dmMesh = mesh.dmMesh();assert(dmMesh);
  topology::Stratum edgeStratum(dmMesh, topology::Stratum::DEPTH, 1);
  const PetscInt eEnd = edgeStratum.end();
  PetscInt eMax;

  err = DMPlexGetHybridBounds(dmMesh, NULL, NULL, &eMax, NULL);PYLITH_CHECK_ERROR(err);

  // Check for fault groups
  PetscBool hasLabel;
  err = DMHasLabel(dmMesh, label(), &hasLabel);PYLITH_CHECK_ERROR(err);
  if (!hasLabel) {
    std::ostringstream msg;
    msg << "Mesh missing group of vertices '" << label() << " defining fault.";
    throw std::runtime_error(msg.str());
  } // if  

  if (strlen(edge()) > 0) {
    PetscBool hasLabel;
    PetscErrorCode err = DMHasLabel(dmMesh, edge(), &hasLabel);PYLITH_CHECK_ERROR(err);
    if (!hasLabel) {
      std::ostringstream msg;
      msg << "Mesh missing group of vertices '" << edge() << "' defining buried edges of fault.";
      throw std::runtime_error(msg.str());
    } // if  
  } // if

  // Verify quadrature scheme is consistent with points collocated
  // with verties. Expect basis functions to be 1.0 at one quadrature
  // point and zero at all others.
  const scalar_array basis = _quadrature->basis();
  const int numBasis = _quadrature->numBasis();
  const int numQuadPts = _quadrature->numQuadPts();
  for (int iQuadPt=0; iQuadPt < numQuadPts; ++iQuadPt) {
    int nonzero = 0;
    const PylithScalar tolerance = 1.0e-6;
    for (int iBasis=0; iBasis < numBasis; ++iBasis) {
      if (fabs(basis[iQuadPt*numBasis+iBasis]) > tolerance) 
        ++nonzero;
    } // for
    if (numBasis != numQuadPts || 1 != nonzero) {
      std::ostringstream msg;
      msg << "Quadrature scheme for fault " << label()
	  << " is incompatible with fault implementation.\n"
	  << "Expected quadrature points collocated with vertices, so that "
	  << "basis functions are \n"
	  << "1.0 at one quadrature point and zero at all other quadrature "
	  << "points.";
      throw std::runtime_error(msg.str());
    } // if
  } // for

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
  const bool includeOnlyCells = true;
  topology::StratumIS cohesiveIS(dmMesh, "material-id", id(), includeOnlyCells);
  const PetscInt* cells = cohesiveIS.points();
  const PetscInt ncells = cohesiveIS.size();

  if (ncells > 0 && eMax < 0) {
    std::ostringstream msg;
    msg << "No hybrid edges found in mesh with cohesive cells for fault '" << label() << "'.";
    throw std::logic_error(msg.str());
  } // if  

  for(PetscInt i = 0; i < ncells; ++i) {
    PetscInt *closure = NULL;
    PetscInt cellNumEdges = 0, closureSize;

    err = DMPlexGetTransitiveClosure(dmMesh, cells[i], PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
    for(PetscInt p = 0; p < closureSize*2; p += 2) {
      const PetscInt point = closure[p];
      if ((point >= eMax) && (point < eEnd)) {
        ++cellNumEdges;
      }
    }
    if (numBasis != cellNumEdges) {
      std::ostringstream msg;
      msg << "Quadrature is incompatible with cell for fault '" << label() << "'. "
	  << "Cell " << cells[i] << " has " << cellNumEdges
	  << " edges but quadrature reference cell has "
	  << numBasis << " edges.";
      throw std::runtime_error(msg.str());
    } // if
    err = DMPlexRestoreTransitiveClosure(dmMesh, cells[i], PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
  } // for

  PYLITH_METHOD_END;
} // verifyConfiguration

// ----------------------------------------------------------------------
// Verify constraints are acceptable.
void
pylith::faults::FaultCohesiveLagrange::checkConstraints(const topology::Field& solution) const
{ // checkConstraints
  PYLITH_METHOD_BEGIN;

  // Check to make sure no vertices connected to the fault are
  // constrained.

  const PetscInt spaceDim = solution.mesh().dimension();
  topology::VecVisitorMesh solutionVisitor(solution);

  const int numVertices = _cohesiveVertices.size();
  PetscInt numConstraints = 0;
  for (int iVertex = 0; iVertex < numVertices; ++iVertex) {
    const int v_negative = _cohesiveVertices[iVertex].negative;
    const int e_lagrange = _cohesiveVertices[iVertex].lagrange;

    if (e_lagrange < 0) { // Skip clamped edges.
      continue;
    } // if

    assert(spaceDim == solutionVisitor.sectionDof(v_negative));
    numConstraints = solutionVisitor.sectionConstraintDof(v_negative);
    if (numConstraints > 0) {
      std::ostringstream msg;
      msg << "Vertex with label '" << v_negative << "' on negative side "
	  << "of fault '" << label() << "' is constrained.\n"
	  << "Fault vertices cannot be constrained.";
      throw std::runtime_error(msg.str());
    } // if

    const int v_positive = _cohesiveVertices[iVertex].positive;
    assert(spaceDim == solutionVisitor.sectionDof(v_positive));
    numConstraints = solutionVisitor.sectionConstraintDof(v_positive);
    if (numConstraints > 0) {
      std::ostringstream msg;
      msg << "Vertex with label '" << v_positive << "' on positive side "
	  << "of fault '" << label() << "' is constrained.\n"
	  << "Fault vertices cannot be constrained.";
      throw std::runtime_error(msg.str());
    } // if
  } // for

  PYLITH_METHOD_END;
} // checkConstraints

// ----------------------------------------------------------------------
// Initialize auxiliary cohesive cell information.
void pylith::faults::FaultCohesiveLagrange::_initializeCohesiveInfo(const topology::Mesh& mesh)
{ // _initializeCohesiveInfo
  PYLITH_METHOD_BEGIN;

  assert(_quadrature);

  const int numBasis = _quadrature->numBasis();
  PetscErrorCode err = 0;

  // Get cohesive cells
  PetscDM dmMesh = mesh.dmMesh();assert(dmMesh);
  const bool includeOnlyCells = true;
  delete _cohesiveIS; _cohesiveIS = new topology::StratumIS(dmMesh, "material-id", id(), includeOnlyCells);assert(_cohesiveIS);
  const PetscInt* cohesiveCells = _cohesiveIS->points();
  const int numCohesiveCells = _cohesiveIS->size();

  // Get domain edges (Lagrange multiplier DOF)
  topology::Stratum edgeStratum(dmMesh, topology::Stratum::DEPTH, 1);
  const PetscInt eEnd = edgeStratum.end();
  PetscInt eMax;
  err = DMPlexGetHybridBounds(dmMesh, NULL, NULL, &eMax, NULL);PYLITH_CHECK_ERROR(err);

  // Get vertices and cells in fault mesh.
  PetscDM faultDMMesh = _faultMesh->dmMesh();assert(faultDMMesh);
  topology::Stratum faultVerticesStratum(faultDMMesh, topology::Stratum::DEPTH, 0);
  const PetscInt fvStart = faultVerticesStratum.begin();
  const PetscInt fvEnd = faultVerticesStratum.end();
  topology::Stratum faultCellsStratum(faultDMMesh, topology::Stratum::HEIGHT, 1);
  const PetscInt fcStart = faultCellsStratum.begin();
  const PetscInt fcEnd = faultCellsStratum.end();
  assert(numCohesiveCells == fcEnd-fcStart);

  topology::SubMeshIS faultMeshIS(*_faultMesh);
  const PetscInt numPoints = faultMeshIS.size();
  const PetscInt* points = faultMeshIS.points();

  _cohesiveToFault.clear();
  typedef std::map<int, int> indexmap_type;
  indexmap_type indexMap;
  _cohesiveVertices.resize(fvEnd-fvStart);

  PetscInt index = 0;
  for (PetscInt iCell = 0; iCell < numCohesiveCells; ++iCell) {
    _cohesiveToFault[cohesiveCells[iCell]] = iCell+fcStart;

    // Get oriented closure
    PetscInt *closure = NULL;
    PetscInt  closureSize, q = 0;
    err = DMPlexGetTransitiveClosure(dmMesh, cohesiveCells[iCell], PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
    for (PetscInt p = 0; p < closureSize*2; p += 2) {
      const PetscInt point = closure[p];
      if ((point >= eMax) && (point < eEnd)) {
        closure[q++] = point;
      } // if
    } // for
    closureSize = q;
    assert(closureSize == numBasis);

    for (int iConstraint = 0; iConstraint < numBasis; ++iConstraint) {
      const PetscInt  e_lagrange = closure[iConstraint];

      const PetscInt *cone = NULL;
      PetscInt coneSize;
      err = DMPlexGetConeSize(dmMesh, e_lagrange, &coneSize);PYLITH_CHECK_ERROR(err);
      assert(coneSize == 2);
      err = DMPlexGetCone(dmMesh, e_lagrange, &cone);PYLITH_CHECK_ERROR(err);
      // Check for clamped vertex
      if (cone[0] == cone[1]) {
        PetscInt v_fault;
        err = PetscFindInt(cone[0], numPoints, points, &v_fault);PYLITH_CHECK_ERROR(err);
        assert(v_fault >= 0);
        if (indexMap.end() == indexMap.find(e_lagrange)) {
          _cohesiveVertices[index].lagrange = -e_lagrange;
          _cohesiveVertices[index].positive = -cone[0];
          _cohesiveVertices[index].negative = -cone[0];
          _cohesiveVertices[index].fault    = -v_fault;
          indexMap[e_lagrange] = index; // add index to map
          ++index;
        } // if
        err = DMSetLabelValue(faultDMMesh, "clamped", v_fault, 1);PYLITH_CHECK_ERROR(err);
        continue;
      }
      const PetscInt v_negative = cone[0];
      const PetscInt v_positive = cone[1];

      PetscInt v_fault;
      err = PetscFindInt(v_negative, numPoints, points, &v_fault);PYLITH_CHECK_ERROR(err);
      assert(v_fault >= 0);
      if (indexMap.end() == indexMap.find(e_lagrange)) {
        _cohesiveVertices[index].lagrange = e_lagrange;
        _cohesiveVertices[index].positive = v_positive;
        _cohesiveVertices[index].negative = v_negative;
        _cohesiveVertices[index].fault = v_fault;
#if 0 // DEBUGGING
        std::cout << "cohesiveVertices[" << index << "]: "
		  << "l: " << e_lagrange
		  << ", p: " << v_positive
		  << ", n: " << v_negative
		  << ", f: " << v_fault
		  << std::endl;
#endif
        indexMap[e_lagrange] = index; // add index to map
        ++index;
      } // if
    } // for
    err = DMPlexRestoreTransitiveClosure(dmMesh, cohesiveCells[iCell], PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
  } // for
  assert(size_t(index) == _cohesiveVertices.size());

  PYLITH_METHOD_END;
} // _initializeCohesiveInfo

// ----------------------------------------------------------------------
// Initialize logger.
void
pylith::faults::FaultCohesiveLagrange::_initializeLogger(void)
{ // initializeLogger
  PYLITH_METHOD_BEGIN;

  delete _logger;
  _logger = new utils::EventLogger;
  assert(_logger);
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

  PYLITH_METHOD_END;
} // initializeLogger

// ----------------------------------------------------------------------
// Transform field from local (fault) coordinate system to
// global coordinate system.
void
pylith::faults::FaultCohesiveLagrange::faultToGlobal(topology::Field* field,
						     const topology::Field& faultOrientation)
{ // faultToGlobal
  PYLITH_METHOD_BEGIN;

  assert(field);

  // Fiber dimension of vector field matches spatial dimension.
  const spatialdata::geocoords::CoordSys* cs = field->mesh().coordsys();assert(cs);
  const int spaceDim = cs->spaceDim();
  scalar_array fieldVertexGlobal(spaceDim);

  // Get fields
  topology::VecVisitorMesh fieldVisitor(*field);
  PetscScalar* fieldArray = fieldVisitor.localArray();

  topology::VecVisitorMesh orientationVisitor(faultOrientation);
  const PetscScalar* orientationArray = orientationVisitor.localArray();

  PetscDM dmMesh = field->mesh().dmMesh();assert(dmMesh);
  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();

  for(PetscInt v = vStart; v < vEnd; ++v) {
    const PetscInt foff = fieldVisitor.sectionOffset(v);
    assert(spaceDim == fieldVisitor.sectionDof(v));

    const PetscInt ooff = orientationVisitor.sectionOffset(v);
    assert(spaceDim*spaceDim == orientationVisitor.sectionDof(v));

    // Rotate from fault to global coordinate system (transpose orientation)
    fieldVertexGlobal = 0.0;
    for(int iDim = 0; iDim < spaceDim; ++iDim) {
      for(int jDim = 0; jDim < spaceDim; ++jDim) {
        fieldVertexGlobal[iDim] += orientationArray[ooff+jDim*spaceDim+iDim] * fieldArray[foff+jDim];
      } // for
    } // for
    for(int iDim = 0; iDim < spaceDim; ++iDim) {
      fieldArray[foff+iDim] = fieldVertexGlobal[iDim];
    } // for
  } // for
  PetscLogFlops((vEnd-vStart) * (2*spaceDim*spaceDim) );
  
#if 0 // DEBUGGING
  field->view("FIELD (GLOBAL)");
#endif

  PYLITH_METHOD_END;
} // faultToGlobal

// ----------------------------------------------------------------------
// Transform field from global coordinate system to local (fault)
// coordinate system.
void
pylith::faults::FaultCohesiveLagrange::globalToFault(topology::Field* field,
						     const topology::Field& faultOrientation)
{ // globalToFault
  PYLITH_METHOD_BEGIN;

  assert(field);

  // Fiber dimension of vector field matches spatial dimension.
  const spatialdata::geocoords::CoordSys* cs = field->mesh().coordsys();assert(cs);
  const int spaceDim = cs->spaceDim();
  scalar_array fieldVertexFault(spaceDim);

  // Get fields
  topology::VecVisitorMesh fieldVisitor(*field);
  PetscScalar* fieldArray = fieldVisitor.localArray();

  topology::VecVisitorMesh orientationVisitor(faultOrientation);
  const PetscScalar* orientationArray = orientationVisitor.localArray();

  PetscDM dmMesh = field->mesh().dmMesh();assert(dmMesh);
  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();

  for(PetscInt v = vStart; v < vEnd; ++v) {
    const PetscInt foff = fieldVisitor.sectionOffset(v);
    assert(spaceDim == fieldVisitor.sectionDof(v));

    const PetscInt ooff = orientationVisitor.sectionOffset(v);
    assert(spaceDim*spaceDim == orientationVisitor.sectionDof(v));

    // Rotate from global coordinate system to fault (orientation)
    fieldVertexFault = 0.0;
    for(int iDim = 0; iDim < spaceDim; ++iDim) {
      for(int jDim = 0; jDim < spaceDim; ++jDim) {
        fieldVertexFault[iDim] += orientationArray[ooff+iDim*spaceDim+jDim] * fieldArray[foff+jDim];
      } // for
    } // for
    for(int iDim = 0; iDim < spaceDim; ++iDim) {
      fieldArray[foff+iDim] = fieldVertexFault[iDim];
    } // for
  } // for
  PetscLogFlops((vEnd-vStart) * (2*spaceDim*spaceDim) );

#if 0 // DEBUGGING
  field->view("FIELD (FAULT)");
#endif

  PYLITH_METHOD_END;
} // faultToGlobal


// ----------------------------------------------------------------------
// Check to see if given vertex is clamped.
bool
pylith::faults::FaultCohesiveLagrange::isClampedVertex(PetscDMLabel clamped,
						       PetscInt vertex)
{ // _isClampedVertex
  PYLITH_METHOD_BEGIN;
  
  bool isClamped = false;
  
  if (clamped) {
    PetscInt value = -1;
    PetscErrorCode err = DMLabelGetValue(clamped, vertex, &value);PYLITH_CHECK_ERROR(err);
    if (value >= 0) {
      isClamped = true;
    } // if
  } // if
  
  PYLITH_METHOD_RETURN(isClamped);
} // _isClampedVertex


// ----------------------------------------------------------------------
// Calculate orientation at fault vertices.
void
pylith::faults::FaultCohesiveLagrange::_calcOrientation(const PylithScalar upDir[3])
{ // _calcOrientation
  PYLITH_METHOD_BEGIN;

  assert(upDir);
  assert(_faultMesh);
  assert(_fields);

  scalar_array up(3);
  for (int i=0; i < 3; ++i) {
    up[i] = upDir[i];
  } // for

  const int cohesiveDim = _faultMesh->dimension();
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int orientationSize = spaceDim * spaceDim;
  const feassemble::CellGeometry& cellGeometry = _quadrature->refGeometry();
  const scalar_array& verticesRef = cellGeometry.vertices();
  const int jacobianSize = (cohesiveDim > 0) ? spaceDim * cohesiveDim : 1;
  scalar_array jacobian(jacobianSize);
  PylithScalar jacobianDet = 0;

  assert(cohesiveDim == _quadrature->cellDim());
  assert(spaceDim == _quadrature->spaceDim());

  // Get vertices and cells in fault mesh.
  PetscDM faultDMMesh = _faultMesh->dmMesh();assert(faultDMMesh);
  topology::Stratum verticesStratum(faultDMMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();

  topology::Stratum cellsStratum(faultDMMesh, topology::Stratum::HEIGHT, 1);
  const PetscInt cStart = cellsStratum.begin();
  const PetscInt cEnd = cellsStratum.end();

  // Containers for orientation information.
  // Allocate orientation field.
  scalar_array orientationVertex(orientationSize);
  _fields->add("orientation", "orientation");
  topology::Field& orientation = _fields->get("orientation");
  const topology::Field& dispRel = _fields->get("relative disp");
  if (spaceDim > 1) orientation.subfieldAdd("strike_dir", spaceDim, topology::Field::VECTOR);
  if (spaceDim > 2) orientation.subfieldAdd("dip_dir", spaceDim, topology::Field::VECTOR);
  orientation.subfieldAdd("normal_dir", spaceDim, topology::Field::VECTOR);
  orientation.subfieldsSetup();
  orientation.newSection(dispRel, orientationSize);
  // Create components for along-strike, up-dip, and normal directions
  if (spaceDim > 1) { 
    orientation.subfieldSetDof("strike_dir", topology::FieldBase::VERTICES_FIELD, spaceDim);
  } // if
  if (spaceDim > 2) {
    orientation.subfieldSetDof("dip_dir", topology::FieldBase::VERTICES_FIELD, spaceDim);
  } // if
  orientation.subfieldSetDof("normal_dir", topology::FieldBase::VERTICES_FIELD, spaceDim);
  orientation.allocate();
  orientation.zeroAll();

  topology::VecVisitorMesh orientationVisitor(orientation);
  PetscScalar* orientationArray = orientationVisitor.localArray();

  // Get section containing coordinates of vertices
  scalar_array coordsCell(numBasis*spaceDim); // :KULDGE: Update numBasis to numCorners after implementing higher order
  topology::CoordsVisitor coordsVisitor(faultDMMesh);

  // Loop over cohesive cells, computing orientation weighted by
  // jacobian at constraint vertices. This involves looping over cells
  // and summing across processors (complete the section) just like a
  // normal FE integration.

  PetscInt *closure = NULL;
  PetscInt closureSize = 0;
  PetscInt cell;
  try {
    for(cell = cStart; cell < cEnd; ++cell) {
      
      // Get orientations at fault cell's vertices.
      coordsVisitor.getClosure(&coordsCell, cell);
      
      PetscErrorCode err = DMPlexGetTransitiveClosure(faultDMMesh, cell, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
      
      // Filter out non-vertices
      PetscInt numVertices = 0;
      for(PetscInt p = 0; p < closureSize*2; p += 2) {
        if ((closure[p] >= vStart) && (closure[p] < vEnd)) {
          closure[numVertices*2]   = closure[p];
          closure[numVertices*2+1] = closure[p+1];
          ++numVertices;
        } // if
      } // for

      for(PetscInt v = 0; v < numVertices; ++v) {
        const PetscInt v_fault = closure[v*2];

        // Compute Jacobian and determinant of Jacobian at vertex
        cellGeometry.jacobian(&jacobian, &jacobianDet, &coordsCell[0], numBasis, spaceDim, &verticesRef[v*cohesiveDim], cohesiveDim);
	
        // Compute orientation
        cellGeometry.orientation(&orientationVertex, jacobian, jacobianDet, up);
	
        // Update orientation
        const PetscInt ooff = orientationVisitor.sectionOffset(v_fault);
	
        for(PetscInt d = 0; d < orientationSize; ++d) {
          orientationArray[ooff+d] += orientationVertex[d];
        } // for
      } // for
      err = DMPlexRestoreTransitiveClosure(faultDMMesh, cell, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
    } // for
  } catch (const std::exception& err) {
    if (closure) {
      DMPlexRestoreTransitiveClosure(faultDMMesh, cell, PETSC_TRUE, &closureSize, &closure);
    } // if
    throw;
  } catch (...) {
    throw;
  } // try/catch
  orientationVisitor.clear();

  //orientation.view("ORIENTATION BEFORE COMPLETE"); // DEBUGGING

  // Assemble orientation information
  orientation.complete();

  // Loop over vertices, make orientation information unit magnitude
  scalar_array vertexDir(orientationSize);
  orientationVisitor.initialize(orientation);
  orientationArray = orientationVisitor.localArray();
  int count = 0;
  for(PetscInt v = vStart; v < vEnd; ++v, ++count) {
    assert(orientationSize == orientationVisitor.sectionDof(v));
    const PetscInt ooff = orientationVisitor.sectionOffset(v);
    for(PetscInt d = 0; d < orientationSize; ++d) {
      orientationVertex[d] = orientationArray[ooff+d];
    } // for
    for (int iDim = 0; iDim < spaceDim; ++iDim) {
      PylithScalar mag = 0;
      for (int jDim = 0, index = iDim * spaceDim; jDim < spaceDim; ++jDim) {
        mag += pow(orientationVertex[index + jDim], 2);
      } // for

      if (mag <= 0.0) {
        std::ostringstream msg;
        msg << "Error calculating fault orientation at fault vertex " << v << ".\n" 
            << "Zero vector in parallel likely indicates inconsistent fault orientation (creation) across processors.\n"
            << "Orientation vector " << iDim << ": (";
        for (int jDim = 0, index = iDim * spaceDim; jDim < spaceDim; ++jDim) {
          msg << " " << orientationVertex[index + jDim];
        } // for
        msg << " )" << std::endl;
        throw std::runtime_error(msg.str());
      } // if

      mag = sqrt(mag);
      for (int jDim = 0, index = iDim * spaceDim; jDim < spaceDim; ++jDim)
        orientationVertex[index + jDim] /= mag;
    } // for
    for(PetscInt d = 0; d < orientationSize; ++d) {
      orientationArray[ooff+d] = orientationVertex[d];
    } // for
  } // for
  orientationVisitor.clear();
  PetscLogFlops(count * orientationSize * 4);

  orientationVisitor.initialize(orientation);
  orientationArray = orientationVisitor.localArray();
  MPI_Comm comm = orientation.mesh().comm();
  if (1 == cohesiveDim) {
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
    int flipLocal = 0;
    if (vEnd > vStart) {
      const PetscInt ooff = orientationVisitor.sectionOffset(vStart);
      assert(orientationSize == orientationVisitor.sectionDof(vStart));
      for (PetscInt d = 0; d < orientationSize; ++d) {
	orientationVertex[d] = orientationArray[ooff+d];
      } // for

      assert(2 == spaceDim);
      const PylithScalar* shearDirVertex = &orientationVertex[0];
      const PylithScalar* normalDirVertex = &orientationVertex[2];
      const PylithScalar shearDirDot = upDir[0] * shearDirVertex[0] + upDir[1] * shearDirVertex[1];
      const PylithScalar normalDirDot = upDir[0] * normalDirVertex[0] + upDir[1] * normalDirVertex[1];
      if (normalDirDot * shearDirDot < 0.0) {
	flipLocal = 1;
      } // if
    } // if
    // Collect flip decisions across all processors
    int flipGlobal = 0;
    MPI_Allreduce(&flipLocal, &flipGlobal, 1, MPI_INT, MPI_SUM, comm);
    
    const int ishear = 0;
    if (flipGlobal > 0) { // flip in any processor wants to flip
      // Flip shear direction
      for(PetscInt v = vStart; v < vEnd; ++v) {
	assert(4 == orientationVisitor.sectionDof(v));
	const PetscInt ooff = orientationVisitor.sectionOffset(v);
        for(PetscInt d = 0; d < orientationSize; ++d) {
          orientationVertex[d] = orientationArray[ooff+d];
        } // for
        for (int iDim = 0; iDim < 2; ++iDim) // flip shear
          orientationVertex[ishear + iDim] *= -1.0;
	
        // Update orientation
        for(PetscInt d = 0; d < orientationSize; ++d) {
          orientationArray[ooff+d] = orientationVertex[d];
        } // for
      } // for
      PetscLogFlops(3 + count * 2);
    } // if

  } else if (2 == cohesiveDim) {
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

    int flipLocal = 0;
    if (vEnd > vStart) {
      const PetscInt ooff = orientationVisitor.sectionOffset(vStart);
      for (PetscInt d = 0; d < orientationSize; ++d) {
	orientationVertex[d] = orientationArray[ooff+d];
      } // for
      
      assert(3 == spaceDim);
      const PylithScalar* dipDirVertex = &orientationVertex[3];
      const PylithScalar* normalDirVertex = &orientationVertex[6];
      const PylithScalar dipDirDot = upDir[0]*dipDirVertex[0] + upDir[1]*dipDirVertex[1] + upDir[2]*dipDirVertex[2];
      const PylithScalar normalDirDot = upDir[0]*normalDirVertex[0] + upDir[1]*normalDirVertex[1] + upDir[2]*normalDirVertex[2];

      if (dipDirDot * normalDirDot < 0.0 || fabs(normalDirVertex[2] + 1.0) < 0.001) {
	// if fault normal is (0,0,+-1) then up-dir dictates reverse
	// motion for case with normal (0,0,1), so we reverse the dip-dir
	// if we have (0,0,-1).
	flipLocal = 1;
      } // if
    } // if

    // Collect flip decisions across all processors
    int flipGlobal = 0;
    MPI_Allreduce(&flipLocal, &flipGlobal, 1, MPI_INT, MPI_SUM, comm);

    const int idip = 3;
    if (flipGlobal > 0) { // flip in any processor wants to flip
      // Flip dip direction
      for(PetscInt v = vStart; v < vEnd; ++v) {
	assert(9 == orientationVisitor.sectionDof(v));
	const PetscInt ooff = orientationVisitor.sectionOffset(v);
        for(PetscInt d = 0; d < orientationSize; ++d) {
          orientationVertex[d] = orientationArray[ooff+d];
        } // for
        for (int iDim = 0; iDim < 3; ++iDim) // flip dip
          orientationVertex[idip + iDim] *= -1.0;
	
        // Update direction
        for(PetscInt d = 0; d < orientationSize; ++d) {
          orientationArray[ooff+d] = orientationVertex[d];
        } // for
      } // for
      PetscLogFlops(5 + count * 3);
    } // if
  } // if

  //orientation.view("ORIENTATION");

  PYLITH_METHOD_END;
} // _calcOrientation

// ----------------------------------------------------------------------
void
pylith::faults::FaultCohesiveLagrange::_calcArea(void)
{ // _calcArea
  PYLITH_METHOD_BEGIN;

  assert(_faultMesh);
  assert(_fields);
  assert(_normalizer);

  // :TODO: Update for higher order?

  // Containers for area information
  const int numBasis = _quadrature->numBasis();
  const int numQuadPts = _quadrature->numQuadPts();
  const int spaceDim = _quadrature->spaceDim();
  const scalar_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == size_t(numQuadPts));

  // Get vertices in fault mesh.
  PetscDM faultDMMesh = _faultMesh->dmMesh();assert(faultDMMesh);
  topology::Stratum verticesStratum(faultDMMesh, topology::Stratum::DEPTH, 0);

  topology::Stratum cellsStratum(faultDMMesh, topology::Stratum::HEIGHT, 1);
  const PetscInt cStart = cellsStratum.begin();
  const PetscInt cEnd = cellsStratum.end();

  // Allocate area field.
  _fields->add("area", "area");
  topology::Field& area = _fields->get("area");
  const topology::Field& dispRel = _fields->get("relative disp");
  area.newSection(dispRel, 1);
  area.allocate();
  area.vectorFieldType(topology::FieldBase::SCALAR);
  const PylithScalar lengthScale = _normalizer->lengthScale();
  area.scale(pow(lengthScale, (spaceDim-1)));
  area.zeroAll();

  topology::VecVisitorMesh areaVisitor(area);
  scalar_array areaCell(numBasis);

  scalar_array coordsCell(numBasis*spaceDim); // :KULDGE: Update numBasis to numCorners after implementing higher order
  topology::CoordsVisitor coordsVisitor(faultDMMesh);

  // Loop over cells in fault mesh, compute area
  for(PetscInt c = cStart; c < cEnd; ++c) {
    areaCell = 0.0;

    // Compute geometry information for current cell
    coordsVisitor.getClosure(&coordsCell, c);
    _quadrature->computeGeometry(&coordsCell[0], coordsCell.size(), c);

    // Get cell geometry information that depends on cell
    const scalar_array& basis = _quadrature->basis();
    const scalar_array& jacobianDet = _quadrature->jacobianDet();

    // Compute area
    for (int iQuad = 0; iQuad < numQuadPts; ++iQuad) {
      const PylithScalar wt = quadWts[iQuad] * jacobianDet[iQuad];
      for (int iBasis = 0; iBasis < numBasis; ++iBasis) {
        const PylithScalar dArea = wt * basis[iQuad * numBasis + iBasis];
        areaCell[iBasis] += dArea;
      } // for
    } // for
    areaVisitor.setClosure(&areaCell[0], areaCell.size(), c, ADD_VALUES);
  } // for
  PetscLogFlops(numQuadPts*(1+numBasis*2)*(cEnd-cStart));

  // Assemble area information
  area.complete();

#if 0 // DEBUGGING
  area.view("AREA");
#endif

  PYLITH_METHOD_END;
} // _calcArea

// ----------------------------------------------------------------------
// Compute change in tractions on fault surface using solution.
void
pylith::faults::FaultCohesiveLagrange::_calcTractionsChange(topology::Field* tractions,
							    const topology::Field& dispT)
{ // _calcTractionsChange
  PYLITH_METHOD_BEGIN;

  assert(tractions);
  assert(_faultMesh);
  assert(_fields);
  assert(_normalizer);

  tractions->label("traction_change");
  tractions->scale(_normalizer->pressureScale());

  // Fiber dimension of tractions matches spatial dimension.
  const int spaceDim = _quadrature->spaceDim();

  // Get fields
  topology::VecVisitorMesh dispTVisitor(dispT);
  const PetscScalar* dispTArray = dispTVisitor.localArray();

  topology::Field& orientation = _fields->get("orientation");
  topology::VecVisitorMesh orientationVisitor(orientation);
  const PetscScalar* orientationArray = orientationVisitor.localArray();

  // Allocate buffer for tractions field (if necessary).
  if (!tractions->localSection()) {
    const topology::Field& dispRel = _fields->get("relative disp");
    tractions->cloneSection(dispRel);
  } // if
  tractions->zeroAll();

  topology::VecVisitorMesh tractionsVisitor(*tractions);
  PetscScalar* tractionsArray = tractionsVisitor.localArray();

  const int numVertices = _cohesiveVertices.size();
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int e_lagrange = _cohesiveVertices[iVertex].lagrange;
    const int v_fault = _cohesiveVertices[iVertex].fault;

    if (e_lagrange < 0) { // skip clamped edges
      continue;
    } // if

    const PetscInt dtoff = dispTVisitor.sectionOffset(e_lagrange);
    assert(spaceDim == dispTVisitor.sectionDof(e_lagrange));

    const PetscInt ooff = orientationVisitor.sectionOffset(v_fault);
    assert(spaceDim*spaceDim == orientationVisitor.sectionDof(v_fault));

    const PetscInt toff = tractionsVisitor.sectionOffset(v_fault);
    assert(spaceDim == tractionsVisitor.sectionDof(v_fault));

    // Rotate from global coordinate system to fault (orientation)
    for(PetscInt iDim = 0; iDim < spaceDim; ++iDim) {
      tractionsArray[toff+iDim] = 0.0;
      for(PetscInt jDim = 0; jDim < spaceDim; ++jDim)
        tractionsArray[toff+iDim] += orientationArray[ooff+iDim*spaceDim+jDim] * dispTArray[dtoff+jDim];
    }
  } // for

  PetscLogFlops(numVertices * (1 + spaceDim) );

#if 0 // DEBUGGING
  tractions->view("TRACTIONS");
#endif

  PYLITH_METHOD_END;
} // _calcTractionsChange

// ----------------------------------------------------------------------
// Allocate buffer for vector field.
void
pylith::faults::FaultCohesiveLagrange::_allocateBufferVectorField(void)
{ // _allocateBufferVectorField
  PYLITH_METHOD_BEGIN;

  assert(_fields);
  if (_fields->hasField("buffer (vector)"))
    PYLITH_METHOD_END;

  // Create vector field; use same shape/chart as relative
  // displacement field.
  assert(_faultMesh);
  _fields->add("buffer (vector)", "buffer");
  topology::Field& buffer = _fields->get("buffer (vector)");
  const topology::Field& dispRel = _fields->get("relative disp");
  buffer.cloneSection(dispRel);
  buffer.zeroAll();
  assert(buffer.vectorFieldType() == topology::FieldBase::VECTOR);

  PYLITH_METHOD_END;
} // _allocateBufferVectorField

// ----------------------------------------------------------------------
// Allocate buffer for scalar field.
void
pylith::faults::FaultCohesiveLagrange::_allocateBufferScalarField(void)
{ // _allocateBufferScalarField
  PYLITH_METHOD_BEGIN;

  assert(_fields);
  if (_fields->hasField("buffer (scalar)"))
    PYLITH_METHOD_END;

  // Create vector field; use same shape/chart as area field.
  assert(_faultMesh);
  _fields->add("buffer (scalar)", "buffer");
  topology::Field& buffer = _fields->get("buffer (scalar)");
  buffer.newSection(topology::FieldBase::VERTICES_FIELD, 1); // :TODO: Update?
  buffer.allocate();
  buffer.vectorFieldType(topology::FieldBase::SCALAR);
  buffer.scale(1.0);
  buffer.zeroAll();
  assert(buffer.vectorFieldType() == topology::FieldBase::SCALAR);

  PYLITH_METHOD_END;
} // _allocateBufferScalarField

// ----------------------------------------------------------------------
//  Get submatrix of Jacobian matrix associated with the negative and
//  positive sides of the fault.
void
pylith::faults::FaultCohesiveLagrange::_getJacobianSubmatrixNP(PetscMat* jacobianSub,
							       std::map<int,int>* indicesMatToSubmat,
							       const topology::Jacobian& jacobian,
							       const topology::SolutionFields& fields)
{ // _getJacobianSubmatrixNP
  PYLITH_METHOD_BEGIN;

  assert(jacobianSub);
  assert(indicesMatToSubmat);

  // Get global order
  PetscSection solutionGlobalSection = fields.solution().globalSection();assert(solutionGlobalSection);

  PetscErrorCode err;

  // Get Jacobian matrix
  const PetscMat jacobianMatrix = jacobian.matrix();
  assert(jacobianMatrix);

  const spatialdata::geocoords::CoordSys* cs = fields.mesh().coordsys();assert(cs);
  const int spaceDim = cs->spaceDim();

  const int numVertices = _cohesiveVertices.size();
  int numIndicesNP = 0;
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int e_lagrange = _cohesiveVertices[iVertex].lagrange;
    if (e_lagrange < 0) { // Ignore clamped edges.
      continue;
    } // if

    // Compute contribution only if Lagrange constraint is local.
    PetscInt goff = 0;
    err = PetscSectionGetOffset(solutionGlobalSection, e_lagrange, &goff);PYLITH_CHECK_ERROR(err);
    if (goff < 0) {
      continue;
    } // if

    numIndicesNP += 2;
  } // for
  int_array indicesNP(numIndicesNP*spaceDim);

  for (int iVertex=0, indexNP=0; iVertex < numVertices; ++iVertex) {
    const int e_lagrange = _cohesiveVertices[iVertex].lagrange;
    const int v_negative = _cohesiveVertices[iVertex].negative;
    const int v_positive = _cohesiveVertices[iVertex].positive;

    if (e_lagrange < 0) { // Ignore clamped edges.
      continue;
    } // if

    // Compute contribution only if Lagrange constraint is local.
    PetscInt gloff = 0;
    err = PetscSectionGetOffset(solutionGlobalSection, e_lagrange, &gloff);PYLITH_CHECK_ERROR(err);
    if (gloff < 0)
      continue;

    PetscInt gnoff = 0;
    err = PetscSectionGetOffset(solutionGlobalSection, v_negative, &gnoff);PYLITH_CHECK_ERROR(err);
    gnoff = gnoff < 0 ? -(gnoff+1) : gnoff;

    PetscInt gpoff = 0;
    err = PetscSectionGetOffset(solutionGlobalSection, v_positive, &gpoff);PYLITH_CHECK_ERROR(err);
    gpoff = gpoff < 0 ? -(gpoff+1) : gpoff;
    
    // Set global order indices
    for(int iDim=0; iDim < spaceDim; ++iDim)
      indicesNP[indexNP*spaceDim+iDim] = gnoff + iDim;
    ++indexNP;
    for(int iDim=0; iDim < spaceDim; ++iDim)
      indicesNP[indexNP*spaceDim+iDim] = gpoff + iDim;
    ++indexNP;
  } // for
  
  // MatGetSubMatrices requires sorted indices
  std::sort(&indicesNP[0], &indicesNP[indicesNP.size()]);  
  
  PetscMat* subMat[1];
  IS indicesIS[1];
  err = ISCreateGeneral(PETSC_COMM_SELF, indicesNP.size(), &indicesNP[0], PETSC_USE_POINTER, &indicesIS[0]);
  PYLITH_CHECK_ERROR(err);
  err = MatGetSubMatrices(jacobianMatrix, 1, indicesIS, indicesIS, MAT_INITIAL_MATRIX, subMat);
  PYLITH_CHECK_ERROR(err);
  err = ISDestroy(&indicesIS[0]); PYLITH_CHECK_ERROR(err);

  *jacobianSub = *subMat[0];
  err = PetscObjectReference((PetscObject) *subMat[0]); 
  PYLITH_CHECK_ERROR(err);
  err = MatDestroyMatrices(1, &subMat[0]); PYLITH_CHECK_ERROR(err);

  // Create map from global indices to local indices (using only the
  // first index as to match the global order.
  indicesMatToSubmat->clear();
  const int indicesNPSize = indicesNP.size();
  for (int i=0; i < indicesNPSize; i+=spaceDim)
    (*indicesMatToSubmat)[indicesNP[i]] = i;

  PYLITH_METHOD_END;
} // _getJacobianSubmatrixNP


// ----------------------------------------------------------------------
// Get cell field associated with integrator.
const pylith::topology::Field&
pylith::faults::FaultCohesiveLagrange::cellField(const char* name,
                                                 const topology::SolutionFields* fields)
{ // cellField
  PYLITH_METHOD_BEGIN;

  if (0 == strcasecmp("partition", name)) {

    PetscDM faultDMMesh = _faultMesh->dmMesh();assert(faultDMMesh);
    topology::Stratum cellsStratum(faultDMMesh, topology::Stratum::HEIGHT, 1);
    const PetscInt cStart = cellsStratum.begin();
    const PetscInt cEnd = cellsStratum.end();

    const int fiberDim = 1;
    _fields->add("partition", "partition", pylith::topology::FieldBase::CELLS_FIELD, fiberDim);
    topology::Field& partition = _fields->get("partition");
    partition.allocate();
    topology::VecVisitorMesh partitionVisitor(partition);
    PetscScalar* partitionArray = partitionVisitor.localArray();

    PetscMPIInt rank;
    PetscErrorCode err = MPI_Comm_rank(_faultMesh->comm(), &rank);PYLITH_CHECK_ERROR(err);
    // Loop over cells in fault mesh, set partition
    for(PetscInt c = cStart; c < cEnd; ++c) {
      const PetscInt off = partitionVisitor.sectionOffset(c);
      assert(1 == partitionVisitor.sectionDof(c));
      partitionArray[off] = rank;
    } // for

    PYLITH_METHOD_RETURN(partition);
  } // if

  // Should not reach this point if requested field was found
  std::ostringstream msg;
  msg << "Request for unknown cell field '" << name << "' for fault '" << label() << ".";
  throw std::runtime_error(msg.str());

  // Satisfy return value
  assert(_fields);
  const topology::Field& buffer = _fields->get("buffer (vector)");
  PYLITH_METHOD_RETURN(buffer);
} // cellField


// End of file 
