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
typedef ALE::ISieveVisitor::IndicesVisitor<RealSection,SieveSubMesh::order_type,PylithInt> IndicesVisitor;

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
					     const PylithScalar upDir[3])
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
  logger.stagePush("FaultFields");

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

  logger.stagePop();

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
  _calcOrientation(upDir);

  // Compute tributary area for each vertex in fault mesh.
  _calcArea();

} // initialize

// ----------------------------------------------------------------------
void
pylith::faults::FaultCohesiveLagrange::splitField(topology::Field<topology::Mesh>* field)
{ // splitField
  assert(0 != field);

  DM             dm       = field->dmMesh();
  PetscSection   section  = field->petscSection();
  const PetscInt spaceDim = field->mesh().dimension();
  PetscInt       numFields, numComp, pStart, pEnd;
  PetscErrorCode err;

  err = PetscSectionGetNumFields(section, &numFields);CHECK_PETSC_ERROR(err);
  // TODO: Does this make sense?
  if (!numFields) return;
  assert(numFields == 2);
  err = PetscSectionGetFieldComponents(section, 0, &numComp);CHECK_PETSC_ERROR(err);
  assert(numComp == spaceDim);
  err = PetscSectionGetFieldComponents(section, 1, &numComp);CHECK_PETSC_ERROR(err);
  assert(numComp == spaceDim);

  const int numVertices = _cohesiveVertices.size();
  for(PetscInt iVertex = 0; iVertex < numVertices; ++iVertex) {
    PetscInt dof;

    const int v_lagrange = _cohesiveVertices[iVertex].lagrange;
    err = PetscSectionGetDof(section, v_lagrange, &dof);CHECK_PETSC_ERROR(err);
    assert(spaceDim == dof);
    err = PetscSectionSetFieldDof(section, v_lagrange, 1, dof);CHECK_PETSC_ERROR(err);
  } // for
  err = PetscSectionGetChart(section, &pStart, &pEnd);CHECK_PETSC_ERROR(err);
  for(PetscInt p = pStart; p < pEnd; ++p) {
    PetscInt dof;

    err = PetscSectionGetFieldDof(section, p, 1, &dof);CHECK_PETSC_ERROR(err);
    if (!dof) {err = PetscSectionSetFieldDof(section, p, 0, spaceDim);CHECK_PETSC_ERROR(err);}
  }
} // splitField

// ----------------------------------------------------------------------
// Integrate contribution of cohesive cells to residual term.
void
pylith::faults::FaultCohesiveLagrange::integrateResidual(
			 const topology::Field<topology::Mesh>& residual,
			 const PylithScalar t,
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

  const int setupEvent = _logger->eventId("FaIR setup");
  const int geometryEvent = _logger->eventId("FaIR geometry");
  const int computeEvent = _logger->eventId("FaIR compute");
  const int restrictEvent = _logger->eventId("FaIR restrict");
  const int updateEvent = _logger->eventId("FaIR update");

  _logger->eventBegin(setupEvent);

  // Get cell geometry information that doesn't depend on cell
  const int spaceDim = _quadrature->spaceDim();

  // Get sections associated with cohesive cells
  DM           residualDM      = residual.dmMesh();
  PetscSection residualSection = residual.petscSection();
  Vec          residualVec     = residual.localVector();
  PetscSection residualGlobalSection;
  PetscScalar *residualArray;
  PetscErrorCode err;
  assert(residualSection);assert(residualVec);
  err = DMGetDefaultGlobalSection(residualDM, &residualGlobalSection);CHECK_PETSC_ERROR(err);

  PetscSection dispTSection = fields->get("disp(t)").petscSection();
  Vec          dispTVec     = fields->get("disp(t)").localVector();
  PetscScalar *dispTArray;
  assert(dispTSection);assert(dispTVec);

  PetscSection dispTIncrSection = fields->get("dispIncr(t->t+dt)").petscSection();
  Vec          dispTIncrVec     = fields->get("dispIncr(t->t+dt)").localVector();
  PetscScalar *dispTIncrArray;
  assert(dispTIncrSection);assert(dispTIncrVec);

  PetscSection dispRelSection = _fields->get("relative disp").petscSection();
  Vec          dispRelVec     = _fields->get("relative disp").localVector();
  PetscScalar *dispRelArray;
  assert(dispRelSection);assert(dispRelVec);

  PetscSection areaSection = _fields->get("area").petscSection();
  Vec          areaVec     = _fields->get("area").localVector();
  PetscScalar *areaArray;
  assert(areaSection);assert(areaVec);

  // Get fault information
  DM dmMesh = fields->mesh().dmMesh();

  _logger->eventEnd(setupEvent);
#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventBegin(computeEvent);
#endif

  // Loop over fault vertices
  const int numVertices = _cohesiveVertices.size();
  err = VecGetArray(dispRelVec, &dispRelArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(areaVec, &areaArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(dispTVec, &dispTArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(dispTIncrVec, &dispTIncrArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(residualVec, &residualArray);CHECK_PETSC_ERROR(err);
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_lagrange = _cohesiveVertices[iVertex].lagrange;
    const int v_fault = _cohesiveVertices[iVertex].fault;
    const int v_negative = _cohesiveVertices[iVertex].negative;
    const int v_positive = _cohesiveVertices[iVertex].positive;
    PetscInt goff;

    // Compute contribution only if Lagrange constraint is local.
    err = PetscSectionGetOffset(residualGlobalSection, v_lagrange, &goff);CHECK_PETSC_ERROR(err);
    if (goff < 0) continue;

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(restrictEvent);
#endif

    // Get relative dislplacement at fault vertex.
    PetscInt drdof, droff;

    err = PetscSectionGetDof(dispRelSection, v_fault, &drdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dispRelSection, v_fault, &droff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == drdof);

    // Get area associated with fault vertex.
    PetscInt adof, aoff;

    err = PetscSectionGetDof(areaSection, v_fault, &adof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(areaSection, v_fault, &aoff);CHECK_PETSC_ERROR(err);
    assert(1 == adof);

    // Get disp(t) at conventional vertices and Lagrange vertex.
    PetscInt dtndof, dtnoff;

    err = PetscSectionGetDof(dispTSection, v_negative, &dtndof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dispTSection, v_negative, &dtnoff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == dtndof);
    PetscInt dtpdof, dtpoff;

    err = PetscSectionGetDof(dispTSection, v_positive, &dtpdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dispTSection, v_positive, &dtpoff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == dtpdof);
    PetscInt dtldof, dtloff;

    err = PetscSectionGetDof(dispTSection, v_lagrange, &dtldof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dispTSection, v_lagrange, &dtloff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == dtldof);

    // Get dispIncr(t->t+dt) at conventional vertices and Lagrange vertex.
    PetscInt dindof, dinoff;

    err = PetscSectionGetDof(dispTIncrSection, v_negative, &dindof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dispTIncrSection, v_negative, &dinoff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == dindof);
    PetscInt dipdof, dipoff;

    err = PetscSectionGetDof(dispTIncrSection, v_positive, &dipdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dispTIncrSection, v_positive, &dipoff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == dipdof);
    PetscInt dildof, diloff;

    err = PetscSectionGetDof(dispTIncrSection, v_lagrange, &dildof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dispTIncrSection, v_lagrange, &diloff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == dildof);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(restrictEvent);
    _logger->eventBegin(computeEvent);
#endif

    // Compute current estimate of displacement at time t+dt using solution increment.
    for(PetscInt d = 0; d < spaceDim; ++d) {
    } // for

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(computeEvent);
    _logger->eventBegin(updateEvent);
#endif

    // Assemble contributions into field
    PetscInt rndof, rnoff, rpdof, rpoff, rldof, rloff;

    err = PetscSectionGetDof(residualSection, v_negative, &rndof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(residualSection, v_negative, &rnoff);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetDof(residualSection, v_positive, &rpdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(residualSection, v_positive, &rpoff);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetDof(residualSection, v_lagrange, &rldof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(residualSection, v_lagrange, &rloff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == rndof);assert(spaceDim == rpdof);assert(spaceDim == rldof);
    for(PetscInt d = 0; d < spaceDim; ++d) {
      residualArray[rnoff+d] +=  areaArray[aoff] * (dispTArray[dtloff+d] + dispTIncrArray[diloff+d]);
      residualArray[rpoff+d] += -areaArray[aoff] * (dispTArray[dtloff+d] + dispTIncrArray[diloff+d]);
      residualArray[rloff+d] += -areaArray[aoff] * (dispTArray[dtpoff+d] + dispTIncrArray[dipoff+d] - dispTArray[dtnoff+d] - dispTIncrArray[dinoff+d] - dispRelArray[droff+d]);
    }

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(updateEvent);
#endif
  } // for
  PetscLogFlops(numVertices*spaceDim*10);
  err = VecRestoreArray(dispRelVec, &dispRelArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(areaVec, &areaArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(dispTVec, &dispTArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(dispTIncrVec, &dispTIncrArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(residualVec, &residualArray);CHECK_PETSC_ERROR(err);

#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventEnd(computeEvent);
#endif
} // integrateResidual

// ----------------------------------------------------------------------
// Compute Jacobian matrix (A) associated with operator.
void
pylith::faults::FaultCohesiveLagrange::integrateJacobian(
				   topology::Jacobian* jacobian,
				   const PylithScalar t,
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
  const int spaceDim = _quadrature->spaceDim();

  // Get sections
  PetscSection areaSection = _fields->get("area").petscSection();
  Vec          areaVec     = _fields->get("area").localVector();
  PetscScalar *areaArray;
  assert(areaSection);assert(areaVec);

  DM           solnDM      = fields->solution().dmMesh();
  PetscSection solnSection = fields->solution().petscSection();
  Vec          solnVec     = fields->solution().localVector();
  PetscSection solnGlobalSection;
  PetscScalar *solnArray;
  PetscErrorCode err;
  assert(solnSection);assert(solnVec);
  err = DMGetDefaultGlobalSection(solnDM, &solnGlobalSection);CHECK_PETSC_ERROR(err);

  // Get fault information
  DM dmMesh = fields->mesh().dmMesh();
  assert(dmMesh);

  // Allocate vectors for vertex values
  scalar_array jacobianVertex(spaceDim*spaceDim);
  int_array indicesL(spaceDim);
  int_array indicesN(spaceDim);
  int_array indicesP(spaceDim);
  int_array indicesRel(spaceDim);
  for (int i=0; i < spaceDim; ++i)
    indicesRel[i] = i;

  // Get sparse matrix
  const PetscMat jacobianMatrix = jacobian->matrix();
  assert(0 != jacobianMatrix);

  _logger->eventEnd(setupEvent);
#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventBegin(computeEvent);
#endif

  const int numVertices = _cohesiveVertices.size();
  err = VecGetArray(areaVec, &areaArray);CHECK_PETSC_ERROR(err);
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_lagrange = _cohesiveVertices[iVertex].lagrange;
    const int v_fault = _cohesiveVertices[iVertex].fault;
    const int v_negative = _cohesiveVertices[iVertex].negative;
    const int v_positive = _cohesiveVertices[iVertex].positive;
    PetscInt gnoff, gpoff, gloff;

    // Compute contribution only if Lagrange constraint is local.
    err = PetscSectionGetOffset(solnGlobalSection, v_lagrange, &gloff);CHECK_PETSC_ERROR(err);
    if (gloff < 0) continue;
    err = PetscSectionGetOffset(solnGlobalSection, v_negative, &gnoff);CHECK_PETSC_ERROR(err);
    gnoff = gnoff < 0 ? -(gnoff+1) : gnoff;
    err = PetscSectionGetOffset(solnGlobalSection, v_positive, &gpoff);CHECK_PETSC_ERROR(err);
    gpoff = gpoff < 0 ? -(gpoff+1) : gpoff;

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(restrictEvent);
#endif

    // Get area associated with fault vertex.
    PetscInt adof, aoff;

    err = PetscSectionGetDof(areaSection, v_fault, &adof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(areaSection, v_fault, &aoff);CHECK_PETSC_ERROR(err);
    assert(1 == adof);

    // Set global order indices
    indicesL = indicesRel + gloff;
    indicesN = indicesRel + gnoff;
    indicesP = indicesRel + gpoff;
    PetscInt cdof;
    err = PetscSectionGetConstraintDof(solnSection, v_negative, &cdof);CHECK_PETSC_ERROR(err);
    assert(0 == cdof);
    err = PetscSectionGetConstraintDof(solnSection, v_positive, &cdof);CHECK_PETSC_ERROR(err);
    assert(0 == cdof);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(restrictEvent);
    _logger->eventBegin(updateEvent);
#endif

    // Set diagonal entries of Jacobian at positive vertex to area
    // associated with vertex.
    for (int iDim=0; iDim < spaceDim; ++iDim)
      jacobianVertex[iDim*spaceDim+iDim] = areaArray[aoff];

    // Values at positive vertex, entry L,P in Jacobian
    MatSetValues(jacobianMatrix,
                 indicesL.size(), &indicesL[0],
                 indicesP.size(), &indicesP[0],
                 &jacobianVertex[0],
                 ADD_VALUES);

    // Values at positive vertex, entry P,L in Jacobian
    MatSetValues(jacobianMatrix,
                 indicesP.size(), &indicesP[0],
                 indicesL.size(), &indicesL[0],
                 &jacobianVertex[0],
                 ADD_VALUES);

    // Values at negative vertex, entry L,N in Jacobian
    jacobianVertex *= -1.0;
    MatSetValues(jacobianMatrix,
                 indicesL.size(), &indicesL[0],
                 indicesN.size(), &indicesN[0],
                 &jacobianVertex[0],
                 ADD_VALUES);

    // Values at negative vertex, entry N,L in Jacobian
    MatSetValues(jacobianMatrix,
                 indicesN.size(), &indicesN[0],
                 indicesL.size(), &indicesL[0],
                 &jacobianVertex[0],
                 ADD_VALUES);

    // Values at Lagrange vertex, entry L,L in Jacobian
    // We must have entries on the diagonal.
    jacobianVertex = 0.0;
    MatSetValues(jacobianMatrix,
                 indicesL.size(), &indicesL[0],
                 indicesL.size(), &indicesL[0],
                 &jacobianVertex[0],
                 ADD_VALUES);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(updateEvent);
#endif

  } // for
  err = VecRestoreArray(areaVec, &areaArray);CHECK_PETSC_ERROR(err);
  PetscLogFlops(numVertices*spaceDim*2);

#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventEnd(computeEvent);
#endif

  _needNewJacobian = false;

#if 0 // DEBUGGING
  sieveMesh->getSendOverlap()->view("Send domain overlap");
  sieveMesh->getRecvOverlap()->view("Receive domain overlap");
#endif
} // integrateJacobian

// ----------------------------------------------------------------------
// Compute Jacobian matrix (A) associated with operator.
void
pylith::faults::FaultCohesiveLagrange::integrateJacobian(
	                  topology::Field<topology::Mesh>* jacobian,
			  const PylithScalar t,
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

  // Add ones to diagonal Jacobian matrix (as field) for
  // convenience. Instead of including the constraints in the Jacobian
  // matrix, we adjust the solution to account for the Lagrange
  // multipliers as part of the solve.

  const PetscInt spaceDim      = _quadrature->spaceDim();
  DM           jacobianDM      = jacobian->dmMesh();
  PetscSection jacobianSection = jacobian->petscSection();
  Vec          jacobianVec     = jacobian->localVector();
  PetscSection jacobianGlobalSection;
  PetscScalar *jacobianArray;
  PetscErrorCode err;
  assert(jacobianSection);assert(jacobianVec);
  err = DMGetDefaultGlobalSection(jacobianDM, &jacobianGlobalSection);CHECK_PETSC_ERROR(err);

  _logger->eventEnd(setupEvent);
#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventBegin(computeEvent);
#endif

  const int numVertices = _cohesiveVertices.size();
  err = VecGetArray(jacobianVec, &jacobianArray);CHECK_PETSC_ERROR(err);
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_lagrange = _cohesiveVertices[iVertex].lagrange;
    PetscInt goff;

    // Compute contribution only if Lagrange constraint is local.
    err = PetscSectionGetOffset(jacobianGlobalSection, v_lagrange, &goff);CHECK_PETSC_ERROR(err);
    if (goff < 0) continue;

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(updateEvent);
#endif
    PetscInt dof, off;

    err = PetscSectionGetDof(jacobianSection, v_lagrange, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(jacobianSection, v_lagrange, &off);CHECK_PETSC_ERROR(err);
    assert(spaceDim == dof);

    for(PetscInt d = 0; d < dof; ++d) {
      jacobianArray[off+d] = 1.0;
    }

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(updateEvent);
#endif
  } // for
  err = VecRestoreArray(jacobianVec, &jacobianArray);CHECK_PETSC_ERROR(err);
  PetscLogFlops(0);

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
  const int restrictEvent = _logger->eventId("FaPr restrict");
  const int updateEvent = _logger->eventId("FaPr update");

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
  for (int i=0; i < spaceDim; ++i)
    indicesRel[i] = i;

  // Get sections
  DM           areaDM      = _fields->get("area").dmMesh();
  PetscSection areaSection = _fields->get("area").petscSection();
  Vec          areaVec     = _fields->get("area").localVector();
  PetscSection areaGlobalSection;
  PetscScalar *areaArray;
  PetscErrorCode err;
  assert(areaSection);assert(areaVec);
  err = DMGetDefaultGlobalSection(areaDM, &areaGlobalSection);CHECK_PETSC_ERROR(err);

  DM           solutionDM      = fields->solution().dmMesh();
  PetscSection solutionSection = fields->solution().petscSection();
  Vec          solutionVec     = fields->solution().localVector();
  PetscSection solutionGlobalSection;
  PetscScalar *solutionArray;
  assert(solutionSection);assert(solutionVec);
  err = DMGetDefaultGlobalSection(solutionDM, &solutionGlobalSection);CHECK_PETSC_ERROR(err);

  _logger->eventEnd(setupEvent);
#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventBegin(computeEvent);
#endif

  PetscMat jacobianNP;
  std::map<int, int> indicesMatToSubmat;
  _getJacobianSubmatrixNP(&jacobianNP, &indicesMatToSubmat, *jacobian, *fields);
  
  const int numVertices = _cohesiveVertices.size();
  err = VecGetArray(areaVec, &areaArray);CHECK_PETSC_ERROR(err);
  for (int iVertex=0, cV = 0; iVertex < numVertices; ++iVertex) {
    const int v_lagrange = _cohesiveVertices[iVertex].lagrange;
    const int v_fault = _cohesiveVertices[iVertex].fault;
    const int v_negative = _cohesiveVertices[iVertex].negative;
    const int v_positive = _cohesiveVertices[iVertex].positive;
    PetscInt goff;

    // Compute contribution only if Lagrange constraint is local.
    err = PetscSectionGetOffset(solutionGlobalSection, v_lagrange, &goff);CHECK_PETSC_ERROR(err);
    if (goff < 0) continue;

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(restrictEvent);
#endif

    // Get area associated with fault vertex.
    PetscInt adof, aoff;

    err = PetscSectionGetDof(areaSection, v_fault, &adof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(areaSection, v_fault, &aoff);CHECK_PETSC_ERROR(err);
    assert(1 == adof);
    PetscInt gnoff, gpoff;

    err = PetscSectionGetOffset(solutionGlobalSection, v_negative, &gnoff);CHECK_PETSC_ERROR(err);
    indicesN = indicesRel + indicesMatToSubmat[gnoff];
    err = MatGetValues(jacobianNP,
                       indicesN.size(), &indicesN[0], indicesN.size(), &indicesN[0],
                       &jacobianVertexN[0]); CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(solutionGlobalSection, v_positive, &gpoff);CHECK_PETSC_ERROR(err);
    indicesP = indicesRel + indicesMatToSubmat[gpoff];
    err = MatGetValues(jacobianNP,
                       indicesP.size(), &indicesP[0], indicesP.size(), &indicesP[0],
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

    // Set global preconditioner index associated with Lagrange constraint vertex.
    PetscInt indexLprecond;
    err = PetscSectionGetOffset(areaGlobalSection, v_lagrange, &indexLprecond);CHECK_PETSC_ERROR(err);
    
    // Set diagonal entries in preconditioned matrix.
    for (int iDim=0; iDim < spaceDim; ++iDim)
      MatSetValue(*precondMatrix,
		  indexLprecond + iDim,
		  indexLprecond + iDim,
		  precondVertexL[iDim],
		  INSERT_VALUES);

#if 0 // DEBUGGING
    std::cout << "1/P_vertex " << *v_lagrange << std::endl;
    for(int iDim = 0; iDim < spaceDim; ++iDim) {
      std::cout << "  " << precondVertexL[iDim] << std::endl;
    } // for
#endif

    
#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(updateEvent);
#endif
  } // for
  err = VecRestoreArray(areaVec, &areaArray);CHECK_PETSC_ERROR(err);
  err = MatDestroy(&jacobianNP);CHECK_PETSC_ERROR(err);
  PetscLogFlops(numVertices*spaceDim*6);

#if !defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(computeEvent);
#endif
} // calcPreconditioner

// ----------------------------------------------------------------------
// Adjust solution from solver with lumped Jacobian to match Lagrange
// multiplier constraints.
void
pylith::faults::FaultCohesiveLagrange::adjustSolnLumped(topology::SolutionFields* const fields,
							const PylithScalar t,
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
  PetscSection areaSection = _fields->get("area").petscSection();
  Vec          areaVec     = _fields->get("area").localVector();
  PetscScalar *areaArray;
  assert(areaSection);assert(areaVec);

  DM           jacobianDM      = jacobian.dmMesh();
  PetscSection jacobianSection = jacobian.petscSection();
  Vec          jacobianVec     = jacobian.localVector();
  PetscSection jacobianGlobalSection;
  PetscScalar *jacobianArray;
  PetscErrorCode err;
  assert(jacobianSection);assert(jacobianVec);
  err = DMGetDefaultGlobalSection(jacobianDM, &jacobianGlobalSection);CHECK_PETSC_ERROR(err);

  PetscSection residualSection = fields->get("residual").petscSection();
  Vec          residualVec     = fields->get("residual").localVector();
  PetscScalar *residualArray;
  assert(residualSection);assert(residualVec);

  PetscSection dispTIncrSection = fields->get("dispIncr(t->t+dt)").petscSection();
  Vec          dispTIncrVec     = fields->get("dispIncr(t->t+dt)").localVector();
  PetscScalar *dispTIncrArray;
  assert(dispTIncrSection);assert(dispTIncrVec);

  PetscSection dispTIncrAdjSection = fields->get("dispIncr adjust").petscSection();
  Vec          dispTIncrAdjVec     = fields->get("dispIncr adjust").localVector();
  PetscScalar *dispTIncrAdjArray;
  assert(dispTIncrAdjSection);assert(dispTIncrAdjVec);

  _logger->eventEnd(setupEvent);

#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventBegin(computeEvent);
#endif

  const int numVertices = _cohesiveVertices.size();
  err = VecGetArray(areaVec, &areaArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(jacobianVec, &jacobianArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(dispTIncrVec, &dispTIncrArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(residualVec, &residualArray);CHECK_PETSC_ERROR(err);
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_lagrange = _cohesiveVertices[iVertex].lagrange;
    const int v_fault = _cohesiveVertices[iVertex].fault;
    const int v_negative = _cohesiveVertices[iVertex].negative;
    const int v_positive = _cohesiveVertices[iVertex].positive;
    PetscInt goff;

    // Compute contribution only if Lagrange constraint is local.
    err = PetscSectionGetOffset(jacobianGlobalSection, v_lagrange, &goff);CHECK_PETSC_ERROR(err);
    if (goff < 0) continue;

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(restrictEvent);
#endif

    // Get residual at cohesive cell's vertices.
    PetscInt rldof, rloff;

    err = PetscSectionGetDof(residualSection, v_lagrange, &rldof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(residualSection, v_lagrange, &rloff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == rldof);

    // Get jacobian at cohesive cell's vertices.
    PetscInt jndof, jnoff, jpdof, jpoff;

    err = PetscSectionGetDof(jacobianSection, v_negative, &jndof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(jacobianSection, v_negative, &jnoff);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetDof(jacobianSection, v_positive, &jpdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(jacobianSection, v_positive, &jpoff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == jndof);assert(spaceDim == jpdof);

    // Get area at fault vertex.
    PetscInt adof, aoff;

    err = PetscSectionGetDof(areaSection, v_fault, &adof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(areaSection, v_fault, &aoff);CHECK_PETSC_ERROR(err);
    assert(1 == adof);
    // TODO assert(areaArray[aoff] > 0.0);

    // Get dispIncr(t) at vertices.
    PetscInt dtndof, dtnoff;

    err = PetscSectionGetDof(dispTIncrSection, v_negative, &dtndof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dispTIncrSection, v_negative, &dtnoff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == dtndof);
    PetscInt dtpdof, dtpoff;

    err = PetscSectionGetDof(dispTIncrSection, v_positive, &dtpdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dispTIncrSection, v_positive, &dtpoff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == dtpdof);
    PetscInt dtldof, dtloff;

    err = PetscSectionGetDof(dispTIncrSection, v_lagrange, &dtldof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dispTIncrSection, v_lagrange, &dtloff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == dtldof);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(restrictEvent);
    _logger->eventBegin(computeEvent);
#endif

    for(PetscInt d = 0; d < spaceDim; ++d) {
      const PylithScalar S = (1.0/jacobianArray[jpoff+d] + 1.0/jacobianArray[jnoff+d]) * areaArray[aoff] * areaArray[aoff];
      // Set Lagrange multiplier value (value from preliminary solve is bogus due to artificial diagonal entry)
      dispTIncrArray[dtloff+d] = 1.0/S * (-residualArray[rloff+d] + areaArray[aoff] * (dispTIncrArray[dtpoff+d] - dispTIncrArray[dtnoff+d]));

      // Adjust displacements to account for Lagrange multiplier values (assumed to be zero in preliminary solve).
      assert(jacobianArray[jnoff+d] > 0.0);
      dispTIncrArray[dtnoff+d] +=  areaArray[aoff] / jacobianArray[jnoff+d]*dispTIncrArray[dtloff+d];

      assert(jacobianArray[jpoff+d] > 0.0);
      dispTIncrArray[dtpoff+d] += -areaArray[aoff] / jacobianArray[jpoff+d]*dispTIncrArray[dtloff+d];
    } // for

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(computeEvent);
    _logger->eventBegin(updateEvent);
#endif

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(updateEvent);
#endif
  } // for
  err = VecRestoreArray(areaVec, &areaArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(jacobianVec, &jacobianArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(dispTIncrVec, &dispTIncrArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(residualVec, &residualArray);CHECK_PETSC_ERROR(err);

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

  DM             dmMesh = mesh.dmMesh();
  PetscBool      hasLabel;
  PetscInt       vStart, vEnd;
  PetscErrorCode err;

  assert(dmMesh);
  err = DMComplexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);
  err = DMComplexHasLabel(dmMesh, label(), &hasLabel);CHECK_PETSC_ERROR(err);
  if (!hasLabel) {
    std::ostringstream msg;
    msg << "Mesh missing group of vertices '" << label() << " for boundary condition.";
    throw std::runtime_error(msg.str());
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
  const int       numCorners = _quadrature->numBasis();
  IS              cellIS;
  const PetscInt *cells;
  PetscInt        numCells;
  err = DMComplexGetStratumIS(dmMesh, "material-id", id(), &cellIS);CHECK_PETSC_ERROR(err);
  assert(cellIS);
  err = ISGetLocalSize(cellIS, &numCells);CHECK_PETSC_ERROR(err);
  err = ISGetIndices(cellIS, &cells);CHECK_PETSC_ERROR(err);
  for(PetscInt c = 0; c < numCells; ++c) {
    PetscInt *closure        = PETSC_NULL;
    PetscInt  cellNumCorners = 0, closureSize;

    err = DMComplexGetTransitiveClosure(dmMesh, cells[c], PETSC_TRUE, &closureSize, &closure);CHECK_PETSC_ERROR(err);
    for(PetscInt p = 0; p < closureSize*2; p += 2) {
      const PetscInt point = closure[p];
      if ((point >= vStart) && (point < vEnd)) {
        ++cellNumCorners;
      }
    }
    if (3 * numCorners != cellNumCorners) {
      std::ostringstream msg;
      msg << "Number of vertices in reference cell (" << numCorners
          << ") is not compatible with number of vertices (" << cellNumCorners
          << ") in cohesive cell " << cells[c] << " for fault '" << label()
          << "'.";
      throw std::runtime_error(msg.str());
    } // if
    err = DMComplexRestoreTransitiveClosure(dmMesh, cells[c], PETSC_TRUE, &closureSize, &closure);CHECK_PETSC_ERROR(err);
  } // for
  err = ISRestoreIndices(cellIS, &cells);CHECK_PETSC_ERROR(err);
} // verifyConfiguration

// ----------------------------------------------------------------------
// Verify constraints are acceptable.
void
pylith::faults::FaultCohesiveLagrange::checkConstraints(const topology::Field<topology::Mesh>& solution) const
{ // checkConstraints Check to make sure no vertices connected to the
  // fault are constrained.

  const PetscInt spaceDim = solution.mesh().dimension();
  PetscSection   section  = solution.petscSection();
  PetscErrorCode err;
  assert(section);

  const int numVertices = _cohesiveVertices.size();
  for(int iVertex = 0; iVertex < numVertices; ++iVertex) {
    const int v_negative = _cohesiveVertices[iVertex].negative;    
    PetscInt fiberDimN, numConstraintsN;

    err = PetscSectionGetDof(section, v_negative, &fiberDimN);CHECK_PETSC_ERROR(err);
    assert(spaceDim == fiberDimN);
    err = PetscSectionGetConstraintDof(section, v_negative, &numConstraintsN);CHECK_PETSC_ERROR(err);
    if (numConstraintsN > 0) {
      std::ostringstream msg;
      msg << "Vertex with label '" << v_negative << "' on negative side "
	  << "of fault '" << label() << "' is constrained.\n"
	  << "Fault vertices cannot be constrained.";
      throw std::runtime_error(msg.str());
    } // if
    const int v_positive = _cohesiveVertices[iVertex].positive;
    PetscInt fiberDimP, numConstraintsP;

    err = PetscSectionGetDof(section, v_positive, &fiberDimP);CHECK_PETSC_ERROR(err);
    assert(spaceDim == fiberDimP);
    err = PetscSectionGetConstraintDof(section, v_positive, &numConstraintsP);CHECK_PETSC_ERROR(err);
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
  DM              dmMesh = mesh.dmMesh();
  IS              cellIS;
  const PetscInt *cells;
  PetscInt        numCells, vStart, vEnd;
  PetscErrorCode  err;

  err = DMComplexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);
  err = DMComplexGetStratumIS(dmMesh, "material-id", id(), &cellIS);CHECK_PETSC_ERROR(err);
  assert(cellIS);
  err = ISGetLocalSize(cellIS, &numCells);CHECK_PETSC_ERROR(err);
  const int numConstraintVert = _quadrature->numBasis();
  const int numCorners = 3 * numConstraintVert; // cohesive cell

  DM              faultDMMesh = _faultMesh->dmMesh();
  IS              subpointMap;
  const PetscInt *points;
  PetscInt        numPoints, fcStart, fcEnd, fvStart, fvEnd;

  assert(faultDMMesh);
  err = DMComplexGetHeightStratum(faultDMMesh, 0, &fcStart, &fcEnd);CHECK_PETSC_ERROR(err);
  err = DMComplexGetDepthStratum(faultDMMesh, 0, &fvStart, &fvEnd);CHECK_PETSC_ERROR(err);
  err = DMComplexGetSubpointMap(faultDMMesh, &subpointMap);CHECK_PETSC_ERROR(err);
  err = ISGetLocalSize(subpointMap, &numPoints);CHECK_PETSC_ERROR(err);
  assert(numCells == fcEnd-fcStart);

  _cohesiveToFault.clear();
  typedef std::map<int, int> indexmap_type;
  indexmap_type indexMap;
  _cohesiveVertices.resize(fvEnd-fvStart);
  PetscInt index = 0;

  err = ISGetIndices(cellIS, &cells);CHECK_PETSC_ERROR(err);
  err = ISGetIndices(subpointMap, &points);CHECK_PETSC_ERROR(err);
  for(PetscInt c = 0; c < numCells; ++c) {
    _cohesiveToFault[cells[c]] = c;

    // Get oriented closure
    PetscInt *closure = PETSC_NULL;
    PetscInt  closureSize, q = 0;

    err = DMComplexGetTransitiveClosure(dmMesh, cells[c], PETSC_TRUE, &closureSize, &closure);CHECK_PETSC_ERROR(err);
    for(PetscInt p = 0; p < closureSize*2; p += 2) {
      const PetscInt point = closure[p];
      if ((point >= vStart) && (point < vEnd)) {
        closure[q++] = point;
      }
    }
    closureSize = q;
    assert(closureSize == numCorners);

    for (int iConstraint = 0; iConstraint < numConstraintVert; ++iConstraint) {
      // normal cohesive vertices i and j and constraint vertex k
      const int indexI = iConstraint;
      const int indexJ = iConstraint + numConstraintVert;
      const int indexK = iConstraint + 2 * numConstraintVert;

      const int v_lagrange = closure[indexK];
      const int v_negative = closure[indexI];
      const int v_positive = closure[indexJ];
      PetscInt  v_fault;

      err = PetscFindInt(v_lagrange, numPoints, points, &v_fault);CHECK_PETSC_ERROR(err);
      assert(v_fault >= 0);
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
    err = DMComplexRestoreTransitiveClosure(dmMesh, cells[c], PETSC_TRUE, &closureSize, &closure);CHECK_PETSC_ERROR(err);
  } // for
  err = ISRestoreIndices(cellIS, &cells);CHECK_PETSC_ERROR(err);
  err = ISRestoreIndices(subpointMap, &points);CHECK_PETSC_ERROR(err);
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
// Transform field from local (fault) coordinate system to
// global coordinate system.
void
pylith::faults::FaultCohesiveLagrange::faultToGlobal(topology::Field<topology::SubMesh>* field,
						     const topology::Field<topology::SubMesh>& faultOrientation)
{ // faultToGlobal
  assert(field);

  // Fiber dimension of vector field matches spatial dimension.
  const spatialdata::geocoords::CoordSys* cs = field->mesh().coordsys();
  assert(cs);
  const int spaceDim = cs->spaceDim();
  scalar_array fieldVertexGlobal(spaceDim);

  // Get sections.
  PetscSection fieldSection = field->petscSection();
  Vec          fieldVec     = field->localVector();
  PetscScalar *fieldArray;
  assert(fieldSection);assert(fieldVec);

  PetscSection orientationSection = faultOrientation.petscSection();
  Vec          orientationVec     = faultOrientation.localVector();
  PetscScalar *orientationArray;
  assert(orientationSection);assert(orientationVec);

  DM             dmMesh = field->mesh().dmMesh();
  PetscInt       vStart, vEnd;
  PetscErrorCode err;

  assert(dmMesh);
  err = DMComplexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);

  err = VecGetArray(fieldVec, &fieldArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(orientationVec, &orientationArray);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt dof, off, odof, ooff;

    err = PetscSectionGetDof(fieldSection, v, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(fieldSection, v, &off);CHECK_PETSC_ERROR(err);
    assert(spaceDim == dof);
    err = PetscSectionGetDof(orientationSection, v, &odof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(orientationSection, v, &ooff);CHECK_PETSC_ERROR(err);
    assert(spaceDim*spaceDim == odof);

    // Rotate from fault to global coordinate system (transpose orientation)
    fieldVertexGlobal = 0.0;
    for(PetscInt d = 0; d < spaceDim; ++d)
      for(PetscInt e = 0; e < spaceDim; ++e)
        fieldVertexGlobal[d] += orientationArray[ooff+e*spaceDim+d] * fieldArray[off+e];
    assert(fieldVertexGlobal.size() == dof);
    for(PetscInt d = 0; d < spaceDim; ++d)
      fieldArray[off+d] = fieldVertexGlobal[d];
  } // for
  PetscLogFlops((vEnd-vStart) * (2*spaceDim*spaceDim) );
  err = VecRestoreArray(fieldVec, &fieldArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(orientationVec, &orientationArray);CHECK_PETSC_ERROR(err);
  
#if 0 // DEBUGGING
  field->view("FIELD (GLOBAL)");
#endif
} // faultToGlobal

// ----------------------------------------------------------------------
// Transform field from global coordinate system to local (fault)
// coordinate system.
void
pylith::faults::FaultCohesiveLagrange::globalToFault(topology::Field<topology::SubMesh>* field,
						     const topology::Field<topology::SubMesh>& faultOrientation)
{ // globalToFault
  assert(field);

  // Fiber dimension of vector field matches spatial dimension.
  const spatialdata::geocoords::CoordSys* cs = field->mesh().coordsys();
  assert(cs);
  const int spaceDim = cs->spaceDim();
  scalar_array fieldVertexFault(spaceDim);

  // Get sections.
  PetscSection fieldSection = field->petscSection();
  Vec          fieldVec     = field->localVector();
  PetscScalar *fieldArray;
  assert(fieldSection);assert(fieldVec);

  PetscSection orientationSection = faultOrientation.petscSection();
  Vec          orientationVec     = faultOrientation.localVector();
  PetscScalar *orientationArray;
  assert(orientationSection);assert(orientationVec);

  DM             dmMesh = field->mesh().dmMesh();
  PetscInt       vStart, vEnd;
  PetscErrorCode err;

  assert(dmMesh);
  err = DMComplexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);

  err = VecGetArray(fieldVec, &fieldArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(orientationVec, &orientationArray);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt dof, off, odof, ooff;

    err = PetscSectionGetDof(fieldSection, v, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(fieldSection, v, &off);CHECK_PETSC_ERROR(err);
    assert(spaceDim == dof);
    err = PetscSectionGetDof(orientationSection, v, &odof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(orientationSection, v, &ooff);CHECK_PETSC_ERROR(err);
    assert(spaceDim*spaceDim == odof);

    // Rotate from global coordinate system to fault (orientation)
    fieldVertexFault = 0.0;
    for(PetscInt d = 0; d < spaceDim; ++d)
      for(PetscInt e = 0; e < spaceDim; ++e)
        fieldVertexFault[d] += orientationArray[ooff+d*spaceDim+e] * fieldArray[off+e];
    assert(fieldVertexFault.size() == dof);
    for(PetscInt d = 0; d < spaceDim; ++d)
      fieldArray[off+d] = fieldVertexFault[d];
  } // for
  PetscLogFlops((vEnd-vStart) * (2*spaceDim*spaceDim) );
  err = VecRestoreArray(fieldVec, &fieldArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(orientationVec, &orientationArray);CHECK_PETSC_ERROR(err);

#if 0 // DEBUGGING
  field->view("FIELD (FAULT)");
#endif
} // faultToGlobal

// ----------------------------------------------------------------------
// Calculate orientation at fault vertices.
void
pylith::faults::FaultCohesiveLagrange::_calcOrientation(const PylithScalar upDir[3])
{ // _calcOrientation
  assert(0 != upDir);
  assert(0 != _faultMesh);
  assert(0 != _fields);

  scalar_array up(3);
  for (int i=0; i < 3; ++i)
    up[i] = upDir[i];

  // Get vertices in fault mesh.
  DM             faultDMMesh = _faultMesh->dmMesh();
  PetscInt       vStart, vEnd, cStart, cEnd;
  PetscErrorCode err;

  assert(faultDMMesh);
  err = DMComplexGetDepthStratum(faultDMMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);
  err = DMComplexGetHeightStratum(faultDMMesh, 0, &cStart, &cEnd);CHECK_PETSC_ERROR(err);

  // Containers for orientation information.
  const int cohesiveDim = _faultMesh->dimension();
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int orientationSize = spaceDim * spaceDim;
  const feassemble::CellGeometry& cellGeometry = _quadrature->refGeometry();
  const scalar_array& verticesRef = cellGeometry.vertices();
  const int jacobianSize = (cohesiveDim > 0) ? spaceDim * cohesiveDim : 1;
  const scalar_array& quadWts = _quadrature->quadWts();
  scalar_array jacobian(jacobianSize);
  PylithScalar jacobianDet = 0;
  scalar_array refCoordsVertex(cohesiveDim);

  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("FaultFields");

  // Allocate orientation field.
  scalar_array orientationVertex(orientationSize);
  _fields->add("orientation", "orientation");
  topology::Field<topology::SubMesh>& orientation   = _fields->get("orientation");
  const topology::Field<topology::SubMesh>& dispRel = _fields->get("relative disp");
  orientation.addField("orientation", cohesiveDim+1);
  orientation.setupFields();
  orientation.newSection(dispRel, orientationSize);
  // Create components for along-strike, up-dip, and normal directions
  orientation.updateDof("orientation", pylith::topology::FieldBase::VERTICES_FIELD, spaceDim);
  orientation.allocate();
  orientation.zero();
  PetscSection orientationSection = orientation.petscSection();
  Vec          orientationVec     = orientation.localVector();
  PetscScalar *orientationArray;

  logger.stagePop();

  // Compute orientation of fault at constraint vertices

  // Get section containing coordinates of vertices
  PetscSection coordSection;
  Vec          coordVec;
  err = DMComplexGetCoordinateSection(faultDMMesh, &coordSection);CHECK_PETSC_ERROR(err);
  err = DMGetCoordinatesLocal(faultDMMesh, &coordVec);CHECK_PETSC_ERROR(err);

  // Set orientation function
  assert(cohesiveDim == _quadrature->cellDim());
  assert(spaceDim == _quadrature->spaceDim());

  // Loop over cohesive cells, computing orientation weighted by
  // jacobian at constraint vertices

  err = VecGetArray(orientationVec, &orientationArray);CHECK_PETSC_ERROR(err);
  for(PetscInt c = cStart; c < cEnd; ++c) {
    PetscInt          *closure = PETSC_NULL;
    PetscInt           closureSize, q = 0;
    const PetscScalar *coords = PETSC_NULL;
    PetscInt           coordsSize;
    scalar_array       coordinatesCell(numBasis * spaceDim);

    // Get orientations at fault cell's vertices.
    err = DMComplexVecGetClosure(faultDMMesh, coordSection, coordVec, c, &coordsSize, &coords);CHECK_PETSC_ERROR(err);
    err = DMComplexGetTransitiveClosure(faultDMMesh, c, PETSC_TRUE, &closureSize, &closure);CHECK_PETSC_ERROR(err);
    for(PetscInt i = 0; i < coordsSize; ++i) {coordinatesCell[i] = coords[i];}

    // Filter out non-vertices
    for(PetscInt p = 0; p < closureSize*2; p += 2) {
      if ((closure[p] >= vStart) && (closure[p] < vEnd)) {
        closure[q*2]   = closure[p];
        closure[q*2+1] = closure[p+1];
        ++q;
      }
    }
    closureSize = q;
    for(PetscInt v = 0; v < closureSize; ++v) {
      // Compute Jacobian and determinant of Jacobian at vertex
      memcpy(&refCoordsVertex[0], &verticesRef[v * cohesiveDim], cohesiveDim * sizeof(PylithScalar));
      cellGeometry.jacobian(&jacobian, &jacobianDet, coordinatesCell, refCoordsVertex);

      // Compute orientation
      cellGeometry.orientation(&orientationVertex, jacobian, jacobianDet, up);

      // Update orientation
      PetscInt off;

      err = PetscSectionGetOffset(orientationSection, closure[v*2], &off);CHECK_PETSC_ERROR(err);
      for(PetscInt d = 0; d < orientationSize; ++d) {
        orientationArray[off+d] += orientationVertex[d];
      }
    } // for
    err = DMComplexVecRestoreClosure(faultDMMesh, coordSection, coordVec, c, &coordsSize, &coords);CHECK_PETSC_ERROR(err);
    err = DMComplexRestoreTransitiveClosure(faultDMMesh, c, PETSC_TRUE, &closureSize, &closure);CHECK_PETSC_ERROR(err);
  } // for
  err = VecRestoreArray(orientationVec, &orientationArray);CHECK_PETSC_ERROR(err);

  //orientation.view("ORIENTATION BEFORE COMPLETE");

  // Assemble orientation information
  orientation.complete();

  // Loop over vertices, make orientation information unit magnitude
  scalar_array vertexDir(orientationSize);
  int count = 0;
  err = VecGetArray(orientationVec, &orientationArray);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v, ++count) {
    PetscInt off;

    err = PetscSectionGetOffset(orientationSection, v, &off);CHECK_PETSC_ERROR(err);
    for(PetscInt d = 0; d < orientationSize; ++d) {
      orientationVertex[d] = orientationArray[off+d];
    }
    for (int iDim = 0; iDim < spaceDim; ++iDim) {
      PylithScalar mag = 0;
      for (int jDim = 0, index = iDim * spaceDim; jDim < spaceDim; ++jDim)
        mag += pow(orientationVertex[index + jDim], 2);
      mag = sqrt(mag);
      assert(mag > 0.0);
      for (int jDim = 0, index = iDim * spaceDim; jDim < spaceDim; ++jDim)
        orientationVertex[index + jDim] /= mag;
    } // for
    for(PetscInt d = 0; d < orientationSize; ++d) {
      orientationArray[off+d] = orientationVertex[d];
    }
  } // for
  err = VecRestoreArray(orientationVec, &orientationArray);CHECK_PETSC_ERROR(err);
  PetscLogFlops(count * orientationSize * 4);

  err = VecGetArray(orientationVec, &orientationArray);CHECK_PETSC_ERROR(err);
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
    PetscInt off;
    int flipLocal = 0;
    if (vEnd > vStart) {
    err = PetscSectionGetOffset(orientationSection, vStart, &off);CHECK_PETSC_ERROR(err);
    for (PetscInt d = 0; d < orientationSize; ++d) {
      orientationVertex[d] = orientationArray[off+d];
    } // for

      assert(2 == spaceDim);
      const PylithScalar* shearDirVertex = &orientationVertex[0];
      const PylithScalar* normalDirVertex = &orientationVertex[2];
      const PylithScalar shearDirDot = 
	upDir[0] * shearDirVertex[0] + upDir[1] * shearDirVertex[1];
      const PylithScalar normalDirDot = 
	upDir[0] * normalDirVertex[0] + upDir[1] * normalDirVertex[1];
      if (normalDirDot * shearDirDot < 0.0) {
	flipLocal = 1;
      } // if
    } // if
    // Collect flip decisions across all processors
    int flipGlobal = 0;
   MPI_Allreduce(&flipLocal, &flipGlobal, 1, MPI_INT, MPI_SUM, ((PetscObject) faultDMMesh)->comm);

    const int ishear = 0;
    const int inormal = 2;
    if (flipGlobal > 0) { // flip in any processor wants to flip
      // Flip shear direction
      for(PetscInt v = vStart; v < vEnd; ++v) {
        PetscInt dof;

        err = PetscSectionGetDof(orientationSection, v, &dof);CHECK_PETSC_ERROR(err);
        err = PetscSectionGetOffset(orientationSection, v, &off);CHECK_PETSC_ERROR(err);
        for(PetscInt d = 0; d < orientationSize; ++d) {
          orientationVertex[d] = orientationArray[off+d];
        }
        assert(4 == dof);
        for (int iDim = 0; iDim < 2; ++iDim) // flip shear
          orientationVertex[ishear + iDim] *= -1.0;
	
        // Update orientation
        for(PetscInt d = 0; d < orientationSize; ++d) {
          orientationArray[off+d] = orientationVertex[d];
        }
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
    PetscInt off;

    int flipLocal = 0;
    if (vEnd > vStart) {
      err = PetscSectionGetOffset(orientationSection, vStart, &off);CHECK_PETSC_ERROR(err);
      for (PetscInt d = 0; d < orientationSize; ++d) {
	orientationVertex[d] = orientationArray[off+d];
      } // for
      
      assert(3 == spaceDim);
      const PylithScalar* dipDirVertex = &orientationVertex[3];
      const PylithScalar* normalDirVertex = &orientationVertex[6];
      const PylithScalar dipDirDot = 
	upDir[0]*dipDirVertex[0] + 
	upDir[1]*dipDirVertex[1] + 
	upDir[2]*dipDirVertex[2];
      const PylithScalar normalDirDot = 
	upDir[0]*normalDirVertex[0] +
	upDir[1]*normalDirVertex[1] +
	upDir[2]*normalDirVertex[2];

      if (dipDirDot * normalDirDot < 0.0 || fabs(normalDirVertex[2] + 1.0) < 0.001) {
	// if fault normal is (0,0,+-1) then up-dir dictates reverse
	// motion for case with normal (0,0,1), so we reverse the dip-dir
	// if we have (0,0,-1).
	flipLocal = 1;
      } // if
    } // if

    // Collect flip decisions across all processors
    int flipGlobal = 0;
    MPI_Allreduce(&flipLocal, &flipGlobal, 1, MPI_INT, MPI_SUM, ((PetscObject) faultDMMesh)->comm);

    const int istrike = 0;
    const int idip = 3;
    const int inormal = 6;
    if (flipGlobal > 0) { // flip in any processor wants to flip
      // Flip dip direction
      for(PetscInt v = vStart; v < vEnd; ++v) {
        PetscInt dof;

        err = PetscSectionGetDof(orientationSection, v, &dof);CHECK_PETSC_ERROR(err);
        err = PetscSectionGetOffset(orientationSection, v, &off);CHECK_PETSC_ERROR(err);
        for(PetscInt d = 0; d < orientationSize; ++d) {
          orientationVertex[d] = orientationArray[off+d];
        }
        assert(9 == dof);
        for (int iDim = 0; iDim < 3; ++iDim) // flip dip
          orientationVertex[idip + iDim] *= -1.0;
	
        // Update direction
        for(PetscInt d = 0; d < orientationSize; ++d) {
          orientationArray[off+d] = orientationVertex[d];
        }
      } // for
      PetscLogFlops(5 + count * 3);
    } // if
  } // if
  err = VecRestoreArray(orientationVec, &orientationArray);CHECK_PETSC_ERROR(err);

  //orientation.view("ORIENTATION");
} // _calcOrientation

// ----------------------------------------------------------------------
void
pylith::faults::FaultCohesiveLagrange::_calcArea(void)
{ // _calcArea
  assert(_faultMesh);
  assert(_fields);
  assert(_normalizer);

  // Containers for area information
  const int cellDim = _quadrature->cellDim();
  const int numBasis = _quadrature->numBasis();
  const int numQuadPts = _quadrature->numQuadPts();
  const int spaceDim = _quadrature->spaceDim();
  const feassemble::CellGeometry& cellGeometry = _quadrature->refGeometry();
  const scalar_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == numQuadPts);
  PylithScalar jacobianDet = 0;
  scalar_array areaCell(numBasis);

  // Get vertices in fault mesh.
  DM             faultDMMesh = _faultMesh->dmMesh();
  PetscInt       vStart, vEnd, cStart, cEnd;
  PetscErrorCode err;

  assert(faultDMMesh);
  err = DMComplexGetDepthStratum(faultDMMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);
  err = DMComplexGetHeightStratum(faultDMMesh, 0, &cStart, &cEnd);CHECK_PETSC_ERROR(err);

  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("FaultFields");

  // Allocate area field.
  _fields->add("area", "area");
  topology::Field<topology::SubMesh>&       area    = _fields->get("area");
  const topology::Field<topology::SubMesh>& dispRel = _fields->get("relative disp");
  area.newSection(dispRel, 1);
  area.allocate();
  area.vectorFieldType(topology::FieldBase::SCALAR);
  const PylithScalar lengthScale = _normalizer->lengthScale();
  area.scale(pow(lengthScale, (spaceDim-1)));
  area.zero();
  PetscSection areaSection = area.petscSection();
  Vec          areaVec     = area.localVector();
  PetscScalar *areaArray;
  assert(areaSection);assert(areaVec);

  logger.stagePop();

  scalar_array coordinatesCell(numBasis*spaceDim);
  PetscSection coordSection;
  Vec          coordVec;
  err = DMComplexGetCoordinateSection(faultDMMesh, &coordSection);CHECK_PETSC_ERROR(err);
  err = DMGetCoordinatesLocal(faultDMMesh, &coordVec);CHECK_PETSC_ERROR(err);

  // Loop over cells in fault mesh, compute area
  for(PetscInt c = cStart; c < cEnd; ++c) {
    areaCell = 0.0;

    // Compute geometry information for current cell
#if defined(PRECOMPUTE_GEOMETRY)
    _quadrature->retrieveGeometry(c);
#else
    const PetscScalar *coords = PETSC_NULL;
    PetscInt           coordsSize;
    err = DMComplexVecGetClosure(faultDMMesh, coordSection, coordVec, c, &coordsSize, &coords);CHECK_PETSC_ERROR(err);
    for(PetscInt i = 0; i < coordsSize; ++i) {coordinatesCell[i] = coords[i];}
    _quadrature->computeGeometry(coordinatesCell, c);
    err = DMComplexVecRestoreClosure(faultDMMesh, coordSection, coordVec, c, &coordsSize, &coords);CHECK_PETSC_ERROR(err);
#endif

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
    err = DMComplexVecSetClosure(faultDMMesh, areaSection, areaVec, c, &areaCell[0], ADD_VALUES);CHECK_PETSC_ERROR(err);
    PetscLogFlops(numQuadPts*(1+numBasis*2));
  } // for

  // Assemble area information
  area.complete();

#if 0 // DEBUGGING
  area.view("AREA");
  faultSieveMesh->getSendOverlap()->view("Send fault overlap");
  faultSieveMesh->getRecvOverlap()->view("Receive fault overlap");
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

  // Get sections.
  PetscSection dispTSection = dispT.petscSection();
  Vec          dispTVec     = dispT.localVector();
  PetscScalar *dispTArray;
  PetscErrorCode err;
  assert(dispTSection);assert(dispTVec);

  PetscSection orientationSection = _fields->get("orientation").petscSection();
  Vec          orientationVec     = _fields->get("orientation").localVector();
  PetscScalar *orientationArray;
  assert(orientationSection);assert(orientationVec);

  // Allocate buffer for tractions field (if necessary).
  PetscSection tractionsSection = tractions->petscSection();
  if (!tractionsSection) {
    ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
    logger.stagePush("FaultFields");

    const topology::Field<topology::SubMesh>& dispRel = _fields->get("relative disp");
    tractions->cloneSection(dispRel);
    tractionsSection = tractions->petscSection();
    logger.stagePop();
  } // if
  Vec          tractionsVec = tractions->localVector();
  PetscScalar *tractionsArray;
  assert(tractionsSection);assert(tractionsVec);
  tractions->zero();

  const PylithScalar pressureScale = tractions->scale();

  const int numVertices = _cohesiveVertices.size();
  err = VecGetArray(dispTVec, &dispTArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(tractionsVec, &tractionsArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(orientationVec, &orientationArray);CHECK_PETSC_ERROR(err);
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_lagrange = _cohesiveVertices[iVertex].lagrange;
    const int v_fault = _cohesiveVertices[iVertex].fault;
    PetscInt dof, off, odof, ooff, tdof, toff;

    err = PetscSectionGetDof(dispTSection, v_lagrange, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dispTSection, v_lagrange, &off);CHECK_PETSC_ERROR(err);
    assert(spaceDim == dof);
    err = PetscSectionGetDof(orientationSection, v_fault, &odof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(orientationSection, v_fault, &ooff);CHECK_PETSC_ERROR(err);
    assert(spaceDim*spaceDim == odof);
    err = PetscSectionGetDof(tractionsSection, v_fault, &tdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(tractionsSection, v_fault, &toff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == tdof);

    // Rotate from global coordinate system to fault (orientation)
    for(PetscInt iDim = 0; iDim < spaceDim; ++iDim) {
      tractionsArray[toff+iDim] = 0.0;
      for(PetscInt jDim = 0; jDim < spaceDim; ++jDim)
        tractionsArray[toff+iDim] += orientationArray[ooff+iDim*spaceDim+jDim] * dispTArray[off+jDim];
    }
  } // for
  err = VecRestoreArray(tractionsVec, &tractionsArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(orientationVec, &orientationArray);CHECK_PETSC_ERROR(err);

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
  logger.stagePush("OutputFields");

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
  logger.stagePush("OutputFields");

  // Create vector field; use same shape/chart as area field.
  assert(0 != _faultMesh);
  _fields->add("buffer (scalar)", "buffer");
  topology::Field<topology::SubMesh>& buffer = _fields->get("buffer (scalar)");
  buffer.newSection(topology::FieldBase::VERTICES_FIELD, 1);
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
  DM           solutionDM      = fields.solution().dmMesh();
  PetscSection solutionSection = fields.solution().petscSection();
  Vec          solutionVec     = fields.solution().localVector();
  PetscSection solutionGlobalSection;
  PetscScalar *solutionArray;
  PetscErrorCode err;
  assert(solutionSection);assert(solutionVec);
  err = DMGetDefaultGlobalSection(solutionDM, &solutionGlobalSection);CHECK_PETSC_ERROR(err);

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
    PetscInt goff;

    // Compute contribution only if Lagrange constraint is local.
    err = PetscSectionGetOffset(solutionGlobalSection, v_lagrange, &goff);CHECK_PETSC_ERROR(err);
    if (goff < 0) continue;
    numIndicesNP += 2;
  } // for
  int_array indicesNP(numIndicesNP*spaceDim);

  for (int iVertex=0, indexNP=0; iVertex < numVertices; ++iVertex) {
    const int v_lagrange = _cohesiveVertices[iVertex].lagrange;
    const int v_negative = _cohesiveVertices[iVertex].negative;
    const int v_positive = _cohesiveVertices[iVertex].positive;
    PetscInt gloff, gnoff, gpoff;

    // Compute contribution only if Lagrange constraint is local.
    err = PetscSectionGetOffset(solutionGlobalSection, v_lagrange, &gloff);CHECK_PETSC_ERROR(err);
    if (gloff < 0) continue;
    err = PetscSectionGetOffset(solutionGlobalSection, v_negative, &gnoff);CHECK_PETSC_ERROR(err);
    gnoff = gnoff < 0 ? -(gnoff+1) : gnoff;
    err = PetscSectionGetOffset(solutionGlobalSection, v_positive, &gpoff);CHECK_PETSC_ERROR(err);
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
// Get cell field associated with integrator.
const pylith::topology::Field<pylith::topology::SubMesh>&
pylith::faults::FaultCohesiveLagrange::cellField(const char* name,
                                                 const topology::SolutionFields* fields)
{ // cellField
  if (0 == strcasecmp("partition", name)) {

    DM             faultDMMesh = _faultMesh->dmMesh();
    PetscInt       cStart, cEnd;
    PetscErrorCode err;

    assert(faultDMMesh);
    err = DMComplexGetHeightStratum(faultDMMesh, 0, &cStart, &cEnd);CHECK_PETSC_ERROR(err);

    ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
    logger.stagePush("OutputFields");

    const int fiberDim = 1;
    _fields->add("partition", "partition", pylith::topology::FieldBase::CELLS_FIELD, fiberDim);
    topology::Field<topology::SubMesh>& partition = _fields->get("partition");
    partition.allocate();
    PetscSection partitionSection = partition.petscSection();
    Vec          partitionVec     = partition.localVector();
    PetscScalar *partitionArray;
    assert(partitionSection);assert(partitionVec);

    MPI_Comm    comm;
    PetscMPIInt rank;
    err = MPI_Comm_rank(((PetscObject) faultDMMesh)->comm, &rank);CHECK_PETSC_ERROR(err);
    // Loop over cells in fault mesh, set partition
    err = VecGetArray(partitionVec, &partitionArray);CHECK_PETSC_ERROR(err);
    for(PetscInt c = cStart; c < cEnd; ++c) {
      PetscInt dof, off;

      err = PetscSectionGetDof(partitionSection, c, &dof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(partitionSection, c, &off);CHECK_PETSC_ERROR(err);
      assert(dof == 1);
      partitionArray[off] = rank;
    } // for
    err = VecRestoreArray(partitionVec, &partitionArray);CHECK_PETSC_ERROR(err);

    logger.stagePop();

    return partition;    

  } // if

  // Should not reach this point if requested field was found
  std::ostringstream msg;
  msg << "Request for unknown cell field '" << name << "' for fault '" << label() << ".";
  throw std::runtime_error(msg.str());

  // Satisfy return values
  assert(0 != _fields);
  const topology::Field<topology::SubMesh>& buffer = _fields->get("buffer (vector)");
  return buffer;
} // cellField


// End of file 
