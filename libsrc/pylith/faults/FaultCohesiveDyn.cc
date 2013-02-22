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

#include "FaultCohesiveDyn.hh" // implementation of object methods

#include "CohesiveTopology.hh" // USES CohesiveTopology
#include "TractPerturbation.hh" // HOLDSA TractPerturbation

#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/feassemble/CellGeometry.hh" // USES CellGeometry
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/SubMesh.hh" // USES SubMesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/Jacobian.hh" // USES Jacobian
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields
#include "pylith/friction/FrictionModel.hh" // USES FrictionModel
#include "pylith/utils/macrodefs.h" // USES CALL_MEMBER_FN
#include "pylith/problems/SolverLinear.hh" // USES SolverLinear

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

// Precomputing geometry significantly increases storage but gives a
// slight speed improvement.
//#define PRECOMPUTE_GEOMETRY

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::SubMesh::SieveMesh SieveSubMesh;

// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::FaultCohesiveDyn::FaultCohesiveDyn(void) :
  _zeroTolerance(1.0e-10),
  _openFreeSurf(true),
  _tractPerturbation(0),
  _friction(0),
  _jacobian(0),
  _ksp(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::faults::FaultCohesiveDyn::~FaultCohesiveDyn(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void pylith::faults::FaultCohesiveDyn::deallocate(void)
{ // deallocate
  FaultCohesiveLagrange::deallocate();

  _tractPerturbation = 0; // :TODO: Use shared pointer
  _friction = 0; // :TODO: Use shared pointer

  delete _jacobian; _jacobian = 0;
  PetscErrorCode err = KSPDestroy(&_ksp);CHECK_PETSC_ERROR(err);
} // deallocate

// ----------------------------------------------------------------------
// Sets the spatial database for the inital tractions
void
pylith::faults::FaultCohesiveDyn::tractPerturbation(TractPerturbation* tract)
{ // tractPerturbation
  _tractPerturbation = tract;
} // tractPerturbation

// ----------------------------------------------------------------------
// Get the friction (constitutive) model.  
void
pylith::faults::FaultCohesiveDyn::frictionModel(friction::FrictionModel* const model)
{ // frictionModel
  _friction = model;
} // frictionModel

// ----------------------------------------------------------------------
// Nondimensional tolerance for detecting near zero values.
void
pylith::faults::FaultCohesiveDyn::zeroTolerance(const PylithScalar value)
{ // zeroTolerance
  if (value < 0.0) {
    std::ostringstream msg;
    msg << "Tolerance (" << value << ") for detecting values near zero for "
      "fault " << label() << " must be nonnegative.";
    throw std::runtime_error(msg.str());
  } // if
  
  _zeroTolerance = value;
} // zeroTolerance

// ----------------------------------------------------------------------
// Set flag used to determine when fault is traction free when it
// opens or it still imposes any initial tractions.
void
pylith::faults::FaultCohesiveDyn::openFreeSurf(const bool value)
{ // openFreeSurf
  _openFreeSurf = value;
} // openFreeSurf

// ----------------------------------------------------------------------
// Initialize fault. Determine orientation and setup boundary
void
pylith::faults::FaultCohesiveDyn::initialize(const topology::Mesh& mesh,
					     const PylithScalar upDir[3])
{ // initialize
  assert(upDir);
  assert(_quadrature);
  assert(_normalizer);

  FaultCohesiveLagrange::initialize(mesh, upDir);

  // Get initial tractions using a spatial database.
  if (_tractPerturbation) {
    const topology::Field<topology::SubMesh>& orientation = _fields->get("orientation");
    _tractPerturbation->initialize(*_faultMesh, orientation, *_normalizer);
  } // if

  // Setup fault constitutive model.
  assert(_friction);
  assert(_faultMesh);
  assert(_fields);
  _friction->normalizer(*_normalizer);
  _friction->initialize(*_faultMesh, _quadrature);

  const spatialdata::geocoords::CoordSys* cs = mesh.coordsys();
  assert(cs);

  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("FaultFields");

  // Create field for relative velocity associated with Lagrange vertex k
  _fields->add("relative velocity", "relative_velocity");
  topology::Field<topology::SubMesh>& velRel = 
    _fields->get("relative velocity");
  topology::Field<topology::SubMesh>& dispRel = _fields->get("relative disp");
  velRel.cloneSection(dispRel);
  velRel.vectorFieldType(topology::FieldBase::VECTOR);
  velRel.scale(_normalizer->lengthScale() / _normalizer->timeScale());

  logger.stagePop();
} // initialize

// ----------------------------------------------------------------------
// Integrate contributions to residual term (r) for operator.
void
pylith::faults::FaultCohesiveDyn::integrateResidual(
			     const topology::Field<topology::Mesh>& residual,
			     const PylithScalar t,
			     topology::SolutionFields* const fields)
{ // integrateResidual
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

  scalar_array tractPerturbVertex(spaceDim);
  PetscSection valuesSection;
  Vec          valuesVec;
  PetscScalar *tractionsArray;
  if (_tractPerturbation) {
    _tractPerturbation->calculate(t);
    
    const topology::Fields<topology::Field<topology::SubMesh> >* params = _tractPerturbation->parameterFields();
    assert(params);
    valuesSection = params->get("value").petscSection();
    valuesVec     = params->get("value").localVector();
    assert(valuesSection);assert(valuesVec);
  } // if

  PetscSection areaSection = _fields->get("area").petscSection();
  Vec          areaVec     = _fields->get("area").localVector();
  PetscScalar *areaArray;
  assert(areaSection);assert(areaVec);

  PetscSection orientationSection = _fields->get("orientation").petscSection();
  Vec          orientationVec     = _fields->get("orientation").localVector();
  PetscScalar *orientationArray;
  assert(orientationSection);assert(orientationVec);

  _logger->eventEnd(setupEvent);
#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventBegin(computeEvent);
#endif

  // Loop over fault vertices
  const int numVertices = _cohesiveVertices.size();
  err = VecGetArray(areaVec, &areaArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(orientationVec, &orientationArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(dispTVec, &dispTArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(dispTIncrVec, &dispTIncrArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(residualVec, &residualArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(valuesVec, &tractionsArray);CHECK_PETSC_ERROR(err);
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

    // Get prescribed traction perturbation at fault vertex.
    if (_tractPerturbation) {
      PetscInt vdof, voff;

      err = PetscSectionGetDof(valuesSection, v_fault, &vdof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(valuesSection, v_fault, &voff);CHECK_PETSC_ERROR(err);
      assert(vdof == spaceDim);
      
      for(PetscInt d = 0; d < spaceDim; ++d) {
        tractPerturbVertex[d] = tractionsArray[voff+d];
      } // for
    } else {
      tractPerturbVertex = 0.0;
    } // if/else

    // Get orientation associated with fault vertex.
    PetscInt odof, ooff;

    err = PetscSectionGetDof(orientationSection, v_fault, &odof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(orientationSection, v_fault, &ooff);CHECK_PETSC_ERROR(err);
    assert(spaceDim*spaceDim == odof);

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

    // Compute slip (in fault coordinates system) from displacements.
    PylithScalar   slipNormal     = 0.0;
    PylithScalar   tractionNormal = 0.0;
    const PetscInt indexN         = spaceDim - 1;
    for(PetscInt d = 0; d < spaceDim; ++d) {
      slipNormal     += orientationArray[ooff+indexN*spaceDim+d] * (dispTArray[dtpoff+d] + dispTIncrArray[dipoff+d] - dispTArray[dtnoff+d] - dispTIncrArray[dinoff+d]);
      tractionNormal += orientationArray[ooff+indexN*spaceDim+d] * (dispTArray[dtloff+d] + dispTIncrArray[diloff+d]);
    } // for
    

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(computeEvent);
    _logger->eventBegin(updateEvent);
#endif
    if (slipNormal < _zeroTolerance || !_openFreeSurf) { 
      // if no opening or flag indicates to still impose initial tractions when fault is open.
      // Assemble contributions into field
      PetscInt rndof, rnoff, rpdof, rpoff;

      err = PetscSectionGetDof(residualSection, v_negative, &rndof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(residualSection, v_negative, &rnoff);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetDof(residualSection, v_positive, &rpdof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(residualSection, v_positive, &rpoff);CHECK_PETSC_ERROR(err);
      assert(spaceDim == rndof);assert(spaceDim == rpdof);
      // Initial (external) tractions oppose (internal) tractions associated with Lagrange multiplier.
      for(PetscInt d = 0; d < spaceDim; ++d) {
        residualArray[rnoff+d] +=  areaArray[aoff] * (dispTArray[dtloff+d] + dispTIncrArray[diloff+d] - tractPerturbVertex[d]);
        residualArray[rpoff+d] += -areaArray[aoff] * (dispTArray[dtloff+d] + dispTIncrArray[diloff+d] - tractPerturbVertex[d]);
      }
    } else { // opening, normal traction should be zero
      if (fabs(tractionNormal) > _zeroTolerance) {
        std::cerr << "WARNING! Fault opening with nonzero traction."
                  << ", v_fault: " << v_fault
                  << ", opening: " << slipNormal
                  << ", normal traction: " << tractionNormal
                  << std::endl;
      } // if
      assert(fabs(tractionNormal) < _zeroTolerance);
    }  // if/else

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(updateEvent);
#endif
  } // for
  PetscLogFlops(numVertices*spaceDim*8);
  err = VecRestoreArray(areaVec, &areaArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(orientationVec, &orientationArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(dispTVec, &dispTArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(dispTIncrVec, &dispTIncrArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(residualVec, &residualArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(valuesVec, &tractionsArray);CHECK_PETSC_ERROR(err);

#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventEnd(computeEvent);
#endif
} // integrateResidual

// ----------------------------------------------------------------------
// Update state variables as needed.
void
pylith::faults::FaultCohesiveDyn::updateStateVars(
				      const PylithScalar t,
				      topology::SolutionFields* const fields)
{ // updateStateVars
  assert(fields);
  assert(_fields);

  _updateRelMotion(*fields);

  const int spaceDim = _quadrature->spaceDim();

  // Allocate arrays for vertex values
  scalar_array tractionTpdtVertex(spaceDim); // Fault coordinate system
  PetscErrorCode err;

  // Get sections
  scalar_array slipVertex(spaceDim);
  PetscSection dispRelSection = _fields->get("relative disp").petscSection();
  Vec          dispRelVec     = _fields->get("relative disp").localVector();
  PetscScalar *dispRelArray;
  assert(dispRelSection);assert(dispRelVec);

  scalar_array slipRateVertex(spaceDim);
  PetscSection velRelSection = _fields->get("relative velocity").petscSection();
  Vec          velRelVec     = _fields->get("relative velocity").localVector();
  PetscScalar *velRelArray;
  assert(velRelSection);assert(velRelVec);

  PetscSection dispTSection = fields->get("disp(t)").petscSection();
  Vec          dispTVec     = fields->get("disp(t)").localVector();
  PetscScalar *dispTArray;
  assert(dispTSection);assert(dispTVec);

  PetscSection dispTIncrSection = fields->get("dispIncr(t->t+dt)").petscSection();
  Vec          dispTIncrVec     = fields->get("dispIncr(t->t+dt)").localVector();
  PetscScalar *dispTIncrArray;
  assert(dispTIncrSection);assert(dispTIncrVec);

  PetscSection orientationSection = _fields->get("orientation").petscSection();
  Vec          orientationVec     = _fields->get("orientation").localVector();
  PetscScalar *orientationArray;
  assert(orientationSection);assert(orientationVec);

  const int numVertices = _cohesiveVertices.size();
  err = VecGetArray(dispRelVec, &dispRelArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(velRelVec, &velRelArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(dispTVec, &dispTArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(dispTIncrVec, &dispTIncrArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(orientationVec, &orientationArray);CHECK_PETSC_ERROR(err);
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_lagrange = _cohesiveVertices[iVertex].lagrange;
    const int v_fault = _cohesiveVertices[iVertex].fault;
    const int v_negative = _cohesiveVertices[iVertex].negative;
    const int v_positive = _cohesiveVertices[iVertex].positive;

    // Get relative displacement
    PetscInt drdof, droff;

    err = PetscSectionGetDof(dispRelSection, v_fault, &drdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dispRelSection, v_fault, &droff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == drdof);

    // Get relative velocity
    PetscInt vrdof, vroff;

    err = PetscSectionGetDof(velRelSection, v_fault, &vrdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(velRelSection, v_fault, &vroff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == vrdof);

    // Get orientation
    PetscInt odof, ooff;

    err = PetscSectionGetDof(orientationSection, v_fault, &odof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(orientationSection, v_fault, &ooff);CHECK_PETSC_ERROR(err);
    assert(spaceDim*spaceDim == odof);

    // Get Lagrange multiplier values from disp(t), and dispIncr(t->t+dt)
    PetscInt dtldof, dtloff;

    err = PetscSectionGetDof(dispTSection, v_lagrange, &dtldof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dispTSection, v_lagrange, &dtloff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == dtldof);
    PetscInt dildof, diloff;

    err = PetscSectionGetDof(dispTIncrSection, v_lagrange, &dildof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dispTIncrSection, v_lagrange, &diloff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == dildof);

    // Compute slip, slip rate, and fault traction (Lagrange
    // multiplier) at time t+dt in fault coordinate system.
    slipVertex = 0.0;
    slipRateVertex = 0.0;
    tractionTpdtVertex = 0.0;
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      for (int jDim=0; jDim < spaceDim; ++jDim) {
        slipVertex[iDim] += orientationArray[ooff+iDim*spaceDim+jDim] * dispRelArray[droff+jDim];
        slipRateVertex[iDim] += orientationArray[ooff+iDim*spaceDim+jDim] * velRelArray[vroff+jDim];
        tractionTpdtVertex[iDim] += orientationArray[ooff+iDim*spaceDim+jDim] * (dispTArray[dtloff+jDim]+dispTIncrArray[diloff+jDim]);
      } // for
    } // for

    // Get friction properties and state variables.
    _friction->retrievePropsStateVars(v_fault);

    // Use fault constitutive model to compute traction associated with
    // friction.
    switch (spaceDim) { // switch
    case 1: { // case 1
      const PylithScalar slipMag = 0.0;
      const PylithScalar slipRateMag = 0.0;
      const PylithScalar tractionNormal = tractionTpdtVertex[0];
      _friction->updateStateVars(t, slipMag, slipRateMag, tractionNormal, v_fault);
      break;
    } // case 1
    case 2: { // case 2
      const PylithScalar slipMag = fabs(slipVertex[0]);
      const PylithScalar slipRateMag = fabs(slipRateVertex[0]);
      const PylithScalar tractionNormal = tractionTpdtVertex[1];
      _friction->updateStateVars(t, slipMag, slipRateMag, tractionNormal, v_fault);
      break;
    } // case 2
    case 3: { // case 3
      const PylithScalar slipMag = 
	sqrt(slipVertex[0]*slipVertex[0] + slipVertex[1]*slipVertex[1]);
      const PylithScalar slipRateMag = 
	sqrt(slipRateVertex[0]*slipRateVertex[0] + 
	     slipRateVertex[1]*slipRateVertex[1]);
      const PylithScalar tractionNormal = tractionTpdtVertex[2];
      _friction->updateStateVars(t, slipMag, slipRateMag, tractionNormal, v_fault);
      break;
    } // case 3
    default:
      assert(0);
      throw std::logic_error("Unknown spatial dimension in FaultCohesiveDyn::updateStateVars().");
    } // switch
  } // for
  err = VecRestoreArray(dispRelVec, &dispRelArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(velRelVec, &velRelArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(dispTVec, &dispTArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(dispTIncrVec, &dispTIncrArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(orientationVec, &orientationArray);CHECK_PETSC_ERROR(err);
} // updateStateVars

// ----------------------------------------------------------------------
// Constrain solution based on friction.
void
pylith::faults::FaultCohesiveDyn::constrainSolnSpace(
				    topology::SolutionFields* const fields,
				    const PylithScalar t,
				    const topology::Jacobian& jacobian)
{ // constrainSolnSpace
  /// Member prototype for _constrainSolnSpaceXD()
  typedef void (pylith::faults::FaultCohesiveDyn::*constrainSolnSpace_fn_type)
    (scalar_array*,
     const PylithScalar,
     const scalar_array&,
     const scalar_array&,
     const scalar_array&,
     const bool);

  assert(fields);
  assert(_quadrature);
  assert(_fields);
  assert(_friction);

  _sensitivitySetup(jacobian);

  // Update time step in friction (can vary).
  _friction->timeStep(_dt);
  const PylithScalar dt = _dt;

  const int spaceDim = _quadrature->spaceDim();
  const int indexN = spaceDim - 1;

  // Allocate arrays for vertex values
  scalar_array tractionTpdtVertex(spaceDim);
  scalar_array dDispRelVertex(spaceDim);
  PetscErrorCode err;

  // Get sections
  scalar_array slipTpdtVertex(spaceDim);
  PetscSection dispRelSection = _fields->get("relative disp").petscSection();
  Vec          dispRelVec     = _fields->get("relative disp").localVector();
  PetscScalar *dispRelArray;
  assert(dispRelSection);assert(dispRelVec);

  scalar_array slipRateVertex(spaceDim);
  PetscSection velRelSection = _fields->get("relative velocity").petscSection();
  Vec          velRelVec     = _fields->get("relative velocity").localVector();
  PetscScalar *velRelArray;
  assert(velRelSection);assert(velRelVec);

  PetscSection orientationSection = _fields->get("orientation").petscSection();
  Vec          orientationVec     = _fields->get("orientation").localVector();
  PetscScalar *orientationArray;
  assert(orientationSection);assert(orientationVec);

  PetscSection dispTSection = fields->get("disp(t)").petscSection();
  Vec          dispTVec     = fields->get("disp(t)").localVector();
  PetscScalar *dispTArray;
  assert(dispTSection);assert(dispTVec);

  scalar_array dDispTIncrVertexN(spaceDim);
  scalar_array dDispTIncrVertexP(spaceDim);
  DM           dispTIncrDM      = fields->get("dispIncr(t->t+dt)").dmMesh();
  PetscSection dispTIncrSection = fields->get("dispIncr(t->t+dt)").petscSection();
  Vec          dispTIncrVec     = fields->get("dispIncr(t->t+dt)").localVector();
  PetscSection dispTIncrGlobalSection;
  PetscScalar *dispTIncrArray;
  assert(dispTIncrSection);assert(dispTIncrVec);
  err = DMGetDefaultGlobalSection(dispTIncrDM, &dispTIncrGlobalSection);CHECK_PETSC_ERROR(err);

  PetscSection dispTIncrAdjSection = fields->get("dispIncr adjust").petscSection();
  Vec          dispTIncrAdjVec     = fields->get("dispIncr adjust").localVector();
  PetscScalar *dispTIncrAdjArray;
  assert(dispTIncrAdjSection);assert(dispTIncrAdjVec);

  scalar_array dTractionTpdtVertex(spaceDim);
  scalar_array dLagrangeTpdtVertex(spaceDim);
  PetscSection sensitivitySection = _fields->get("sensitivity dLagrange").petscSection();
  Vec          sensitivityVec     = _fields->get("sensitivity dLagrange").localVector();
  PetscScalar *sensitivityArray;
  assert(sensitivitySection);assert(sensitivityVec);

  constrainSolnSpace_fn_type constrainSolnSpaceFn;
  switch (spaceDim) { // switch
  case 1:
    constrainSolnSpaceFn = 
      &pylith::faults::FaultCohesiveDyn::_constrainSolnSpace1D;
    break;
  case 2: 
    constrainSolnSpaceFn = 
      &pylith::faults::FaultCohesiveDyn::_constrainSolnSpace2D;
    break;
  case 3:
    constrainSolnSpaceFn = 
      &pylith::faults::FaultCohesiveDyn::_constrainSolnSpace3D;
    break;
  default :
    assert(0);
    throw std::logic_error("Unknown spatial dimension in "
			   "FaultCohesiveDyn::constrainSolnSpace().");
  } // switch


#if 0 // DEBUGGING
  dispRelSection->view("BEFORE RELATIVE DISPLACEMENT");
  dispIncrSection->view("BEFORE DISP INCR (t->t+dt)");
#endif

  const int numVertices = _cohesiveVertices.size();
  err = VecGetArray(dispTVec, &dispTArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(dispTIncrVec, &dispTIncrArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(dispTIncrAdjVec, &dispTIncrAdjArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(orientationVec, &orientationArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(sensitivityVec, &sensitivityArray);CHECK_PETSC_ERROR(err);
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_lagrange = _cohesiveVertices[iVertex].lagrange;
    const int v_fault = _cohesiveVertices[iVertex].fault;
    const int v_negative = _cohesiveVertices[iVertex].negative;
    const int v_positive = _cohesiveVertices[iVertex].positive;

    // Get displacement values
    PetscInt dtndof, dtnoff;

    err = PetscSectionGetDof(dispTSection, v_negative, &dtndof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dispTSection, v_negative, &dtnoff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == dtndof);
    PetscInt dtpdof, dtpoff;

    err = PetscSectionGetDof(dispTSection, v_positive, &dtpdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dispTSection, v_positive, &dtpoff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == dtpdof);

    // Get displacement increment values.
    PetscInt dindof, dinoff;

    err = PetscSectionGetDof(dispTIncrSection, v_negative, &dindof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dispTIncrSection, v_negative, &dinoff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == dindof);
    PetscInt dipdof, dipoff;

    err = PetscSectionGetDof(dispTIncrSection, v_positive, &dipdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dispTIncrSection, v_positive, &dipoff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == dipdof);

    // Get orientation
    PetscInt odof, ooff;

    err = PetscSectionGetDof(orientationSection, v_fault, &odof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(orientationSection, v_fault, &ooff);CHECK_PETSC_ERROR(err);
    assert(spaceDim*spaceDim == odof);

    // Get Lagrange multiplier values from disp(t), and dispIncr(t->t+dt)
    PetscInt dtldof, dtloff;

    err = PetscSectionGetDof(dispTSection, v_lagrange, &dtldof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dispTSection, v_lagrange, &dtloff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == dtldof);
    PetscInt dildof, diloff;

    err = PetscSectionGetDof(dispTIncrSection, v_lagrange, &dildof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dispTIncrSection, v_lagrange, &diloff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == dildof);

    // Step 1: Prevent nonphysical trial solutions. The product of the
    // normal traction and normal slip must be nonnegative (forbid
    // interpenetration with tension or opening with compression).
    
    // Compute slip, slip rate, and Lagrange multiplier at time t+dt
    // in fault coordinate system.
    slipTpdtVertex = 0.0;
    slipRateVertex = 0.0;
    tractionTpdtVertex = 0.0;
    for(PetscInt d = 0; d < spaceDim; ++d) {
      for(PetscInt e = 0; e < spaceDim; ++e) {
        slipTpdtVertex[d] += orientationArray[ooff+d*spaceDim+e] * (dispTArray[dtpoff+e] + dispTIncrArray[dipoff+e] - dispTArray[dtnoff+e] - dispTIncrArray[dinoff+e]);
        slipRateVertex[d] += orientationArray[ooff+d*spaceDim+e] * (dispTIncrArray[dipoff+e] - dispTIncrArray[dinoff+e]) / dt;
        tractionTpdtVertex[d] += orientationArray[ooff+d*spaceDim+e] * (dispTArray[dtloff+e] + dispTIncrArray[diloff+e]);
      } // for
      if (fabs(slipRateVertex[d]) < _zeroTolerance) {
        slipRateVertex[d] = 0.0;
      } // if
    } // for
    if (fabs(slipTpdtVertex[indexN]) < _zeroTolerance) {
      slipTpdtVertex[indexN] = 0.0;
    } // if

    PylithScalar dSlipVertexNormal = 0.0;
    PylithScalar dTractionTpdtVertexNormal = 0.0;
    if (slipTpdtVertex[indexN]*tractionTpdtVertex[indexN] < 0.0) {
#if 0 // DEBUGGING
      std::cout << "STEP 1 CORRECTING NONPHYSICAL SLIP/TRACTIONS"
		<< ", v_fault: " << v_fault
		<< ", slipNormal: " << slipTpdtVertex[indexN]
		<< ", tractionNormal: " << tractionTpdtVertex[indexN]
		<< std::endl;
#endif
      // Don't know what behavior is appropriate so set smaller of
      // traction and slip to zero (should be appropriate if problem
      // is nondimensionalized correctly).
      if (fabs(slipTpdtVertex[indexN]) > fabs(tractionTpdtVertex[indexN])) {
        // slip is bigger, so force normal traction back to zero
        dTractionTpdtVertexNormal = -tractionTpdtVertex[indexN];
        tractionTpdtVertex[indexN] = 0.0;
      } else {
        // traction is bigger, so force slip back to zero
        dSlipVertexNormal = -slipTpdtVertex[indexN];
        slipTpdtVertex[indexN] = 0.0;
      } // if/else
    } // if
    if (slipTpdtVertex[indexN] < 0.0) {
#if 0 // DEBUGGING
      std::cout << "STEP 1 CORRECTING INTERPENETRATION"
		<< ", v_fault: " << v_fault
		<< ", slipNormal: " << slipTpdtVertex[indexN]
		<< ", tractionNormal: " << tractionTpdtVertex[indexN]
		<< std::endl;
#endif
      dSlipVertexNormal = -slipTpdtVertex[indexN];
      slipTpdtVertex[indexN] = 0.0;
    } // if

    // Step 2: Apply friction criterion to trial solution to get
    // change in Lagrange multiplier (dTractionTpdtVertex) in fault
    // coordinate system.

    // Get friction properties and state variables.
    _friction->retrievePropsStateVars(v_fault);

    // Use fault constitutive model to compute traction associated with
    // friction.
    dTractionTpdtVertex = 0.0;
    const bool iterating = true; // Iterating to get friction
    CALL_MEMBER_FN(*this, constrainSolnSpaceFn)(&dTractionTpdtVertex, t, slipTpdtVertex, slipRateVertex, tractionTpdtVertex, iterating);

    // Rotate increment in traction back to global coordinate system.
    dLagrangeTpdtVertex = 0.0;
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      for (int jDim=0; jDim < spaceDim; ++jDim) {
        dLagrangeTpdtVertex[iDim] += orientationArray[ooff+jDim*spaceDim+iDim] * dTractionTpdtVertex[jDim];
      } // for

      // Add in potential contribution from adjusting Lagrange
      // multiplier for fault normal DOF of trial solution in Step 1.
      dLagrangeTpdtVertex[iDim] += orientationArray[ooff+indexN*spaceDim+iDim] * dTractionTpdtVertexNormal;
    } // for

#if 0 // debugging
    std::cout << "v_fault: " << v_fault;
    std::cout << ", slipVertex: ";
    for (int iDim=0; iDim < spaceDim; ++iDim)
      std::cout << "  " << slipTpdtVertex[iDim];
    std::cout << ",  slipRateVertex: ";
    for (int iDim=0; iDim < spaceDim; ++iDim)
      std::cout << "  " << slipRateVertex[iDim];
    std::cout << ",  tractionTpdtVertex: ";
    for (int iDim=0; iDim < spaceDim; ++iDim)
      std::cout << "  " << tractionTpdtVertex[iDim];
    std::cout << ",  dTractionTpdtVertex: ";
    for (int iDim=0; iDim < spaceDim; ++iDim)
      std::cout << "  " << dTractionTpdtVertex[iDim];
    std::cout << std::endl;
#endif
     
    // Set change in Lagrange multiplier
    PetscInt sdof, soff;

    err = PetscSectionGetDof(sensitivitySection, v_fault, &sdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(sensitivitySection, v_fault, &soff);CHECK_PETSC_ERROR(err);
    assert(dLagrangeTpdtVertex.size() == sdof);
    for(PetscInt d = 0; d < sdof; ++d) {
      sensitivityArray[soff+d] = dLagrangeTpdtVertex[d];
    }

#if 0 // UNNECESSARY?
    // Update displacement in trial solution (if necessary) so that it
    // conforms to physical constraints.
    if (0.0 != dSlipVertexNormal) {
      // Compute relative displacement from slip.
      dDispRelVertex = 0.0;
      for (int iDim=0; iDim < spaceDim; ++iDim) {
        dDispRelVertex[iDim] += orientationArray[ooff+indexN*spaceDim+iDim] * dSlipVertexNormal;

        dDispTIncrVertexN[iDim] = -0.5*dDispRelVertex[iDim];
        dDispTIncrVertexP[iDim] = +0.5*dDispRelVertex[iDim];
      } // for

      // Update displacement field
      assert(dDispTIncrVertexN.size() == dispIncrSection->getFiberDimension(v_negative));
      dispIncrAdjSection->updateAddPoint(v_negative, &dDispTIncrVertexN[0]);
      
      assert(dDispTIncrVertexP.size() == dispIncrSection->getFiberDimension(v_positive));
      dispIncrAdjSection->updateAddPoint(v_positive, &dDispTIncrVertexP[0]);    
    } // if
#endif

  } // for
  err = VecRestoreArray(dispTVec, &dispTArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(dispTIncrVec, &dispTIncrArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(dispTIncrAdjVec, &dispTIncrAdjArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(orientationVec, &orientationArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(sensitivityVec, &sensitivityArray);CHECK_PETSC_ERROR(err);

  // Step 3: Calculate change in displacement field corresponding to
  // change in Lagrange multipliers imposed by friction criterion.

  // Solve sensitivity problem for negative side of the fault.
  bool negativeSideFlag = true;
  _sensitivityUpdateJacobian(negativeSideFlag, jacobian, *fields);
  _sensitivityReformResidual(negativeSideFlag);
  _sensitivitySolve();
  _sensitivityUpdateSoln(negativeSideFlag);

  // Solve sensitivity problem for positive side of the fault.
  negativeSideFlag = false;
  _sensitivityUpdateJacobian(negativeSideFlag, jacobian, *fields);
  _sensitivityReformResidual(negativeSideFlag);
  _sensitivitySolve();
  _sensitivityUpdateSoln(negativeSideFlag);

  // Step 4: Update Lagrange multipliers and displacement fields based
  // on changes imposed by friction criterion in Step 2 (change in
  // Lagrange multipliers) and Step 3 (slip associated with change in
  // Lagrange multipliers).
  //
  // Use line search to find best update. This improves convergence
  // because it accounts for feedback between the fault constitutive
  // model and the deformation. We also search in log space because
  // some fault constitutive models depend on the log of slip rate.

  const PylithScalar residualTol = _zeroTolerance; // L2 misfit in tractions
  const int maxIter = 16;
  PylithScalar logAlphaL = log10(_zeroTolerance); // minimum step
  PylithScalar logAlphaR = log10(1.0); // maximum step
  PylithScalar logAlphaM = 0.5*(logAlphaL + logAlphaR);
  PylithScalar logAlphaML = 0.5*(logAlphaL + logAlphaM);
  PylithScalar logAlphaMR = 0.5*(logAlphaM + logAlphaR);
  PylithScalar residualL = _constrainSolnSpaceNorm(pow(10.0, logAlphaL), t, fields);
  PylithScalar residualML = _constrainSolnSpaceNorm(pow(10.0, logAlphaML), t, fields);
  PylithScalar residualM = _constrainSolnSpaceNorm(pow(10.0, logAlphaM), t, fields);
  PylithScalar residualMR = _constrainSolnSpaceNorm(pow(10.0, logAlphaMR), t, fields);
  PylithScalar residualR = _constrainSolnSpaceNorm(pow(10.0, logAlphaR), t, fields);
  for (int iter=0; iter < maxIter; ++iter) {
    if (residualM < residualTol || residualR < residualTol)
      // if residual is very small, we prefer the full step
      break;

#if 0
    const int rank = _faultMesh->sieveMesh()->commRank();
    std::cout << "["<<rank<<"] alphaL: " << pow(10.0, logAlphaL)
	      << ", residuaL: " << residualL
	      << ", alphaM: " << pow(10.0, logAlphaM)
	      << ", residualM: " << residualM
	      << ", alphaR: " << pow(10.0, logAlphaR)
	      << ", residualR: " << residualR
	      << std::endl;
#endif

    if (residualL < residualML && residualL < residualM && residualL < residualMR && residualL < residualR) {
      logAlphaL = logAlphaL;
      logAlphaR = logAlphaM;
      residualL = residualL;
      residualR = residualM;
      residualM = residualML;
    } else if (residualML <= residualL  && residualML < residualM && residualML < residualMR && residualML < residualR) {
      logAlphaL = logAlphaL;
      logAlphaR = logAlphaM;
      residualL = residualL;
      residualR = residualM;
      residualM = residualML;
    } else if (residualM <= residualL  && residualM <= residualML && residualM < residualMR && residualM < residualR) {
      logAlphaL = logAlphaML;
      logAlphaR = logAlphaMR;
      residualL = residualML;
      residualR = residualMR;
      residualM = residualM;
    } else if (residualMR <= residualL  && residualMR <= residualML && residualMR <= residualM && residualMR < residualR) {
      logAlphaL = logAlphaM;
      logAlphaR = logAlphaR;
      residualL = residualM;
      residualR = residualR;
      residualM = residualMR;
    } else if (residualR <= residualL  && residualR <= residualML && residualR <= residualM && residualR <= residualMR) {
      logAlphaL = logAlphaM;
      logAlphaR = logAlphaR;
      residualL = residualM;
      residualR = residualR;
      residualM = residualMR;
    } else {
      assert(0);
      throw std::logic_error("Unknown case in constrain solution space "
			     "line search.");
    } // if/else
    logAlphaM = (logAlphaL + logAlphaR) / 2.0;
    logAlphaML = (logAlphaL + logAlphaM) / 2.0;
    logAlphaMR = (logAlphaM + logAlphaR) / 2.0;

    residualML = _constrainSolnSpaceNorm(pow(10.0, logAlphaML), t, fields);
    residualMR = _constrainSolnSpaceNorm(pow(10.0, logAlphaMR), t, fields);

  } // for
  // Account for possibility that end points have lowest residual
  if (residualR <= residualM || residualR < residualTol) {
    logAlphaM = logAlphaR;
    residualM = residualR;
  } else if (residualL < residualM) {
    logAlphaM = logAlphaL;
    residualM = residualL;
  } // if/else
  const PylithScalar alpha = pow(10.0, logAlphaM); // alphaM is our best guess
#if 0 // DEBUGGING
  std::cout << "ALPHA: " << alpha
	    << ", residual: " << residualM
	    << std::endl;
#endif

  scalar_array slipTVertex(spaceDim);
  scalar_array dSlipTpdtVertex(spaceDim);
  scalar_array dispRelVertex(spaceDim);
  PetscSection sensitivityRelSection = _fields->get("sensitivity relative disp").petscSection();
  Vec          sensitivityRelVec     = _fields->get("sensitivity relative disp").localVector();
  PetscScalar *sensitivityRelArray;
  assert(sensitivityRelSection);assert(sensitivityRelVec);

  err = VecGetArray(dispTVec, &dispTArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(dispTIncrVec, &dispTIncrArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(dispTIncrAdjVec, &dispTIncrAdjArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(orientationVec, &orientationArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(sensitivityVec, &sensitivityArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(sensitivityRelVec, &sensitivityRelArray);CHECK_PETSC_ERROR(err);
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_fault = _cohesiveVertices[iVertex].fault;
    const int v_lagrange = _cohesiveVertices[iVertex].lagrange;
    const int v_negative = _cohesiveVertices[iVertex].negative;
    const int v_positive = _cohesiveVertices[iVertex].positive;

    // Get change in Lagrange multiplier computed from friction criterion.
    PetscInt sdof, soff;

    err = PetscSectionGetDof(sensitivitySection, v_fault, &sdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(sensitivitySection, v_fault, &soff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == sdof);

    // Get change in relative displacement from sensitivity solve.
    PetscInt srdof, sroff;

    err = PetscSectionGetDof(sensitivityRelSection, v_fault, &srdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(sensitivityRelSection, v_fault, &sroff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == srdof);

    // Get current relative displacement for updating.
    PetscInt drdof, droff;

    err = PetscSectionGetDof(dispRelSection, v_fault, &drdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dispRelSection, v_fault, &droff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == drdof);

    // Get orientation.
    PetscInt odof, ooff;

    err = PetscSectionGetDof(orientationSection, v_fault, &odof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(orientationSection, v_fault, &ooff);CHECK_PETSC_ERROR(err);
    assert(spaceDim*spaceDim == odof);

    // Get displacement.
    PetscInt dtndof, dtnoff;

    err = PetscSectionGetDof(dispTSection, v_negative, &dtndof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dispTSection, v_negative, &dtnoff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == dtndof);
    PetscInt dtpdof, dtpoff;

    err = PetscSectionGetDof(dispTSection, v_positive, &dtpdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dispTSection, v_positive, &dtpoff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == dtpdof);

    // Get displacement increment (trial solution).
    PetscInt dindof, dinoff;

    err = PetscSectionGetDof(dispTIncrSection, v_negative, &dindof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dispTIncrSection, v_negative, &dinoff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == dindof);
    PetscInt dipdof, dipoff;

    err = PetscSectionGetDof(dispTIncrSection, v_positive, &dipdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dispTIncrSection, v_positive, &dipoff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == dipdof);

    // Get Lagrange multiplier at time t
    PetscInt dtldof, dtloff;

    err = PetscSectionGetDof(dispTSection, v_lagrange, &dtldof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dispTSection, v_lagrange, &dtloff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == dtldof);

    // Get Lagrange multiplier increment (trial solution)
    PetscInt dildof, diloff;

    err = PetscSectionGetDof(dispTIncrSection, v_lagrange, &dildof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dispTIncrSection, v_lagrange, &diloff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == dildof);

    // Scale perturbation in relative displacements and change in
    // Lagrange multipliers by alpha using only shear components.
    slipTVertex = 0.0;
    slipTpdtVertex = 0.0;
    dSlipTpdtVertex = 0.0;
    tractionTpdtVertex = 0.0;
    dTractionTpdtVertex = 0.0;
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      for (int jDim=0; jDim < spaceDim; ++jDim) {
        slipTVertex[iDim] += orientationArray[ooff+iDim*spaceDim+jDim] * (dispTArray[dtpoff+jDim] - dispTArray[dtnoff+jDim]);
        slipTpdtVertex[iDim] += orientationArray[ooff+iDim*spaceDim+jDim] * (dispTArray[dtpoff+jDim] - dispTArray[dtnoff+jDim] + dispTIncrArray[dipoff+jDim] - dispTIncrArray[dinoff+jDim]);
        dSlipTpdtVertex[iDim] += orientationArray[ooff+iDim*spaceDim+jDim] * alpha*sensitivityRelArray[sroff+jDim];
        tractionTpdtVertex[iDim] += orientationArray[ooff+iDim*spaceDim+jDim] * (dispTArray[dtloff+jDim] + dispTIncrArray[diloff+jDim]);
        dTractionTpdtVertex[iDim] += orientationArray[ooff+iDim*spaceDim+jDim] * alpha*sensitivityArray[soff+jDim];
      } // for
    } // for

    // FIRST, correct nonphysical trial solutions.
    // Order of steps 5a-5c is important!
    if ((slipTpdtVertex[indexN] + dSlipTpdtVertex[indexN]) * (tractionTpdtVertex[indexN] + dTractionTpdtVertex[indexN])	< 0.0) {
      // Step 5a: Prevent nonphysical trial solutions. The product of the
      // normal traction and normal slip must be nonnegative (forbid
      // interpenetration with tension or opening with compression).
      
#if 0 // DEBUGGING
      std::cout << "STEP 5a CORRECTING NONPHYSICAL SLIP/TRACTIONS"
		<< ", v_fault: " << v_fault
		<< ", slipNormal: " << slipTpdtVertex[indexN] + dSlipTpdtVertex[indexN]
		<< ", tractionNormal: " << tractionTpdtVertex[indexN] + dTractionTpdtVertex[indexN]
		<< std::endl;
#endif
      // Don't know what behavior is appropriate so set smaller of
      // traction and slip to zero (should be appropriate if problem
      // is nondimensionalized correctly).
      if (fabs(slipTpdtVertex[indexN] + dSlipTpdtVertex[indexN]) > fabs(tractionTpdtVertex[indexN] + dTractionTpdtVertex[indexN])) {
        // slip is bigger, so force normal traction back to zero
        dTractionTpdtVertex[indexN] = -tractionTpdtVertex[indexN];
      } else {
        // traction is bigger, so force slip back to zero
        dSlipTpdtVertex[indexN] = -slipTpdtVertex[indexN];
      } // if/else

    } else if (slipTpdtVertex[indexN] + dSlipTpdtVertex[indexN] > _zeroTolerance) {
      // Step 5b: Insure fault traction is zero when opening (if alpha=1
      // this should be enforced already, but will not be properly
      // enforced when alpha < 1).
      for (int iDim=0; iDim < spaceDim; ++iDim) {
        dTractionTpdtVertex[iDim] = -tractionTpdtVertex[iDim];
      } // for
      
    } else if (slipTpdtVertex[indexN] + dSlipTpdtVertex[indexN] < 0.0) {
      // Step 5c: Prevent interpenetration.
#if 0 // DEBUGGING
      std::cout << "STEP 5b CORRECTING INTERPENETATION"
		<< ", v_fault: " << v_fault
		<< ", slipNormal: " << slipTpdtVertex[indexN] + dSlipTpdtVertex[indexN]
		<< std::endl;
#endif
      dSlipTpdtVertex[indexN] = -slipTpdtVertex[indexN];
      
    } // if/else      

#if 0 // UNNECESSARY????
    // Step 5d: Prevent over-correction in shear slip resulting in backslip
    PylithScalar slipDot = 0.0;
    PylithScalar tractionSlipDot = 0.0;
    for (int iDim=0; iDim < indexN; ++iDim)  {
      // Compute dot product between slip and increment in slip (want positive)
      slipDot += (slipTpdtVertex[iDim] - slipTVertex[iDim]) * (slipTpdtVertex[iDim] + dSlipTpdtVertex[iDim] - slipTVertex[iDim]);
      // Compute dot product of traction and slip
      tractionSlipDot += (tractionTpdtVertex[iDim] + dTractionTpdtVertex[iDim]) * (slipTpdtVertex[iDim] + dSlipTpdtVertex[iDim]);
    } // for
    if (slipDot < 0.0 && sqrt(fabs(slipDot)) > _zeroTolerance && tractionSlipDot < 0.0) {
#if 0 // DEBUGGING
      std::cout << "STEP 5d CORRECTING BACKSLIP"
		<< ", v_fault: " << v_fault
		<< ", slipDot: " << slipDot
		<< ", tractionSlipDot: " << tractionSlipDot
		<< std::endl;
#endif
      // Correct backslip, use bisection as guess      
      for (int iDim=0; iDim < indexN; ++iDim) {
        dTractionTpdtVertex[iDim] *= 0.5;
        dSlipTpdtVertex[iDim] = 0.5*(slipTVertex[iDim] - slipTpdtVertex[iDim]);
      } // for
    } // if/else
#endif
    
    // Update current estimate of slip from t to t+dt.
    slipTpdtVertex += dSlipTpdtVertex;
    
    // Compute relative displacement from slip.
    dispRelVertex = 0.0;
    dDispRelVertex = 0.0;
    dLagrangeTpdtVertex = 0.0;
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      for (int jDim=0; jDim < spaceDim; ++jDim) {
        dispRelVertex[iDim] += orientationArray[ooff+jDim*spaceDim+iDim] * slipTpdtVertex[jDim];
        dDispRelVertex[iDim] += orientationArray[ooff+jDim*spaceDim+iDim] * dSlipTpdtVertex[jDim];
        dLagrangeTpdtVertex[iDim] += orientationArray[ooff+jDim*spaceDim+iDim] * dTractionTpdtVertex[jDim];
      } // for

      dDispTIncrVertexN[iDim] = -0.5*dDispRelVertex[iDim];
      dDispTIncrVertexP[iDim] = +0.5*dDispRelVertex[iDim];
    } // for

#if 0 // debugging
    std::cout << "v_fault: " << v_fault;
    std::cout << ", tractionTpdtVertex: ";
    for (int iDim=0; iDim < spaceDim; ++iDim)
      std::cout << "  " << tractionTpdtVertex[iDim];
    std::cout << ", dTractionTpdtVertex: ";
    for (int iDim=0; iDim < spaceDim; ++iDim)
      std::cout << "  " << dTractionTpdtVertex[iDim];
    std::cout << ", slipTpdtVertex: ";
    for (int iDim=0; iDim < spaceDim; ++iDim)
      std::cout << "  " << slipTpdtVertex[iDim]-dSlipTpdtVertex[iDim];
    std::cout << ",  dSlipTpdtVertex: ";
    for (int iDim=0; iDim < spaceDim; ++iDim)
      std::cout << "  " << dSlipTpdtVertex[iDim];
    std::cout << std::endl;
#endif

    // Compute contribution to adjusting solution only if Lagrange
    // constraint is local (the adjustment is assembled across processors).
    PetscInt goff;
    err = PetscSectionGetOffset(dispTIncrGlobalSection, v_lagrange, &goff);CHECK_PETSC_ERROR(err);
    if (goff >= 0) {
      // Update Lagrange multiplier increment.
      PetscInt dialdof, dialoff;

      err = PetscSectionGetDof(dispTIncrAdjSection, v_lagrange, &dialdof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(dispTIncrAdjSection, v_lagrange, &dialoff);CHECK_PETSC_ERROR(err);
      assert(spaceDim == dialdof);
      for(PetscInt d = 0; d < dialdof; ++d) {
        dispTIncrAdjArray[dialoff+d] += dLagrangeTpdtVertex[d];
      }
      // Update displacement field
      PetscInt diandof, dianoff;

      err = PetscSectionGetDof(dispTIncrAdjSection, v_negative, &diandof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(dispTIncrAdjSection, v_negative, &dianoff);CHECK_PETSC_ERROR(err);
      assert(spaceDim == diandof);
      for(PetscInt d = 0; d < diandof; ++d) {
        dispTIncrAdjArray[dianoff+d] += dDispTIncrVertexN[d];
      }
      PetscInt diapdof, diapoff;

      err = PetscSectionGetDof(dispTIncrAdjSection, v_positive, &diapdof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(dispTIncrAdjSection, v_positive, &diapoff);CHECK_PETSC_ERROR(err);
      assert(spaceDim == diapdof);
      for(PetscInt d = 0; d < diapdof; ++d) {
        dispTIncrAdjArray[diapoff+d] += dDispTIncrVertexP[d];
      }
    } // if

  } // for
  err = VecRestoreArray(dispTVec, &dispTArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(dispTIncrVec, &dispTIncrArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(dispTIncrAdjVec, &dispTIncrAdjArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(orientationVec, &orientationArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(sensitivityVec, &sensitivityArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(sensitivityRelVec, &sensitivityRelArray);CHECK_PETSC_ERROR(err);

#if 0 // DEBUGGING
  //dLagrangeTpdtSection->view("AFTER dLagrange");
  dispIncrAdjSection->view("AFTER DISP INCR adjust");
  dispIncrSection->view("AFTER DISP INCR");
#endif
} // constrainSolnSpace

// ----------------------------------------------------------------------
// Adjust solution from solver with lumped Jacobian to match Lagrange
// multiplier constraints.
void
pylith::faults::FaultCohesiveDyn::adjustSolnLumped(
			 topology::SolutionFields* const fields,
			 const PylithScalar t,
			 const topology::Field<topology::Mesh>& jacobian)
{ // adjustSolnLumped
  /// Member prototype for _constrainSolnSpaceXD()
  typedef void (pylith::faults::FaultCohesiveDyn::*constrainSolnSpace_fn_type)
    (scalar_array*,
     const PylithScalar,
     const scalar_array&,
     const scalar_array&,
     const scalar_array&,
     const bool);

  assert(fields);
  assert(_quadrature);

  // Cohesive cells with conventional vertices i and j, and constraint
  // vertex k require three adjustments to the solution:
  //
  //   * DOF k: Compute increment in Lagrange multipliers
  //            dl_k = S^{-1} (-C_ki (A_i^{-1} r_i - C_kj A_j^{-1} r_j + u_i - u_j) - d_k)
  //            S = C_ki (A_i^{-1} + A_j^{-1}) C_ki^T
  //
  //   * Adjust Lagrange multipliers to match friction criterion
  //
  //   * DOF k: Adjust displacement increment (solution) to create slip
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

  // Allocate arrays for vertex values
  scalar_array tractionTpdtVertex(spaceDim);
  scalar_array lagrangeTpdtVertex(spaceDim);
  scalar_array dTractionTpdtVertex(spaceDim);
  scalar_array dLagrangeTpdtVertex(spaceDim);
  PetscErrorCode err;

  // Update time step in friction (can vary).
  _friction->timeStep(_dt);

  // Get section information
  scalar_array slipVertex(spaceDim);
  PetscSection dispRelSection = _fields->get("relative disp").petscSection();
  Vec          dispRelVec     = _fields->get("relative disp").localVector();
  PetscScalar *dispRelArray;
  assert(dispRelSection);assert(dispRelVec);

  scalar_array slipRateVertex(spaceDim);
  PetscSection velRelSection = _fields->get("relative velocity").petscSection();
  Vec          velRelVec     = _fields->get("relative velocity").localVector();
  PetscScalar *velRelArray;
  assert(velRelSection);assert(velRelVec);

  PetscSection orientationSection = _fields->get("orientation").petscSection();
  Vec          orientationVec     = _fields->get("orientation").localVector();
  PetscScalar *orientationArray;
  assert(orientationSection);assert(orientationVec);

  PetscSection areaSection = _fields->get("area").petscSection();
  Vec          areaVec     = _fields->get("area").localVector();
  PetscScalar *areaArray;
  assert(areaSection);assert(areaVec);

  PetscSection dispTSection = fields->get("disp(t)").petscSection();
  Vec          dispTVec     = fields->get("disp(t)").localVector();
  PetscScalar *dispTArray;
  assert(dispTSection);assert(dispTVec);

  scalar_array dispIncrVertexN(spaceDim);
  scalar_array dispIncrVertexP(spaceDim);
  scalar_array lagrangeTIncrVertex(spaceDim);
  PetscSection dispTIncrSection = fields->get("dispIncr(t->t+dt)").petscSection();
  Vec          dispTIncrVec     = fields->get("dispIncr(t->t+dt)").localVector();
  PetscScalar *dispTIncrArray;
  assert(dispTIncrSection);assert(dispTIncrVec);

  PetscSection dispTIncrAdjSection = fields->get("dispIncr adjust").petscSection();
  Vec          dispTIncrAdjVec     = fields->get("dispIncr adjust").localVector();
  PetscScalar *dispTIncrAdjArray;
  assert(dispTIncrAdjSection);assert(dispTIncrAdjVec);

  PetscSection jacobianSection = jacobian.petscSection();
  Vec          jacobianVec     = jacobian.localVector();
  PetscScalar *jacobianArray;
  assert(jacobianSection);assert(jacobianVec);

  DM           residualDM      = fields->get("residual").dmMesh();
  PetscSection residualSection = fields->get("residual").petscSection();
  Vec          residualVec     = fields->get("residual").localVector();
  PetscSection residualGlobalSection;
  PetscScalar *residualArray;
  assert(residualSection);assert(residualVec);
  err = DMGetDefaultGlobalSection(residualDM, &residualGlobalSection);CHECK_PETSC_ERROR(err);

  constrainSolnSpace_fn_type constrainSolnSpaceFn;
  switch (spaceDim) { // switch
  case 1:
    constrainSolnSpaceFn = 
      &pylith::faults::FaultCohesiveDyn::_constrainSolnSpace1D;
    break;
  case 2: 
    constrainSolnSpaceFn = 
      &pylith::faults::FaultCohesiveDyn::_constrainSolnSpace2D;
    break;
  case 3:
    constrainSolnSpaceFn = 
      &pylith::faults::FaultCohesiveDyn::_constrainSolnSpace3D;
    break;
  default :
    assert(0);
    throw std::logic_error("Unknown spatial dimension in "
			   "FaultCohesiveDyn::adjustSolnLumped.");
  } // switch

  _logger->eventEnd(setupEvent);

#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventBegin(computeEvent);
#endif

  const int numVertices = _cohesiveVertices.size();
  err = VecGetArray(dispRelVec, &dispRelArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(velRelVec, &velRelArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(dispTVec, &dispTArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(dispTIncrVec, &dispTIncrArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(dispTIncrAdjVec, &dispTIncrAdjArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(areaVec, &areaArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(orientationVec, &orientationArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(residualVec, &residualArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(jacobianVec, &jacobianArray);CHECK_PETSC_ERROR(err);
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_lagrange = _cohesiveVertices[iVertex].lagrange;
    const int v_fault = _cohesiveVertices[iVertex].fault;
    const int v_negative = _cohesiveVertices[iVertex].negative;
    const int v_positive = _cohesiveVertices[iVertex].positive;

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
    PetscScalar area;

    err = PetscSectionGetDof(areaSection, v_fault, &adof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(areaSection, v_fault, &aoff);CHECK_PETSC_ERROR(err);
    assert(1 == adof);
    area = areaArray[aoff];
    assert(area > 0.0);

    // Get disp(t) at Lagrange vertex.
    PetscInt dtldof, dtloff;

    err = PetscSectionGetDof(dispTSection, v_lagrange, &dtldof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dispTSection, v_lagrange, &dtloff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == dtldof);

    // Get dispIncr(t) at cohesive cell's vertices.
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

    // Get relative displacement at fault vertex.
    PetscInt drdof, droff;

    err = PetscSectionGetDof(dispRelSection, v_fault, &drdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dispRelSection, v_fault, &droff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == drdof);

    // Get relative velocity at fault vertex.
    PetscInt vrdof, vroff;

    err = PetscSectionGetDof(velRelSection, v_fault, &vrdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(velRelSection, v_fault, &vroff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == vrdof);
    
    // Get fault orientation at fault vertex.
    PetscInt odof, ooff;

    err = PetscSectionGetDof(orientationSection, v_fault, &odof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(orientationSection, v_fault, &ooff);CHECK_PETSC_ERROR(err);
    assert(spaceDim*spaceDim == odof);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(restrictEvent);
    _logger->eventBegin(computeEvent);
#endif

    // Adjust solution as in prescribed rupture, updating the Lagrange
    // multipliers and the corresponding displacment increments.
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      assert(jacobianArray[jpoff+iDim] > 0.0);
      assert(jacobianArray[jnoff+iDim] > 0.0);
      const PylithScalar S = (1.0/jacobianArray[jpoff+iDim] + 1.0/jacobianArray[jnoff+iDim]) * area*area;
      assert(S > 0.0);
      lagrangeTIncrVertex[iDim] = 1.0/S * (-residualArray[rloff+iDim] + area * (dispTIncrArray[dipoff+iDim] - dispTIncrArray[dinoff+iDim]));
      dispIncrVertexN[iDim] =  area / jacobianArray[jnoff+iDim]*lagrangeTIncrVertex[iDim];
      dispIncrVertexP[iDim] = -area / jacobianArray[jpoff+iDim]*lagrangeTIncrVertex[iDim];
    } // for

    // Compute slip, slip rate, and Lagrange multiplier at time t+dt
    // in fault coordinate system.
    slipVertex = 0.0;
    slipRateVertex = 0.0;
    tractionTpdtVertex = 0.0;
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      for (int jDim=0; jDim < spaceDim; ++jDim) {
        slipVertex[iDim] += orientationArray[ooff+iDim*spaceDim+jDim] * dispRelArray[droff+jDim];
        slipRateVertex[iDim] += orientationArray[ooff+iDim*spaceDim+jDim] * velRelArray[vroff+jDim];
        tractionTpdtVertex[iDim] += orientationArray[ooff+iDim*spaceDim+jDim] * (dispTArray[dtloff+jDim] + lagrangeTIncrVertex[jDim]);
      } // for
    } // for
    
    // Get friction properties and state variables.
    _friction->retrievePropsStateVars(v_fault);

    // Use fault constitutive model to compute traction associated with
    // friction.
    dTractionTpdtVertex = 0.0;
    const bool iterating = false; // No iteration for friction in lumped soln
    CALL_MEMBER_FN(*this, constrainSolnSpaceFn)(&dTractionTpdtVertex, t, slipVertex, slipRateVertex, tractionTpdtVertex, iterating);

    // Rotate traction back to global coordinate system.
    dLagrangeTpdtVertex = 0.0;
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      for (int jDim=0; jDim < spaceDim; ++jDim) {
        dLagrangeTpdtVertex[iDim] += orientationArray[ooff+jDim*spaceDim+iDim] * dTractionTpdtVertex[jDim];
      } // for
    } // for

#if 0 // debugging
    std::cout << "dispIncrP: ";
    for (int iDim=0; iDim < spaceDim; ++iDim)
      std::cout << "  " << dispIncrVertexP[iDim];
    std::cout << ", dispIncrN: "; 
    for (int iDim=0; iDim < spaceDim; ++iDim)
      std::cout << "  " << dispIncrVertexN[iDim];
    std::cout << ", slipVertex: ";
    for (int iDim=0; iDim < spaceDim; ++iDim)
      std::cout << "  " << slipVertex[iDim];
    std::cout << ", slipRateVertex: ";
    for (int iDim=0; iDim < spaceDim; ++iDim)
      std::cout << "  " << slipRateVertex[iDim];
    std::cout << ", orientationVertex: ";
    for (int iDim=0; iDim < spaceDim*spaceDim; ++iDim)
      std::cout << "  " << orientationArray[ooff+iDim];
    std::cout << ", tractionVertex: ";
    for (int iDim=0; iDim < spaceDim; ++iDim)
      std::cout << "  " << tractionTpdtVertex[iDim];
    std::cout << ", lagrangeTVertex: ";
    for (int iDim=0; iDim < spaceDim; ++iDim)
      std::cout << "  " << lagrangeTVertex[iDim];
    std::cout << ", lagrangeTIncrVertex: ";
    for (int iDim=0; iDim < spaceDim; ++iDim)
      std::cout << "  " << lagrangeTIncrVertex[iDim];
    std::cout << ", dTractionTpdtVertex: ";
    for (int iDim=0; iDim < spaceDim; ++iDim)
      std::cout << "  " << dTractionTpdtVertex[iDim];
    std::cout << ", dLagrangeTpdtVertex: ";
    for (int iDim=0; iDim < spaceDim; ++iDim)
      std::cout << "  " << dLagrangeTpdtVertex[iDim];
    std::cout << std::endl;
#endif

    // Compute change in displacement.
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      assert(jacobianArray[jpoff+iDim] > 0.0);
      assert(jacobianArray[jnoff+iDim] > 0.0);

      dispIncrVertexN[iDim] += area * dLagrangeTpdtVertex[iDim] / jacobianArray[jnoff+iDim];
      dispIncrVertexP[iDim] -= area * dLagrangeTpdtVertex[iDim] / jacobianArray[jpoff+iDim];

      // Update increment in Lagrange multiplier.
      lagrangeTIncrVertex[iDim] += dLagrangeTpdtVertex[iDim];
    } // for

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(computeEvent);
    _logger->eventBegin(updateEvent);
#endif

    // Compute contribution to adjusting solution only if Lagrange
    // constraint is local (the adjustment is assembled across processors).
    PetscInt goff;
    err = PetscSectionGetOffset(residualGlobalSection, v_lagrange, &goff);CHECK_PETSC_ERROR(err);
    if (goff >= 0) {
      // Adjust displacements to account for Lagrange multiplier values
      // (assumed to be zero in preliminary solve).
      // Update displacement field
      PetscInt diandof, dianoff;

      err = PetscSectionGetDof(dispTIncrAdjSection, v_negative, &diandof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(dispTIncrAdjSection, v_negative, &dianoff);CHECK_PETSC_ERROR(err);
      assert(spaceDim == diandof);
      for(PetscInt d = 0; d < diandof; ++d) {
        dispTIncrAdjArray[dianoff+d] += dispIncrVertexN[d];
      }
      PetscInt diapdof, diapoff;

      err = PetscSectionGetDof(dispTIncrAdjSection, v_positive, &diapdof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(dispTIncrAdjSection, v_positive, &diapoff);CHECK_PETSC_ERROR(err);
      assert(spaceDim == diapdof);
      for(PetscInt d = 0; d < diapdof; ++d) {
        dispTIncrAdjArray[diapoff+d] += dispIncrVertexP[d];
      }
    } // if

    // The Lagrange multiplier and relative displacement are NOT
    // assembled across processors, so update even if Lagrange vertex
    // is not local.

    // Set Lagrange multiplier value. Value from preliminary solve is
    // bogus due to artificial diagonal entry in Jacobian of 1.0.
    for(PetscInt d = 0; d < dildof; ++d) {
      dispTIncrArray[diloff+d] = lagrangeTIncrVertex[d];
    }

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(updateEvent);
#endif
  } // for
  err = VecRestoreArray(dispRelVec, &dispRelArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(velRelVec, &velRelArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(dispTVec, &dispTArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(dispTIncrVec, &dispTIncrArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(dispTIncrAdjVec, &dispTIncrAdjArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(areaVec, &areaArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(orientationVec, &orientationArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(residualVec, &residualArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(jacobianVec, &jacobianArray);CHECK_PETSC_ERROR(err);
  PetscLogFlops(numVertices*spaceDim*(17 + // adjust solve
                                      9 + // updates
                                      spaceDim*9));

#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventEnd(computeEvent);
#endif

#if 0 // DEBUGGING
  //dLagrangeTpdtSection->view("AFTER dLagrange");
  //dispIncrSection->view("AFTER DISP INCR (t->t+dt)");
#endif
} // adjustSolnLumped

// ----------------------------------------------------------------------
// Get vertex field associated with integrator.
const pylith::topology::Field<pylith::topology::SubMesh>&
pylith::faults::FaultCohesiveDyn::vertexField(const char* name,
					      const topology::SolutionFields* fields)
{ // vertexField
  assert(_faultMesh);
  assert(_quadrature);
  assert(_normalizer);
  assert(_fields);
  assert(_friction);

  const int cohesiveDim = _faultMesh->dimension();
  const int spaceDim = _quadrature->spaceDim();

  const topology::Field<topology::SubMesh>& orientation = _fields->get("orientation");

  PylithScalar scale = 0.0;
  int fiberDim = 0;
  if (0 == strcasecmp("slip", name)) {
    const topology::Field<topology::SubMesh>& dispRel = _fields->get("relative disp");
    _allocateBufferVectorField();
    topology::Field<topology::SubMesh>& buffer =  _fields->get("buffer (vector)");
    buffer.copy(dispRel);
    buffer.label("slip");
    FaultCohesiveLagrange::globalToFault(&buffer, orientation);
    return buffer;
  } else if (0 == strcasecmp("slip_rate", name)) {
    const topology::Field<topology::SubMesh>& velRel = _fields->get("relative velocity");
    _allocateBufferVectorField();
    topology::Field<topology::SubMesh>& buffer = _fields->get("buffer (vector)");
    buffer.copy(velRel);
    buffer.label("slip_rate");
    FaultCohesiveLagrange::globalToFault(&buffer, orientation);
    return buffer;
  } else if (cohesiveDim > 0 && 0 == strcasecmp("strike_dir", name)) {
    PetscSection orientationSection = _fields->get("orientation").petscSection();
    Vec          orientationVec     = _fields->get("orientation").localVector();
    assert(orientationSection);assert(orientationVec);
    _allocateBufferVectorField();
    topology::Field<topology::SubMesh>& buffer = _fields->get("buffer (vector)");
    buffer.copy(orientationSection, 0, PETSC_DETERMINE, orientationVec);
    buffer.label("strike_dir");
    buffer.scale(1.0);
    return buffer;

  } else if (2 == cohesiveDim && 0 == strcasecmp("dip_dir", name)) {
    PetscSection orientationSection = _fields->get("orientation").petscSection();
    Vec          orientationVec     = _fields->get("orientation").localVector();
    assert(orientationSection);assert(orientationVec);
    _allocateBufferVectorField();
    topology::Field<topology::SubMesh>& buffer = _fields->get("buffer (vector)");
    buffer.copy(orientationSection, 1, PETSC_DETERMINE, orientationVec);
    buffer.label("dip_dir");
    buffer.scale(1.0);
    return buffer;
  } else if (0 == strcasecmp("normal_dir", name)) {
    PetscSection orientationSection = _fields->get("orientation").petscSection();
    Vec          orientationVec     = _fields->get("orientation").localVector();
    assert(orientationSection);assert(orientationVec);
    _allocateBufferVectorField();
    topology::Field<topology::SubMesh>& buffer = _fields->get("buffer (vector)");
    buffer.copy(orientationSection, cohesiveDim, PETSC_DETERMINE, orientationVec);
    buffer.label("normal_dir");
    buffer.scale(1.0);
    return buffer;
  } else if (0 == strcasecmp("traction", name)) {
    assert(fields);
    const topology::Field<topology::Mesh>& dispT = fields->get("disp(t)");
    _allocateBufferVectorField();
    topology::Field<topology::SubMesh>& buffer = _fields->get("buffer (vector)");
    _calcTractions(&buffer, dispT);
    return buffer;
  } else if (_friction->hasPropStateVar(name)) {
    return _friction->getField(name);
  } else if (_tractPerturbation && _tractPerturbation->hasParameter(name)) {
    const topology::Field<topology::SubMesh>& param = _tractPerturbation->vertexField(name, fields);
    if (param.vectorFieldType() == topology::FieldBase::VECTOR) {
      _allocateBufferVectorField();
      topology::Field<topology::SubMesh>& buffer = _fields->get("buffer (vector)");
      buffer.copy(param);
      FaultCohesiveLagrange::globalToFault(&buffer, orientation);
      return buffer;
    } else {
      return param;
    } // if/else

  } else {
    std::ostringstream msg;
    msg << "Request for unknown vertex field '" << name << "' for fault '"
        << label() << "'.";
    throw std::runtime_error(msg.str());
  } // else

  // Should never get here.
  throw std::logic_error("Unknown field in FaultCohesiveDyn::vertexField().");

  // Satisfy return values
  assert(_fields);
  const topology::Field<topology::SubMesh>& buffer = _fields->get("buffer (vector)");

  return buffer;
} // vertexField

// ----------------------------------------------------------------------
// Compute tractions on fault surface using solution.
void
pylith::faults::FaultCohesiveDyn::_calcTractions(
    topology::Field<topology::SubMesh>* tractions,
    const topology::Field<topology::Mesh>& dispT)
{ // _calcTractions
  assert(tractions);
  assert(_faultMesh);
  assert(_fields);
  assert(_normalizer);

  // Fiber dimension of tractions matches spatial dimension.
  const int spaceDim = _quadrature->spaceDim();
  PetscErrorCode err;

  // Get sections.
  PetscSection dispTSection = dispT.petscSection();
  Vec          dispTVec     = dispT.localVector();
  PetscScalar *dispTArray;
  assert(dispTSection);assert(dispTVec);

  PetscSection orientationSection = _fields->get("orientation").petscSection();
  Vec          orientationVec     = _fields->get("orientation").localVector();
  PetscScalar *orientationArray;
  assert(orientationSection);assert(orientationVec);

  // Allocate buffer for tractions field (if necessary).
  PetscSection tractionsSection = tractions->petscSection();
  Vec          tractionsVec;
  PetscScalar *tractionsArray;
  if (!tractionsSection) {
    ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
    logger.stagePush("FaultFields");

    const topology::Field<topology::SubMesh>& dispRel = _fields->get("relative disp");
    tractions->cloneSection(dispRel);

    logger.stagePop();
  } // if
  const PylithScalar pressureScale = _normalizer->pressureScale();
  tractions->label("traction");
  tractions->scale(pressureScale);
  tractions->zero();
  tractionsSection = tractions->petscSection();
  tractionsVec     = tractions->localVector();

  const int numVertices = _cohesiveVertices.size();
  err = VecGetArray(dispTVec, &dispTArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(orientationVec, &orientationArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(tractionsVec, &tractionsArray);CHECK_PETSC_ERROR(err);
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_lagrange = _cohesiveVertices[iVertex].lagrange;
    const int v_fault = _cohesiveVertices[iVertex].fault;
    PetscInt dtldof, dtloff;

    err = PetscSectionGetDof(dispTSection, v_lagrange, &dtldof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dispTSection, v_lagrange, &dtloff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == dtldof);
    PetscInt odof, ooff;

    err = PetscSectionGetDof(orientationSection, v_fault, &odof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(orientationSection, v_fault, &ooff);CHECK_PETSC_ERROR(err);
    assert(spaceDim*spaceDim == odof);
    PetscInt tdof, toff;

    err = PetscSectionGetDof(tractionsSection, v_fault, &tdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(tractionsSection, v_fault, &toff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == tdof);

    // Rotate tractions to fault coordinate system.
    for(PetscInt d = 0; d < spaceDim; ++d) {
      tractionsArray[toff+d] = 0.0;
      for(PetscInt e = 0; e < spaceDim; ++e) {
        tractionsArray[toff+d] += orientationArray[ooff+d*spaceDim+e] * dispTArray[dtloff+e];
      } // for
    } // for
  } // for

  PetscLogFlops(numVertices * (1 + spaceDim) );

#if 0 // DEBUGGING
  tractions->view("TRACTIONS");
#endif

} // _calcTractions

// ----------------------------------------------------------------------
// Update relative displacement and velocity (slip and slip rate)
// associated with Lagrange vertex k corresponding to diffential
// velocity between conventional vertices i and j.
void
pylith::faults::FaultCohesiveDyn::_updateRelMotion(const topology::SolutionFields& fields)
{ // _updateRelMotion
  assert(_fields);

  const int spaceDim = _quadrature->spaceDim();
  PetscErrorCode err;

  // Get section information
  PetscSection dispTSection = fields.get("disp(t)").petscSection();
  Vec          dispTVec     = fields.get("disp(t)").localVector();
  PetscScalar *dispTArray;
  assert(dispTSection);assert(dispTVec);

  PetscSection dispTIncrSection = fields.get("dispIncr(t->t+dt)").petscSection();
  Vec          dispTIncrVec     = fields.get("dispIncr(t->t+dt)").localVector();
  PetscScalar *dispTIncrArray;
  assert(dispTIncrSection);assert(dispTIncrVec);

  PetscSection dispRelSection = _fields->get("relative disp").petscSection();
  Vec          dispRelVec     = _fields->get("relative disp").localVector();
  PetscScalar *dispRelArray;
  assert(dispRelSection);assert(dispRelVec);

  PetscSection velocitySection = fields.get("velocity(t)").petscSection();
  Vec          velocityVec     = fields.get("velocity(t)").localVector();
  PetscScalar *velocityArray;
  assert(velocitySection);assert(velocityVec);

  PetscSection velRelSection = _fields->get("relative velocity").petscSection();
  Vec          velRelVec     = _fields->get("relative velocity").localVector();
  PetscScalar *velRelArray;
  assert(velRelSection);assert(velRelVec);

  const int numVertices = _cohesiveVertices.size();
  err = VecGetArray(dispTVec, &dispTArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(dispTIncrVec, &dispTIncrArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(dispRelVec, &dispRelArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(velRelVec, &velRelArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(velocityVec, &velocityArray);CHECK_PETSC_ERROR(err);
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_fault = _cohesiveVertices[iVertex].fault;
    const int v_negative = _cohesiveVertices[iVertex].negative;
    const int v_positive = _cohesiveVertices[iVertex].positive;

    // Get displacement values
    PetscInt dtndof, dtnoff;

    err = PetscSectionGetDof(dispTSection, v_negative, &dtndof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dispTSection, v_negative, &dtnoff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == dtndof);
    PetscInt dtpdof, dtpoff;

    err = PetscSectionGetDof(dispTSection, v_positive, &dtpdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dispTSection, v_positive, &dtpoff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == dtpdof);
    PetscInt dindof, dinoff;

    err = PetscSectionGetDof(dispTIncrSection, v_negative, &dindof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dispTIncrSection, v_negative, &dinoff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == dindof);
    PetscInt dipdof, dipoff;

    err = PetscSectionGetDof(dispTIncrSection, v_positive, &dipdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dispTIncrSection, v_positive, &dipoff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == dipdof);

    // Update relative displacement field.
    PetscInt drdof, droff;

    err = PetscSectionGetDof(dispRelSection, v_fault, &drdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dispRelSection, v_fault, &droff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == drdof);
    for(PetscInt d = 0; d < spaceDim; ++d) {
      const PylithScalar value = dispTArray[dtpoff+d] + dispTIncrArray[dipoff+d] - dispTArray[dtnoff+d] - dispTIncrArray[dinoff+d];
      dispRelArray[droff+d] = fabs(value) > _zeroTolerance ? value : 0.0;
    } // for

    // Get velocity values
    PetscInt vndof, vnoff;

    err = PetscSectionGetDof(velocitySection, v_negative, &vndof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(velocitySection, v_negative, &vnoff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == vndof);
    PetscInt vpdof, vpoff;

    err = PetscSectionGetDof(velocitySection, v_positive, &vpdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(velocitySection, v_positive, &vpoff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == vpdof);

    // Update relative velocity field.
    PetscInt vrdof, vroff;

    err = PetscSectionGetDof(velRelSection, v_fault, &vrdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(velRelSection, v_fault, &vroff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == vrdof);
    for(PetscInt d = 0; d < spaceDim; ++d) {
      const PylithScalar value = velocityArray[vpoff+d] - velocityArray[vnoff+d];
      velRelArray[vroff+d] = fabs(value) > _zeroTolerance ? value : 0.0;
    } // for
  } // for
  err = VecRestoreArray(dispTVec, &dispTArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(dispTIncrVec, &dispTIncrArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(dispRelVec, &dispRelArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(velRelVec, &velRelArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(velocityVec, &velocityArray);CHECK_PETSC_ERROR(err);

  PetscLogFlops(numVertices*spaceDim*spaceDim*4);
} // _updateRelMotion

// ----------------------------------------------------------------------
// Setup sensitivity problem to compute change in slip given change in Lagrange multipliers.
void
pylith::faults::FaultCohesiveDyn::_sensitivitySetup(const topology::Jacobian& jacobian)
{ // _sensitivitySetup
  assert(_fields);
  assert(_quadrature);

  const int spaceDim = _quadrature->spaceDim();

  // Setup fields involved in sensitivity solve.
  if (!_fields->hasField("sensitivity solution")) {
    _fields->add("sensitivity solution", "sensitivity_soln");
    topology::Field<topology::SubMesh>& solution = _fields->get("sensitivity solution");
    const topology::Field<topology::SubMesh>& dispRel = _fields->get("relative disp");
    solution.cloneSection(dispRel);
    solution.createScatter(solution.mesh());
  } // if
  const topology::Field<topology::SubMesh>& solution = _fields->get("sensitivity solution");

  if (!_fields->hasField("sensitivity residual")) {
    _fields->add("sensitivity residual", "sensitivity_residual");
    topology::Field<topology::SubMesh>& residual = _fields->get("sensitivity residual");
    residual.cloneSection(solution);
    residual.createScatter(solution.mesh());
  } // if

  if (!_fields->hasField("sensitivity relative disp")) {
    _fields->add("sensitivity relative disp", "sensitivity_relative_disp");
    topology::Field<topology::SubMesh>& dispRel = _fields->get("sensitivity relative disp");
    dispRel.cloneSection(solution);
  } // if
  topology::Field<topology::SubMesh>& dispRel = _fields->get("sensitivity relative disp");
  dispRel.zero();

  if (!_fields->hasField("sensitivity dLagrange")) {
    _fields->add("sensitivity dLagrange", "sensitivity_dlagrange");
    topology::Field<topology::SubMesh>& dLagrange = _fields->get("sensitivity dLagrange");
    dLagrange.cloneSection(solution);
  } // if
  topology::Field<topology::SubMesh>& dLagrange = _fields->get("sensitivity dLagrange");
  dLagrange.zero();

  // Setup Jacobian sparse matrix for sensitivity solve.
  if (0 == _jacobian)
    _jacobian = new topology::Jacobian(solution, jacobian.matrixType());
  assert(_jacobian);
  _jacobian->zero();

  // Setup PETSc KSP linear solver.
  if (0 == _ksp) {
    PetscErrorCode err = 0;
    err = KSPCreate(_faultMesh->comm(), &_ksp);CHECK_PETSC_ERROR(err);
    err = KSPSetInitialGuessNonzero(_ksp, PETSC_FALSE);CHECK_PETSC_ERROR(err);
    PylithScalar rtol = 0.0;
    PylithScalar atol = 0.0;
    PylithScalar dtol = 0.0;
    int maxIters = 0;
    err = KSPGetTolerances(_ksp, &rtol, &atol, &dtol, &maxIters);CHECK_PETSC_ERROR(err);
    rtol = 1.0e-3*_zeroTolerance;
    atol = 1.0e-5*_zeroTolerance;
    err = KSPSetTolerances(_ksp, rtol, atol, dtol, maxIters);CHECK_PETSC_ERROR(err);

    PC pc;
    err = KSPGetPC(_ksp, &pc);CHECK_PETSC_ERROR(err);
    err = PCSetType(pc, PCJACOBI);CHECK_PETSC_ERROR(err);
    err = KSPSetType(_ksp, KSPGMRES);CHECK_PETSC_ERROR(err);

    err = KSPAppendOptionsPrefix(_ksp, "friction_");CHECK_PETSC_ERROR(err);
    err = KSPSetFromOptions(_ksp);CHECK_PETSC_ERROR(err);
  } // if
} // _sensitivitySetup

// ----------------------------------------------------------------------
// Update the Jacobian values for the sensitivity solve.
void
pylith::faults::FaultCohesiveDyn::_sensitivityUpdateJacobian(const bool negativeSide,
                                                             const topology::Jacobian& jacobian,
                                                             const topology::SolutionFields& fields)
{ // _sensitivityUpdateJacobian
  assert(_quadrature);
  assert(_fields);

  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int subnrows = numBasis*spaceDim;
  const int submatrixSize = subnrows * subnrows;

  // Get solution field
  const topology::Field<topology::Mesh>& solutionDomain = fields.solution();
  DM           solutionDomainDM      = solutionDomain.dmMesh();
  PetscSection solutionDomainSection = solutionDomain.petscSection();
  Vec          solutionDomainVec     = solutionDomain.localVector();
  PetscSection solutionDomainGlobalSection;
  PetscScalar *solutionDomainArray;
  PetscErrorCode err;
  assert(solutionDomainSection);assert(solutionDomainVec);
  err = DMGetDefaultGlobalSection(solutionDomainDM, &solutionDomainGlobalSection);CHECK_PETSC_ERROR(err);

  // Get cohesive cells
  DM              dmMesh = fields.mesh().dmMesh();
  IS              cellsCohesiveIS;
  const PetscInt *cellsCohesive;
  PetscInt        numCohesiveCells, vStart, vEnd;

  assert(dmMesh);
  err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);
  err = DMPlexGetStratumIS(dmMesh, "material-id", id(), &cellsCohesiveIS);CHECK_PETSC_ERROR(err);
  err = ISGetLocalSize(cellsCohesiveIS, &numCohesiveCells);CHECK_PETSC_ERROR(err);
  err = ISGetIndices(cellsCohesiveIS, &cellsCohesive);CHECK_PETSC_ERROR(err);

  // Visitor for Jacobian matrix associated with domain.
  scalar_array jacobianSubCell(submatrixSize);
  const PetscMat jacobianDomainMatrix = jacobian.matrix();
  assert(jacobianDomainMatrix);

  // Get fault Sieve mesh
  DM faultDMMesh = _faultMesh->dmMesh();
  assert(faultDMMesh);

  // Get sensitivity solution field
  DM           solutionFaultDM      = _fields->get("sensitivity solution").dmMesh();
  PetscSection solutionFaultSection = _fields->get("sensitivity solution").petscSection();
  Vec          solutionFaultVec     = _fields->get("sensitivity solution").localVector();
  PetscSection solutionFaultGlobalSection;
  PetscScalar *solutionFaultArray;
  assert(solutionFaultSection);assert(solutionFaultVec);
  err = DMGetDefaultGlobalSection(solutionFaultDM, &solutionFaultGlobalSection);CHECK_PETSC_ERROR(err);

  // Visitor for Jacobian matrix associated with fault.
  assert(_jacobian);
  const PetscMat jacobianFaultMatrix = _jacobian->matrix();
  assert(jacobianFaultMatrix);

  const int iCone = (negativeSide) ? 0 : 1;

  IS *cellsIS = (numCohesiveCells > 0) ? new IS[numCohesiveCells] : 0;
  int_array indicesGlobal(subnrows);
  int_array indicesLocal(numCohesiveCells*subnrows);
  int_array indicesPerm(subnrows);
  for(PetscInt c = 0; c < numCohesiveCells; ++c) {
    // Get cone for cohesive cell
    PetscInt          *closure = PETSC_NULL;
    PetscInt           closureSize, q = 0;
    err = DMPlexGetTransitiveClosure(dmMesh, cellsCohesive[c], PETSC_TRUE, &closureSize, &closure);CHECK_PETSC_ERROR(err);
    // Filter out non-vertices
    for(PetscInt p = 0; p < closureSize*2; p += 2) {
      if ((closure[p] >= vStart) && (closure[p] < vEnd)) {
        closure[q] = closure[p];
        ++q;
      }
    }
    closureSize = q;
    assert(closureSize == 3*numBasis);

    // Get indices
    for(int iBasis = 0; iBasis < numBasis; ++iBasis) {
      // negative side of the fault: iCone=0
      // positive side of the fault: iCone=1
      const int v_domain = closure[iCone*numBasis+iBasis];
      PetscInt goff;

      err = PetscSectionGetOffset(solutionDomainGlobalSection, v_domain, &goff);CHECK_PETSC_ERROR(err);
      for(int iDim = 0, iB = iBasis*spaceDim, gind = goff < 0 ? -(goff+1) : goff; iDim < spaceDim; ++iDim) {
        indicesGlobal[iB+iDim] = gind + iDim;
      } // for
    } // for
    err = DMPlexRestoreTransitiveClosure(dmMesh, cellsCohesive[c], PETSC_TRUE, &closureSize, &closure);CHECK_PETSC_ERROR(err);

    for (int i=0; i < subnrows; ++i) {
      indicesPerm[i]  = i;
    } // for
    err = PetscSortIntWithArray(indicesGlobal.size(), &indicesGlobal[0], &indicesPerm[0]);CHECK_PETSC_ERROR(err);

    for (int i=0; i < subnrows; ++i) {
      indicesLocal[c*subnrows+indicesPerm[i]] = i;
    } // for
    cellsIS[c] = PETSC_NULL;
    err = ISCreateGeneral(PETSC_COMM_SELF, indicesGlobal.size(), &indicesGlobal[0], PETSC_COPY_VALUES, &cellsIS[c]);CHECK_PETSC_ERROR(err);
  } // for

  PetscMat* submatrices = 0;
  err = MatGetSubMatrices(jacobianDomainMatrix, numCohesiveCells, cellsIS, cellsIS, MAT_INITIAL_MATRIX, &submatrices);CHECK_PETSC_ERROR(err);

  for(PetscInt c = 0; c < numCohesiveCells; ++c) {
    // Get values for submatrix associated with cohesive cell
    jacobianSubCell = 0.0;
    err = MatGetValues(submatrices[c], subnrows, &indicesLocal[c*subnrows], subnrows, &indicesLocal[c*subnrows],
                       &jacobianSubCell[0]);CHECK_PETSC_ERROR_MSG(err, "Restrict from PETSc Mat failed.");

    // Insert cell contribution into PETSc Matrix
    PetscInt c_fault = _cohesiveToFault[cellsCohesive[c]];
    err = DMPlexMatSetClosure(faultDMMesh, solutionFaultSection, solutionFaultGlobalSection,  jacobianFaultMatrix, c_fault, &jacobianSubCell[0], INSERT_VALUES);
    CHECK_PETSC_ERROR_MSG(err, "Update to PETSc Mat failed.");

    // Destory IS for cohesiveCell
    err = ISDestroy(&cellsIS[c]);CHECK_PETSC_ERROR(err);
  } // for

  err = MatDestroyMatrices(numCohesiveCells, &submatrices);CHECK_PETSC_ERROR(err);
  delete[] cellsIS; cellsIS = 0;
  err = ISRestoreIndices(cellsCohesiveIS, &cellsCohesive);CHECK_PETSC_ERROR(err);
  err = ISDestroy(&cellsCohesiveIS);CHECK_PETSC_ERROR(err);

  _jacobian->assemble("final_assembly");

#if 0 // DEBUGGING
  std::cout << "DOMAIN JACOBIAN" << std::endl;
  jacobian.view();
  std::cout << "SENSITIVITY JACOBIAN" << std::endl;
  _jacobian->view();
#endif
} // _sensitivityUpdateJacobian

// ----------------------------------------------------------------------
// Reform residual for sensitivity problem.
void
pylith::faults::FaultCohesiveDyn::_sensitivityReformResidual(const bool negativeSide)
{ // _sensitivityReformResidual
  /** Compute residual -L^T dLagrange
   *
   * Note: We need all entries for L, even those on other processors,
   * so we compute L rather than extract entries from the Jacobian.
   */

  const PylithScalar signFault = (negativeSide) ?  1.0 : -1.0;

  // Get cell information
  const int numQuadPts = _quadrature->numQuadPts();
  const scalar_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == numQuadPts);
  const int spaceDim = _quadrature->spaceDim();
  const int numBasis = _quadrature->numBasis();


  scalar_array basisProducts(numBasis*numBasis);

  // Get fault cell information
  DM             faultDMMesh = _faultMesh->dmMesh();
  PetscInt       cStart, cEnd;
  PetscErrorCode err;

  assert(faultDMMesh);
  err = DMPlexGetHeightStratum(faultDMMesh, 0, &cStart, &cEnd);CHECK_PETSC_ERROR(err);
  const int numCells = cEnd-cStart;

  // Get sections
  scalar_array coordinatesCell(numBasis*spaceDim);
  PetscSection coordSection;
  Vec          coordVec;
  err = DMPlexGetCoordinateSection(faultDMMesh, &coordSection);CHECK_PETSC_ERROR(err);
  err = DMGetCoordinatesLocal(faultDMMesh, &coordVec);CHECK_PETSC_ERROR(err);

  PetscSection dLagrangeSection = _fields->get("sensitivity dLagrange").petscSection();
  Vec          dLagrangeVec     = _fields->get("sensitivity dLagrange").localVector();
  assert(dLagrangeSection);assert(dLagrangeVec);

  scalar_array residualCell(numBasis*spaceDim);
  topology::Field<topology::SubMesh>& residual = _fields->get("sensitivity residual");
  PetscSection residualSection = residual.petscSection();
  Vec          residualVec     = residual.localVector();
  assert(residualSection);assert(residualVec);
  residual.zero();

  // Loop over cells
  for(PetscInt c = cStart; c < cEnd; ++c) {
    // Compute geometry
    const PetscScalar *coords = PETSC_NULL;
    PetscInt           coordsSize;
    err = DMPlexVecGetClosure(faultDMMesh, coordSection, coordVec, c, &coordsSize, &coords);CHECK_PETSC_ERROR(err);
    for(PetscInt i = 0; i < coordsSize; ++i) {coordinatesCell[i] = coords[i];}
    _quadrature->computeGeometry(coordinatesCell, c);
    err = DMPlexVecRestoreClosure(faultDMMesh, coordSection, coordVec, c, &coordsSize, &coords);CHECK_PETSC_ERROR(err);

    // Restrict input fields to cell
    const PetscScalar *dLagrangeArray = PETSC_NULL;
    PetscInt           dLagrangeSize;
    err = DMPlexVecGetClosure(faultDMMesh, dLagrangeSection, dLagrangeVec, c, &dLagrangeSize, &dLagrangeArray);CHECK_PETSC_ERROR(err);

    // Get cell geometry information that depends on cell
    const scalar_array& basis = _quadrature->basis();
    const scalar_array& jacobianDet = _quadrature->jacobianDet();

    // Compute product of basis functions.
    // Want values summed over quadrature points
    basisProducts = 0.0;
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      const PylithScalar wt = quadWts[iQuad] * jacobianDet[iQuad];

      for(int iBasis = 0, iQ=iQuad*numBasis; iBasis < numBasis; ++iBasis) {
        const PylithScalar valI = wt*basis[iQ+iBasis];
	
        for(int jBasis = 0; jBasis < numBasis; ++jBasis) {
	  
          basisProducts[iBasis*numBasis+jBasis] += valI*basis[iQ+jBasis];
        } // for
      } // for
    } // for

    residualCell = 0.0;
    
    for(int iBasis = 0; iBasis < numBasis; ++iBasis) {
      for(int jBasis = 0; jBasis < numBasis; ++jBasis) {
        const PylithScalar l = signFault * basisProducts[iBasis*numBasis+jBasis];
        for(PetscInt d = 0; d < spaceDim; ++d) {
          residualCell[iBasis*spaceDim+d] += l * dLagrangeArray[jBasis*spaceDim+d];
        } // for
      } // for
    } // for
    err = DMPlexVecRestoreClosure(faultDMMesh, dLagrangeSection, dLagrangeVec, c, &dLagrangeSize, &dLagrangeArray);CHECK_PETSC_ERROR(err);

    // Assemble cell contribution into field
    err = DMPlexVecSetClosure(faultDMMesh, residualSection, residualVec, c, &residualCell[0], ADD_VALUES);CHECK_PETSC_ERROR(err);
  } // for
} // _sensitivityReformResidual

// ----------------------------------------------------------------------
// Solve sensitivity problem.
void
pylith::faults::FaultCohesiveDyn::_sensitivitySolve(void)
{ // _sensitivitySolve
  assert(_fields);
  assert(_jacobian);
  assert(_ksp);

  topology::Field<topology::SubMesh>& residual = _fields->get("sensitivity residual");
  topology::Field<topology::SubMesh>& solution = _fields->get("sensitivity solution");

  // Assemble residual over processors.
  residual.complete();

  // Update PetscVector view of field.
  residual.scatterSectionToVector();

  PetscErrorCode err = 0;
  const PetscMat jacobianMat = _jacobian->matrix();
  err = KSPSetOperators(_ksp, jacobianMat, jacobianMat,
    DIFFERENT_NONZERO_PATTERN); CHECK_PETSC_ERROR(err);

  const PetscVec residualVec = residual.globalVector();
  const PetscVec solutionVec = solution.globalVector();
  err = KSPSolve(_ksp, residualVec, solutionVec); CHECK_PETSC_ERROR(err);

  // Update section view of field.
  solution.scatterVectorToSection();

#if 0 // DEBUGGING
  residual.view("SENSITIVITY RESIDUAL");
  solution.view("SENSITIVITY SOLUTION");
#endif
} // _sensitivitySolve

// ----------------------------------------------------------------------
// Update the relative displacement field values based on the
// sensitivity solve.
void
pylith::faults::FaultCohesiveDyn::_sensitivityUpdateSoln(const bool negativeSide)
{ // _sensitivityUpdateSoln
  assert(_fields);
  assert(_quadrature);

  const int spaceDim = _quadrature->spaceDim();
  PetscErrorCode err;

  scalar_array dispVertex(spaceDim);
  PetscSection solutionSection = _fields->get("sensitivity solution").petscSection();
  Vec          solutionVec     = _fields->get("sensitivity solution").localVector();
  PetscScalar *solutionArray;
  assert(solutionSection);assert(solutionVec);
  PetscSection dispRelSection = _fields->get("sensitivity relative disp").petscSection();
  Vec          dispRelVec     = _fields->get("sensitivity relative disp").localVector();
  PetscScalar *dispRelArray;
  assert(dispRelSection);assert(dispRelVec);
  PetscSection dLagrangeTpdtSection = _fields->get("sensitivity dLagrange").petscSection();
  Vec          dLagrangeTpdtVec     = _fields->get("sensitivity dLagrange").localVector();
  PetscScalar *dLagrangeTpdtArray;
  assert(dLagrangeTpdtSection);assert(dLagrangeTpdtVec);

  const PylithScalar sign = (negativeSide) ? -1.0 : 1.0;

  const int numVertices = _cohesiveVertices.size();
  err = VecGetArray(dispRelVec, &dispRelArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(solutionVec, &solutionArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(dLagrangeTpdtVec, &dLagrangeTpdtArray);CHECK_PETSC_ERROR(err);
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_fault = _cohesiveVertices[iVertex].fault;
    PetscInt sdof, soff;

    err = PetscSectionGetDof(solutionSection, v_fault, &sdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(solutionSection, v_fault, &soff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == sdof);

    // If no change in the Lagrange multiplier computed from friction criterion, there are no updates, so continue.
    PetscInt dldof, dloff;

    err = PetscSectionGetDof(dLagrangeTpdtSection, v_fault, &dldof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dLagrangeTpdtSection, v_fault, &dloff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == dldof);
    PylithScalar dLagrangeTpdtVertexMag = 0.0;
    for(PetscInt d = 0; d < spaceDim; ++d) {
      dLagrangeTpdtVertexMag += dLagrangeTpdtArray[dloff+d]*dLagrangeTpdtArray[dloff+d];
    } // for
    if (0.0 == dLagrangeTpdtVertexMag) continue;

    // Update relative displacements associated with sensitivity solve solution
    PetscInt drdof, droff;

    err = PetscSectionGetDof(dispRelSection, v_fault, &drdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dispRelSection, v_fault, &droff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == drdof);
    for(PetscInt d = 0; d < drdof; ++d) {
      dispRelArray[droff+d] += sign*solutionArray[soff+d];
    }
  } // for
  err = VecRestoreArray(dispRelVec, &dispRelArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(solutionVec, &solutionArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(dLagrangeTpdtVec, &dLagrangeTpdtArray);CHECK_PETSC_ERROR(err);
} // _sensitivityUpdateSoln


// ----------------------------------------------------------------------
// Compute norm of residual associated with matching fault
// constitutive model using update from sensitivity solve. We use
// this in a line search to find a good update (required because
// fault constitutive model may have a complex nonlinear feedback
// with deformation).
PylithScalar
pylith::faults::FaultCohesiveDyn::_constrainSolnSpaceNorm(const PylithScalar alpha,
							  const PylithScalar t,
							  topology::SolutionFields* const fields)
{ // _constrainSolnSpaceNorm
  /// Member prototype for _constrainSolnSpaceXD()
  typedef void (pylith::faults::FaultCohesiveDyn::*constrainSolnSpace_fn_type)
    (scalar_array*,
     const PylithScalar,
     const scalar_array&,
     const scalar_array&,
     const scalar_array&,
     const bool);

  // Update time step in friction (can vary).
  _friction->timeStep(_dt);
  const PylithScalar dt = _dt;

  const int spaceDim = _quadrature->spaceDim();
  const int indexN = spaceDim - 1;
  PetscErrorCode err;

  constrainSolnSpace_fn_type constrainSolnSpaceFn;
  switch (spaceDim) { // switch
  case 1:
    constrainSolnSpaceFn = 
      &pylith::faults::FaultCohesiveDyn::_constrainSolnSpace1D;
    break;
  case 2: 
    constrainSolnSpaceFn = 
      &pylith::faults::FaultCohesiveDyn::_constrainSolnSpace2D;
    break;
  case 3:
    constrainSolnSpaceFn = 
      &pylith::faults::FaultCohesiveDyn::_constrainSolnSpace3D;
    break;
  default :
    assert(0);
    throw std::logic_error("Unknown spatial dimension in "
			   "FaultCohesiveDyn::constrainSolnSpace().");
  } // switch

  // Get sections
  scalar_array slipTpdtVertex(spaceDim); // fault coordinates
  scalar_array slipRateVertex(spaceDim); // fault coordinates
  scalar_array tractionTpdtVertex(spaceDim); // fault coordinates
  scalar_array tractionMisfitVertex(spaceDim); // fault coordinates

  PetscSection orientationSection = _fields->get("orientation").petscSection();
  Vec          orientationVec     = _fields->get("orientation").localVector();
  PetscScalar *orientationArray;

  PetscSection dLagrangeTpdtSection = _fields->get("sensitivity dLagrange").petscSection();
  Vec          dLagrangeTpdtVec     = _fields->get("sensitivity dLagrange").localVector();
  PetscScalar *dLagrangeTpdtArray;
  assert(dLagrangeTpdtSection);assert(dLagrangeTpdtVec);

  PetscSection sensDispRelSection = _fields->get("sensitivity relative disp").petscSection();
  Vec          sensDispRelVec     = _fields->get("sensitivity relative disp").localVector();
  PetscScalar *sensDispRelArray;

  PetscSection dispTSection = fields->get("disp(t)").petscSection();
  Vec          dispTVec     = fields->get("disp(t)").localVector();
  PetscScalar *dispTArray;
  assert(dispTSection);assert(dispTVec);

  DM           dispTIncrDM      = fields->get("dispIncr(t->t+dt)").dmMesh();
  PetscSection dispTIncrSection = fields->get("dispIncr(t->t+dt)").petscSection();
  Vec          dispTIncrVec     = fields->get("dispIncr(t->t+dt)").localVector();
  PetscSection dispTIncrGlobalSection;
  PetscScalar *dispTIncrArray;
  assert(dispTIncrSection);assert(dispTIncrVec);
  err = DMGetDefaultGlobalSection(dispTIncrDM, &dispTIncrGlobalSection);CHECK_PETSC_ERROR(err);

  bool isOpening = false;
  PylithScalar norm2 = 0.0;
  int numVertices = _cohesiveVertices.size();
  err = VecGetArray(orientationVec, &orientationArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(dispTVec, &dispTArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(dispTIncrVec, &dispTIncrArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(sensDispRelVec, &sensDispRelArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(dLagrangeTpdtVec, &dLagrangeTpdtArray);CHECK_PETSC_ERROR(err);
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_lagrange = _cohesiveVertices[iVertex].lagrange;
    const int v_fault = _cohesiveVertices[iVertex].fault;
    const int v_negative = _cohesiveVertices[iVertex].negative;
    const int v_positive = _cohesiveVertices[iVertex].positive;
    PetscInt goff;

    // Compute contribution only if Lagrange constraint is local.
    err = PetscSectionGetOffset(dispTIncrGlobalSection, v_lagrange, &goff);CHECK_PETSC_ERROR(err);
    if (goff < 0) continue;

    // Get displacement values
    PetscInt dtndof, dtnoff;

    err = PetscSectionGetDof(dispTSection, v_negative, &dtndof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dispTSection, v_negative, &dtnoff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == dtndof);
    PetscInt dtpdof, dtpoff;

    err = PetscSectionGetDof(dispTSection, v_positive, &dtpdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dispTSection, v_positive, &dtpoff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == dtpdof);

    // Get displacement increment values.
    PetscInt dindof, dinoff;

    err = PetscSectionGetDof(dispTIncrSection, v_negative, &dindof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dispTIncrSection, v_negative, &dinoff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == dindof);
    PetscInt dipdof, dipoff;

    err = PetscSectionGetDof(dispTIncrSection, v_positive, &dipdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dispTIncrSection, v_positive, &dipoff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == dipdof);

    // Get orientation
    PetscInt odof, ooff;

    err = PetscSectionGetDof(orientationSection, v_fault, &odof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(orientationSection, v_fault, &ooff);CHECK_PETSC_ERROR(err);
    assert(spaceDim*spaceDim == odof);

    // Get change in relative displacement from sensitivity solve.
    PetscInt sdrdof, sdroff;

    err = PetscSectionGetDof(sensDispRelSection, v_fault, &sdrdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(sensDispRelSection, v_fault, &sdroff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == sdrdof);

    // Get Lagrange multiplier values from disp(t), and dispIncr(t->t+dt)
    PetscInt dtldof, dtloff;

    err = PetscSectionGetDof(dispTSection, v_lagrange, &dtldof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dispTSection, v_lagrange, &dtloff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == dtldof);
    PetscInt dildof, diloff;

    err = PetscSectionGetDof(dispTIncrSection, v_lagrange, &dildof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dispTIncrSection, v_lagrange, &diloff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == dildof);
    PetscInt sdldof, sdloff;

    err = PetscSectionGetDof(dLagrangeTpdtSection, v_fault, &sdldof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dLagrangeTpdtSection, v_fault, &sdloff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == sdldof);

    // Compute slip, slip rate, and traction at time t+dt as part of
    // line search.
    slipTpdtVertex = 0.0;
    slipRateVertex = 0.0;
    tractionTpdtVertex = 0.0;
    for(PetscInt d = 0; d < spaceDim; ++d) {
      for(PetscInt e = 0; e < spaceDim; ++e) {
        slipTpdtVertex[d] += orientationArray[ooff+d*spaceDim+e] *
          (dispTArray[dtpoff+e] + dispTIncrArray[dipoff+e] - dispTArray[dtnoff+e] - dispTIncrArray[dinoff+e] + alpha*sensDispRelArray[sdroff+e]);
        slipRateVertex[d] += orientationArray[ooff+d*spaceDim+e] * (dispTIncrArray[dipoff+e] - dispTIncrArray[dinoff+e] + alpha*sensDispRelArray[sdroff+e]) / dt;
        tractionTpdtVertex[d] += orientationArray[ooff+d*spaceDim+e] * (dispTArray[dtloff+e] + dispTIncrArray[diloff+e] + alpha*dLagrangeTpdtArray[sdloff+e]);
      } // for
      if (fabs(slipRateVertex[d]) < _zeroTolerance) {
        slipRateVertex[d] = 0.0;
      } // if
    } // for
    if (fabs(slipTpdtVertex[indexN]) < _zeroTolerance) {
      slipTpdtVertex[indexN] = 0.0;
    } // if

    // FIRST, correct nonphysical trial solutions.
    // Order of steps a-c is important!

    if (slipTpdtVertex[indexN]*tractionTpdtVertex[indexN] < 0.0) {
      // Step a: Prevent nonphysical trial solutions. The product of the
      // normal traction and normal slip must be nonnegative (forbid
      // interpenetration with tension or opening with compression).
      
      // Don't know what behavior is appropriate so set smaller of
      // traction and slip to zero (should be appropriate if problem
      // is nondimensionalized correctly).
      if (fabs(slipTpdtVertex[indexN]) > fabs(tractionTpdtVertex[indexN])) {
        // fault opening is bigger, so force normal traction back to zero
        tractionTpdtVertex[indexN] = 0.0;
      } else {
        // traction is bigger, so force fault opening back to zero
        slipTpdtVertex[indexN] = 0.0;
      } // if/else

    } else if (slipTpdtVertex[indexN] > _zeroTolerance) {
      // Step b: Ensure fault traction is zero when opening (if
      // alpha=1 this should be enforced already, but will not be
      // properly enforced when alpha < 1).
      
      for(PetscInt d = 0; d < spaceDim; ++d) {
        tractionTpdtVertex[d] = 0.0;
      } // for
    } else if (slipTpdtVertex[indexN] < 0.0) {
      // Step c: Prevent interpenetration.

      slipTpdtVertex[indexN] = 0.0;
    } // if

    if (slipTpdtVertex[indexN] > _zeroTolerance) {
      isOpening = true;
    } // if

    // Apply friction criterion to trial solution to get change in
    // Lagrange multiplier (dLagrangeTpdtVertex) in fault coordinate
    // system.
    
    // Get friction properties and state variables.
    _friction->retrievePropsStateVars(v_fault);
    
    // Use fault constitutive model to compute traction associated with
    // friction.
    tractionMisfitVertex = 0.0;
    const bool iterating = true; // Iterating to get friction
    CALL_MEMBER_FN(*this, constrainSolnSpaceFn)(&tractionMisfitVertex, t,
                                                slipTpdtVertex, slipRateVertex, tractionTpdtVertex,
                                                iterating);

#if 0 // DEBUGGING
    std::cout << "alpha: " << alpha
	      << ", v_fault: " << v_fault;
    std::cout << ", misfit:";
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      std::cout << " " << tractionMisfitVertex[iDim];
    } // for
    std::cout << ", slip:";
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      std::cout << " " << slipTpdtVertex[iDim];
    } // for
    std::cout << ", traction:";
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      std::cout << " " << tractionTpdtVertex[iDim];
    } // for
    std::cout << ", dDispRel:";
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      std::cout << " " << sensDispRelArray[sdroff+iDim];
    } // for
    std::cout << std::endl;
#endif

    for(PetscInt d = 0; d < spaceDim; ++d) {
      norm2 += tractionMisfitVertex[d]*tractionMisfitVertex[d];
    } // for
  } // for
  err = VecRestoreArray(orientationVec, &orientationArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(dispTVec, &dispTArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(dispTIncrVec, &dispTIncrArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(sensDispRelVec, &sensDispRelArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(dLagrangeTpdtVec, &dLagrangeTpdtArray);CHECK_PETSC_ERROR(err);

  if (isOpening && alpha < 1.0) {
    norm2 = PYLITH_MAXFLOAT;
  } // if

  PetscScalar norm2Total = 0.0;
  PetscInt numVerticesTotal = 0;
  err = MPI_Allreduce(&norm2, &norm2Total, 1, MPIU_SCALAR, MPI_SUM, fields->mesh().comm());
  err = MPI_Allreduce(&numVertices, &numVerticesTotal, 1, MPIU_INT, MPI_SUM, fields->mesh().comm());

  assert(numVerticesTotal > 0);
  return sqrt(norm2Total) / numVerticesTotal;
} // _constrainSolnSpaceNorm


// ----------------------------------------------------------------------
// Constrain solution space in 1-D.
void
pylith::faults::FaultCohesiveDyn::_constrainSolnSpace1D(scalar_array* dTractionTpdt,
	 const PylithScalar t,
         const scalar_array& slip,
         const scalar_array& sliprate,
	 const scalar_array& tractionTpdt,
	 const bool iterating)
{ // _constrainSolnSpace1D
  assert(dTractionTpdt);

  if (fabs(slip[0]) < _zeroTolerance) {
    // if compression, then no changes to solution
  } else {
    // if tension, then traction is zero.
    
    const PylithScalar dlp = -tractionTpdt[0];
    (*dTractionTpdt)[0] = dlp;
  } // else
  
  PetscLogFlops(2);
} // _constrainSolnSpace1D

// ----------------------------------------------------------------------
// Constrain solution space in 2-D.
void
pylith::faults::FaultCohesiveDyn::_constrainSolnSpace2D(scalar_array* dTractionTpdt,
	 const PylithScalar t,
         const scalar_array& slip,
         const scalar_array& slipRate,
	 const scalar_array& tractionTpdt,
	 const bool iterating)
{ // _constrainSolnSpace2D
  assert(dTractionTpdt);

  const PylithScalar slipMag = fabs(slip[0]);
  const PylithScalar slipRateMag = fabs(slipRate[0]);

  const PylithScalar tractionNormal = tractionTpdt[1];
  const PylithScalar tractionShearMag = fabs(tractionTpdt[0]);

  if (fabs(slip[1]) < _zeroTolerance && tractionNormal < -_zeroTolerance) {
    // if in compression and no opening
    const PylithScalar frictionStress = 
      _friction->calcFriction(t, slipMag, slipRateMag, tractionNormal);
    if (tractionShearMag > frictionStress || (iterating && slipRateMag > 0.0)) {
      // traction is limited by friction, so have sliding OR
      // friction exceeds traction due to overshoot in slip

      if (tractionShearMag > 0.0) {
	// Update traction increment based on value required to stick
	// versus friction
	const PylithScalar dlp = -(tractionShearMag - frictionStress) *
	  tractionTpdt[0] / tractionShearMag;
	(*dTractionTpdt)[0] = dlp;
      } else {
	// No shear stress and no friction.
      } // if/else
    } else {
      // friction exceeds value necessary to stick
      // no changes to solution
      if (iterating) {
	assert(0.0 == slipRateMag);
      } // if
    } // if/else
  } else {
    // if in tension, then traction is zero.
    (*dTractionTpdt)[0] = -tractionTpdt[0];
    (*dTractionTpdt)[1] = -tractionTpdt[1];
  } // else

  PetscLogFlops(8);
} // _constrainSolnSpace2D

// ----------------------------------------------------------------------
// Constrain solution space in 3-D.
void
pylith::faults::FaultCohesiveDyn::_constrainSolnSpace3D(scalar_array* dTractionTpdt,
	 const PylithScalar t,
         const scalar_array& slip,
         const scalar_array& slipRate,
	 const scalar_array& tractionTpdt,
	 const bool iterating)
{ // _constrainSolnSpace3D
  assert(dTractionTpdt);

  const PylithScalar slipShearMag = sqrt(slip[0] * slip[0] +
             slip[1] * slip[1]);
  PylithScalar slipRateMag = sqrt(slipRate[0]*slipRate[0] + 
            slipRate[1]*slipRate[1]);
  
  const PylithScalar tractionNormal = tractionTpdt[2];
  const PylithScalar tractionShearMag = 
    sqrt(tractionTpdt[0] * tractionTpdt[0] +
	 tractionTpdt[1] * tractionTpdt[1]);
  
  if (fabs(slip[2]) < _zeroTolerance && tractionNormal < -_zeroTolerance) {
    // if in compression and no opening
    const PylithScalar frictionStress = 
      _friction->calcFriction(t, slipShearMag, slipRateMag, tractionNormal);

    if (tractionShearMag > frictionStress || (iterating && slipRateMag > 0.0)) {
      // traction is limited by friction, so have sliding OR
      // friction exceeds traction due to overshoot in slip
      
      if (tractionShearMag > 0.0) {
	// Update traction increment based on value required to stick
	// versus friction
	const PylithScalar dlp = -(tractionShearMag - frictionStress) * 
	  tractionTpdt[0] / tractionShearMag;
	const PylithScalar dlq = -(tractionShearMag - frictionStress) * 
	  tractionTpdt[1] / tractionShearMag;
	
	(*dTractionTpdt)[0] = dlp;
	(*dTractionTpdt)[1] = dlq;
      } else {
	// No shear stress and no friction.
      } // if/else	
      
    } else {
      // else friction exceeds value necessary, so stick
      // no changes to solution
      if (iterating) {
	assert(0.0 == slipRateMag);
      } // if
    } // if/else
  } else {
    // if in tension, then traction is zero.
    (*dTractionTpdt)[0] = -tractionTpdt[0];
    (*dTractionTpdt)[1] = -tractionTpdt[1];
    (*dTractionTpdt)[2] = -tractionTpdt[2];
  } // else

  PetscLogFlops(22);
} // _constrainSolnSpace3D


// End of file 
