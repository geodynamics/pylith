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

#include "FaultCohesiveDyn.hh" // implementation of object methods

#include "CohesiveTopology.hh" // USES CohesiveTopology
#include "TractPerturbation.hh" // HOLDSA TractPerturbation

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/Jacobian.hh" // USES Jacobian
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/topology/Stratum.hh" // USES Stratum, StratumIS

#include "pylith/friction/FrictionModel.hh" // USES FrictionModel
#include "pylith/problems/SolverLinear.hh" // USES SolverLinear

#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/feassemble/CellGeometry.hh" // USES CellGeometry

#include "pylith/utils/EventLogger.hh" // USES EventLogger
#include "pylith/utils/constdefs.h" // USES PYLITH_MAXFLOAT
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

#include <iostream> // TEMPORARY

//#define DETAILED_EVENT_LOGGING

// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::FaultCohesiveDyn::FaultCohesiveDyn(void) :
  _zeroTolerance(1.0e-10),
  _tractPerturbation(0),
  _friction(0),
  _jacobian(0),
  _ksp(0),
  _openFreeSurf(true)
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
  PYLITH_METHOD_BEGIN;

  FaultCohesiveLagrange::deallocate();

  _tractPerturbation = 0; // :TODO: Use shared pointer
  _friction = 0; // :TODO: Use shared pointer

  delete _jacobian; _jacobian = 0;
  PetscErrorCode err = KSPDestroy(&_ksp);PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_END;
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
  PYLITH_METHOD_BEGIN;

  assert(upDir);
  assert(_quadrature);
  assert(_normalizer);

  FaultCohesiveLagrange::initialize(mesh, upDir);

  // Get initial tractions using a spatial database.
  if (_tractPerturbation) {
    const topology::Field& orientation = _fields->get("orientation");
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

  // Create field for relative velocity associated with Lagrange vertex k
  _fields->add("relative velocity", "relative_velocity");
  topology::Field& velRel = _fields->get("relative velocity");
  topology::Field& dispRel = _fields->get("relative disp");
  velRel.cloneSection(dispRel);
  velRel.vectorFieldType(topology::FieldBase::VECTOR);
  velRel.scale(_normalizer->lengthScale() / _normalizer->timeScale());

  PYLITH_METHOD_END;
} // initialize

// ----------------------------------------------------------------------
// Integrate contributions to residual term (r) for operator.
void
pylith::faults::FaultCohesiveDyn::integrateResidual(const topology::Field& residual,
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
  PetscSection residualGlobalSection = residual.globalSection();assert(residualGlobalSection);

  topology::VecVisitorMesh residualVisitor(residual);
  PetscScalar* residualArray = residualVisitor.localArray();

  topology::Field& dispT = fields->get("disp(t)");
  topology::VecVisitorMesh dispTVisitor(dispT);
  const PetscScalar* dispTArray = dispTVisitor.localArray();

  topology::Field& dispTIncr = fields->get("dispIncr(t->t+dt)");
  topology::VecVisitorMesh dispTIncrVisitor(dispTIncr);
  const PetscScalar* dispTIncrArray = dispTIncrVisitor.localArray();

  scalar_array tractPerturbVertex(spaceDim);
  topology::VecVisitorMesh* tractionsVisitor = 0;
  PetscScalar *tractionsArray = NULL;
  if (_tractPerturbation) {
    _tractPerturbation->calculate(t);
    
    const topology::Fields* params = _tractPerturbation->parameterFields();assert(params);
    const topology::Field& tractions = params->get("value");

    tractionsVisitor = new topology::VecVisitorMesh(tractions);
    tractionsArray = tractionsVisitor->localArray();
  } // if

  topology::Field& area = _fields->get("area");
  topology::VecVisitorMesh areaVisitor(area);
  const PetscScalar* areaArray = areaVisitor.localArray();

  topology::Field& orientation = _fields->get("orientation");
  topology::VecVisitorMesh orientationVisitor(orientation);
  const PetscScalar* orientationArray = orientationVisitor.localArray();

  _logger->eventEnd(setupEvent);
#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventBegin(computeEvent);
#endif

  // Loop over fault vertices
  PetscErrorCode err = 0;
  const int numVertices = _cohesiveVertices.size();
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int e_lagrange = _cohesiveVertices[iVertex].lagrange;
    const int v_fault = _cohesiveVertices[iVertex].fault;
    const int v_negative = _cohesiveVertices[iVertex].negative;
    const int v_positive = _cohesiveVertices[iVertex].positive;

    // Skip clamped vertices
    if (e_lagrange < 0) {
      continue;
    } // if

    // Compute contribution only if Lagrange constraint is local.
    PetscInt goff = 0;
    err = PetscSectionGetOffset(residualGlobalSection, e_lagrange, &goff);PYLITH_CHECK_ERROR(err);
    if (goff < 0)
      continue;

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(restrictEvent);
#endif

    // Get prescribed traction perturbation at fault vertex.
    if (_tractPerturbation) {
      const PetscInt toff = tractionsVisitor->sectionOffset(v_fault);
      assert(spaceDim == tractionsVisitor->sectionDof(v_fault));
      for(PetscInt d = 0; d < spaceDim; ++d) {
        tractPerturbVertex[d] = tractionsArray[toff+d];
      } // for
    } else {
      tractPerturbVertex = 0.0;
    } // if/else

    // Get orientation associated with fault vertex.
    const PetscInt ooff = orientationVisitor.sectionOffset(v_fault);
    assert(spaceDim*spaceDim == orientationVisitor.sectionDof(v_fault));

    // Get area associated with fault vertex.
    const PetscInt aoff = areaVisitor.sectionOffset(v_fault);
    assert(1 == areaVisitor.sectionDof(v_fault));

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

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(restrictEvent);
    _logger->eventBegin(computeEvent);
#endif

    // Compute slip (in fault coordinates system) from displacements.
    PylithScalar slipNormal = 0.0;
    PylithScalar tractionNormal = 0.0;
    const PetscInt indexN = spaceDim - 1;
    for(PetscInt d = 0; d < spaceDim; ++d) {
      slipNormal += orientationArray[ooff+indexN*spaceDim+d] * (dispTArray[dtpoff+d] + dispTIncrArray[dipoff+d] - dispTArray[dtnoff+d] - dispTIncrArray[dinoff+d]);
      tractionNormal += orientationArray[ooff+indexN*spaceDim+d] * (dispTArray[dtloff+d] + dispTIncrArray[diloff+d]);
    } // for
    
#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(computeEvent);
    _logger->eventBegin(updateEvent);
#endif
    if (slipNormal < _zeroTolerance || !_openFreeSurf) { 
      // if no opening or flag indicates to still impose initial tractions when fault is open.
      // Assemble contributions into field
      const PetscInt rnoff = residualVisitor.sectionOffset(v_negative);
      assert(spaceDim == residualVisitor.sectionDof(v_negative));

      const PetscInt rpoff = residualVisitor.sectionOffset(v_positive);
      assert(spaceDim == residualVisitor.sectionDof(v_positive));

      // Initial (external) tractions oppose (internal) tractions associated with Lagrange multiplier.
      for(PetscInt d = 0; d < spaceDim; ++d) {
        residualArray[rnoff+d] +=  areaArray[aoff] * (dispTArray[dtloff+d] + dispTIncrArray[diloff+d] - tractPerturbVertex[d]);
        residualArray[rpoff+d] += -areaArray[aoff] * (dispTArray[dtloff+d] + dispTIncrArray[diloff+d] - tractPerturbVertex[d]);
      } // for
    } else { // opening, normal traction should be zero
      std::ostringstream msg;
      if (fabs(tractionNormal) > _zeroTolerance) {
        msg << "WARNING! Fault opening with nonzero traction."
                  << ", v_fault: " << v_fault
                  << ", opening: " << slipNormal
                  << ", normal traction: " << tractionNormal
                  << std::endl;
        throw std::runtime_error(msg.str());
      } // if
    }  // if/else

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(updateEvent);
#endif
  } // for
  PetscLogFlops(numVertices*spaceDim*8);
  delete tractionsVisitor; tractionsVisitor = 0;

#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventEnd(computeEvent);
#endif

  PYLITH_METHOD_END;
} // integrateResidual

// ----------------------------------------------------------------------
// Update state variables as needed.
void
pylith::faults::FaultCohesiveDyn::updateStateVars(const PylithScalar t,
						  topology::SolutionFields* const fields)
{ // updateStateVars
  PYLITH_METHOD_BEGIN;

  assert(fields);
  assert(_fields);

  _updateRelMotion(*fields);

  const int spaceDim = _quadrature->spaceDim();

  // Allocate arrays for vertex values
  scalar_array tractionTpdtVertex(spaceDim); // Fault coordinate system

  // Get fields.
  topology::Field& dispT = fields->get("disp(t)");
  topology::VecVisitorMesh dispTVisitor(dispT);
  const PetscScalar* dispTArray = dispTVisitor.localArray();

  topology::Field& dispTIncr = fields->get("dispIncr(t->t+dt)");
  topology::VecVisitorMesh dispTIncrVisitor(dispTIncr);
  const PetscScalar* dispTIncrArray = dispTIncrVisitor.localArray();

  scalar_array slipVertex(spaceDim);
  topology::Field& dispRel = _fields->get("relative disp");
  topology::VecVisitorMesh dispRelVisitor(dispRel);
  const PetscScalar* dispRelArray = dispRelVisitor.localArray();

  scalar_array slipRateVertex(spaceDim);
  topology::Field& velRel = _fields->get("relative velocity");
  topology::VecVisitorMesh velRelVisitor(velRel);
  const PetscScalar* velRelArray = velRelVisitor.localArray();

  topology::Field& orientation = _fields->get("orientation");
  topology::VecVisitorMesh orientationVisitor(orientation);
  const PetscScalar* orientationArray = orientationVisitor.localArray();

  const int numVertices = _cohesiveVertices.size();
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int e_lagrange = _cohesiveVertices[iVertex].lagrange;
    const int v_fault = _cohesiveVertices[iVertex].fault;

    // Skip clamped vertices
    if (e_lagrange < 0) {
      continue;
    } // if

    // Get relative displacement
    const PetscInt droff = dispRelVisitor.sectionOffset(v_fault);
    assert(spaceDim == dispRelVisitor.sectionDof(v_fault));

    // Get relative velocity
    const PetscInt vroff = velRelVisitor.sectionOffset(v_fault);
    assert(spaceDim == velRelVisitor.sectionDof(v_fault));

    // Get orientation
    const PetscInt ooff = orientationVisitor.sectionOffset(v_fault);
    assert(spaceDim*spaceDim == orientationVisitor.sectionDof(v_fault));

    // Get Lagrange multiplier values from disp(t), and dispIncr(t->t+dt)
    const PetscInt dtloff = dispTVisitor.sectionOffset(e_lagrange);
    assert(spaceDim == dispTVisitor.sectionDof(e_lagrange));

    const PetscInt diloff = dispTIncrVisitor.sectionOffset(e_lagrange);
    assert(spaceDim == dispTIncrVisitor.sectionDof(e_lagrange));

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

  PYLITH_METHOD_END;
} // updateStateVars

// ----------------------------------------------------------------------
// Constrain solution based on friction.
void
pylith::faults::FaultCohesiveDyn::constrainSolnSpace(topology::SolutionFields* const fields,
						     const PylithScalar t,
						     const topology::Jacobian& jacobian)
{ // constrainSolnSpace
  PYLITH_METHOD_BEGIN;

  /// Member prototype for _constrainSolnSpaceXD()
  typedef void (pylith::faults::FaultCohesiveDyn::*constrainSolnSpace_fn_type)
    (scalar_array*,
     const PylithScalar,
     const scalar_array&,
     const scalar_array&,
     const scalar_array&,
     const PylithScalar,
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

  // Get sections
  scalar_array slipTpdtVertex(spaceDim);
  scalar_array slipRateVertex(spaceDim);
  topology::VecVisitorMesh dispRelVisitor(_fields->get("relative disp"));

  topology::VecVisitorMesh orientationVisitor(_fields->get("orientation"));
  const PetscScalar* orientationArray = orientationVisitor.localArray();

  topology::VecVisitorMesh dispTVisitor(fields->get("disp(t)"));
  const PetscScalar* dispTArray = dispTVisitor.localArray();

  scalar_array dDispTIncrVertexN(spaceDim);
  scalar_array dDispTIncrVertexP(spaceDim);
  topology::VecVisitorMesh dispTIncrVisitor(fields->get("dispIncr(t->t+dt)"));
  const PetscScalar* dispTIncrArray = dispTIncrVisitor.localArray();

  PetscSection dispTIncrGlobalSection = fields->get("dispIncr(t->t+dt)").globalSection();assert(dispTIncrGlobalSection);

  topology::VecVisitorMesh dispTIncrAdjVisitor(fields->get("dispIncr adjust"));
  PetscScalar* dispTIncrAdjArray = dispTIncrAdjVisitor.localArray();

  scalar_array dTractionTpdtVertex(spaceDim);
  scalar_array dLagrangeTpdtVertex(spaceDim);
  topology::VecVisitorMesh dLagrangeVisitor(_fields->get("sensitivity dLagrange"));
  PetscScalar* dLagrangeArray = dLagrangeVisitor.localArray();

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

  const int numVertices = _cohesiveVertices.size();
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int e_lagrange = _cohesiveVertices[iVertex].lagrange;
    const int v_fault = _cohesiveVertices[iVertex].fault;
    const int v_negative = _cohesiveVertices[iVertex].negative;
    const int v_positive = _cohesiveVertices[iVertex].positive;

    // Skip clamped vertices
    if (e_lagrange < 0) {
      continue;
    } // if

    // Get displacement values
    const PetscInt dtnoff = dispTVisitor.sectionOffset(v_negative);
    assert(spaceDim == dispTVisitor.sectionDof(v_negative));

    const PetscInt dtpoff = dispTVisitor.sectionOffset(v_positive);
    assert(spaceDim == dispTVisitor.sectionDof(v_positive));

    const PetscInt dtloff = dispTVisitor.sectionOffset(e_lagrange);
    assert(spaceDim == dispTVisitor.sectionDof(e_lagrange));

    // Get displacement increment values.
    const PetscInt dinoff = dispTIncrVisitor.sectionOffset(v_negative);
    assert(spaceDim == dispTIncrVisitor.sectionDof(v_negative));

    const PetscInt dipoff = dispTIncrVisitor.sectionOffset(v_positive);
    assert(spaceDim == dispTIncrVisitor.sectionDof(v_positive));

    const PetscInt diloff = dispTIncrVisitor.sectionOffset(e_lagrange);
    assert(spaceDim == dispTIncrVisitor.sectionDof(e_lagrange));

    // Get orientation
    const PetscInt ooff = orientationVisitor.sectionOffset(v_fault);
    assert(spaceDim*spaceDim == orientationVisitor.sectionDof(v_fault));

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
    const PylithScalar jacobianShearVertex = 0.0;
    const bool iterating = true; // Iterating to get friction
    CALL_MEMBER_FN(*this, constrainSolnSpaceFn)(&dTractionTpdtVertex, t, slipTpdtVertex, slipRateVertex, tractionTpdtVertex, jacobianShearVertex, iterating);

    // Rotate increment in traction back to global coordinate system.
    dLagrangeTpdtVertex = 0.0;
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      for (int jDim=0; jDim < spaceDim; ++jDim) {
        dLagrangeTpdtVertex[iDim] += orientationArray[ooff+jDim*spaceDim+iDim] * dTractionTpdtVertex[jDim];
      } // for

      // :TODO: BRAD - add stuff here for updating slip

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
    const PetscInt soff = dLagrangeVisitor.sectionOffset(v_fault);
    assert(spaceDim == dLagrangeVisitor.sectionDof(v_fault));
    for(PetscInt d = 0; d < spaceDim; ++d) {
      dLagrangeArray[soff+d] = dLagrangeTpdtVertex[d];
    } // for
    
  } // for
  dispTIncrAdjVisitor.clear();
  dLagrangeVisitor.clear();

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

#if 0 // DEBUGGING
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

  topology::VecVisitorMesh sensDispRelVisitor(_fields->get("sensitivity relative disp"));
  PetscScalar* sensDispRelArray = sensDispRelVisitor.localArray();

  dispTIncrAdjVisitor.initialize(fields->get("dispIncr adjust"));
  dispTIncrAdjArray = dispTIncrAdjVisitor.localArray();

  dLagrangeVisitor.initialize(_fields->get("sensitivity dLagrange"));
  dLagrangeArray = dLagrangeVisitor.localArray();

  PetscErrorCode err = 0;
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_fault = _cohesiveVertices[iVertex].fault;
    const int e_lagrange = _cohesiveVertices[iVertex].lagrange;
    const int v_negative = _cohesiveVertices[iVertex].negative;
    const int v_positive = _cohesiveVertices[iVertex].positive;

    // Skip clamped vertices
    if (e_lagrange < 0) {
      continue;
    } // if

    // Get change in Lagrange multiplier computed from friction criterion.
    const PetscInt soff = dLagrangeVisitor.sectionOffset(v_fault);
    assert(spaceDim == dLagrangeVisitor.sectionDof(v_fault));

    // Get change in relative displacement from sensitivity solve.
    const PetscInt sroff = sensDispRelVisitor.sectionOffset(v_fault);
    assert(spaceDim == sensDispRelVisitor.sectionDof(v_fault));

    // Get current relative displacement for updating.
    assert(spaceDim == dispRelVisitor.sectionDof(v_fault));

    // Get orientation.
    const PetscInt ooff = orientationVisitor.sectionOffset(v_fault);
    assert(spaceDim*spaceDim == orientationVisitor.sectionDof(v_fault));

    // Get displacement.
    const PetscInt dtnoff = dispTVisitor.sectionOffset(v_negative);
    assert(spaceDim == dispTVisitor.sectionDof(v_negative));

    const PetscInt dtpoff = dispTVisitor.sectionOffset(v_positive);
    assert(spaceDim == dispTVisitor.sectionDof(v_positive));

    const PetscInt dtloff = dispTVisitor.sectionOffset(e_lagrange);
    assert(spaceDim == dispTVisitor.sectionDof(e_lagrange));

    // Get displacement increment (trial solution).
    const PetscInt dinoff = dispTIncrVisitor.sectionOffset(v_negative);
    assert(spaceDim == dispTIncrVisitor.sectionDof(v_negative));

    const PetscInt dipoff = dispTIncrVisitor.sectionOffset(v_positive);
    assert(spaceDim == dispTIncrVisitor.sectionDof(v_positive));

    const PetscInt diloff = dispTIncrVisitor.sectionOffset(e_lagrange);
    assert(spaceDim == dispTIncrVisitor.sectionDof(e_lagrange));

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
        dSlipTpdtVertex[iDim] += orientationArray[ooff+iDim*spaceDim+jDim] * alpha*sensDispRelArray[sroff+jDim];
        tractionTpdtVertex[iDim] += orientationArray[ooff+iDim*spaceDim+jDim] * (dispTArray[dtloff+jDim] + dispTIncrArray[diloff+jDim]);
        dTractionTpdtVertex[iDim] += orientationArray[ooff+iDim*spaceDim+jDim] * alpha*dLagrangeArray[soff+jDim];
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
    err = PetscSectionGetOffset(dispTIncrGlobalSection, e_lagrange, &goff);PYLITH_CHECK_ERROR(err);
    if (goff >= 0) {
      // Get offsets in displacement increment adjustment.
      const PetscInt dialoff = dispTIncrAdjVisitor.sectionOffset(e_lagrange);
      assert(spaceDim == dispTIncrAdjVisitor.sectionDof(e_lagrange));

      const PetscInt dianoff = dispTIncrAdjVisitor.sectionOffset(v_negative);
      assert(spaceDim == dispTIncrAdjVisitor.sectionDof(v_negative));

      const PetscInt diapoff = dispTIncrAdjVisitor.sectionOffset(v_positive);
      assert(spaceDim == dispTIncrAdjVisitor.sectionDof(v_positive));

      // Update Lagrange multiplier increment.
      for(PetscInt d = 0; d < spaceDim; ++d) {
        dispTIncrAdjArray[dialoff+d] += dLagrangeTpdtVertex[d];
        dispTIncrAdjArray[dianoff+d] += dDispTIncrVertexN[d];
        dispTIncrAdjArray[diapoff+d] += dDispTIncrVertexP[d];
      } // for
    } // if
  } // for

  PYLITH_METHOD_END;
} // constrainSolnSpace

// ----------------------------------------------------------------------
// Adjust solution from solver with lumped Jacobian to match Lagrange
// multiplier constraints.
void
pylith::faults::FaultCohesiveDyn::adjustSolnLumped(topology::SolutionFields* const fields,
						   const PylithScalar t,
						   const topology::Field& jacobian)
{ // adjustSolnLumped
  PYLITH_METHOD_BEGIN;

  /// Member prototype for _constrainSolnSpaceXD()
  typedef void (pylith::faults::FaultCohesiveDyn::*constrainSolnSpace_fn_type)
    (scalar_array*,
     const PylithScalar,
     const scalar_array&,
     const scalar_array&,
     const scalar_array&,
     const PylithScalar,
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
  const int computeEvent = _logger->eventId("FaAS compute");
#if defined(DETAILED_EVENT_LOGGING)
  const int geometryEvent = _logger->eventId("FaAS geometry");
  const int restrictEvent = _logger->eventId("FaAS restrict");
  const int updateEvent = _logger->eventId("FaAS update");
#endif

  _logger->eventBegin(setupEvent);

  // Get cell information and setup storage for cell data
  const int spaceDim = _quadrature->spaceDim();

  // Allocate arrays for vertex values
  scalar_array tractionTpdtVertex(spaceDim);
  scalar_array lagrangeTpdtVertex(spaceDim);
  scalar_array dTractionTpdtVertex(spaceDim);
  scalar_array dLagrangeTpdtVertex(spaceDim);

  // Update time step in friction (can vary).
  _friction->timeStep(_dt);

  // Get section information
  scalar_array slipVertex(spaceDim);
  scalar_array dispRelVertex(spaceDim);
  topology::VecVisitorMesh dispRelVisitor(_fields->get("relative disp"));
  PetscScalar* dispRelArray = dispRelVisitor.localArray();

  scalar_array slipRateVertex(spaceDim);

  topology::VecVisitorMesh areaVisitor(_fields->get("area"));
  const PetscScalar* areaArray = areaVisitor.localArray();

  topology::VecVisitorMesh orientationVisitor(_fields->get("orientation"));
  const PetscScalar* orientationArray = orientationVisitor.localArray();

  topology::VecVisitorMesh dispTVisitor(fields->get("disp(t)"));
  const PetscScalar* dispTArray = dispTVisitor.localArray();

  scalar_array dispIncrVertexN(spaceDim);
  scalar_array dispIncrVertexP(spaceDim);
  scalar_array lagrangeTIncrVertex(spaceDim);
  topology::VecVisitorMesh dispTIncrVisitor(fields->get("dispIncr(t->t+dt)"));
  PetscScalar* dispTIncrArray = dispTIncrVisitor.localArray();

  topology::VecVisitorMesh dispTIncrAdjVisitor(fields->get("dispIncr adjust"));
  PetscScalar* dispTIncrAdjArray = dispTIncrAdjVisitor.localArray();

  topology::VecVisitorMesh jacobianVisitor(jacobian);
  const PetscScalar* jacobianArray = jacobianVisitor.localArray();

  topology::VecVisitorMesh residualVisitor(fields->get("residual"));
  const PetscScalar* residualArray = residualVisitor.localArray();

  PetscSection solnGlobalSection = fields->get("dispIncr(t->t+dt)").globalSection();assert(solnGlobalSection);

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
    throw std::logic_error("Unknown spatial dimension in FaultCohesiveDyn::adjustSolnLumped.");
  } // switch

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

    // Skip clamped vertices
    if (e_lagrange < 0) {
      continue;
    } // if

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(restrictEvent);
#endif

    // Get residual at cohesive cell's vertices.
    const PetscInt rloff = residualVisitor.sectionOffset(e_lagrange);
    assert(spaceDim == residualVisitor.sectionDof(e_lagrange));

    // Get jacobian at cohesive cell's vertices.
    const PetscInt jnoff = jacobianVisitor.sectionOffset(v_negative);
    assert(spaceDim == jacobianVisitor.sectionDof(v_negative));

    const PetscInt jpoff = jacobianVisitor.sectionOffset(v_positive);
    assert(spaceDim == jacobianVisitor.sectionDof(v_positive));

    // Get disp(t) at Lagrange vertex.
    const PetscInt dtloff = dispTVisitor.sectionOffset(e_lagrange);
    assert(spaceDim == dispTVisitor.sectionDof(e_lagrange));

    // Get dispIncr(t) at cohesive cell's vertices.
    const PetscInt dinoff = dispTIncrVisitor.sectionOffset(v_negative);
    assert(spaceDim == dispTIncrVisitor.sectionDof(v_negative));

    const PetscInt dipoff = dispTIncrVisitor.sectionOffset(v_positive);
    assert(spaceDim == dispTIncrVisitor.sectionDof(v_positive));

    const PetscInt diloff = dispTIncrVisitor.sectionOffset(e_lagrange);
    assert(spaceDim == dispTIncrVisitor.sectionDof(e_lagrange));

    // Get relative displacement at fault vertex.
    const PetscInt droff = dispRelVisitor.sectionOffset(v_fault);
    assert(spaceDim == dispRelVisitor.sectionDof(v_fault));

    // Get area at fault vertex.
    const PetscInt aoff = areaVisitor.sectionOffset(v_fault);
    assert(1 == areaVisitor.sectionDof(v_fault));
    const PetscScalar areaVertex = areaArray[aoff];
    assert(areaVertex > 0.0);

    // Get fault orientation at fault vertex.
    const PetscInt ooff = orientationVisitor.sectionOffset(v_fault);
    assert(spaceDim*spaceDim == orientationVisitor.sectionDof(v_fault));

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(restrictEvent);
    _logger->eventBegin(computeEvent);
#endif

    // Adjust solution as in prescribed rupture, updating the Lagrange
    // multipliers and the corresponding displacment increments.
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      assert(jacobianArray[jpoff+iDim] > 0.0);
      assert(jacobianArray[jnoff+iDim] > 0.0);
      const PylithScalar S = (1.0/jacobianArray[jpoff+iDim] + 1.0/jacobianArray[jnoff+iDim]) * areaVertex*areaVertex;
      assert(S > 0.0);
      lagrangeTIncrVertex[iDim] = 1.0/S * (-residualArray[rloff+iDim] + areaVertex * (dispTIncrArray[dipoff+iDim] - dispTIncrArray[dinoff+iDim]));
      dispIncrVertexN[iDim] =  areaVertex / jacobianArray[jnoff+iDim]*lagrangeTIncrVertex[iDim];
      dispIncrVertexP[iDim] = -areaVertex / jacobianArray[jpoff+iDim]*lagrangeTIncrVertex[iDim];
    } // for

    // Compute slip, slip rate, and Lagrange multiplier at time t+dt
    // in fault coordinate system.
    slipVertex = 0.0;
    slipRateVertex = 0.0;
    tractionTpdtVertex = 0.0;
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      for (int jDim=0; jDim < spaceDim; ++jDim) {
        slipVertex[iDim] += orientationArray[ooff+iDim*spaceDim+jDim] * dispRelArray[droff+jDim];
        tractionTpdtVertex[iDim] += orientationArray[ooff+iDim*spaceDim+jDim] * (dispTArray[dtloff+jDim] + lagrangeTIncrVertex[jDim]);
      } // for
    } // for
    // Jacobian is diagonal and isotropic, so it is invariant with
    // respect to rotation and contains one unique term.
    const PylithScalar jacobianShearVertex = -1.0 / (areaVertex * (1.0 / jacobianArray[jnoff+0] + 1.0 / jacobianArray[jpoff+0]));
    
    // Get friction properties and state variables.
    _friction->retrievePropsStateVars(v_fault);

    // Use fault constitutive model to compute traction associated with
    // friction.
    dTractionTpdtVertex = 0.0;

    const bool iterating = false; // No iteration for friction in lumped soln
    CALL_MEMBER_FN(*this, constrainSolnSpaceFn)(&dTractionTpdtVertex, t, slipVertex, slipRateVertex, tractionTpdtVertex, jacobianShearVertex, iterating);

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

      dispIncrVertexN[iDim] += areaVertex * dLagrangeTpdtVertex[iDim] / jacobianArray[jnoff+iDim];
      dispIncrVertexP[iDim] -= areaVertex * dLagrangeTpdtVertex[iDim] / jacobianArray[jpoff+iDim];

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
    err = PetscSectionGetOffset(solnGlobalSection, e_lagrange, &goff);PYLITH_CHECK_ERROR(err);
    if (goff >= 0) {
      const PetscInt dianoff = dispTIncrAdjVisitor.sectionOffset(v_negative);
      assert(spaceDim == dispTIncrAdjVisitor.sectionDof(v_negative));

      const PetscInt diapoff = dispTIncrAdjVisitor.sectionOffset(v_positive);
      assert(spaceDim == dispTIncrAdjVisitor.sectionDof(v_positive));

      // Adjust displacements to account for Lagrange multiplier values
      // (assumed to be zero in preliminary solve).
      // Update displacement field
      for(PetscInt d = 0; d < spaceDim; ++d) {
        dispTIncrAdjArray[dianoff+d] += dispIncrVertexN[d];
        dispTIncrAdjArray[diapoff+d] += dispIncrVertexP[d];
      } // for
    } // if

    // The Lagrange multiplier and relative displacement are NOT
    // assembled across processors, so update even if Lagrange vertex
    // is not local.

    // Set Lagrange multiplier value. Value from preliminary solve is
    // bogus due to artificial diagonal entry in Jacobian of 1.0.
    for(PetscInt d = 0; d < spaceDim; ++d) {
      dispTIncrArray[diloff+d] = lagrangeTIncrVertex[d];
    } // for

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(updateEvent);
#endif
  } // for
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

  PYLITH_METHOD_END;
} // adjustSolnLumped

// ----------------------------------------------------------------------
// Get vertex field associated with integrator.
const pylith::topology::Field&
pylith::faults::FaultCohesiveDyn::vertexField(const char* name,
					      const topology::SolutionFields* fields)
{ // vertexField
  PYLITH_METHOD_BEGIN;

  assert(_faultMesh);
  assert(_quadrature);
  assert(_normalizer);
  assert(_fields);
  assert(_friction);

  const int cohesiveDim = _faultMesh->dimension();

  const topology::Field& orientation = _fields->get("orientation");

  if (0 == strcasecmp("slip", name)) {
    const topology::Field& dispRel = _fields->get("relative disp");
    _allocateBufferVectorField();
    topology::Field& buffer =  _fields->get("buffer (vector)");
    buffer.copy(dispRel);
    buffer.label("slip");
    FaultCohesiveLagrange::globalToFault(&buffer, orientation);
    PYLITH_METHOD_RETURN(buffer);

  } else if (0 == strcasecmp("slip_rate", name)) {
    const topology::Field& velRel = _fields->get("relative velocity");
    _allocateBufferVectorField();
    topology::Field& buffer = _fields->get("buffer (vector)");
    buffer.copy(velRel);
    buffer.label("slip_rate");
    FaultCohesiveLagrange::globalToFault(&buffer, orientation);
    PYLITH_METHOD_RETURN(buffer);

  } else if (cohesiveDim > 0 && 0 == strcasecmp("strike_dir", name)) {
    _allocateBufferVectorField();
    topology::Field& buffer = _fields->get("buffer (vector)");
    buffer.copySubfield(orientation, "strike_dir");
    PYLITH_METHOD_RETURN(buffer);

  } else if (2 == cohesiveDim && 0 == strcasecmp("dip_dir", name)) {
    _allocateBufferVectorField();
    topology::Field& buffer = _fields->get("buffer (vector)");
    buffer.copySubfield(orientation, "dip_dir");
    PYLITH_METHOD_RETURN(buffer);

  } else if (0 == strcasecmp("normal_dir", name)) {
    _allocateBufferVectorField();
    topology::Field& buffer = _fields->get("buffer (vector)");
    buffer.copySubfield(orientation, "normal_dir");
    PYLITH_METHOD_RETURN(buffer);

  } else if (0 == strcasecmp("traction", name)) {
    assert(fields);
    const topology::Field& dispT = fields->get("disp(t)");
    _allocateBufferVectorField();
    topology::Field& buffer = _fields->get("buffer (vector)");
    _calcTractions(&buffer, dispT);
    PYLITH_METHOD_RETURN(buffer);

  } else if (_friction->hasPropStateVar(name)) {
    PYLITH_METHOD_RETURN(_friction->getField(name));

  } else if (_tractPerturbation && _tractPerturbation->hasParameter(name)) {
    const topology::Field& param = _tractPerturbation->vertexField(name, fields);
    if (param.vectorFieldType() == topology::FieldBase::VECTOR) {
      _allocateBufferVectorField();
      topology::Field& buffer = _fields->get("buffer (vector)");
      buffer.copy(param);
      FaultCohesiveLagrange::globalToFault(&buffer, orientation);
      PYLITH_METHOD_RETURN(buffer);
    } else {
      PYLITH_METHOD_RETURN(param);
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
  const topology::Field& buffer = _fields->get("buffer (vector)");

  PYLITH_METHOD_RETURN(buffer);
} // vertexField

// ----------------------------------------------------------------------
// Compute tractions on fault surface using solution.
void
pylith::faults::FaultCohesiveDyn::_calcTractions(topology::Field* tractions,
						 const topology::Field& dispT)
{ // _calcTractions
  PYLITH_METHOD_BEGIN;

  assert(tractions);
  assert(_faultMesh);
  assert(_fields);
  assert(_normalizer);

  // Fiber dimension of tractions matches spatial dimension.
  const int spaceDim = _quadrature->spaceDim();

  // Get fields.
  topology::VecVisitorMesh dispTVisitor(dispT);
  const PetscScalar* dispTArray = dispTVisitor.localArray();

  topology::VecVisitorMesh orientationVisitor(_fields->get("orientation"));
  const PetscScalar* orientationArray = orientationVisitor.localArray();

  // Allocate buffer for tractions field (if necessary).
  if (!tractions->localSection()) {
    const topology::Field& dispRel = _fields->get("relative disp");
    tractions->cloneSection(dispRel);
  } // if
  const PylithScalar pressureScale = _normalizer->pressureScale();
  tractions->label("traction");
  tractions->scale(pressureScale);
  tractions->zeroAll();

  topology::VecVisitorMesh tractionsVisitor(*tractions);
  PetscScalar* tractionsArray = tractionsVisitor.localArray();

  const int numVertices = _cohesiveVertices.size();
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int e_lagrange = _cohesiveVertices[iVertex].lagrange;
    const int v_fault = _cohesiveVertices[iVertex].fault;

    // Skip clamped vertices
    if (e_lagrange < 0) {
      continue;
    } // if

    const PetscInt dtloff = dispTVisitor.sectionOffset(e_lagrange);
    assert(spaceDim == dispTVisitor.sectionDof(e_lagrange));

    const PetscInt ooff = orientationVisitor.sectionOffset(v_fault);
    assert(spaceDim*spaceDim == orientationVisitor.sectionDof(v_fault));

    const PetscInt toff = tractionsVisitor.sectionOffset(v_fault);
    assert(spaceDim == tractionsVisitor.sectionDof(v_fault));

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

  PYLITH_METHOD_END;
} // _calcTractions

// ----------------------------------------------------------------------
// Update relative displacement and velocity (slip and slip rate)
// associated with Lagrange vertex k corresponding to diffential
// velocity between conventional vertices i and j.
void
pylith::faults::FaultCohesiveDyn::_updateRelMotion(const topology::SolutionFields& fields)
{ // _updateRelMotion
  PYLITH_METHOD_BEGIN;

  assert(_fields);

  const int spaceDim = _quadrature->spaceDim();

  // Get section information
  topology::VecVisitorMesh dispTVisitor(fields.get("disp(t)"));
  const PetscScalar* dispTArray = dispTVisitor.localArray();

  topology::VecVisitorMesh dispTIncrVisitor(fields.get("dispIncr(t->t+dt)"));
  const PetscScalar* dispTIncrArray = dispTIncrVisitor.localArray();

  topology::VecVisitorMesh velocityVisitor(fields.get("velocity(t)"));
  const PetscScalar* velocityArray = velocityVisitor.localArray();

  topology::VecVisitorMesh dispRelVisitor(_fields->get("relative disp"));
  PetscScalar* dispRelArray = dispRelVisitor.localArray();

  topology::VecVisitorMesh velRelVisitor(_fields->get("relative velocity"));
  PetscScalar* velRelArray = velRelVisitor.localArray();

  const int numVertices = _cohesiveVertices.size();
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_fault = _cohesiveVertices[iVertex].fault;
    const int v_negative = _cohesiveVertices[iVertex].negative;
    const int v_positive = _cohesiveVertices[iVertex].positive;

    // Skip clamped vertices
    if (v_fault < 0) {
      continue;
    } // if

    // Get displacement offsets.
    const PetscInt dtnoff = dispTVisitor.sectionOffset(v_negative);
    assert(spaceDim == dispTVisitor.sectionDof(v_negative));
    
    const PetscInt dtpoff = dispTVisitor.sectionOffset(v_positive);
    assert(spaceDim == dispTVisitor.sectionDof(v_positive));

    const PetscInt dinoff = dispTIncrVisitor.sectionOffset(v_negative);
    assert(spaceDim == dispTIncrVisitor.sectionDof(v_negative));

    const PetscInt dipoff = dispTIncrVisitor.sectionOffset(v_positive);
    assert(spaceDim == dispTIncrVisitor.sectionDof(v_positive));

    // Get velocity offsets.
    const PetscInt vnoff = velocityVisitor.sectionOffset(v_negative);
    assert(spaceDim == velocityVisitor.sectionDof(v_negative));

    const PetscInt vpoff = velocityVisitor.sectionOffset(v_positive);
    assert(spaceDim == velocityVisitor.sectionDof(v_positive));

    // Relative displacement/velocity offsets.
    const PetscInt droff = dispRelVisitor.sectionOffset(v_fault);
    assert(spaceDim == dispRelVisitor.sectionDof(v_fault));

    const PetscInt vroff = velRelVisitor.sectionOffset(v_fault);
    assert(spaceDim == velRelVisitor.sectionDof(v_fault));

    for(PetscInt d = 0; d < spaceDim; ++d) {
      const PylithScalar dispValue = dispTArray[dtpoff+d] + dispTIncrArray[dipoff+d] - dispTArray[dtnoff+d] - dispTIncrArray[dinoff+d];
      dispRelArray[droff+d] = fabs(dispValue) > _zeroTolerance ? dispValue : 0.0;

      const PylithScalar velValue = velocityArray[vpoff+d] - velocityArray[vnoff+d];
      velRelArray[vroff+d] = fabs(velValue) > _zeroTolerance ? velValue : 0.0;
    } // for

  } // for
  PetscLogFlops(numVertices*spaceDim*spaceDim*4);

  PYLITH_METHOD_END;
} // _updateRelMotion

// ----------------------------------------------------------------------
// Setup sensitivity problem to compute change in slip given change in Lagrange multipliers.
void
pylith::faults::FaultCohesiveDyn::_sensitivitySetup(const topology::Jacobian& jacobian)
{ // _sensitivitySetup
  PYLITH_METHOD_BEGIN;

  assert(_fields);
  assert(_quadrature);

  // Setup fields involved in sensitivity solve.
  if (!_fields->hasField("sensitivity solution")) {
    _fields->add("sensitivity solution", "sensitivity_soln");
    topology::Field& solution = _fields->get("sensitivity solution");
    const topology::Field& dispRel = _fields->get("relative disp");
    solution.cloneSection(dispRel);
    solution.createScatter(solution.mesh());
  } // if
  const topology::Field& solution = _fields->get("sensitivity solution");

  if (!_fields->hasField("sensitivity residual")) {
    _fields->add("sensitivity residual", "sensitivity_residual");
    topology::Field& residual = _fields->get("sensitivity residual");
    residual.cloneSection(solution);
    residual.createScatter(solution.mesh());
  } // if

  if (!_fields->hasField("sensitivity relative disp")) {
    _fields->add("sensitivity relative disp", "sensitivity_relative_disp");
    topology::Field& dispRel = _fields->get("sensitivity relative disp");
    dispRel.cloneSection(solution);
  } // if
  topology::Field& dispRel = _fields->get("sensitivity relative disp");
  dispRel.zeroAll();

  if (!_fields->hasField("sensitivity dLagrange")) {
    _fields->add("sensitivity dLagrange", "sensitivity_dlagrange");
    topology::Field& dLagrange = _fields->get("sensitivity dLagrange");
    dLagrange.cloneSection(solution);
    topology::VecVisitorMesh::optimizeClosure(dLagrange);
  } // if
  topology::Field& dLagrange = _fields->get("sensitivity dLagrange");
  dLagrange.zeroAll();

  // Setup Jacobian sparse matrix for sensitivity solve.
  if (!_jacobian) {
    _jacobian = new topology::Jacobian(solution, jacobian.matrixType());
  } // if
  assert(_jacobian);
  _jacobian->zero();

  // Setup PETSc KSP linear solver.
  if (!_ksp) {
    PetscErrorCode err = 0;
    err = KSPCreate(_faultMesh->comm(), &_ksp);PYLITH_CHECK_ERROR(err);
    err = KSPSetInitialGuessNonzero(_ksp, PETSC_FALSE);PYLITH_CHECK_ERROR(err);
    PylithScalar rtol = 0.0;
    PylithScalar atol = 0.0;
    PylithScalar dtol = 0.0;
    int maxIters = 0;
    err = KSPGetTolerances(_ksp, &rtol, &atol, &dtol, &maxIters);PYLITH_CHECK_ERROR(err);
    rtol = 1.0e-3*_zeroTolerance;
    atol = 1.0e-5*_zeroTolerance;
    err = KSPSetTolerances(_ksp, rtol, atol, dtol, maxIters);PYLITH_CHECK_ERROR(err);

    PC pc;
    err = KSPGetPC(_ksp, &pc);PYLITH_CHECK_ERROR(err);
    err = PCSetType(pc, PCJACOBI);PYLITH_CHECK_ERROR(err);
    err = KSPSetType(_ksp, KSPGMRES);PYLITH_CHECK_ERROR(err);

    err = KSPAppendOptionsPrefix(_ksp, "friction_");PYLITH_CHECK_ERROR(err);
    err = KSPSetFromOptions(_ksp);PYLITH_CHECK_ERROR(err);
  } // if

  PYLITH_METHOD_END;
} // _sensitivitySetup

// ----------------------------------------------------------------------
// Update the Jacobian values for the sensitivity solve.
void
pylith::faults::FaultCohesiveDyn::_sensitivityUpdateJacobian(const bool negativeSide,
                                                             const topology::Jacobian& jacobian,
                                                             const topology::SolutionFields& fields)
{ // _sensitivityUpdateJacobian
  PYLITH_METHOD_BEGIN;

  assert(_quadrature);
  assert(_fields);

  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int subnrows = numBasis*spaceDim;
  const int submatrixSize = subnrows * subnrows;

  PetscErrorCode err = 0;

  // Get solution field
  const topology::Field& solutionDomain = fields.solution();
  PetscSection solutionDomainSection = solutionDomain.localSection();assert(solutionDomainSection);
  PetscVec solutionDomainVec = solutionDomain.localVector();assert(solutionDomainVec);
  PetscSection solutionDomainGlobalSection = solutionDomain.globalSection();assert(solutionDomainGlobalSection);

  // Get cohesive cells
  PetscDM dmMesh = fields.mesh().dmMesh();assert(dmMesh);
  assert(_cohesiveIS);
  const PetscInt *cellsCohesive = _cohesiveIS->points();
  const PetscInt numCohesiveCells = _cohesiveIS->size();

  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();

  // Visitor for Jacobian matrix associated with domain.
  scalar_array jacobianSubCell(submatrixSize);
  const PetscMat jacobianDomainMatrix = jacobian.matrix();assert(jacobianDomainMatrix);

  // Get fault mesh
  PetscDM faultDMMesh = _faultMesh->dmMesh();assert(faultDMMesh);

  // Get sensitivity solution field
  PetscSection solutionFaultSection = _fields->get("sensitivity solution").localSection();assert(solutionFaultSection);
  PetscVec solutionFaultVec = _fields->get("sensitivity solution").localVector();assert(solutionFaultVec);
  PetscSection solutionFaultGlobalSection = _fields->get("sensitivity solution").globalSection();assert(solutionFaultGlobalSection);
  assert(_jacobian);
  const PetscMat jacobianFaultMatrix = _jacobian->matrix();assert(jacobianFaultMatrix);

  const int iCone = (negativeSide) ? 0 : 1;

  PetscIS* cellsIS = (numCohesiveCells > 0) ? new PetscIS[numCohesiveCells] : 0;
  int_array indicesGlobal(subnrows);
  int_array indicesLocal(numCohesiveCells*subnrows);
  int_array indicesPerm(subnrows);
  for (PetscInt c = 0; c < numCohesiveCells; ++c) {
    // Get cone for cohesive cell
    const PetscInt *cone;
    PetscInt        coneSize;
    PetscInt       *closureA = NULL, *closureB = NULL;
    PetscInt        closureSizeA, closureSizeB, q;

    err = DMPlexGetCone(dmMesh, cellsCohesive[c], &cone);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetConeSize(dmMesh, cellsCohesive[c], &coneSize);PYLITH_CHECK_ERROR(err);
    assert(coneSize >= 4);
    err = DMPlexGetTransitiveClosure(dmMesh, cone[0], PETSC_TRUE, &closureSizeA, &closureA);PYLITH_CHECK_ERROR(err);
    // Filter out non-vertices
    q = 0;
    for(PetscInt p = 0; p < closureSizeA*2; p += 2) {
      if ((closureA[p] >= vStart) && (closureA[p] < vEnd)) {
        closureA[q] = closureA[p];
        ++q;
      } // if
    } // for
    closureSizeA = q;
    err = DMPlexGetTransitiveClosure(dmMesh, cone[1], PETSC_TRUE, &closureSizeB, &closureB);PYLITH_CHECK_ERROR(err);
    // Filter out non-vertices
    q = 0;
    for(PetscInt p = 0; p < closureSizeB*2; p += 2) {
      if ((closureB[p] >= vStart) && (closureB[p] < vEnd)) {
        closureB[q] = closureB[p];
        ++q;
      } // if
    } // for
    closureSizeB = q;
    assert(closureSizeA == numBasis);
    assert(closureSizeB == numBasis);

    // Get indices
    for (int iBasis = 0; iBasis < numBasis; ++iBasis) {
      // negative side of the fault: iCone=0
      // positive side of the fault: iCone=1
      const int v_domain = iCone ? closureB[iBasis] :  closureA[iBasis];
      PetscInt goff;

      err = PetscSectionGetOffset(solutionDomainGlobalSection, v_domain, &goff);PYLITH_CHECK_ERROR(err);
      for (int iDim = 0, iB = iBasis*spaceDim, gind = goff < 0 ? -(goff+1) : goff; iDim < spaceDim; ++iDim) {
        indicesGlobal[iB+iDim] = gind + iDim;
      } // for

    } // for
    err = DMPlexRestoreTransitiveClosure(dmMesh, cone[0], PETSC_TRUE, &closureSizeA, &closureA);PYLITH_CHECK_ERROR(err);
    err = DMPlexRestoreTransitiveClosure(dmMesh, cone[1], PETSC_TRUE, &closureSizeB, &closureB);PYLITH_CHECK_ERROR(err);

    for (int i=0; i < subnrows; ++i) {
      indicesPerm[i]  = i;
    } // for
    err = PetscSortIntWithArray(indicesGlobal.size(), &indicesGlobal[0], &indicesPerm[0]);PYLITH_CHECK_ERROR(err);

    for (int i=0; i < subnrows; ++i) {
      indicesLocal[c*subnrows+indicesPerm[i]] = i;
    } // for
    cellsIS[c] = NULL;
    err = ISCreateGeneral(PETSC_COMM_SELF, indicesGlobal.size(), &indicesGlobal[0], PETSC_COPY_VALUES, &cellsIS[c]);PYLITH_CHECK_ERROR(err);

  } // for

  PetscMat* submatrices = NULL;
  err = MatGetSubMatrices(jacobianDomainMatrix, numCohesiveCells, cellsIS, cellsIS, MAT_INITIAL_MATRIX, &submatrices);PYLITH_CHECK_ERROR(err);

  for (PetscInt c = 0; c < numCohesiveCells; ++c) {
    // Get values for submatrix associated with cohesive cell
    jacobianSubCell = 0.0;
    err = MatGetValues(submatrices[c], subnrows, &indicesLocal[c*subnrows], subnrows, &indicesLocal[c*subnrows],
                       &jacobianSubCell[0]);PYLITH_CHECK_ERROR_MSG(err, "Restrict from PETSc Mat failed.");

    // Insert cell contribution into PETSc Matrix
    PetscInt c_fault = _cohesiveToFault[cellsCohesive[c]];

    err = DMPlexMatSetClosure(faultDMMesh, solutionFaultSection, solutionFaultGlobalSection,  jacobianFaultMatrix, c_fault, &jacobianSubCell[0], INSERT_VALUES);PYLITH_CHECK_ERROR_MSG(err, "Update to PETSc Mat failed.");

    // Destory IS for cohesiveCell
    err = ISDestroy(&cellsIS[c]);PYLITH_CHECK_ERROR(err);
  } // for

  err = MatDestroyMatrices(numCohesiveCells, &submatrices);PYLITH_CHECK_ERROR(err);
  delete[] cellsIS; cellsIS = 0;

  _jacobian->assemble("final_assembly");

#if 0 // DEBUGGING
  //std::cout << "DOMAIN JACOBIAN" << std::endl;
  //jacobian.view();
  std::cout << "SENSITIVITY JACOBIAN" << std::endl;
  _jacobian->view();
#endif

  PYLITH_METHOD_END;
} // _sensitivityUpdateJacobian

// ----------------------------------------------------------------------
// Reform residual for sensitivity problem.
void
pylith::faults::FaultCohesiveDyn::_sensitivityReformResidual(const bool negativeSide)
{ // _sensitivityReformResidual
  PYLITH_METHOD_BEGIN;

  /** Compute residual -L^T dLagrange
   *
   * Note: We need all entries for L, even those on other processors,
   * so we compute L rather than extract entries from the Jacobian.
   */

  const PylithScalar signFault = (negativeSide) ?  1.0 : -1.0;

  // Get cell information
  const size_t numQuadPts = _quadrature->numQuadPts();
  const scalar_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == numQuadPts);
  const int spaceDim = _quadrature->spaceDim();
  const int numBasis = _quadrature->numBasis();

  scalar_array basisProducts(numBasis*numBasis);

  // Get fault cell information
  PetscDM faultDMMesh = _faultMesh->dmMesh();assert(faultDMMesh);
  topology::Stratum cellsStratum(faultDMMesh, topology::Stratum::HEIGHT, 1);
  const PetscInt cStart = cellsStratum.begin();
  const PetscInt cEnd = cellsStratum.end();

  // Get sections
  scalar_array coordsCell(numBasis*spaceDim); // :KULDGE: Update numBasis to numCorners after implementing higher order
  topology::CoordsVisitor coordsVisitor(faultDMMesh);

  scalar_array dLagrangeCell(numBasis*spaceDim);
  topology::VecVisitorMesh dLagrangeVisitor(_fields->get("sensitivity dLagrange"));

  scalar_array residualCell(numBasis*spaceDim);
  topology::Field& residual = _fields->get("sensitivity residual");
  topology::VecVisitorMesh residualVisitor(residual);
  residual.zeroAll();

  // Loop over cells
  for(PetscInt c = cStart; c < cEnd; ++c) {
    // Compute geometry
    coordsVisitor.getClosure(&coordsCell, c);
    _quadrature->computeGeometry(&coordsCell[0], coordsCell.size(), c);

    // Restrict input fields to cell
    dLagrangeVisitor.getClosure(&dLagrangeCell, c);

    // Get cell geometry information that depends on cell
    const scalar_array& basis = _quadrature->basis();
    const scalar_array& jacobianDet = _quadrature->jacobianDet();

    // Compute product of basis functions.
    // Want values summed over quadrature points
    basisProducts = 0.0;
    for (size_t iQuad=0; iQuad < numQuadPts; ++iQuad) {
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
          residualCell[iBasis*spaceDim+d] += l * dLagrangeCell[jBasis*spaceDim+d];
        } // for
      } // for
    } // for

    // Assemble cell contribution into field
    residualVisitor.setClosure(&residualCell[0], residualCell.size(), c, ADD_VALUES);
  } // for

  PYLITH_METHOD_END;
} // _sensitivityReformResidual

// ----------------------------------------------------------------------
// Solve sensitivity problem.
void
pylith::faults::FaultCohesiveDyn::_sensitivitySolve(void)
{ // _sensitivitySolve
  PYLITH_METHOD_BEGIN;

  assert(_fields);
  assert(_jacobian);
  assert(_ksp);

  topology::Field& residual = _fields->get("sensitivity residual");
  topology::Field& solution = _fields->get("sensitivity solution");

  // Assemble residual over processors.
  residual.complete();

  // Update PetscVector view of field.
  residual.scatterLocalToGlobal();

  PetscErrorCode err = 0;
  const PetscMat jacobianMat = _jacobian->matrix();
  err = KSPSetOperators(_ksp, jacobianMat, jacobianMat);PYLITH_CHECK_ERROR(err);

  const PetscVec residualVec = residual.globalVector();
  const PetscVec solutionVec = solution.globalVector();
  err = KSPSolve(_ksp, residualVec, solutionVec);PYLITH_CHECK_ERROR(err);

  // Update section view of field.
  solution.scatterGlobalToLocal();

#if 0 // DEBUGGING
  residual.view("SENSITIVITY RESIDUAL");
  solution.view("SENSITIVITY SOLUTION");
#endif

  PYLITH_METHOD_END;
} // _sensitivitySolve

// ----------------------------------------------------------------------
// Update the relative displacement field values based on the
// sensitivity solve.
void
pylith::faults::FaultCohesiveDyn::_sensitivityUpdateSoln(const bool negativeSide)
{ // _sensitivityUpdateSoln
  PYLITH_METHOD_BEGIN;

  assert(_fields);
  assert(_quadrature);

  const int spaceDim = _quadrature->spaceDim();

  topology::VecVisitorMesh solutionVisitor(_fields->get("sensitivity solution"));
  const PetscScalar* solutionArray = solutionVisitor.localArray();

  topology::VecVisitorMesh dispRelVisitor(_fields->get("sensitivity relative disp"));
  PetscScalar* dispRelArray = dispRelVisitor.localArray();

  topology::VecVisitorMesh dLagrangeVisitor(_fields->get("sensitivity dLagrange"));
  PetscScalar* dLagrangeArray = dLagrangeVisitor.localArray();

  const PylithScalar sign = (negativeSide) ? -1.0 : 1.0;
  const int numVertices = _cohesiveVertices.size();
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_fault = _cohesiveVertices[iVertex].fault;

    // Skip clamped vertices
    if (v_fault < 0) {
      continue;
    } // if

    const PetscInt dloff = dLagrangeVisitor.sectionOffset(v_fault);
    assert(spaceDim == dLagrangeVisitor.sectionDof(v_fault));

    // If no change in the Lagrange multiplier computed from friction criterion, there are no updates, so continue.
    PylithScalar dLagrangeVertexMag = 0.0;
    for(PetscInt d = 0; d < spaceDim; ++d) {
      dLagrangeVertexMag += dLagrangeArray[dloff+d]*dLagrangeArray[dloff+d];
    } // for
    if (0.0 == dLagrangeVertexMag) {
      continue;
    } // if

    const PetscInt soff = solutionVisitor.sectionOffset(v_fault);
    assert(spaceDim == solutionVisitor.sectionDof(v_fault));

    const PetscInt droff = dispRelVisitor.sectionOffset(v_fault);
    assert(spaceDim == dispRelVisitor.sectionDof(v_fault));

    // Update relative displacements associated with sensitivity solve solution
    for(PetscInt d = 0; d < spaceDim; ++d) {
      dispRelArray[droff+d] += sign*solutionArray[soff+d];
    } // for
  } // for

  PYLITH_METHOD_END;
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
  PYLITH_METHOD_BEGIN;

  /// Member prototype for _constrainSolnSpaceXD()
  typedef void (pylith::faults::FaultCohesiveDyn::*constrainSolnSpace_fn_type)
    (scalar_array*,
     const PylithScalar,
     const scalar_array&,
     const scalar_array&,
     const scalar_array&,
     const PylithScalar,
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

  // Get fields
  scalar_array slipTpdtVertex(spaceDim); // fault coordinates
  scalar_array slipRateVertex(spaceDim); // fault coordinates
  scalar_array tractionTpdtVertex(spaceDim); // fault coordinates
  scalar_array tractionMisfitVertex(spaceDim); // fault coordinates

  topology::VecVisitorMesh orientationVisitor(_fields->get("orientation"));
  const PetscScalar* orientationArray = orientationVisitor.localArray();

  topology::VecVisitorMesh dLagrangeVisitor(_fields->get("sensitivity dLagrange"));
  const PetscScalar* dLagrangeArray = dLagrangeVisitor.localArray();

  topology::VecVisitorMesh sensDispRelVisitor(_fields->get("sensitivity relative disp"));
  const PetscScalar* sensDispRelArray = sensDispRelVisitor.localArray();

  topology::VecVisitorMesh dispTVisitor(fields->get("disp(t)"));
  const PetscScalar* dispTArray = dispTVisitor.localArray();

  topology::Field& dispTIncr = fields->get("dispIncr(t->t+dt)");
  topology::VecVisitorMesh dispTIncrVisitor(dispTIncr);
  const PetscScalar* dispTIncrArray = dispTIncrVisitor.localArray();
  PetscSection dispTIncrGlobalSection = dispTIncr.globalSection();assert(dispTIncrGlobalSection);

  bool isOpening = false;
  PylithScalar norm2 = 0.0;
  int numVertices = _cohesiveVertices.size();
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int e_lagrange = _cohesiveVertices[iVertex].lagrange;
    const int v_fault = _cohesiveVertices[iVertex].fault;
    const int v_negative = _cohesiveVertices[iVertex].negative;
    const int v_positive = _cohesiveVertices[iVertex].positive;

    // Skip clamped vertices
    if (e_lagrange < 0) {
      continue;
    } // if

    // Compute contribution only if Lagrange constraint is local.
    PetscInt goff;
    err = PetscSectionGetOffset(dispTIncrGlobalSection, e_lagrange, &goff);PYLITH_CHECK_ERROR(err);
    if (goff < 0) {
      continue;
    } // if

    // Get displacement values
    const PetscInt dtnoff = dispTVisitor.sectionOffset(v_negative);
    assert(spaceDim == dispTVisitor.sectionDof(v_negative));
    
    const PetscInt dtpoff = dispTVisitor.sectionOffset(v_positive);
    assert(spaceDim == dispTVisitor.sectionDof(v_positive));
    
    const PetscInt dtloff = dispTVisitor.sectionOffset(e_lagrange);
    assert(spaceDim == dispTVisitor.sectionDof(e_lagrange));
    
    // Get displacement increment values.
    const PetscInt dinoff = dispTIncrVisitor.sectionOffset(v_negative);
    assert(spaceDim == dispTIncrVisitor.sectionDof(v_negative));
    
    const PetscInt dipoff = dispTIncrVisitor.sectionOffset(v_positive);
    assert(spaceDim == dispTIncrVisitor.sectionDof(v_positive));
    
    const PetscInt diloff = dispTIncrVisitor.sectionOffset(e_lagrange);
    assert(spaceDim == dispTIncrVisitor.sectionDof(e_lagrange));
    
    // Get orientation
    const PetscInt ooff = orientationVisitor.sectionOffset(v_fault);
    assert(spaceDim*spaceDim == orientationVisitor.sectionDof(v_fault));
    
    // Get change in relative displacement from sensitivity solve.
    const PetscInt sdroff = sensDispRelVisitor.sectionOffset(v_fault);
    assert(spaceDim == sensDispRelVisitor.sectionDof(v_fault));

    const PetscInt sdloff = dLagrangeVisitor.sectionOffset(v_fault);
    assert(spaceDim == dLagrangeVisitor.sectionDof(v_fault));

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
        tractionTpdtVertex[d] += orientationArray[ooff+d*spaceDim+e] * (dispTArray[dtloff+e] + dispTIncrArray[diloff+e] + alpha*dLagrangeArray[sdloff+e]);
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
    const PylithScalar jacobianShearVertex = 0.0;
    const bool iterating = true; // Iterating to get friction
    CALL_MEMBER_FN(*this, constrainSolnSpaceFn)(&tractionMisfitVertex, t,
                                                slipTpdtVertex, slipRateVertex, tractionTpdtVertex, jacobianShearVertex,
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

  if (isOpening && alpha < 1.0) {
    norm2 = PYLITH_MAXFLOAT;
  } // if

  PetscScalar norm2Total = 0.0;
  PetscInt numVerticesTotal = 0;
  err = MPI_Allreduce(&norm2, &norm2Total, 1, MPIU_SCALAR, MPI_SUM, fields->mesh().comm());
  err = MPI_Allreduce(&numVertices, &numVerticesTotal, 1, MPIU_INT, MPI_SUM, fields->mesh().comm());

  assert(numVerticesTotal > 0);
  PYLITH_METHOD_RETURN(sqrt(norm2Total) / numVerticesTotal);
} // _constrainSolnSpaceNorm


// ----------------------------------------------------------------------
// Constrain solution space in 1-D.
void
pylith::faults::FaultCohesiveDyn::_constrainSolnSpace1D(scalar_array* dTractionTpdt,
							const PylithScalar t,
							const scalar_array& slip,
							const scalar_array& sliprate,
							const scalar_array& tractionTpdt,
							const PylithScalar jacobianShear,
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
							const PylithScalar jacobianShear,
							const bool iterating)
{ // _constrainSolnSpace2D
  assert(dTractionTpdt);

  PylithScalar slipMag = fabs(slip[0]);
  const PylithScalar slipRateMag = fabs(slipRate[0]);

  const PylithScalar tractionNormal = tractionTpdt[1];
  const PylithScalar tractionShearMag = fabs(tractionTpdt[0]);
  
  if (fabs(slip[1]) < _zeroTolerance && tractionNormal < -_zeroTolerance) {
    // if in compression and no opening
    PylithScalar frictionStress = _friction->calcFriction(t, slipMag, slipRateMag, tractionNormal);

    if (tractionShearMag > frictionStress || (iterating && slipRateMag > 0.0)) {
      // traction is limited by friction, so have sliding OR
      // friction exceeds traction due to overshoot in slip

      if (tractionShearMag > 0.0) {
#if 1 // New Newton stuff
	if (0.0 != jacobianShear) {
	  assert(jacobianShear < 0.0);
	  // Use Newton to get better update
	  const int maxiter = 32;
	  PylithScalar slipMagCur = slipMag;
	  PylithScalar slipRateMagCur = slipRateMag;
	  PylithScalar tractionShearMagCur = tractionShearMag;
	  const PylithScalar slipMag0 = fabs(slip[0] - slipRate[0] * _dt);
	  for (int iter=0; iter < maxiter; ++iter) {
	    const PylithScalar frictionDeriv = _friction->calcFrictionDeriv(t, slipMagCur, slipRateMagCur, tractionNormal);
	    slipMag = slipMagCur;
	    if (slipMag > 0.0) {
	      // Use Newton (in log slip space) to get better update in slip & traction.
	      // D_{i+1} = exp(ln(D_i) - (T-T_f)/(D_i * (jacobian - frictionDeriv))
	      slipMagCur = exp(log(slipMag) - (tractionShearMagCur - frictionStress) / (slipMag * (jacobianShear - frictionDeriv)));
	    } else {
	      // Use Newton (in linear slip space) to get better update in slip & traction.
	      // D_{i+1} = D_i - (T-T_f)/(jacobian - frictionDeriv)
	      slipMagCur = slipMag - (tractionShearMagCur - frictionStress) / (jacobianShear - frictionDeriv);
	    } // if
	    tractionShearMagCur += (slipMagCur - slipMag) * jacobianShear;
	    slipRateMagCur = (slipMagCur - slipMag0) / _dt;
	    frictionStress = _friction->calcFriction(t, slipMagCur, slipRateMagCur, tractionNormal);
	    if (fabs(tractionShearMagCur - frictionStress) < _zeroTolerance) {
	      break;
	    } // if
	  } // for
	} // if
#endif

	// Update traction increment based on value required to stick
	// versus friction
	const PylithScalar dlp = -(tractionShearMag - frictionStress) * tractionTpdt[0] / tractionShearMag;
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
							const PylithScalar jacobianShear,
							const bool iterating)
{ // _constrainSolnSpace3D
  assert(dTractionTpdt);

  PylithScalar slipMag = sqrt(slip[0] * slip[0] + slip[1] * slip[1]);
  const PylithScalar slipRateMag = sqrt(slipRate[0]*slipRate[0] + slipRate[1]*slipRate[1]);
  
  const PylithScalar tractionNormal = tractionTpdt[2];
  const PylithScalar tractionShearMag = sqrt(tractionTpdt[0] * tractionTpdt[0] + tractionTpdt[1] * tractionTpdt[1]);
  
  if (fabs(slip[2]) < _zeroTolerance && tractionNormal < -_zeroTolerance) {
    // if in compression and no opening
    PylithScalar frictionStress = _friction->calcFriction(t, slipMag, slipRateMag, tractionNormal);

    if (tractionShearMag > frictionStress || (iterating && slipRateMag > 0.0)) {
      // traction is limited by friction, so have sliding OR
      // friction exceeds traction due to overshoot in slip
      
      if (tractionShearMag > 0.0) {
#if 1 // New Newton stuff
	if (0.0 != jacobianShear) {
	  assert(jacobianShear < 0.0);
	  // Use Newton to get better update
	  const int maxiter = 32;
	  PylithScalar slipMagCur = slipMag;
	  PylithScalar slipRateMagCur = slipRateMag;
	  PylithScalar tractionShearMagCur = tractionShearMag;
	  const PylithScalar slipMag0 = sqrt(pow(slip[0]-slipRate[0]*_dt, 2) + pow(slip[1]-slipRate[1]*_dt, 2));
	  for (int iter=0; iter < maxiter; ++iter) {
	    const PylithScalar frictionDeriv = _friction->calcFrictionDeriv(t, slipMagCur, slipRateMagCur, tractionNormal);
	    slipMag = slipMagCur;
	    if (slipMag > 0.0) {
	      // Use Newton (in log slip space) to get better update in slip & traction.
	      // D_{i+1} = exp(ln(D_i) - (T-T_f)/(D_i * (jacobian - frictionDeriv))
	      slipMagCur = exp(log(slipMag) - (tractionShearMagCur - frictionStress) / (slipMag * (jacobianShear - frictionDeriv)));
	    } else {
	      // Use Newton (in linear slip space) to get better update in slip & traction.
	      // D_{i+1} = D_i - (T-T_f)/(jacobian - frictionDeriv)
	      slipMagCur = slipMag - (tractionShearMagCur - frictionStress) / (jacobianShear - frictionDeriv);
	    } // if
	    tractionShearMagCur += (slipMagCur - slipMag) * jacobianShear;
	    slipRateMagCur = (slipMagCur - slipMag0) / _dt;
	    frictionStress = _friction->calcFriction(t, slipMagCur, slipRateMagCur, tractionNormal);
	    if (fabs(tractionShearMagCur - frictionStress) < _zeroTolerance) {
	      break;
	    } // if
	  } // for
	} // if
#endif

	// Update traction increment based on value required to stick
	// versus friction
	const PylithScalar dlp = -(tractionShearMag - frictionStress) * tractionTpdt[0] / tractionShearMag;
	const PylithScalar dlq = -(tractionShearMag - frictionStress) * tractionTpdt[1] / tractionShearMag;
	
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
