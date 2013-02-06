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

#include "FaultCohesiveDynKin.hh" // implementation of object methods

#include "CohesiveTopology.hh" // USES CohesiveTopology
#include "TractPerturbation.hh" // HOLDSA TractPerturbation
#include "DKSelector.hh" // USES DKSelector
#include "EqKinSrc.hh" // USES EqKinSrc

#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/feassemble/CellGeometry.hh" // USES CellGeometry
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/SubMesh.hh" // USES SubMesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/FieldsNew.hh" // USES FieldsNew
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
typedef pylith::topology::Mesh::RealSection RealSection;
typedef pylith::topology::Mesh::RealUniformSection RealUniformSection;
typedef pylith::topology::SubMesh::SieveMesh SieveSubMesh;

typedef pylith::topology::Field<pylith::topology::SubMesh>::RestrictVisitor RestrictVisitor;
typedef pylith::topology::Field<pylith::topology::SubMesh>::UpdateAddVisitor UpdateAddVisitor;
typedef ALE::ISieveVisitor::IndicesVisitor<RealSection,SieveSubMesh::order_type,PylithInt> IndicesVisitor;

// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::FaultCohesiveDynKin::FaultCohesiveDynKin(void) :
  _zeroTolerance(1.0e-10),
  _openFreeSurf(true),
  _dkSelector(0),
  _tractPerturbation(0),
  _friction(0),
  _jacobian(0),
  _ksp(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::faults::FaultCohesiveDynKin::~FaultCohesiveDynKin(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void pylith::faults::FaultCohesiveDynKin::deallocate(void)
{ // deallocate
  FaultCohesiveLagrange::deallocate();

  _tractPerturbation = 0; // :TODO: Use shared pointer
  _dkSelector = 0; // :TODO: Used shared pointer
  _friction = 0; // :TODO: Use shared pointer

  delete _jacobian; _jacobian = 0;
  PetscErrorCode err = KSPDestroy(&_ksp);CHECK_PETSC_ERROR(err);
} // deallocate

// ----------------------------------------------------------------------
// Sets the spatial database for the inital tractions
void
pylith::faults::FaultCohesiveDynKin::tractPerturbation(TractPerturbation* tract)
{ // tractPerturbation
  _tractPerturbation = tract;
} // tractPerturbation

// ----------------------------------------------------------------------
// Sets the spatial database for the dynamic-kinematic selector
void
pylith::faults::FaultCohesiveDynKin::dkSelector(DKSelector* dksel)
{ // DKSelector
  _dkSelector = dksel;
} // dkSelector

// ----------------------------------------------------------------------
// Get the friction (constitutive) model.  
void
pylith::faults::FaultCohesiveDynKin::frictionModel(friction::FrictionModel* const model)
{ // frictionModel
  _friction = model;
} // frictionModel

// ----------------------------------------------------------------------
// Set kinematic earthquake source.
void
pylith::faults::FaultCohesiveDynKin::eqsrcs(const char* const * names,
                                         const int numNames,
                                         EqKinSrc** sources,
                                         const int numSources)
{ // eqsrcs
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
// Nondimensional tolerance for detecting near zero values.
void
pylith::faults::FaultCohesiveDynKin::zeroTolerance(const PylithScalar value)
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
pylith::faults::FaultCohesiveDynKin::openFreeSurf(const bool value)
{ // openFreeSurf
  _openFreeSurf = value;
} // openFreeSurf

// ----------------------------------------------------------------------
// Initialize fault. Determine orientation and setup boundary
void
pylith::faults::FaultCohesiveDynKin::initialize(const topology::Mesh& mesh,
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

  // Get the DynKin Selector using a spatial database.
  if (_dkSelector) {
    _dkSelector->initialize(*_faultMesh, *_normalizer);
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

  // Initialize the eqSrcs
  const srcs_type::const_iterator srcsEnd = _eqSrcs.end();
  for (srcs_type::iterator s_iter = _eqSrcs.begin(); s_iter != srcsEnd; ++s_iter) {
    EqKinSrc* src = s_iter->second;
    assert(0 != src);
    src->initialize(*_faultMesh, *_normalizer);
  } // for

  // Get the DynKin Selector using a spatial database. 
  if (_dkSelector) {
     _dkSelector->initialize(*_faultMesh, *_normalizer);
     _fields->add("Dynamic Kinematic Selector","dynamic_kinematic_selector");
  } // if should be useless
     
  logger.stagePop();
} // initialize

// ----------------------------------------------------------------------
// Integrate contributions to residual term (r) for operator.
void
pylith::faults::FaultCohesiveDynKin::integrateResidual(
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
  scalar_array residualVertexN(spaceDim);
  scalar_array residualVertexP(spaceDim);
  scalar_array residualVertexL(spaceDim);
  const ALE::Obj<RealSection>& residualSection = residual.section();
  assert(!residualSection.isNull());

  const ALE::Obj<RealSection>& dispTSection = fields->get("disp(t)").section();
  assert(!dispTSection.isNull());

  const ALE::Obj<RealSection>& dispTIncrSection = 
    fields->get("dispIncr(t->t+dt)").section();
  assert(!dispTIncrSection.isNull());

  scalar_array dispTpdtVertexN(spaceDim);
  scalar_array dispTpdtVertexP(spaceDim);
  scalar_array dispTpdtVertexL(spaceDim);

  scalar_array tractPerturbVertex(spaceDim);
  ALE::Obj<RealUniformSection> tractPerturbSection;
  int tractPerturbIndex = 0;
  int tractPerturbFiberDim = 0;
  if (_tractPerturbation) {
    _tractPerturbation->calculate(t);
    
    const topology::FieldsNew<topology::SubMesh>* params = _tractPerturbation->parameterFields();
    assert(params);
    tractPerturbSection = params->section();
    assert(!tractPerturbSection.isNull());

    tractPerturbFiberDim = params->fiberDim();
    tractPerturbIndex = params->sectionIndex("value");
    const int valueFiberDim = params->sectionFiberDim("value");
    assert(valueFiberDim == spaceDim);
  } // if

  const ALE::Obj<RealSection>& dispRelSection =
    _fields->get("relative disp").section();
  assert(!dispRelSection.isNull());

  const ALE::Obj<RealSection>& areaSection = _fields->get("area").section();
  assert(!areaSection.isNull());

  const ALE::Obj<RealSection>& orientationSection = 
    _fields->get("orientation").section();
  assert(!orientationSection.isNull());

  // Get fault information
  const ALE::Obj<SieveMesh>& sieveMesh = fields->mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::order_type>& globalOrder =
      sieveMesh->getFactory()->getGlobalOrder(sieveMesh, "default",
					      residualSection);
  assert(!globalOrder.isNull());

  // Get the dkSelector
  //if (_dkSelector) {
  topology::Field<topology::SubMesh>& dk = _fields->get("Dynamic Kinematic Selector");
  _dkSelector->dk(&dk);
  const ALE::Obj<RealSection>& dkSelSection = dk.section();
  assert(!dkSelSection.isNull());
  //}
  //else { // should never be here
  //  std::ostringstream msg;                                           
  //  msg << "No Dynamic Kinematic Selector available.";         
  //  throw std::runtime_error(msg.str());
  //} 

  _logger->eventEnd(setupEvent);
#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventBegin(computeEvent);
#endif

  // Loop over fault vertices
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

    // Common needs
    
    // Get area associated with fault vertex.
    assert(1 == areaSection->getFiberDimension(v_fault));
    assert(areaSection->restrictPoint(v_fault));
    const PylithScalar areaVertex = *areaSection->restrictPoint(v_fault);

    // Get disp(t) at conventional vertices and Lagrange vertex.
    assert(spaceDim == dispTSection->getFiberDimension(v_negative));
    const PylithScalar* dispTVertexN = dispTSection->restrictPoint(v_negative);
    assert(dispTVertexN);

    assert(spaceDim == dispTSection->getFiberDimension(v_positive));
    const PylithScalar* dispTVertexP = dispTSection->restrictPoint(v_positive);
    assert(dispTVertexP);

    assert(spaceDim == dispTSection->getFiberDimension(v_lagrange));
    const PylithScalar* dispTVertexL = dispTSection->restrictPoint(v_lagrange);
    assert(dispTVertexL);    

    // Get dispIncr(t->t+dt) at conventional vertices and Lagrange vertex.
    assert(spaceDim == dispTIncrSection->getFiberDimension(v_negative));
    const PylithScalar* dispTIncrVertexN =
      dispTIncrSection->restrictPoint(v_negative);
    assert(dispTIncrVertexN);

    assert(spaceDim == dispTIncrSection->getFiberDimension(v_positive));
    const PylithScalar* dispTIncrVertexP =
      dispTIncrSection->restrictPoint(v_positive);
    assert(dispTIncrVertexP);

    assert(spaceDim == dispTIncrSection->getFiberDimension(v_lagrange));
    const PylithScalar* dispTIncrVertexL =
      dispTIncrSection->restrictPoint(v_lagrange);
    assert(dispTIncrVertexL);

    if (dkSelSection[iVertex] > 0.5) { // Kinematic Case

      // Get relative dislplacement at fault vertex.
      assert(spaceDim == dispRelSection->getFiberDimension(v_fault));
      const PylithScalar* dispRelVertex = dispRelSection->restrictPoint(v_fault);
      assert(dispRelVertex);

    } else { // Dynamic Case

      // Get prescribed traction perturbation at fault vertex.
      if (_tractPerturbation) {
        assert(tractPerturbFiberDim == tractPerturbSection->getFiberDimension(v_fault));
        const PylithScalar* paramsVertex = tractPerturbSection->restrictPoint(v_fault);
        assert(paramsVertex);
  
        for (int iDim=0; iDim < spaceDim; ++iDim) {
          tractPerturbVertex[iDim] = paramsVertex[tractPerturbIndex+iDim];
        } // for
      } else {
        tractPerturbVertex = 0.0;
      } // if/else

      // Get orientation associated with fault vertex.
      assert(spaceDim*spaceDim == orientationSection->getFiberDimension(v_fault));
      const PylithScalar* orientationVertex = orientationSection->restrictPoint(v_fault);
      assert(orientationVertex);

    } // if/else

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(restrictEvent);
    _logger->eventBegin(computeEvent);
#endif

    // Compute current estimate of displacement at time t+dt using
    // solution increment.
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      dispTpdtVertexN[iDim] = dispTVertexN[iDim] + dispTIncrVertexN[iDim];
      dispTpdtVertexP[iDim] = dispTVertexP[iDim] + dispTIncrVertexP[iDim];
      dispTpdtVertexL[iDim] = dispTVertexL[iDim] + dispTIncrVertexL[iDim];
    } // for
    
    if (dkSelSection[iVertex] > 0.5) { // Kinematic Case

      residualVertexN = areaVertex * dispTpdtVertexL;
      residualVertexP = -residualVertexN;
  
      residualVertexL = 0.0;
      for (int iDim=0; iDim < spaceDim; ++iDim) {
        residualVertexL[iDim] = -areaVertex * 
          (dispTpdtVertexP[iDim] - dispTpdtVertexN[iDim] - dispRelVertex[iDim]);
      } // for

    } else { // Dynamic Case

      // Compute slip (in fault coordinates system) from displacements.
      PylithScalar slipNormal = 0.0;
      PylithScalar tractionNormal = 0.0;
      const int indexN = spaceDim - 1;
      for (int jDim=0; jDim < spaceDim; ++jDim) {
        slipNormal += orientationVertex[indexN*spaceDim+jDim] * 
        	(dispTpdtVertexP[jDim] - dispTpdtVertexN[jDim]);
        tractionNormal += 
        	orientationVertex[indexN*spaceDim+jDim] * dispTpdtVertexL[jDim];
      } // for
    
      residualVertexN = 0.0;
      residualVertexL = 0.0;
      if (slipNormal < _zeroTolerance || !_openFreeSurf) { 
        // if no opening or flag indicates to still impose initial
        // tractions when fault is open.
        //
        // Initial (external) tractions oppose (internal) tractions
        // associated with Lagrange multiplier.
        residualVertexN = areaVertex * (dispTpdtVertexL - tractPerturbVertex);

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
      residualVertexP = -residualVertexN;

    } // if/else

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(computeEvent);
    _logger->eventBegin(updateEvent);
#endif

    // Assemble contributions into field
    assert(residualVertexN.size() == 
	   residualSection->getFiberDimension(v_negative));
    residualSection->updateAddPoint(v_negative, &residualVertexN[0]);

    assert(residualVertexP.size() == 
	   residualSection->getFiberDimension(v_positive));
    residualSection->updateAddPoint(v_positive, &residualVertexP[0]);

    if (dkSelSection[iVertex] > 0.5) { // Kinematic Case

      assert(residualVertexL.size() == 
            residualSection->getFiberDimension(v_lagrange));
      residualSection->updateAddPoint(v_lagrange, &residualVertexL[0]);

    } // if

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(updateEvent);
#endif
  } // for
  PetscLogFlops(numVertices*spaceDim*8);

#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventEnd(computeEvent);
#endif
} // integrateResidual

// ----------------------------------------------------------------------
// Update state variables as needed.
void
pylith::faults::FaultCohesiveDynKin::updateStateVars(
				      const PylithScalar t,
				      topology::SolutionFields* const fields)
{ // updateStateVars
  assert(fields);
  assert(_fields);

  _updateRelMotion(*fields);

  const int spaceDim = _quadrature->spaceDim();

  // Allocate arrays for vertex values
  scalar_array tractionTpdtVertex(spaceDim); // Fault coordinate system

  // Get sections
  scalar_array slipVertex(spaceDim);
  const ALE::Obj<RealSection>& dispRelSection = 
    _fields->get("relative disp").section();
  assert(!dispRelSection.isNull());

  scalar_array slipRateVertex(spaceDim);
  const ALE::Obj<RealSection>& velRelSection =
      _fields->get("relative velocity").section();
  assert(!velRelSection.isNull());

  const ALE::Obj<RealSection>& dispTSection = fields->get("disp(t)").section();
  assert(!dispTSection.isNull());

  const ALE::Obj<RealSection>& dispTIncrSection =
      fields->get("dispIncr(t->t+dt)").section();
  assert(!dispTIncrSection.isNull());

  const ALE::Obj<RealSection>& orientationSection =
      _fields->get("orientation").section();
  assert(!orientationSection.isNull());

  const int numVertices = _cohesiveVertices.size();
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_lagrange = _cohesiveVertices[iVertex].lagrange;
    const int v_fault = _cohesiveVertices[iVertex].fault;
    const int v_negative = _cohesiveVertices[iVertex].negative;
    const int v_positive = _cohesiveVertices[iVertex].positive;

    // Get relative displacement
    assert(spaceDim == dispRelSection->getFiberDimension(v_fault));
    const PylithScalar* dispRelVertex = dispRelSection->restrictPoint(v_fault);
    assert(dispRelVertex);

    // Get relative velocity
    assert(spaceDim == velRelSection->getFiberDimension(v_fault));
    const PylithScalar* velRelVertex = velRelSection->restrictPoint(v_fault);
    assert(velRelVertex);

    // Get orientation
    assert(spaceDim*spaceDim == orientationSection->getFiberDimension(v_fault));
    const PylithScalar* orientationVertex = 
      orientationSection->restrictPoint(v_fault);

    // Get Lagrange multiplier values from disp(t), and dispIncr(t->t+dt)
    assert(spaceDim == dispTSection->getFiberDimension(v_lagrange));
    const PylithScalar* lagrangeTVertex = dispTSection->restrictPoint(v_lagrange);
    assert(spaceDim == dispTIncrSection->getFiberDimension(v_lagrange));
    const PylithScalar* lagrangeTIncrVertex = 
      dispTIncrSection->restrictPoint(v_lagrange);

    // Compute slip, slip rate, and fault traction (Lagrange
    // multiplier) at time t+dt in fault coordinate system.
    slipVertex = 0.0;
    slipRateVertex = 0.0;
    tractionTpdtVertex = 0.0;
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      for (int jDim=0; jDim < spaceDim; ++jDim) {
	slipVertex[iDim] += orientationVertex[iDim*spaceDim+jDim] *
	  dispRelVertex[jDim];
	slipRateVertex[iDim] += orientationVertex[iDim*spaceDim+jDim] *
	  velRelVertex[jDim];
	tractionTpdtVertex[iDim] += orientationVertex[iDim*spaceDim+jDim] *
	  (lagrangeTVertex[jDim]+lagrangeTIncrVertex[jDim]);
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
      _friction->updateStateVars(t, slipMag, slipRateMag, tractionNormal, 
				 v_fault);
      break;
    } // case 1
    case 2: { // case 2
      const PylithScalar slipMag = fabs(slipVertex[0]);
      const PylithScalar slipRateMag = fabs(slipRateVertex[0]);
      const PylithScalar tractionNormal = tractionTpdtVertex[1];
      _friction->updateStateVars(t, slipMag, slipRateMag, tractionNormal, 
				 v_fault);
      break;
    } // case 2
    case 3: { // case 3
      const PylithScalar slipMag = 
	sqrt(slipVertex[0]*slipVertex[0] + slipVertex[1]*slipVertex[1]);
      const PylithScalar slipRateMag = 
	sqrt(slipRateVertex[0]*slipRateVertex[0] + 
	     slipRateVertex[1]*slipRateVertex[1]);
      const PylithScalar tractionNormal = tractionTpdtVertex[2];
      _friction->updateStateVars(t, slipMag, slipRateMag, tractionNormal, 
				 v_fault);
      break;
    } // case 3
    default:
      assert(0);
      throw std::logic_error("Unknown spatial dimension in "
			     "FaultCohesiveDynKin::updateStateVars().");
    } // switch
  } // for
} // updateStateVars

// ----------------------------------------------------------------------
// Constrain solution based on friction.
void
pylith::faults::FaultCohesiveDynKin::constrainSolnSpace(
				    topology::SolutionFields* const fields,
				    const PylithScalar t,
				    const topology::Jacobian& jacobian)
{ // constrainSolnSpace
  /// Member prototype for _constrainSolnSpaceXD()
  typedef void (pylith::faults::FaultCohesiveDynKin::*constrainSolnSpace_fn_type)
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

  // Get the dkSelector
  if (_dkSelector) {
    topology::Field<topology::SubMesh>& dk = _fields->get("Dynamic Kinematic Selector");
    _dkSelector->dk(&dk);
    const ALE::Obj<RealSection>& dkSelSection = dk.section();
    assert(!dkSelSection.isNull());
  } else { // should never be here
    std::ostringstream msg;
    msg << "No Dynamic Kinematic Selector available.";
    throw std::runtime_error(msg.str());
  }

  scalar_array slipTpdtVertex(spaceDim);
  const ALE::Obj<RealSection>& dispRelSection = 
    _fields->get("relative disp").section();
  assert(!dispRelSection.isNull());

  scalar_array slipRateVertex(spaceDim);
  const ALE::Obj<RealSection>& velRelSection =
      _fields->get("relative velocity").section();
  assert(!velRelSection.isNull());

  const ALE::Obj<RealSection>& orientationSection =
      _fields->get("orientation").section();
  assert(!orientationSection.isNull());

  const ALE::Obj<RealSection>& dispTSection = fields->get("disp(t)").section();
  assert(!dispTSection.isNull());

  scalar_array dDispTIncrVertexN(spaceDim);
  scalar_array dDispTIncrVertexP(spaceDim);
  const ALE::Obj<RealSection>& dispIncrSection =
      fields->get("dispIncr(t->t+dt)").section();
  assert(!dispIncrSection.isNull());

  const ALE::Obj<RealSection>& dispIncrAdjSection =
      fields->get("dispIncr adjust").section();
  assert(!dispIncrAdjSection.isNull());

  scalar_array dTractionTpdtVertex(spaceDim);
  scalar_array dLagrangeTpdtVertex(spaceDim);
  const ALE::Obj<RealSection>& dLagrangeTpdtSection =
      _fields->get("sensitivity dLagrange").section();
  assert(!dLagrangeTpdtSection.isNull());

  constrainSolnSpace_fn_type constrainSolnSpaceFn;
  switch (spaceDim) { // switch
  case 1:
    constrainSolnSpaceFn = 
      &pylith::faults::FaultCohesiveDynKin::_constrainSolnSpace1D;
    break;
  case 2: 
    constrainSolnSpaceFn = 
      &pylith::faults::FaultCohesiveDynKin::_constrainSolnSpace2D;
    break;
  case 3:
    constrainSolnSpaceFn = 
      &pylith::faults::FaultCohesiveDynKin::_constrainSolnSpace3D;
    break;
  default :
    assert(0);
    throw std::logic_error("Unknown spatial dimension in "
			   "FaultCohesiveDynKin::constrainSolnSpace().");
  } // switch


#if 0 // DEBUGGING
  dispRelSection->view("BEFORE RELATIVE DISPLACEMENT");
  dispIncrSection->view("BEFORE DISP INCR (t->t+dt)");
#endif

  const int numVertices = _cohesiveVertices.size();
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    if (dkSelSection[iVertex] < 0.5) // Dynamic Case
      continue;

    const int v_lagrange = _cohesiveVertices[iVertex].lagrange;
    const int v_fault = _cohesiveVertices[iVertex].fault;
    const int v_negative = _cohesiveVertices[iVertex].negative;
    const int v_positive = _cohesiveVertices[iVertex].positive;

    // Get displacement values
    assert(spaceDim == dispTSection->getFiberDimension(v_negative));
    const PylithScalar* dispTVertexN = dispTSection->restrictPoint(v_negative);
    assert(dispTVertexN);

    assert(spaceDim == dispTSection->getFiberDimension(v_positive));
    const PylithScalar* dispTVertexP = dispTSection->restrictPoint(v_positive);
    assert(dispTVertexP);

    // Get displacement increment values.
    assert(spaceDim == dispIncrSection->getFiberDimension(v_negative));
    const PylithScalar* dispIncrVertexN = dispIncrSection->restrictPoint(v_negative);
    assert(dispIncrVertexN);

    assert(spaceDim == dispIncrSection->getFiberDimension(v_positive));
    const PylithScalar* dispIncrVertexP = dispIncrSection->restrictPoint(v_positive);
    assert(dispIncrVertexP);

    // Get orientation
    assert(spaceDim*spaceDim == orientationSection->getFiberDimension(v_fault));
    const PylithScalar* orientationVertex = orientationSection->restrictPoint(v_fault);

    // Get Lagrange multiplier values from disp(t), and dispIncr(t->t+dt)
    assert(spaceDim == dispTSection->getFiberDimension(v_lagrange));
    const PylithScalar* lagrangeTVertex = dispTSection->restrictPoint(v_lagrange);
    assert(lagrangeTVertex);

    assert(spaceDim == dispIncrSection->getFiberDimension(v_lagrange));
    const PylithScalar* lagrangeTIncrVertex = dispIncrSection->restrictPoint(v_lagrange);
    assert(lagrangeTIncrVertex);

    // Step 1: Prevent nonphysical trial solutions. The product of the
    // normal traction and normal slip must be nonnegative (forbid
    // interpenetration with tension or opening with compression).
    
    // Compute slip, slip rate, and Lagrange multiplier at time t+dt
    // in fault coordinate system.
    slipTpdtVertex = 0.0;
    slipRateVertex = 0.0;
    tractionTpdtVertex = 0.0;
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      for (int jDim=0; jDim < spaceDim; ++jDim) {
	slipTpdtVertex[iDim] += orientationVertex[iDim*spaceDim+jDim] *
	  (dispTVertexP[jDim] + dispIncrVertexP[jDim]
	   - dispTVertexN[jDim] - dispIncrVertexN[jDim]);
	slipRateVertex[iDim] += orientationVertex[iDim*spaceDim+jDim] *
	  (dispIncrVertexP[jDim] - dispIncrVertexN[jDim]) / dt;
	tractionTpdtVertex[iDim] += orientationVertex[iDim*spaceDim+jDim] *
	  (lagrangeTVertex[jDim] + lagrangeTIncrVertex[jDim]);
      } // for
      if (fabs(slipRateVertex[iDim]) < _zeroTolerance) {
	slipRateVertex[iDim] = 0.0;
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
    CALL_MEMBER_FN(*this,
		   constrainSolnSpaceFn)(&dTractionTpdtVertex,
					 t, slipTpdtVertex, slipRateVertex,
					 tractionTpdtVertex,
					 iterating);

    // Rotate increment in traction back to global coordinate system.
    dLagrangeTpdtVertex = 0.0;
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      for (int jDim=0; jDim < spaceDim; ++jDim) {
	dLagrangeTpdtVertex[iDim] += 
	  orientationVertex[jDim*spaceDim+iDim] * dTractionTpdtVertex[jDim];
      } // for

      // Add in potential contribution from adjusting Lagrange
      // multiplier for fault normal DOF of trial solution in Step 1.
      dLagrangeTpdtVertex[iDim] += 
	orientationVertex[indexN*spaceDim+iDim] * dTractionTpdtVertexNormal;
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
    assert(dLagrangeTpdtVertex.size() ==
        dLagrangeTpdtSection->getFiberDimension(v_fault));
    dLagrangeTpdtSection->updatePoint(v_fault, &dLagrangeTpdtVertex[0]);

#if 0 // UNNECESSARY?
    // Update displacement in trial solution (if necessary) so that it
    // conforms to physical constraints.
    if (0.0 != dSlipVertexNormal) {
      // Compute relative displacement from slip.
      dDispRelVertex = 0.0;
      for (int iDim=0; iDim < spaceDim; ++iDim) {
	dDispRelVertex[iDim] += 
	  orientationVertex[indexN*spaceDim+iDim] * dSlipVertexNormal;

      dDispTIncrVertexN[iDim] = -0.5*dDispRelVertex[iDim];
      dDispTIncrVertexP[iDim] = +0.5*dDispRelVertex[iDim];
      } // for

      // Update displacement field
      assert(dDispTIncrVertexN.size() ==
	     dispIncrSection->getFiberDimension(v_negative));
      dispIncrAdjSection->updateAddPoint(v_negative, &dDispTIncrVertexN[0]);
      
      assert(dDispTIncrVertexP.size() ==
	     dispIncrSection->getFiberDimension(v_positive));
      dispIncrAdjSection->updateAddPoint(v_positive, &dDispTIncrVertexP[0]);    
    } // if
#endif

  } // for

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
  const ALE::Obj<RealSection>& sensDispRelSection = _fields->get("sensitivity relative disp").section();
  assert(!sensDispRelSection.isNull());

  const ALE::Obj<SieveMesh>& sieveMesh = fields->mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::order_type>& globalOrder =
    sieveMesh->getFactory()->getGlobalOrder(sieveMesh, "default", 
					    dispIncrSection);
  assert(!globalOrder.isNull());

  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    if (dkSelSection[iVertex] < 0.5 ) // Dynamic Case
      continue;

    const int v_fault = _cohesiveVertices[iVertex].fault;
    const int v_lagrange = _cohesiveVertices[iVertex].lagrange;
    const int v_negative = _cohesiveVertices[iVertex].negative;
    const int v_positive = _cohesiveVertices[iVertex].positive;

    // Get change in Lagrange multiplier computed from friction criterion.
    dLagrangeTpdtSection->restrictPoint(v_fault, &dLagrangeTpdtVertex[0],
					dLagrangeTpdtVertex.size());

    // Get change in relative displacement from sensitivity solve.
    assert(spaceDim == sensDispRelSection->getFiberDimension(v_fault));
    const PylithScalar* sensDispRelVertex = 
      sensDispRelSection->restrictPoint(v_fault);
    assert(sensDispRelVertex);

    // Get current relative displacement for updating.
    dispRelSection->restrictPoint(v_fault, &dispRelVertex[0],
				  dispRelVertex.size());

    // Get orientation.
    assert(spaceDim*spaceDim == orientationSection->getFiberDimension(v_fault));
    const PylithScalar* orientationVertex = 
      orientationSection->restrictPoint(v_fault);
    assert(orientationVertex);

    // Get displacement.
    assert(spaceDim == dispTSection->getFiberDimension(v_negative));
    const PylithScalar* dispTVertexN = dispTSection->restrictPoint(v_negative);
    assert(dispTVertexN);

    assert(spaceDim == dispTSection->getFiberDimension(v_positive));
    const PylithScalar* dispTVertexP = dispTSection->restrictPoint(v_positive);
    assert(dispTVertexP);

    // Get displacement increment (trial solution).
    assert(spaceDim == dispIncrSection->getFiberDimension(v_negative));
    const PylithScalar* dispIncrVertexN = 
      dispIncrSection->restrictPoint(v_negative);
    assert(dispIncrVertexN);

    assert(spaceDim == dispIncrSection->getFiberDimension(v_positive));
    const PylithScalar* dispIncrVertexP = 
      dispIncrSection->restrictPoint(v_positive);
    assert(dispIncrVertexP);

    // Get Lagrange multiplier at time t
    assert(spaceDim == dispTSection->getFiberDimension(v_lagrange));
    const PylithScalar* lagrangeTVertex = dispTSection->restrictPoint(v_lagrange);
    assert(lagrangeTVertex);

    // Get Lagrange multiplier increment (trial solution)
    assert(spaceDim == dispIncrSection->getFiberDimension(v_lagrange));
    const PylithScalar* lagrangeTIncrVertex = 
      dispIncrSection->restrictPoint(v_lagrange);
    assert(lagrangeTIncrVertex);

    // Scale perturbation in relative displacements and change in
    // Lagrange multipliers by alpha using only shear components.
    slipTVertex = 0.0;
    slipTpdtVertex = 0.0;
    dSlipTpdtVertex = 0.0;
    tractionTpdtVertex = 0.0;
    dTractionTpdtVertex = 0.0;
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      for (int jDim=0; jDim < spaceDim; ++jDim) {
	slipTVertex[iDim] += orientationVertex[iDim*spaceDim+jDim] * 
	  (dispTVertexP[jDim] - dispTVertexN[jDim]);
	slipTpdtVertex[iDim] += orientationVertex[iDim*spaceDim+jDim] * 
	  (dispTVertexP[jDim] - dispTVertexN[jDim] +
	   dispIncrVertexP[jDim] - dispIncrVertexN[jDim]);
	dSlipTpdtVertex[iDim] += orientationVertex[iDim*spaceDim+jDim] * 
	  alpha*sensDispRelVertex[jDim];
	tractionTpdtVertex[iDim] += orientationVertex[iDim*spaceDim+jDim] *
	  (lagrangeTVertex[jDim] + lagrangeTIncrVertex[jDim]);
	dTractionTpdtVertex[iDim] += orientationVertex[iDim*spaceDim+jDim] * 
	  alpha*dLagrangeTpdtVertex[jDim];
      } // for
    } // for

    // FIRST, correct nonphysical trial solutions.
    // Order of steps 5a-5c is important!
    if ((slipTpdtVertex[indexN] + dSlipTpdtVertex[indexN]) * 
	(tractionTpdtVertex[indexN] + dTractionTpdtVertex[indexN])
	< 0.0) {
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
      if (fabs(slipTpdtVertex[indexN] + dSlipTpdtVertex[indexN]) > 
	  fabs(tractionTpdtVertex[indexN] + dTractionTpdtVertex[indexN])) {
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
      slipDot += 
	(slipTpdtVertex[iDim] - slipTVertex[iDim]) * 
	(slipTpdtVertex[iDim] + dSlipTpdtVertex[iDim] - slipTVertex[iDim]);
      // Compute dot product of traction and slip
      tractionSlipDot += (tractionTpdtVertex[iDim] + dTractionTpdtVertex[iDim])
	* (slipTpdtVertex[iDim] + dSlipTpdtVertex[iDim]);
    } // for
    if (slipDot < 0.0 &&
	sqrt(fabs(slipDot)) > _zeroTolerance && 
	tractionSlipDot < 0.0) {
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
	dispRelVertex[iDim] += orientationVertex[jDim*spaceDim+iDim] *
	  slipTpdtVertex[jDim];
	dDispRelVertex[iDim] += orientationVertex[jDim*spaceDim+iDim] *
	  dSlipTpdtVertex[jDim];
	dLagrangeTpdtVertex[iDim] += orientationVertex[jDim*spaceDim+iDim] * 
	  dTractionTpdtVertex[jDim];
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
    if (globalOrder->isLocal(v_lagrange)) {

      // Update Lagrange multiplier increment.
      assert(dLagrangeTpdtVertex.size() ==
	     dispIncrSection->getFiberDimension(v_lagrange));
      dispIncrAdjSection->updateAddPoint(v_lagrange, &dLagrangeTpdtVertex[0]);

      // Update displacement field
      assert(dDispTIncrVertexN.size() ==
	     dispIncrSection->getFiberDimension(v_negative));
      dispIncrAdjSection->updateAddPoint(v_negative, &dDispTIncrVertexN[0]);
      
      assert(dDispTIncrVertexP.size() ==
	     dispIncrSection->getFiberDimension(v_positive));
      dispIncrAdjSection->updateAddPoint(v_positive, &dDispTIncrVertexP[0]);
    } // if

  } // for

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
pylith::faults::FaultCohesiveDynKin::adjustSolnLumped(
			 topology::SolutionFields* const fields,
			 const PylithScalar t,
			 const topology::Field<topology::Mesh>& jacobian)
{ // adjustSolnLumped
  /// Member prototype for _constrainSolnSpaceXD()
  typedef void (pylith::faults::FaultCohesiveDynKin::*constrainSolnSpace_fn_type)
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

  // Update time step in friction (can vary).
  _friction->timeStep(_dt);

  // Get section information
  scalar_array dispRelVertex(spaceDim);
  scalar_array slipVertex(spaceDim);
  const ALE::Obj<RealSection>& dispRelSection = 
    _fields->get("relative disp").section();
  assert(!dispRelSection.isNull());

  scalar_array slipRateVertex(spaceDim);
  const ALE::Obj<RealSection>& velRelSection =
      _fields->get("relative velocity").section();
  assert(!velRelSection.isNull());

  const ALE::Obj<RealSection>& orientationSection =
      _fields->get("orientation").section();
  assert(!orientationSection.isNull());

  const ALE::Obj<RealSection>& areaSection = _fields->get("area").section();
  assert(!areaSection.isNull());

  const ALE::Obj<RealSection>& dispTSection = fields->get("disp(t)").section();
  assert(!dispTSection.isNull());

  scalar_array dispIncrVertexN(spaceDim);
  scalar_array dispIncrVertexP(spaceDim);
  scalar_array lagrangeTIncrVertex(spaceDim);
  const ALE::Obj<RealSection>& dispIncrSection =
      fields->get("dispIncr(t->t+dt)").section();
  assert(!dispIncrSection.isNull());

  const ALE::Obj<RealSection>& dispIncrAdjSection = fields->get(
    "dispIncr adjust").section();
  assert(!dispIncrAdjSection.isNull());

  const ALE::Obj<RealSection>& jacobianSection = jacobian.section();
  assert(!jacobianSection.isNull());

  const ALE::Obj<RealSection>& residualSection =
      fields->get("residual").section();

  const ALE::Obj<SieveMesh>& sieveMesh = fields->mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::order_type>& globalOrder =
    sieveMesh->getFactory()->getGlobalOrder(sieveMesh, "default", 
					    jacobianSection);
  assert(!globalOrder.isNull());

  constrainSolnSpace_fn_type constrainSolnSpaceFn;
  switch (spaceDim) { // switch
  case 1:
    constrainSolnSpaceFn = 
      &pylith::faults::FaultCohesiveDynKin::_constrainSolnSpace1D;
    break;
  case 2: 
    constrainSolnSpaceFn = 
      &pylith::faults::FaultCohesiveDynKin::_constrainSolnSpace2D;
    break;
  case 3:
    constrainSolnSpaceFn = 
      &pylith::faults::FaultCohesiveDynKin::_constrainSolnSpace3D;
    break;
  default :
    assert(0);
    throw std::logic_error("Unknown spatial dimension in "
			   "FaultCohesiveDynKin::adjustSolnLumped.");
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

    // Get residual at cohesive cell's vertices.
    assert(spaceDim == residualSection->getFiberDimension(v_lagrange));
    const PylithScalar* residualVertexL = residualSection->restrictPoint(v_lagrange);
    assert(residualVertexL);

    // Get jacobian at cohesive cell's vertices.
    assert(spaceDim == jacobianSection->getFiberDimension(v_negative));
    const PylithScalar* jacobianVertexN = jacobianSection->restrictPoint(v_negative);
    assert(jacobianVertexN);

    assert(spaceDim == jacobianSection->getFiberDimension(v_positive));
    const PylithScalar* jacobianVertexP = jacobianSection->restrictPoint(v_positive);
    assert(jacobianVertexP);

    // Get area at fault vertex.
    assert(1 == areaSection->getFiberDimension(v_fault));
    assert(areaSection->restrictPoint(v_fault));
    const PylithScalar areaVertex = *areaSection->restrictPoint(v_fault);
    assert(areaVertex > 0.0);

    // Get disp(t) at Lagrange vertex.
    assert(spaceDim == dispTSection->getFiberDimension(v_lagrange));
    const PylithScalar* lagrangeTVertex = dispTSection->restrictPoint(v_lagrange);
    assert(lagrangeTVertex);

    // Get dispIncr(t) at cohesive cell's vertices.
    dispIncrSection->restrictPoint(v_negative, &dispIncrVertexN[0],
				    dispIncrVertexN.size());
    dispIncrSection->restrictPoint(v_positive, &dispIncrVertexP[0],
				    dispIncrVertexP.size());
    dispIncrSection->restrictPoint(v_lagrange, &lagrangeTIncrVertex[0],
				    lagrangeTIncrVertex.size());

    // Get relative displacement at fault vertex.
    dispRelSection->restrictPoint(v_fault, &dispRelVertex[0], 
				  dispRelVertex.size());

    // Get relative velocity at fault vertex.
    assert(spaceDim == velRelSection->getFiberDimension(v_fault));
    const PylithScalar* velRelVertex = velRelSection->restrictPoint(v_fault);
    assert(velRelVertex);
    
    // Get fault orientation at fault vertex.
    assert(spaceDim*spaceDim == orientationSection->getFiberDimension(v_fault));
    const PylithScalar* orientationVertex = 
      orientationSection->restrictPoint(v_fault);
    assert(orientationVertex);
    

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(restrictEvent);
    _logger->eventBegin(computeEvent);
#endif

    // Adjust solution as in prescribed rupture, updating the Lagrange
    // multipliers and the corresponding displacment increments.
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      assert(jacobianVertexP[iDim] > 0.0);
      assert(jacobianVertexN[iDim] > 0.0);
      const PylithScalar S = (1.0/jacobianVertexP[iDim] + 1.0/jacobianVertexN[iDim]) *
	areaVertex * areaVertex;
      assert(S > 0.0);
      lagrangeTIncrVertex[iDim] = 1.0/S * 
	(-residualVertexL[iDim] +
	 areaVertex * (dispIncrVertexP[iDim] - dispIncrVertexN[iDim]));

      assert(jacobianVertexN[iDim] > 0.0);
      dispIncrVertexN[iDim] = 
	+areaVertex / jacobianVertexN[iDim]*lagrangeTIncrVertex[iDim];

      assert(jacobianVertexP[iDim] > 0.0);
      dispIncrVertexP[iDim] = 
	-areaVertex / jacobianVertexP[iDim]*lagrangeTIncrVertex[iDim];

    } // for

    // Compute slip, slip rate, and Lagrange multiplier at time t+dt
    // in fault coordinate system.
    slipVertex = 0.0;
    slipRateVertex = 0.0;
    tractionTpdtVertex = 0.0;
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      for (int jDim=0; jDim < spaceDim; ++jDim) {
	slipVertex[iDim] += orientationVertex[iDim*spaceDim+jDim] *
	  dispRelVertex[jDim];
	slipRateVertex[iDim] += orientationVertex[iDim*spaceDim+jDim] *
	  velRelVertex[jDim];
	tractionTpdtVertex[iDim] += orientationVertex[iDim*spaceDim+jDim] *
	  (lagrangeTVertex[jDim] + lagrangeTIncrVertex[jDim]);
      } // for
    } // for
    
    // Get friction properties and state variables.
    _friction->retrievePropsStateVars(v_fault);

    // Use fault constitutive model to compute traction associated with
    // friction.
    dTractionTpdtVertex = 0.0;
    const bool iterating = false; // No iteration for friction in lumped soln
    CALL_MEMBER_FN(*this,
		   constrainSolnSpaceFn)(&dTractionTpdtVertex,
					 t, slipVertex, slipRateVertex,
					 tractionTpdtVertex,
					 iterating);

    // Rotate traction back to global coordinate system.
    dLagrangeTpdtVertex = 0.0;
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      for (int jDim=0; jDim < spaceDim; ++jDim) {
	dLagrangeTpdtVertex[iDim] += 
	  orientationVertex[jDim*spaceDim+iDim] * dTractionTpdtVertex[jDim];
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
      std::cout << "  " << orientationVertex[iDim];
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
      assert(jacobianVertexP[iDim] > 0.0);
      assert(jacobianVertexN[iDim] > 0.0);

      dispIncrVertexN[iDim] += 
	areaVertex * dLagrangeTpdtVertex[iDim] / jacobianVertexN[iDim];
      dispIncrVertexP[iDim] -= 
	areaVertex * dLagrangeTpdtVertex[iDim] / jacobianVertexP[iDim];

      // Update increment in Lagrange multiplier.
      lagrangeTIncrVertex[iDim] += dLagrangeTpdtVertex[iDim];
    } // for

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(computeEvent);
    _logger->eventBegin(updateEvent);
#endif

    // Compute contribution to adjusting solution only if Lagrange
    // constraint is local (the adjustment is assembled across processors).
    if (globalOrder->isLocal(v_lagrange)) {
      // Adjust displacements to account for Lagrange multiplier values
      // (assumed to be zero in preliminary solve).
      assert(dispIncrVertexN.size() == 
	     dispIncrAdjSection->getFiberDimension(v_negative));
      dispIncrAdjSection->updateAddPoint(v_negative, &dispIncrVertexN[0]);
      
      assert(dispIncrVertexP.size() == 
	     dispIncrAdjSection->getFiberDimension(v_positive));
      dispIncrAdjSection->updateAddPoint(v_positive, &dispIncrVertexP[0]);
    } // if

    // The Lagrange multiplier and relative displacement are NOT
    // assembled across processors, so update even if Lagrange vertex
    // is not local.

    // Set Lagrange multiplier value. Value from preliminary solve is
    // bogus due to artificial diagonal entry in Jacobian of 1.0.
    assert(lagrangeTIncrVertex.size() == 
	   dispIncrSection->getFiberDimension(v_lagrange));
    dispIncrSection->updatePoint(v_lagrange, &lagrangeTIncrVertex[0]);

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
} // adjustSolnLumped

// ----------------------------------------------------------------------
// Get vertex field associated with integrator.
const pylith::topology::Field<pylith::topology::SubMesh>&
pylith::faults::FaultCohesiveDynKin::vertexField(const char* name,
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

  const int slipStrLen = strlen("final_slip");
  const int timeStrLen = strlen("slip_time");

  PylithScalar scale = 0.0;
  int fiberDim = 0;
  if (0 == strcasecmp("slip", name)) {
    const topology::Field<topology::SubMesh>& dispRel = 
      _fields->get("relative disp");
    _allocateBufferVectorField();
    topology::Field<topology::SubMesh>& buffer =
        _fields->get("buffer (vector)");
    buffer.copy(dispRel);
    buffer.label("slip");
    FaultCohesiveLagrange::globalToFault(&buffer, orientation);
    return buffer;

  } else if (0 == strcasecmp("slip_rate", name)) {
    const topology::Field<topology::SubMesh>& velRel = 
      _fields->get("relative velocity");
    _allocateBufferVectorField();
    topology::Field<topology::SubMesh>& buffer =
        _fields->get("buffer (vector)");
    buffer.copy(velRel);
    buffer.label("slip_rate");
    FaultCohesiveLagrange::globalToFault(&buffer, orientation);
    return buffer;

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
    buffer.label("strike_dir");
    buffer.scale(1.0);
    buffer.copy(dirSection);
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
    buffer.label("dip_dir");
    buffer.scale(1.0);
    buffer.copy(dirSection);
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
    buffer.label("normal_dir");
    buffer.scale(1.0);
    buffer.copy(dirSection);
    return buffer;

  } else if (0 == strcasecmp("traction", name)) {
    assert(fields);
    const topology::Field<topology::Mesh>& dispT = fields->get("disp(t)");
    _allocateBufferVectorField();
    topology::Field<topology::SubMesh>& buffer =
        _fields->get("buffer (vector)");
    _calcTractions(&buffer, dispT);
    return buffer;

  } else if (_friction->hasPropStateVar(name)) {
    return _friction->getField(name);

  } else if (_tractPerturbation && _tractPerturbation->hasParameter(name)) {
    const topology::Field<topology::SubMesh>& param = _tractPerturbation->vertexField(name, fields);
    if (param.vectorFieldType() == topology::FieldBase::VECTOR) {
      _allocateBufferVectorField();
      topology::Field<topology::SubMesh>& buffer =
        _fields->get("buffer (vector)");
      buffer.copy(param);
      FaultCohesiveLagrange::globalToFault(&buffer, orientation);
      return buffer;
    } else {
      return param;
    } // if/else

  } else if (0 == strncasecmp("final_slip_X", name, slipStrLen)) {
    const std::string value = std::string(name).substr(slipStrLen + 1);

    const srcs_type::const_iterator s_iter = _eqSrcs.find(value);
    assert(s_iter != _eqSrcs.end());

    // Need to append name of rupture to final slip label. Because
    // Field is const, we use a buffer.
    _allocateBufferVectorField();
    topology::Field<topology::SubMesh>& buffer =
        _fields->get("buffer (vector)");
    buffer.copy(s_iter->second->finalSlip());
    assert(value.length() > 0);
    const std::string& label = (_eqSrcs.size() > 1) ?
      std::string("final_slip_") + std::string(value) : "final_slip";
    buffer.label(label.c_str());

    return buffer;

  } else if (0 == strncasecmp("slip_time_X", name, timeStrLen)) {
    const std::string value = std::string(name).substr(timeStrLen + 1);
    const srcs_type::const_iterator s_iter = _eqSrcs.find(value);
    assert(s_iter != _eqSrcs.end());

    // Need to append name of rupture to final slip label. Because
    // Field is const, we use a buffer.
    _allocateBufferScalarField();
    topology::Field<topology::SubMesh>& buffer =
        _fields->get("buffer (scalar)");
    buffer.copy(s_iter->second->slipTime());
    assert(value.length() > 0);
    const std::string& label = (_eqSrcs.size() > 1) ?
      std::string("slip_time_") + std::string(value) : "slip_time";
    buffer.label(label.c_str());

    return buffer;

  } else {
    std::ostringstream msg;
    msg << "Request for unknown vertex field '" << name << "' for fault '"
        << label() << "'.";
    throw std::runtime_error(msg.str());
  } // else

  // Should never get here.
  throw std::logic_error("Unknown field in FaultCohesiveDynKin::vertexField().");

  // Satisfy return values
  assert(_fields);
  const topology::Field<topology::SubMesh>& buffer = _fields->get(
    "buffer (vector)");

  return buffer;
} // vertexField

// ----------------------------------------------------------------------
// Compute tractions on fault surface using solution.
void
pylith::faults::FaultCohesiveDynKin::_calcTractions(
    topology::Field<topology::SubMesh>* tractions,
    const topology::Field<topology::Mesh>& dispT)
{ // _calcTractions
  assert(tractions);
  assert(_faultMesh);
  assert(_fields);
  assert(_normalizer);

  // Fiber dimension of tractions matches spatial dimension.
  const int spaceDim = _quadrature->spaceDim();
  scalar_array tractionsVertex(spaceDim);

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
    logger.stagePush("FaultFields");

    const topology::Field<topology::SubMesh>& dispRel = 
      _fields->get("relative disp");
    tractions->cloneSection(dispRel);

    logger.stagePop();
  } // if
  const PylithScalar pressureScale = _normalizer->pressureScale();
  tractions->label("traction");
  tractions->scale(pressureScale);
  tractions->zero();

  const int numVertices = _cohesiveVertices.size();
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_lagrange = _cohesiveVertices[iVertex].lagrange;
    const int v_fault = _cohesiveVertices[iVertex].fault;

    assert(spaceDim == dispTSection->getFiberDimension(v_lagrange));
    const PylithScalar* dispTVertex = dispTSection->restrictPoint(v_lagrange);
    assert(dispTVertex);

    assert(spaceDim*spaceDim == 
	   orientationSection->getFiberDimension(v_fault));
    const PylithScalar* orientationVertex = 
      orientationSection->restrictPoint(v_fault);
    assert(orientationVertex);

    // Rotate tractions to fault coordinate system.
    tractionsVertex = 0.0;
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      for (int jDim=0; jDim < spaceDim; ++jDim) {
	tractionsVertex[iDim] += orientationVertex[iDim*spaceDim+jDim] *
	  dispTVertex[jDim];
      } // for
    } // for

    assert(tractionsVertex.size() == 
	   tractionsSection->getFiberDimension(v_fault));
    tractionsSection->updatePoint(v_fault, &tractionsVertex[0]);
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
pylith::faults::FaultCohesiveDynKin::_updateRelMotion(const topology::SolutionFields& fields)
{ // _updateRelMotion
  assert(_fields);

  const int spaceDim = _quadrature->spaceDim();

  // Get section information
  const ALE::Obj<RealSection>& dispTSection =
    fields.get("disp(t)").section();
  assert(!dispTSection.isNull());

  const ALE::Obj<RealSection>& dispIncrSection =
    fields.get("dispIncr(t->t+dt)").section();
  assert(!dispIncrSection.isNull());

  scalar_array dispRelVertex(spaceDim);
  const ALE::Obj<RealSection>& dispRelSection =
    _fields->get("relative disp").section();
  assert(!dispRelSection.isNull());

  const ALE::Obj<RealSection>& velocitySection =
      fields.get("velocity(t)").section();
  assert(!velocitySection.isNull());

  scalar_array velRelVertex(spaceDim);
  const ALE::Obj<RealSection>& velRelSection =
      _fields->get("relative velocity").section();
  assert(!velRelSection.isNull());

  const int numVertices = _cohesiveVertices.size();
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_fault = _cohesiveVertices[iVertex].fault;
    const int v_negative = _cohesiveVertices[iVertex].negative;
    const int v_positive = _cohesiveVertices[iVertex].positive;

    // Get displacement values
    assert(spaceDim == dispTSection->getFiberDimension(v_negative));
    const PylithScalar* dispTVertexN = dispTSection->restrictPoint(v_negative);
    assert(dispTVertexN);

    assert(spaceDim == dispTSection->getFiberDimension(v_positive));
    const PylithScalar* dispTVertexP = dispTSection->restrictPoint(v_positive);
    assert(dispTVertexP);

    assert(spaceDim == dispIncrSection->getFiberDimension(v_negative));
    const PylithScalar* dispIncrVertexN = 
      dispIncrSection->restrictPoint(v_negative);
    assert(dispIncrVertexN);

    assert(spaceDim == dispIncrSection->getFiberDimension(v_positive));
    const PylithScalar* dispIncrVertexP = 
      dispIncrSection->restrictPoint(v_positive);
    assert(dispIncrVertexP);

    // Compute relative displacememt
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      const PylithScalar value = 
	dispTVertexP[iDim] + dispIncrVertexP[iDim] 
	- dispTVertexN[iDim] -  dispIncrVertexN[iDim];
      dispRelVertex[iDim] = fabs(value) > _zeroTolerance ? value : 0.0;
    } // for

    // Update relative displacement field.
    assert(dispRelVertex.size() == 
	   dispRelSection->getFiberDimension(v_fault));
    dispRelSection->updatePoint(v_fault, &dispRelVertex[0]);

    // Get velocity values
    assert(spaceDim == velocitySection->getFiberDimension(v_negative));
    const PylithScalar* velocityVertexN = velocitySection->restrictPoint(v_negative);
    assert(velocityVertexN);

    assert(spaceDim == velocitySection->getFiberDimension(v_positive));
    const PylithScalar* velocityVertexP = velocitySection->restrictPoint(v_positive);
    assert(velocityVertexP);

    // Compute relative velocity
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      const PylithScalar value = velocityVertexP[iDim] - velocityVertexN[iDim];
      velRelVertex[iDim] = fabs(value) > _zeroTolerance ? value : 0.0;
    } // for

    // Update relative velocity field.
    assert(velRelVertex.size() == 
	   velRelSection->getFiberDimension(v_fault));
    velRelSection->updatePoint(v_fault, &velRelVertex[0]);
  } // for

  PetscLogFlops(numVertices*spaceDim*spaceDim*4);
} // _updateRelMotion

// ----------------------------------------------------------------------
// Setup sensitivity problem to compute change in slip given change in Lagrange multipliers.
void
pylith::faults::FaultCohesiveDynKin::_sensitivitySetup(const topology::Jacobian& jacobian)
{ // _sensitivitySetup
  assert(_fields);
  assert(_quadrature);

  const int spaceDim = _quadrature->spaceDim();

  // Setup fields involved in sensitivity solve.
  if (!_fields->hasField("sensitivity solution")) {
    _fields->add("sensitivity solution", "sensitivity_soln");
    topology::Field<topology::SubMesh>& solution =
        _fields->get("sensitivity solution");
    const topology::Field<topology::SubMesh>& dispRel =
        _fields->get("relative disp");
    solution.cloneSection(dispRel);
    solution.createScatter(solution.mesh());
  } // if
  const topology::Field<topology::SubMesh>& solution =
      _fields->get("sensitivity solution");

  if (!_fields->hasField("sensitivity residual")) {
    _fields->add("sensitivity residual", "sensitivity_residual");
    topology::Field<topology::SubMesh>& residual =
        _fields->get("sensitivity residual");
    residual.cloneSection(solution);
    residual.createScatter(solution.mesh());
  } // if

  if (!_fields->hasField("sensitivity relative disp")) {
    _fields->add("sensitivity relative disp", "sensitivity_relative_disp");
    topology::Field<topology::SubMesh>& dispRel =
        _fields->get("sensitivity relative disp");
    dispRel.cloneSection(solution);
  } // if
  topology::Field<topology::SubMesh>& dispRel =
    _fields->get("sensitivity relative disp");
  dispRel.zero();

  if (!_fields->hasField("sensitivity dLagrange")) {
    _fields->add("sensitivity dLagrange", "sensitivity_dlagrange");
    topology::Field<topology::SubMesh>& dLagrange =
        _fields->get("sensitivity dLagrange");
    dLagrange.cloneSection(solution);
  } // if
  topology::Field<topology::SubMesh>& dLagrange =
    _fields->get("sensitivity dLagrange");
  dLagrange.zero();

  // Setup Jacobian sparse matrix for sensitivity solve.
  if (0 == _jacobian)
    _jacobian = new topology::Jacobian(solution, jacobian.matrixType());
  assert(_jacobian);
  _jacobian->zero();

  // Setup PETSc KSP linear solver.
  if (0 == _ksp) {
    PetscErrorCode err = 0;
    err = KSPCreate(_faultMesh->comm(), &_ksp); CHECK_PETSC_ERROR(err);
    err = KSPSetInitialGuessNonzero(_ksp, PETSC_FALSE); CHECK_PETSC_ERROR(err);
    PylithScalar rtol = 0.0;
    PylithScalar atol = 0.0;
    PylithScalar dtol = 0.0;
    int maxIters = 0;
    err = KSPGetTolerances(_ksp, &rtol, &atol, &dtol, &maxIters); 
    CHECK_PETSC_ERROR(err);
    rtol = 1.0e-3*_zeroTolerance;
    atol = 1.0e-5*_zeroTolerance;
    err = KSPSetTolerances(_ksp, rtol, atol, dtol, maxIters);
    CHECK_PETSC_ERROR(err);

    PC pc;
    err = KSPGetPC(_ksp, &pc); CHECK_PETSC_ERROR(err);
    err = PCSetType(pc, PCJACOBI); CHECK_PETSC_ERROR(err);
    err = KSPSetType(_ksp, KSPGMRES); CHECK_PETSC_ERROR(err);

    err = KSPAppendOptionsPrefix(_ksp, "friction_");
    err = KSPSetFromOptions(_ksp); CHECK_PETSC_ERROR(err);
  } // if
} // _sensitivitySetup

// ----------------------------------------------------------------------
// Update the Jacobian values for the sensitivity solve.
void
pylith::faults::FaultCohesiveDynKin::_sensitivityUpdateJacobian(const bool negativeSide,
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
  const ALE::Obj<RealSection>& solutionDomainSection = solutionDomain.section();
  assert(!solutionDomainSection.isNull());

  // Get cohesive cells
  const ALE::Obj<SieveMesh>& sieveMesh = fields.mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& cellsCohesive = sieveMesh->getLabelStratum("material-id", id());
  assert(!cellsCohesive.isNull());
  const SieveMesh::label_sequence::iterator cellsCohesiveBegin = cellsCohesive->begin();
  const SieveMesh::label_sequence::iterator cellsCohesiveEnd = cellsCohesive->end();

  // Visitor for Jacobian matrix associated with domain.
  scalar_array jacobianSubCell(submatrixSize);
  const PetscMat jacobianDomainMatrix = jacobian.matrix();
  assert(jacobianDomainMatrix);
  const ALE::Obj<SieveMesh::order_type>& globalOrderDomain =
    sieveMesh->getFactory()->getGlobalOrder(sieveMesh, "default", solutionDomainSection);
  assert(!globalOrderDomain.isNull());
  const ALE::Obj<SieveMesh::sieve_type>& sieve = sieveMesh->getSieve();
  assert(!sieve.isNull());
  const int closureSize = int(pow(sieve->getMaxConeSize(), sieveMesh->depth()));
  assert(closureSize >= 0);
  ALE::ISieveVisitor::NConeRetriever<SieveMesh::sieve_type> ncV(*sieve, closureSize);

  // Get fault Sieve mesh
  const ALE::Obj<SieveSubMesh>& faultSieveMesh = _faultMesh->sieveMesh();
  assert(!faultSieveMesh.isNull());

  // Get sensitivity solution field
  const ALE::Obj<RealSection>& solutionFaultSection =
    _fields->get("sensitivity solution").section();
  assert(!solutionFaultSection.isNull());

  // Visitor for Jacobian matrix associated with fault.
  assert(_jacobian);
  const PetscMat jacobianFaultMatrix = _jacobian->matrix();
  assert(jacobianFaultMatrix);
  const ALE::Obj<SieveSubMesh::order_type>& globalOrderFault =
    faultSieveMesh->getFactory()->getGlobalOrder(faultSieveMesh, "default", solutionFaultSection);
  assert(!globalOrderFault.isNull());
  // We would need to request unique points here if we had an interpolated mesh
  IndicesVisitor jacobianFaultVisitor(*solutionFaultSection,
				      *globalOrderFault, closureSize*spaceDim);

  const int iCone = (negativeSide) ? 0 : 1;

  const int numCohesiveCells = cellsCohesive->size();
  IS* cellsIS = (numCohesiveCells > 0) ? new IS[numCohesiveCells] : 0;
  int_array indicesGlobal(subnrows);
  int_array indicesLocal(numCohesiveCells*subnrows);
  int_array indicesPerm(subnrows);

  PetscErrorCode err = 0;
  int iCohesiveCell = 0;
  for (SieveMesh::label_sequence::iterator c_iter=cellsCohesiveBegin;
       c_iter != cellsCohesiveEnd;
       ++c_iter, ++iCohesiveCell) {
    // Get cone for cohesive cell
    ncV.clear();
    ALE::ISieveTraversal<SieveMesh::sieve_type>::orientedClosure(*sieve, *c_iter, ncV);
    const int coneSize = ncV.getSize();
    assert(coneSize == 3*numBasis);
    const SieveMesh::point_type* cohesiveCone = ncV.getPoints();
    assert(cohesiveCone);

    // Get indices
    for (int iBasis = 0; iBasis < numBasis; ++iBasis) {
      // negative side of the fault: iCone=0
      // positive side of the fault: iCone=1
      const int v_domain = cohesiveCone[iCone*numBasis+iBasis];

      for (int iDim=0, iB=iBasis*spaceDim; iDim < spaceDim; ++iDim) {
	indicesGlobal[iB+iDim] = globalOrderDomain->getIndex(v_domain) + iDim;
      } // for
    } // for

    for (int i=0; i < subnrows; ++i) {
      indicesPerm[i]  = i;
    } // for
    err = PetscSortIntWithArray(indicesGlobal.size(), &indicesGlobal[0], &indicesPerm[0]);CHECK_PETSC_ERROR(err);

    for (int i=0; i < subnrows; ++i) {
      indicesLocal[iCohesiveCell*subnrows+indicesPerm[i]] = i;
    } // for
    cellsIS[iCohesiveCell] = PETSC_NULL;
    err = ISCreateGeneral(PETSC_COMM_SELF, indicesGlobal.size(), &indicesGlobal[0], PETSC_COPY_VALUES, &cellsIS[iCohesiveCell]);CHECK_PETSC_ERROR(err);

  } // for

  PetscMat* submatrices = 0;
  err = MatGetSubMatrices(jacobianDomainMatrix, numCohesiveCells, cellsIS, cellsIS, MAT_INITIAL_MATRIX, &submatrices);CHECK_PETSC_ERROR(err);

  iCohesiveCell = 0;
  for (SieveMesh::label_sequence::iterator c_iter=cellsCohesiveBegin;
       c_iter != cellsCohesiveEnd;
       ++c_iter, ++iCohesiveCell) {
    // Get values for submatrix associated with cohesive cell
    jacobianSubCell = 0.0;
    err = MatGetValues(submatrices[iCohesiveCell], 
		       subnrows, &indicesLocal[iCohesiveCell*subnrows],
		       subnrows, &indicesLocal[iCohesiveCell*subnrows],
		       &jacobianSubCell[0]);CHECK_PETSC_ERROR_MSG(err, "Restrict from PETSc Mat failed.");

    // Insert cell contribution into PETSc Matrix
    const SieveMesh::point_type c_fault = _cohesiveToFault[*c_iter];
    jacobianFaultVisitor.clear();
    err = updateOperator(jacobianFaultMatrix, *faultSieveMesh->getSieve(),
			 jacobianFaultVisitor, c_fault,
			 &jacobianSubCell[0], INSERT_VALUES);
    CHECK_PETSC_ERROR_MSG(err, "Update to PETSc Mat failed.");

    // Destory IS for cohesiveCell
    err = ISDestroy(&cellsIS[iCohesiveCell]);CHECK_PETSC_ERROR(err);
  } // for

  err = MatDestroyMatrices(numCohesiveCells, &submatrices);CHECK_PETSC_ERROR(err);
  delete[] cellsIS; cellsIS = 0;

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
pylith::faults::FaultCohesiveDynKin::_sensitivityReformResidual(const bool negativeSide)
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
  const ALE::Obj<SieveMesh>& faultSieveMesh = _faultMesh->sieveMesh();
  assert(!faultSieveMesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& cells =
    faultSieveMesh->heightStratum(0);
  assert(!cells.isNull());
  const SieveSubMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveSubMesh::label_sequence::iterator cellsEnd = cells->end();
  const int numCells = cells->size();

  // Get sections
  scalar_array coordinatesCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& coordinates = 
    faultSieveMesh->getRealSection("coordinates");
  assert(!coordinates.isNull());
  RestrictVisitor coordsVisitor(*coordinates, 
				coordinatesCell.size(), &coordinatesCell[0]);

  scalar_array dLagrangeCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& dLagrangeSection = 
    _fields->get("sensitivity dLagrange").section();
  assert(!dLagrangeSection.isNull());
  RestrictVisitor dLagrangeVisitor(*dLagrangeSection, 
				   dLagrangeCell.size(), &dLagrangeCell[0]);

  scalar_array residualCell(numBasis*spaceDim);
  topology::Field<topology::SubMesh>& residual =
      _fields->get("sensitivity residual");
  const ALE::Obj<RealSection>& residualSection = residual.section();
  UpdateAddVisitor residualVisitor(*residualSection, &residualCell[0]);

  residual.zero();

  // Loop over cells
  for (SieveSubMesh::label_sequence::iterator c_iter=cellsBegin;
       c_iter != cellsEnd;
       ++c_iter) {
    // Compute geometry
    coordsVisitor.clear();
    faultSieveMesh->restrictClosure(*c_iter, coordsVisitor);
    _quadrature->computeGeometry(coordinatesCell, *c_iter);

    // Restrict input fields to cell
    dLagrangeVisitor.clear();
    faultSieveMesh->restrictClosure(*c_iter, dLagrangeVisitor);

    // Get cell geometry information that depends on cell
    const scalar_array& basis = _quadrature->basis();
    const scalar_array& jacobianDet = _quadrature->jacobianDet();

    // Compute product of basis functions.
    // Want values summed over quadrature points
    basisProducts = 0.0;
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      const PylithScalar wt = quadWts[iQuad] * jacobianDet[iQuad];

      for (int iBasis=0, iQ=iQuad*numBasis; iBasis < numBasis; ++iBasis) {
        const PylithScalar valI = wt*basis[iQ+iBasis];
	
        for (int jBasis=0; jBasis < numBasis; ++jBasis) {
	  
	  basisProducts[iBasis*numBasis+jBasis] += valI*basis[iQ+jBasis];
	} // for
      } // for
    } // for

    residualCell = 0.0;
    
    for (int iBasis=0; iBasis < numBasis; ++iBasis) {
      for (int jBasis=0; jBasis < numBasis; ++jBasis) {
	const PylithScalar l = signFault * basisProducts[iBasis*numBasis+jBasis];
	for (int iDim=0; iDim < spaceDim; ++iDim) {
	  residualCell[iBasis*spaceDim+iDim] += 
	    l * dLagrangeCell[jBasis*spaceDim+iDim];
	} // for
      } // for
    } // for

    // Assemble cell contribution into field
    residualVisitor.clear();
    faultSieveMesh->updateClosure(*c_iter, residualVisitor);    
  } // for
} // _sensitivityReformResidual

// ----------------------------------------------------------------------
// Solve sensitivity problem.
void
pylith::faults::FaultCohesiveDynKin::_sensitivitySolve(void)
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

  const PetscVec residualVec = residual.vector();
  const PetscVec solutionVec = solution.vector();
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
pylith::faults::FaultCohesiveDynKin::_sensitivityUpdateSoln(const bool negativeSide)
{ // _sensitivityUpdateSoln
  assert(_fields);
  assert(_quadrature);

  const int spaceDim = _quadrature->spaceDim();

  scalar_array dispVertex(spaceDim);
  const ALE::Obj<RealSection>& solutionSection =
      _fields->get("sensitivity solution").section();
  assert(!solutionSection.isNull());
  const ALE::Obj<RealSection>& dispRelSection =
    _fields->get("sensitivity relative disp").section();
  assert(!dispRelSection.isNull());
  const ALE::Obj<RealSection>& dLagrangeTpdtSection =
      _fields->get("sensitivity dLagrange").section();
  assert(!dLagrangeTpdtSection.isNull());

  const PylithScalar sign = (negativeSide) ? -1.0 : 1.0;

  const int numVertices = _cohesiveVertices.size();
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_fault = _cohesiveVertices[iVertex].fault;

    solutionSection->restrictPoint(v_fault, &dispVertex[0], dispVertex.size());

    // If no change in the Lagrange multiplier computed from friction
    // criterion, there are no updates, so continue.
    assert(spaceDim == dLagrangeTpdtSection->getFiberDimension(v_fault));
    const PylithScalar* dLagrangeTpdtVertex = dLagrangeTpdtSection->restrictPoint(v_fault);
    assert(dLagrangeTpdtVertex);
    PylithScalar dLagrangeTpdtVertexMag = 0.0;
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      dLagrangeTpdtVertexMag += dLagrangeTpdtVertex[iDim]*dLagrangeTpdtVertex[iDim];
    } // for
    if (0.0 == dLagrangeTpdtVertexMag)
      continue;

    // Update relative displacements associated with sensitivity solve
    // solution
    dispVertex *= sign;

    assert(dispVertex.size() == dispRelSection->getFiberDimension(v_fault));
    dispRelSection->updateAddPoint(v_fault, &dispVertex[0]);
  } // for
} // _sensitivityUpdateSoln


// ----------------------------------------------------------------------
// Compute norm of residual associated with matching fault
// constitutive model using update from sensitivity solve. We use
// this in a line search to find a good update (required because
// fault constitutive model may have a complex nonlinear feedback
// with deformation).
PylithScalar
pylith::faults::FaultCohesiveDynKin::_constrainSolnSpaceNorm(const PylithScalar alpha,
							  const PylithScalar t,
							  topology::SolutionFields* const fields)
{ // _constrainSolnSpaceNorm
  /// Member prototype for _constrainSolnSpaceXD()
  typedef void (pylith::faults::FaultCohesiveDynKin::*constrainSolnSpace_fn_type)
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

  constrainSolnSpace_fn_type constrainSolnSpaceFn;
  switch (spaceDim) { // switch
  case 1:
    constrainSolnSpaceFn = 
      &pylith::faults::FaultCohesiveDynKin::_constrainSolnSpace1D;
    break;
  case 2: 
    constrainSolnSpaceFn = 
      &pylith::faults::FaultCohesiveDynKin::_constrainSolnSpace2D;
    break;
  case 3:
    constrainSolnSpaceFn = 
      &pylith::faults::FaultCohesiveDynKin::_constrainSolnSpace3D;
    break;
  default :
    assert(0);
    throw std::logic_error("Unknown spatial dimension in "
			   "FaultCohesiveDynKin::constrainSolnSpace().");
  } // switch

  // Get sections
  scalar_array slipTpdtVertex(spaceDim); // fault coordinates
  scalar_array slipRateVertex(spaceDim); // fault coordinates
  scalar_array tractionTpdtVertex(spaceDim); // fault coordinates
  scalar_array tractionMisfitVertex(spaceDim); // fault coordinates

  const ALE::Obj<RealSection>& orientationSection = _fields->get("orientation").section();
  assert(!orientationSection.isNull());

  const ALE::Obj<RealSection>& dLagrangeTpdtSection = _fields->get("sensitivity dLagrange").section();
  assert(!dLagrangeTpdtSection.isNull());

  const ALE::Obj<RealSection>& sensDispRelSection = _fields->get("sensitivity relative disp").section();
  assert(!sensDispRelSection.isNull());

  const ALE::Obj<RealSection>& dispTSection = fields->get("disp(t)").section();
  assert(!dispTSection.isNull());

  const ALE::Obj<RealSection>& dispIncrSection = fields->get("dispIncr(t->t+dt)").section();
  assert(!dispIncrSection.isNull());


  // Get fault information
  const ALE::Obj<SieveMesh>& sieveMesh = fields->mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::order_type>& globalOrder =
      sieveMesh->getFactory()->getGlobalOrder(sieveMesh, "default", dispIncrSection);
  assert(!globalOrder.isNull());

  bool isOpening = false;
  PylithScalar norm2 = 0.0;
  int numVertices = _cohesiveVertices.size();
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_lagrange = _cohesiveVertices[iVertex].lagrange;
    const int v_fault = _cohesiveVertices[iVertex].fault;
    const int v_negative = _cohesiveVertices[iVertex].negative;
    const int v_positive = _cohesiveVertices[iVertex].positive;

    // Compute contribution only if Lagrange constraint is local.
    if (!globalOrder->isLocal(v_lagrange))
      continue;

    // Get displacement values
    assert(spaceDim == dispTSection->getFiberDimension(v_negative));
    const PylithScalar* dispTVertexN = dispTSection->restrictPoint(v_negative);
    assert(dispTVertexN);

    assert(spaceDim == dispTSection->getFiberDimension(v_positive));
    const PylithScalar* dispTVertexP = dispTSection->restrictPoint(v_positive);
    assert(dispTVertexP);

    // Get displacement increment values.
    assert(spaceDim == dispIncrSection->getFiberDimension(v_negative));
    const PylithScalar* dispIncrVertexN = dispIncrSection->restrictPoint(v_negative);
    assert(dispIncrVertexN);

    assert(spaceDim == dispIncrSection->getFiberDimension(v_positive));
    const PylithScalar* dispIncrVertexP = dispIncrSection->restrictPoint(v_positive);
    assert(dispIncrVertexP);

    // Get orientation
    assert(spaceDim*spaceDim == orientationSection->getFiberDimension(v_fault));
    const PylithScalar* orientationVertex = orientationSection->restrictPoint(v_fault);

    // Get change in relative displacement from sensitivity solve.
    assert(spaceDim == sensDispRelSection->getFiberDimension(v_fault));
    const PylithScalar* dDispRelVertex = sensDispRelSection->restrictPoint(v_fault);
    assert(dDispRelVertex);

    // Get Lagrange multiplier values from disp(t), and dispIncr(t->t+dt)
    assert(spaceDim == dispTSection->getFiberDimension(v_lagrange));
    const PylithScalar* lagrangeTVertex = dispTSection->restrictPoint(v_lagrange);
    assert(lagrangeTVertex);

    assert(spaceDim == dispIncrSection->getFiberDimension(v_lagrange));
    const PylithScalar* lagrangeTIncrVertex = dispIncrSection->restrictPoint(v_lagrange);
    assert(lagrangeTIncrVertex);

    assert(spaceDim == dLagrangeTpdtSection->getFiberDimension(v_fault));
    const PylithScalar* dLagrangeTpdtVertex = dLagrangeTpdtSection->restrictPoint(v_fault);
    assert(dLagrangeTpdtVertex);

    // Compute slip, slip rate, and traction at time t+dt as part of
    // line search.
    slipTpdtVertex = 0.0;
    slipRateVertex = 0.0;
    tractionTpdtVertex = 0.0;
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      for (int jDim=0; jDim < spaceDim; ++jDim) {
	slipTpdtVertex[iDim] += orientationVertex[iDim*spaceDim+jDim] *
	  (dispTVertexP[jDim] + dispIncrVertexP[jDim]
	   - dispTVertexN[jDim] - dispIncrVertexN[jDim] +
	   alpha*dDispRelVertex[jDim]);
	slipRateVertex[iDim] += orientationVertex[iDim*spaceDim+jDim] *
	  (dispIncrVertexP[jDim] - dispIncrVertexN[jDim] +
	   alpha*dDispRelVertex[jDim]) / dt;
	tractionTpdtVertex[iDim] += orientationVertex[iDim*spaceDim+jDim] *
	  (lagrangeTVertex[jDim] + lagrangeTIncrVertex[jDim] +
	   alpha*dLagrangeTpdtVertex[jDim]);
      } // for
      if (fabs(slipRateVertex[iDim]) < _zeroTolerance) {
	slipRateVertex[iDim] = 0.0;
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
      // Step b: Insure fault traction is zero when opening (if
      // alpha=1 this should be enforced already, but will not be
      // properly enforced when alpha < 1).
      
      for (int iDim=0; iDim < spaceDim; ++iDim) {
	tractionTpdtVertex[iDim] = 0.0;
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
    CALL_MEMBER_FN(*this,
		   constrainSolnSpaceFn)(&tractionMisfitVertex, t,
					 slipTpdtVertex, slipRateVertex,
					 tractionTpdtVertex,
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
      std::cout << " " << dDispRelVertex[iDim];
    } // for
    std::cout << std::endl;
#endif

    for (int iDim=0; iDim < spaceDim; ++iDim) {
      norm2 += tractionMisfitVertex[iDim]*tractionMisfitVertex[iDim];
    } // for
  } // for

  if (isOpening && alpha < 1.0) {
    norm2 = PYLITH_MAXFLOAT;
  } // if

  PylithScalar norm2Total = 0.0;
  int numVerticesTotal = 0;
  if (sizeof(PylithScalar) == 8) {
    MPI_Allreduce(&norm2, &norm2Total, 1, MPI_DOUBLE, MPI_SUM, fields->mesh().comm());
  } else {
    MPI_Allreduce(&norm2, &norm2Total, 1, MPI_FLOAT, MPI_SUM, fields->mesh().comm());
  } // if/else
  MPI_Allreduce(&numVertices, &numVerticesTotal, 1, MPI_INT, MPI_SUM, fields->mesh().comm());

  assert(numVerticesTotal > 0);
  return sqrt(norm2Total) / numVerticesTotal;
} // _constrainSolnSpaceNorm


// ----------------------------------------------------------------------
// Constrain solution space in 1-D.
void
pylith::faults::FaultCohesiveDynKin::_constrainSolnSpace1D(scalar_array* dTractionTpdt,
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
pylith::faults::FaultCohesiveDynKin::_constrainSolnSpace2D(scalar_array* dTractionTpdt,
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
pylith::faults::FaultCohesiveDynKin::_constrainSolnSpace3D(scalar_array* dTractionTpdt,
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
