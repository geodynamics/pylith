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
// Copyright (c) 2010 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "FaultCohesiveDyn.hh" // implementation of object methods

#include "CohesiveTopology.hh" // USES CohesiveTopology

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
typedef pylith::topology::Mesh::RealSection RealSection;
typedef pylith::topology::SubMesh::SieveMesh SieveSubMesh;
typedef pylith::topology::Mesh::RestrictVisitor RestrictVisitor;

// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::FaultCohesiveDyn::FaultCohesiveDyn(void) :
  _dbInitialTract(0),
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

  _dbInitialTract = 0; // :TODO: Use shared pointer
  _friction = 0; // :TODO: Use shared pointer

  delete _jacobian; _jacobian = 0;
  if (0 != _ksp) {
    PetscErrorCode err = KSPDestroy(&_ksp); _ksp = 0;
    CHECK_PETSC_ERROR(err);
  } // if
} // deallocate

// ----------------------------------------------------------------------
// Sets the spatial database for the inital tractions
void
pylith::faults::FaultCohesiveDyn::dbInitialTract(spatialdata::spatialdb::SpatialDB* db)
{ // dbInitial
  _dbInitialTract = db;
} // dbInitial

// ----------------------------------------------------------------------
// Get the friction (constitutive) model.  
void
pylith::faults::FaultCohesiveDyn::frictionModel(friction::FrictionModel* const model)
{ // frictionModel
  _friction = model;
} // frictionModel

// ----------------------------------------------------------------------
// Initialize fault. Determine orientation and setup boundary
void
pylith::faults::FaultCohesiveDyn::initialize(const topology::Mesh& mesh,
					      const double upDir[3])
{ // initialize
  assert(0 != upDir);
  assert(0 != _quadrature);
  assert(0 != _normalizer);

  FaultCohesiveLagrange::initialize(mesh, upDir);

  // Get initial tractions using a spatial database.
  _setupInitialTractions();

  // Setup fault constitutive model.
  assert(0 != _friction);
  assert(0 != _faultMesh);
  assert(0 != _fields);
  _friction->normalizer(*_normalizer);
  _friction->initialize(*_faultMesh, _quadrature, _fields->get("area"));

  const spatialdata::geocoords::CoordSys* cs = mesh.coordsys();
  assert(0 != cs);


  // Create field for slip rate associated with Lagrange vertex k
  _fields->add("slip rate", "slip_rate");
  topology::Field<topology::SubMesh>& slipRate = _fields->get("slip rate");
  topology::Field<topology::SubMesh>& slip = _fields->get("slip");
  slipRate.cloneSection(slip);
  slipRate.vectorFieldType(topology::FieldBase::VECTOR);

  //logger.stagePop();
} // initialize

// ----------------------------------------------------------------------
// Integrate contributions to residual term (r) for operator.
void
pylith::faults::FaultCohesiveDyn::integrateResidual(
			     const topology::Field<topology::Mesh>& residual,
			     const double t,
			     topology::SolutionFields* const fields)
{ // integrateResidual
  assert(0 != fields);
  assert(0 != _fields);

  // Initial fault tractions have been assembled, so they do not need
  // assembling across processors.

  FaultCohesiveLagrange::integrateResidual(residual, t, fields);

  // No contribution if no initial tractions are specified.
  if (0 == _dbInitialTract)
    return;

  const int spaceDim = _quadrature->spaceDim();

  // Get sections
  double_array forcesInitialVertex(spaceDim);
  const ALE::Obj<RealSection>& forcesInitialSection = 
    _fields->get("initial forces").section();
  assert(!forcesInitialSection.isNull());

  double_array residualVertex(spaceDim);
  const ALE::Obj<RealSection>& residualSection = residual.section();
  assert(!residualSection.isNull());

  const ALE::Obj<RealSection>& slipSection = _fields->get("slip").section();
  assert(!slipSection.isNull());

  const ALE::Obj<SieveMesh>& sieveMesh = fields->mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::order_type>& globalOrder =
    sieveMesh->getFactory()->getGlobalOrder(sieveMesh, "default", 
					    residualSection);
  assert(!globalOrder.isNull());

  const int numVertices = _cohesiveVertices.size();
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_fault = _cohesiveVertices[iVertex].fault;
    const int v_lagrange = _cohesiveVertices[iVertex].lagrange;
    const int v_negative = _cohesiveVertices[iVertex].negative;
    const int v_positive = _cohesiveVertices[iVertex].positive;

    // Compute contribution only if Lagrange constraint is local.
    if (!globalOrder->isLocal(v_lagrange))
      continue;

    // Get initial forces at fault vertex. Forces are in the global
    // coordinate system so no rotation is necessary.
    forcesInitialSection->restrictPoint(v_fault, 
					&forcesInitialVertex[0],
					forcesInitialVertex.size());

    assert(spaceDim == slipSection->getFiberDimension(v_fault));
    const double* slipVertex = slipSection->restrictPoint(v_fault);
    assert(0 != slipVertex);
    
    // only apply initial tractions if there is no opening
    if (0.0 == slipVertex[spaceDim-1]) {
      residualVertex = forcesInitialVertex;

      assert(residualVertex.size() == 
	     residualSection->getFiberDimension(v_positive));
      residualSection->updateAddPoint(v_positive, &residualVertex[0]);
      
      residualVertex *= -1.0;
      assert(residualVertex.size() == 
	     residualSection->getFiberDimension(v_negative));
      residualSection->updateAddPoint(v_negative, &residualVertex[0]);
    } // if
  } // for

  PetscLogFlops(numVertices*spaceDim);
} // integrateResidual

// ----------------------------------------------------------------------
// Update state variables as needed.
void
pylith::faults::FaultCohesiveDyn::updateStateVars(
				      const double t,
				      topology::SolutionFields* const fields)
{ // updateStateVars
  assert(0 != fields);
  assert(0 != _fields);

  _updateSlipRate(*fields);

  const int spaceDim = _quadrature->spaceDim();

  // Allocate arrays for vertex values
  double_array tractionTVertex(spaceDim);
  double_array tractionTpdtVertex(spaceDim);
  double_array slipTpdtVertex(spaceDim);
  double_array lagrangeTpdtVertex(spaceDim);

  // Get sections
  double_array slipVertex(spaceDim);
  const ALE::Obj<RealSection>& slipSection = _fields->get("slip").section();
  assert(!slipSection.isNull());

  double_array slipRateVertex(spaceDim);
  const ALE::Obj<RealSection>& slipRateSection =
      _fields->get("slip rate").section();
  assert(!slipRateSection.isNull());

  const ALE::Obj<RealSection>& areaSection = _fields->get("area").section();
  assert(!areaSection.isNull());

  double_array lagrangeTVertex(spaceDim);
  const ALE::Obj<RealSection>& dispTSection = fields->get("disp(t)").section();
  assert(!dispTSection.isNull());

  double_array lagrangeTIncrVertex(spaceDim);
  const ALE::Obj<RealSection>& dispTIncrSection =
      fields->get("dispIncr(t->t+dt)").section();
  assert(!dispTIncrSection.isNull());

  const int numVertices = _cohesiveVertices.size();
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_lagrange = _cohesiveVertices[iVertex].lagrange;
    const int v_fault = _cohesiveVertices[iVertex].fault;
    const int v_negative = _cohesiveVertices[iVertex].negative;
    const int v_positive = _cohesiveVertices[iVertex].positive;

    // Get slip
    slipSection->restrictPoint(v_fault, &slipVertex[0], slipVertex.size());

    // Get slip rate
    slipRateSection->restrictPoint(v_fault, &slipRateVertex[0],
      slipRateVertex.size());

    // Get total fault area asssociated with vertex (assembled over all cells)
    const double* areaVertex = areaSection->restrictPoint(v_fault);
    assert(0 != areaVertex);
    assert(1 == areaSection->getFiberDimension(v_fault));

    // Get Lagrange multiplier values from disp(t), and dispIncr(t->t+dt)
    dispTSection->restrictPoint(v_lagrange, &lagrangeTVertex[0],
      lagrangeTVertex.size());
    dispTIncrSection->restrictPoint(v_lagrange, &lagrangeTIncrVertex[0],
      lagrangeTIncrVertex.size());

    // Compute Lagrange multiplier at time t+dt
    lagrangeTpdtVertex = lagrangeTVertex + lagrangeTIncrVertex;

    // :KLUDGE: Solution at Lagrange constraint vertices is the
    // Lagrange multiplier value, which is currently the force.
    // Compute traction by dividing force by area
    assert(*areaVertex > 0);
    tractionTVertex = lagrangeTVertex / (*areaVertex);
      tractionTpdtVertex = lagrangeTpdtVertex / (*areaVertex);

    // Get friction properties and state variables.
    _friction->retrievePropsAndVars(v_fault);

    // Use fault constitutive model to compute traction associated with
    // friction.
    switch (spaceDim) { // switch
    case 1: { // case 1
      const double slipMag = 0.0;
      const double slipRateMag = 0.0;
      const double tractionNormal = tractionTpdtVertex[0];
      _friction->updateStateVars(slipMag, slipRateMag, tractionNormal, v_fault);
      break;
    } // case 1
    case 2: { // case 2
      const double slipMag = fabs(slipVertex[0]);
      const double slipRateMag = fabs(slipRateVertex[0]);
      const double tractionNormal = tractionTpdtVertex[1];
      _friction->updateStateVars(slipMag, slipRateMag, tractionNormal, v_fault);
      break;
    } // case 2
    case 3: { // case 3
      const double slipMag = 
	sqrt(slipVertex[0]*slipVertex[0] + slipVertex[1]*slipVertex[1]);
      const double slipRateMag = 
	sqrt(slipRateVertex[0]*slipRateVertex[0] + 
	     slipRateVertex[1]*slipRateVertex[1]);
      const double tractionNormal = tractionTpdtVertex[2];
      _friction->updateStateVars(slipMag, slipRateMag, tractionNormal, v_fault);
      break;
    } // case 3
    default:
      assert(0);
      throw std::logic_error("Unknown spatial dimension in "
			     "FaultCohesiveDyn::updateStateVars().");
    } // switch
  } // for
} // updateStateVars

// ----------------------------------------------------------------------
// Constrain solution based on friction.
void
pylith::faults::FaultCohesiveDyn::constrainSolnSpace(
				    topology::SolutionFields* const fields,
				    const double t,
				    const topology::Jacobian& jacobian)
{ // constrainSolnSpace
  /// Member prototype for _constrainSolnSpaceXD()
  typedef void (pylith::faults::FaultCohesiveDyn::*constrainSolnSpace_fn_type)
    (double_array*,
     const double_array&,
     const double_array&,
     const double_array&,
     const double);

  assert(0 != fields);
  assert(0 != _quadrature);
  assert(0 != _fields);
  assert(0 != _friction);

  _updateSlipRate(*fields);
  _sensitivitySetup(jacobian);

  // Update time step in friction (can vary).
  _friction->timeStep(_dt);

  const int spaceDim = _quadrature->spaceDim();

  // Allocate arrays for vertex values
  double_array tractionTVertex(spaceDim);
  double_array tractionTpdtVertex(spaceDim);
  double_array slipTpdtVertex(spaceDim);
  double_array lagrangeTpdtVertex(spaceDim);

  // Get sections
  double_array slipVertex(spaceDim);
  double_array dSlipVertex(spaceDim);
  const ALE::Obj<RealSection>& slipSection = _fields->get("slip").section();
  assert(!slipSection.isNull());

  double_array slipRateVertex(spaceDim);
  const ALE::Obj<RealSection>& slipRateSection =
      _fields->get("slip rate").section();
  assert(!slipRateSection.isNull());

  const ALE::Obj<RealSection>& areaSection = _fields->get("area").section();
  assert(!areaSection.isNull());

  double_array orientationVertex(spaceDim * spaceDim);
  const ALE::Obj<RealSection>& orientationSection =
      _fields->get("orientation").section();
  assert(!orientationSection.isNull());

  double_array lagrangeTVertex(spaceDim);
  double_array dispTVertexN(spaceDim);
  double_array dispTVertexP(spaceDim);
  const ALE::Obj<RealSection>& dispTSection = fields->get("disp(t)").section();
  assert(!dispTSection.isNull());

  double_array lagrangeTIncrVertex(spaceDim);
  double_array dispTIncrVertexN(spaceDim);
  double_array dispTIncrVertexP(spaceDim);
  const ALE::Obj<RealSection>& dispTIncrSection =
      fields->get("dispIncr(t->t+dt)").section();
  assert(!dispTIncrSection.isNull());

  double_array dLagrangeTpdtVertex(spaceDim);
  const ALE::Obj<RealSection>& dLagrangeTpdtSection =
      _fields->get("sensitivity dLagrange").section();
  assert(!dLagrangeTpdtSection.isNull());

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
  slipSection->view("SLIP");
  areaSection->view("AREA");
  dispTSection->view("DISP (t)");
  dispTIncrSection->view("DISP INCR (t->t+dt)");
#endif

  const int numVertices = _cohesiveVertices.size();
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_lagrange = _cohesiveVertices[iVertex].lagrange;
    const int v_fault = _cohesiveVertices[iVertex].fault;

    // Get slip
    slipSection->restrictPoint(v_fault, &slipVertex[0], slipVertex.size());

    // Get slip rate
    slipRateSection->restrictPoint(v_fault, &slipRateVertex[0],
      slipRateVertex.size());

    // Get total fault area asssociated with vertex (assembled over all cells)
    const double* areaVertex = areaSection->restrictPoint(v_fault);
    assert(0 != areaVertex);
    assert(1 == areaSection->getFiberDimension(v_fault));

    // Get Lagrange multiplier values from disp(t), and dispIncr(t->t+dt)
    dispTSection->restrictPoint(v_lagrange, &lagrangeTVertex[0],
      lagrangeTVertex.size());
    dispTIncrSection->restrictPoint(v_lagrange, &lagrangeTIncrVertex[0],
      lagrangeTIncrVertex.size());

    // Compute Lagrange multiplier at time t+dt
    lagrangeTpdtVertex = lagrangeTVertex + lagrangeTIncrVertex;
    dLagrangeTpdtVertex = 0.0;

    // :KLUDGE: Solution at Lagrange constraint vertices is the
    // Lagrange multiplier value, which is currently the force.
    // Compute traction by dividing force by area
    assert(*areaVertex > 0);
    tractionTVertex = lagrangeTVertex / (*areaVertex);
    tractionTpdtVertex = lagrangeTpdtVertex / (*areaVertex);

    // Get friction properties and state variables.
    _friction->retrievePropsAndVars(v_fault);

    // Use fault constitutive model to compute traction associated with
    // friction.
    CALL_MEMBER_FN(*this,
		   constrainSolnSpaceFn)(&dLagrangeTpdtVertex,
					 slipVertex, slipRateVertex,
					 tractionTpdtVertex, *areaVertex);

    assert(dLagrangeTpdtVertex.size() ==
        dLagrangeTpdtSection->getFiberDimension(v_fault));
    dLagrangeTpdtSection->updatePoint(v_fault, &dLagrangeTpdtVertex[0]);
  } // for

  // Solve sensitivity problem for negative side of the fault.
  bool negativeSide = true;
  _sensitivityUpdateJacobian(negativeSide, jacobian, *fields);
  _sensitivityReformResidual(negativeSide);
  _sensitivitySolve();
  _sensitivityUpdateSoln(negativeSide);

  // Solve sensitivity problem for positive side of the fault.
  negativeSide = false;
  _sensitivityUpdateJacobian(negativeSide, jacobian, *fields);
  _sensitivityReformResidual(negativeSide);
  _sensitivitySolve();
  _sensitivityUpdateSoln(negativeSide);

  // Update slip field based on solution of sensitivity problem and
  // increment in Lagrange multipliers.
  const ALE::Obj<RealSection>& dispRelSection =
    _fields->get("sensitivity dispRel").section();
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_fault = _cohesiveVertices[iVertex].fault;
    const int v_lagrange = _cohesiveVertices[iVertex].lagrange;
    const int v_negative = _cohesiveVertices[iVertex].negative;
    const int v_positive = _cohesiveVertices[iVertex].positive;

    // Get fault orientation
    orientationSection->restrictPoint(v_fault, &orientationVertex[0],
    orientationVertex.size());

    // Get change in relative displacement.
    const double* dispRelVertex = dispRelSection->restrictPoint(v_fault);
    assert(0 != dispRelVertex);
    assert(spaceDim == dispRelSection->getFiberDimension(v_fault));

    // Get Lagrange multiplier at time t
    dispTSection->restrictPoint(v_lagrange, &lagrangeTVertex[0],
				lagrangeTVertex.size());

    // Get Lagrange multiplier increment at time t
    dispTIncrSection->restrictPoint(v_lagrange, &lagrangeTIncrVertex[0],
				    lagrangeTIncrVertex.size());

    // Get change in Lagrange multiplier.
    dLagrangeTpdtSection->restrictPoint(v_fault, &dLagrangeTpdtVertex[0],
					dLagrangeTpdtVertex.size());
    // Compute change in slip.
    dSlipVertex = 0.0;
    for (int iDim = 0; iDim < spaceDim; ++iDim)
      for (int kDim = 0; kDim < spaceDim; ++kDim)
        dSlipVertex[iDim] += 
	  orientationVertex[iDim*spaceDim+kDim] * dispRelVertex[kDim];

    // Do not allow fault interpenetration and set fault opening to
    // zero if fault is under compression.
    const int indexN = spaceDim - 1;
    if (dSlipVertex[indexN] < 0.0)
      dSlipVertex[indexN] = 0.0;
    const double lagrangeTpdtNormal = lagrangeTVertex[indexN] + 
      lagrangeTIncrVertex[indexN] + dLagrangeTpdtVertex[indexN];
    if (lagrangeTpdtNormal < 0.0)
      dSlipVertex[indexN] = 0.0;

    // Set change in slip.
    assert(dSlipVertex.size() ==
        slipSection->getFiberDimension(v_fault));
    slipSection->updateAddPoint(v_fault, &dSlipVertex[0]);
    
    // Update Lagrange multiplier increment.
    assert(dLagrangeTpdtVertex.size() ==
	   dispTIncrSection->getFiberDimension(v_lagrange));
    dispTIncrSection->updateAddPoint(v_lagrange, &dLagrangeTpdtVertex[0]);

    // Compute change in displacement field.
    dispTIncrVertexN = 0.0;
    for (int iDim = 0; iDim < spaceDim; ++iDim)
      for (int kDim = 0; kDim < spaceDim; ++kDim)
        dispTIncrVertexN[iDim] += 
	  orientationVertex[kDim*spaceDim+iDim] * dSlipVertex[kDim];
    
    // Update displacement field
    dispTIncrVertexN *= -0.5;
    assert(dispTIncrVertexN.size() ==
	   dispTIncrSection->getFiberDimension(v_negative));
    dispTIncrSection->updateAddPoint(v_negative, &dispTIncrVertexN[0]);
    
    dispTIncrVertexP = -dispTIncrVertexN;
    assert(dispTIncrVertexP.size() ==
	   dispTIncrSection->getFiberDimension(v_positive));
    dispTIncrSection->updateAddPoint(v_positive, &dispTIncrVertexP[0]);
    
  } // for

#if 0 // DEBUGGING
  dLagrangeTpdtSection->view("AFTER dLagrange");
  dispTIncrSection->view("AFTER DISP INCR (t->t+dt)");
  slipSection->view("AFTER SLIP");
#endif
} // constrainSolnSpace

// ----------------------------------------------------------------------
// Adjust solution from solver with lumped Jacobian to match Lagrange
// multiplier constraints.
void
pylith::faults::FaultCohesiveDyn::adjustSolnLumped(
			 topology::SolutionFields* const fields,
			 const topology::Field<topology::Mesh>& jacobian)
{ // adjustSolnLumped
  /// Member prototype for _constrainSolnSpaceXD()
  typedef void (pylith::faults::FaultCohesiveDyn::*constrainSolnSpace_fn_type)
    (double_array*,
     const double_array&,
     const double_array&,
     const double_array&,
     const double);

  /// Member prototype for _sensitivitySolveLumpedXD()
  typedef void (pylith::faults::FaultCohesiveDyn::*sensitivitySolveLumped_fn_type)
    (double_array*,
     const double_array&,
     const double_array&,
     const double_array&);

  /// Member prototype for _adjustSolnLumpedXD()
  typedef void (pylith::faults::FaultCohesiveDyn::*adjustSolnLumped_fn_type)
    (double_array*, double_array*, double_array*,
     const double_array&, const double_array&,
     const double_array&, const double_array&,
     const double_array&, const double_array&,
     const double_array&, const double_array&);


  assert(0 != fields);
  assert(0 != _quadrature);

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
  const int orientationSize = spaceDim * spaceDim;

  // Allocate arrays for vertex values
  double_array tractionTVertex(spaceDim);
  double_array tractionTpdtVertex(spaceDim);
  double_array slipTpdtVertex(spaceDim);
  double_array lagrangeTpdtVertex(spaceDim);
  double_array dLagrangeTpdtVertex(spaceDim);

  // Update time step in friction (can vary).
  _friction->timeStep(_dt);

  // Get section information
  double_array slipVertex(spaceDim);
  const ALE::Obj<RealSection>& slipSection = _fields->get("slip").section();
  assert(!slipSection.isNull());

  double_array slipRateVertex(spaceDim);
  const ALE::Obj<RealSection>& slipRateSection =
      _fields->get("slip rate").section();
  assert(!slipRateSection.isNull());

  const ALE::Obj<RealSection>& areaSection = _fields->get("area").section();
  assert(!areaSection.isNull());

  double_array orientationVertex(orientationSize);
  const ALE::Obj<RealSection>& orientationSection =
      _fields->get("orientation").section();
  assert(!orientationSection.isNull());

  double_array dispTVertexN(spaceDim);
  double_array dispTVertexP(spaceDim);
  double_array lagrangeTVertex(spaceDim);
  const ALE::Obj<RealSection>& dispTSection = fields->get("disp(t)").section();
  assert(!dispTSection.isNull());

  double_array dispTIncrVertexN(spaceDim);
  double_array dispTIncrVertexP(spaceDim);
  double_array lagrangeTIncrVertex(spaceDim);
  const ALE::Obj<RealSection>& dispTIncrSection =
      fields->get("dispIncr(t->t+dt)").section();
  assert(!dispTIncrSection.isNull());

  const ALE::Obj<RealSection>& dispTIncrAdjSection = fields->get(
    "dispIncr adjust").section();
  assert(!dispTIncrAdjSection.isNull());

  double_array jacobianVertexN(spaceDim);
  double_array jacobianVertexP(spaceDim);
  const ALE::Obj<RealSection>& jacobianSection = jacobian.section();
  assert(!jacobianSection.isNull());

  double_array residualVertexN(spaceDim);
  double_array residualVertexP(spaceDim);
  const ALE::Obj<RealSection>& residualSection =
      fields->get("residual").section();

  const ALE::Obj<SieveMesh>& sieveMesh = fields->mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::order_type>& globalOrder =
    sieveMesh->getFactory()->getGlobalOrder(sieveMesh, "default", 
					    jacobianSection);
  assert(!globalOrder.isNull());

  adjustSolnLumped_fn_type adjustSolnLumpedFn;
  constrainSolnSpace_fn_type constrainSolnSpaceFn;
  sensitivitySolveLumped_fn_type sensitivitySolveLumpedFn;
  switch (spaceDim) { // switch
  case 1:
    adjustSolnLumpedFn = 
      &pylith::faults::FaultCohesiveDyn::_adjustSolnLumped1D;
    constrainSolnSpaceFn = 
      &pylith::faults::FaultCohesiveDyn::_constrainSolnSpace1D;
    sensitivitySolveLumpedFn =
      &pylith::faults::FaultCohesiveDyn::_sensitivitySolveLumped1D;
    break;
  case 2: 
    adjustSolnLumpedFn = 
      &pylith::faults::FaultCohesiveDyn::_adjustSolnLumped2D;
    constrainSolnSpaceFn = 
      &pylith::faults::FaultCohesiveDyn::_constrainSolnSpace2D;
    sensitivitySolveLumpedFn =
      &pylith::faults::FaultCohesiveDyn::_sensitivitySolveLumped2D;
    break;
  case 3:
    adjustSolnLumpedFn = 
      &pylith::faults::FaultCohesiveDyn::_adjustSolnLumped3D;
    constrainSolnSpaceFn = 
      &pylith::faults::FaultCohesiveDyn::_constrainSolnSpace3D;
    sensitivitySolveLumpedFn =
      &pylith::faults::FaultCohesiveDyn::_sensitivitySolveLumped3D;
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
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_lagrange = _cohesiveVertices[iVertex].lagrange;
    const int v_fault = _cohesiveVertices[iVertex].fault;
    const int v_negative = _cohesiveVertices[iVertex].negative;
    const int v_positive = _cohesiveVertices[iVertex].positive;

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(restrictEvent);
#endif

    // Get slip
    slipSection->restrictPoint(v_fault, &slipVertex[0], slipVertex.size());

    // Get slip rate
    slipRateSection->restrictPoint(v_fault, &slipRateVertex[0],
				   slipRateVertex.size());
    
    // Get total fault area asssociated with vertex (assembled over all cells)
    const double* areaVertex = areaSection->restrictPoint(v_fault);
    assert(0 != areaVertex);
    assert(1 == areaSection->getFiberDimension(v_fault));
    
    // Get fault orientation
    orientationSection->restrictPoint(v_fault, &orientationVertex[0],
				      orientationVertex.size());
    
    // Get Jacobian at vertices on positive and negative sides of the fault.
    jacobianSection->restrictPoint(v_negative, &jacobianVertexN[0],
				   jacobianVertexN.size());
    jacobianSection->restrictPoint(v_positive, &jacobianVertexP[0],
				   jacobianVertexP.size());
    
    // Get residual at cohesive cell's vertices.
    residualSection->restrictPoint(v_negative, &residualVertexN[0], 
				   residualVertexN.size());
    residualSection->restrictPoint(v_positive, &residualVertexP[0], 
				   residualVertexP.size());

    // Get disp(t) at cohesive cell's vertices.
    dispTSection->restrictPoint(v_negative, &dispTVertexN[0], 
				dispTVertexN.size());
    dispTSection->restrictPoint(v_positive, &dispTVertexP[0], 
				dispTVertexP.size());
    
    // Get Lagrange multiplier values from disp(t), and dispIncr(t->t+dt)
    dispTSection->restrictPoint(v_lagrange, &lagrangeTVertex[0],
				lagrangeTVertex.size());
    dispTIncrSection->restrictPoint(v_lagrange, &lagrangeTIncrVertex[0],
				    lagrangeTIncrVertex.size());
    

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(restrictEvent);
    _logger->eventBegin(computeEvent);
#endif

    CALL_MEMBER_FN(*this,
		   adjustSolnLumpedFn)(&lagrangeTIncrVertex,
				       &dispTIncrVertexN,
				       &dispTIncrVertexP,
				       slipVertex,
				       orientationVertex,
				       dispTVertexN,
				       dispTVertexP,
				       residualVertexN,
				       residualVertexP,
				       jacobianVertexN,
				       jacobianVertexP);

    
    // Compute Lagrange multiplier at time t+dt
    lagrangeTpdtVertex = lagrangeTVertex + lagrangeTIncrVertex;
    dLagrangeTpdtVertex = 0.0;
    
    // :KLUDGE: Solution at Lagrange constraint vertices is the
    // Lagrange multiplier value, which is currently the force.
    // Compute traction by dividing force by area
    assert(*areaVertex > 0);
    tractionTVertex = lagrangeTVertex / (*areaVertex);
    tractionTpdtVertex = lagrangeTpdtVertex / (*areaVertex);
    
    // Get friction properties and state variables.
    _friction->retrievePropsAndVars(v_fault);

    CALL_MEMBER_FN(*this,
		   constrainSolnSpaceFn)(&dLagrangeTpdtVertex,
					 slipVertex, slipRateVertex,
					 tractionTpdtVertex, *areaVertex);
    CALL_MEMBER_FN(*this,
       sensitivitySolveLumpedFn)(&slipVertex,
           dLagrangeTpdtVertex, jacobianVertexN, jacobianVertexP);

    lagrangeTIncrVertex += dLagrangeTpdtVertex;

    // :TODO: Refactor this into sensitivitySolveLumpedXD().
    switch (spaceDim) { // switch
    case 1: {
      assert(jacobianVertexN[0] > 0.0);
      assert(jacobianVertexP[0] > 0.0);

      const double dlp = lagrangeTIncrVertex[0];

      // Update displacements at negative vertex
      dispTIncrVertexN[0] = +1.0 / jacobianVertexN[0] * dlp;
  
      // Update displacements at positive vertex
      dispTIncrVertexP[0] = -1.0 / jacobianVertexP[0] * dlp;
  
      break;
    } // case 1
    case 2: {
      assert(jacobianVertexN[0] > 0.0);
      assert(jacobianVertexN[1] > 0.0);
      assert(jacobianVertexP[0] > 0.0);
      assert(jacobianVertexP[1] > 0.0);

      // Check to make sure Jacobian is same at all DOF for
      // vertices i and j (means S is diagonal with equal enties).
      assert(jacobianVertexN[0] == jacobianVertexN[1]);
      assert(jacobianVertexP[0] == jacobianVertexP[1]);

      const double Cpx = orientationVertex[0];
      const double Cpy = orientationVertex[1];
      const double Cqx = orientationVertex[2];
      const double Cqy = orientationVertex[3];

      const double dlp = lagrangeTIncrVertex[0];
      const double dlq = lagrangeTIncrVertex[1];

      const double dlx = Cpx * dlp + Cqx * dlq;
      const double dly = Cpy * dlp + Cqy * dlq;
  
      // Update displacements at negative vertex.
      dispTIncrVertexN[0] = dlx / jacobianVertexN[0];
      dispTIncrVertexN[1] = dly / jacobianVertexN[0];
  
      // Update displacements at positive vertex.
      dispTIncrVertexP[0] = -dlx / jacobianVertexP[0];
      dispTIncrVertexP[1] = -dly / jacobianVertexP[0];

      break;
    } // case 2
    case 3: {
      assert(jacobianVertexN[0] > 0.0);
      assert(jacobianVertexN[1] > 0.0);
      assert(jacobianVertexN[2] > 0.0);
      assert(jacobianVertexP[0] > 0.0);
      assert(jacobianVertexP[1] > 0.0);
      assert(jacobianVertexP[2] > 0.0);

      // Check to make sure Jacobian is same at all DOF for
      // vertices i and j (means S is diagonal with equal enties).
      assert(jacobianVertexN[0] == jacobianVertexN[1] && 
	     jacobianVertexN[0] == jacobianVertexN[2]);
      assert(jacobianVertexP[0] == jacobianVertexP[1] && 
	     jacobianVertexP[0] == jacobianVertexP[2]);

      const double Cpx = orientationVertex[0];
      const double Cpy = orientationVertex[1];
      const double Cpz = orientationVertex[2];
      const double Cqx = orientationVertex[3];
      const double Cqy = orientationVertex[4];
      const double Cqz = orientationVertex[5];
      const double Crx = orientationVertex[6];
      const double Cry = orientationVertex[7];
      const double Crz = orientationVertex[8];

      const double dlp = lagrangeTIncrVertex[0];
      const double dlq = lagrangeTIncrVertex[1];
      const double dlr = lagrangeTIncrVertex[2];

      const double dlx = Cpx * dlp + Cqx * dlq + Crx * dlr;
      const double dly = Cpy * dlp + Cqy * dlq + Cry * dlr;
      const double dlz = Cpz * dlp + Cqz * dlq + Crz * dlr;

      // Update displacements at negative vertex.
      dispTIncrVertexN[0] = dlx / jacobianVertexN[0];
      dispTIncrVertexN[1] = dly / jacobianVertexN[1];
      dispTIncrVertexN[2] = dlz / jacobianVertexN[2];

      // Update displacements at positive vertex.
      dispTIncrVertexP[0] = -dlx / jacobianVertexP[0];
      dispTIncrVertexP[1] = -dly / jacobianVertexP[1];
      dispTIncrVertexP[2] = -dlz / jacobianVertexP[2];

      break;
    } // case 3
    default:
      assert(0);
      throw std::logic_error("Unknown spatial dimension in "
			     "FaultCohesiveDyn::adjustSolnLumped().");
    } // switch

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(computeEvent);
    _logger->eventBegin(updateEvent);
#endif

    // Compute contribution to adjusting solution only if Lagrange
    // constraint is local (the adjustment is assembled across processors).
    if (globalOrder->isLocal(v_lagrange)) {
      // Adjust displacements to account for Lagrange multiplier values
      // (assumed to be zero in perliminary solve).
      assert(dispTIncrVertexN.size() == 
	     dispTIncrAdjSection->getFiberDimension(v_negative));
      dispTIncrAdjSection->updateAddPoint(v_negative, &dispTIncrVertexN[0]);
      
      assert(dispTIncrVertexP.size() == 
	     dispTIncrAdjSection->getFiberDimension(v_positive));
      dispTIncrAdjSection->updateAddPoint(v_positive, &dispTIncrVertexP[0]);
    } // if

    // The Lagrange multiplier and slip are NOT assembled across processors.

    // Set Lagrange multiplier value. Value from preliminary solve is
    // bogus due to artificial diagonal entry of 1.0.
    assert(lagrangeTIncrVertex.size() == 
	   dispTIncrSection->getFiberDimension(v_lagrange));
    dispTIncrSection->updatePoint(v_lagrange, &lagrangeTIncrVertex[0]);

    // Update the slip estimate based on adjustment to the Lagrange
    // multiplier values.
    assert(slipVertex.size() ==
        slipSection->getFiberDimension(v_fault));
    slipSection->updatePoint(v_fault, &slipVertex[0]);
#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(updateEvent);
#endif
  } // for

#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventEnd(computeEvent);
#endif
} // adjustSolnLumped

// ----------------------------------------------------------------------
// Get vertex field associated with integrator.
const pylith::topology::Field<pylith::topology::SubMesh>&
pylith::faults::FaultCohesiveDyn::vertexField(const char* name,
                                               const topology::SolutionFields* fields)
{ // vertexField
  assert(0 != _faultMesh);
  assert(0 != _quadrature);
  assert(0 != _normalizer);
  assert(0 != _fields);
  assert(0 != _friction);

  const int cohesiveDim = _faultMesh->dimension();
  const int spaceDim = _quadrature->spaceDim();

  double scale = 0.0;
  int fiberDim = 0;
  if (0 == strcasecmp("slip", name)) {
    const topology::Field<topology::SubMesh>& slip = _fields->get("slip");
    return slip;

  } else if (0 == strcasecmp("slip_rate", name)) {
    const topology::Field<topology::SubMesh>& slipRate =
      _fields->get("slip rate");
    return slipRate;

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

  } else if (0 == strcasecmp("initial_traction", name)) {
    assert(0 != _dbInitialTract);
    _allocateBufferVectorField();
    topology::Field<topology::SubMesh>& buffer =
        _fields->get("buffer (vector)");
    _getInitialTractions(&buffer);
    return buffer;

  } else if (0 == strcasecmp("traction", name)) {
    assert(0 != fields);
    const topology::Field<topology::Mesh>& dispT = fields->get("disp(t)");
    _allocateBufferVectorField();
    topology::Field<topology::SubMesh>& buffer =
        _fields->get("buffer (vector)");
    _calcTractions(&buffer, dispT);
    return buffer;

  } else if (_friction->hasProperty(name) || _friction->hasStateVar(name)) {
    assert(0 != _fields);
    if (!_fields->hasField("buffer (other)"))
      _fields->add("buffer (other)", "buffer");
    topology::Field<topology::SubMesh>& buffer = _fields->get("buffer (other)");
    _friction->getField(&buffer, name);
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
pylith::faults::FaultCohesiveDyn::cellField(const char* name,
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
void
pylith::faults::FaultCohesiveDyn::_setupInitialTractions(void)
{ // _setupInitialTractions
  assert(0 != _normalizer);
  assert(0 != _quadrature);

  // If no initial tractions specified, leave method
  if (0 == _dbInitialTract)
    return;

  assert(0 != _normalizer);
  const double pressureScale = _normalizer->pressureScale();
  const double lengthScale = _normalizer->lengthScale();

  // Get quadrature information
  const int numQuadPts = _quadrature->numQuadPts();
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const double_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == numQuadPts);

  double_array quadPtsGlobal(numQuadPts*spaceDim);

  // Create section to hold initial tractions.
  _fields->add("initial forces", "initial_forces");
  topology::Field<topology::SubMesh>& forcesInitial = 
    _fields->get("initial forces");
  topology::Field<topology::SubMesh>& slip = _fields->get("slip");
  forcesInitial.cloneSection(slip);
  forcesInitial.scale(pressureScale);
  const ALE::Obj<RealSection>& forcesInitialSection = forcesInitial.section();
  assert(!forcesInitialSection.isNull());
  double_array forcesInitialCell(numBasis*spaceDim);
  double_array tractionQuadPt(spaceDim);
  topology::Mesh::UpdateAddVisitor forcesInitialVisitor(*forcesInitialSection,
        &forcesInitialCell[0]);

  assert(0 != _dbInitialTract);
  _dbInitialTract->open();
  switch (spaceDim) { // switch
  case 1: {
    const char* valueNames[] = { "traction-normal" };
    _dbInitialTract->queryVals(valueNames, 1);
    break;
  } // case 1
  case 2: {
    const char* valueNames[] = { "traction-shear", "traction-normal" };
    _dbInitialTract->queryVals(valueNames, 2);
    break;
  } // case 2
  case 3: {
    const char* valueNames[] = { "traction-shear-leftlateral",
				 "traction-shear-updip", "traction-normal" };
    _dbInitialTract->queryVals(valueNames, 3);
    break;
  } // case 3
  default:
    std::cerr << "Bad spatial dimension '" << spaceDim << "'." << std::endl;
    assert(0);
    throw std::logic_error("Bad spatial dimension in Neumann.");
  } // switch
  
  // Get cells associated with fault
  const ALE::Obj<SieveSubMesh>& faultSieveMesh = _faultMesh->sieveMesh();
  assert(!faultSieveMesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& cells = 
    faultSieveMesh->heightStratum(0);
  assert(!cells.isNull());
  const SieveSubMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveSubMesh::label_sequence::iterator cellsEnd = cells->end();

  const spatialdata::geocoords::CoordSys* cs = _faultMesh->coordsys();
  assert(0 != cs);

#if !defined(PRECOMPUTE_GEOMETRY)
  double_array coordinatesCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& coordinates =
    faultSieveMesh->getRealSection("coordinates");
  RestrictVisitor coordsVisitor(*coordinates,
        coordinatesCell.size(), &coordinatesCell[0]);
#endif

  for (SieveSubMesh::label_sequence::iterator c_iter=cellsBegin;
       c_iter != cellsEnd;
       ++c_iter) {
    // Compute geometry information for current cell
#if defined(PRECOMPUTE_GEOMETRY)
    _quadrature->retrieveGeometry(*c_iter);
#else
    coordsVisitor.clear();
    faultSieveMesh->restrictClosure(*c_iter, coordsVisitor);
    _quadrature->computeGeometry(coordinatesCell, *c_iter);
#endif

    const double_array& quadPtsNonDim = _quadrature->quadPts();
    quadPtsGlobal = quadPtsNonDim;
    _normalizer->dimensionalize(&quadPtsGlobal[0], quadPtsGlobal.size(),
        lengthScale);
    forcesInitialCell = 0.0;

    // Loop over quadrature points in cell and query database
    for (int iQuadPt=0, index=0;
        iQuadPt < numQuadPts;
        ++iQuadPt, index+=spaceDim) {

      tractionQuadPt = 0.0;
      int err = _dbInitialTract->query(&tractionQuadPt[0], spaceDim,
          &quadPtsGlobal[index], spaceDim, cs);
      if (err) {
        std::ostringstream msg;
        msg << "Could not find parameters for physical properties at \n" << "(";
        for (int i = 0; i < spaceDim; ++i)
          msg << "  " << quadPtsGlobal[index + i];
        msg << ") in friction model " << label() << "\n"
            << "using spatial database '" << _dbInitialTract->label() << "'.";
        throw std::runtime_error(msg.str());
      } // if
      tractionQuadPt /= pressureScale;

      // Get cell geometry information that depends on cell
      const double_array& basis = _quadrature->basis();
      const double_array& jacobianDet = _quadrature->jacobianDet();

      // Integrate tractions over cell.
      const double wt = quadWts[iQuadPt] * jacobianDet[iQuadPt];
      for (int iBasis=0; iBasis < numBasis; ++iBasis) {
	const double valI = wt*basis[iQuadPt*numBasis+iBasis];
	for (int jBasis=0; jBasis < numBasis; ++jBasis) {
	  const double valIJ = valI * basis[iQuadPt*numBasis+jBasis];
	  for (int iDim=0; iDim < spaceDim; ++iDim)
	    forcesInitialCell[iBasis*spaceDim+iDim] += 
	      tractionQuadPt[iDim] * valIJ;
	} // for
      } // for
    } // for
    // Assemble cell contribution into field
    forcesInitialVisitor.clear();
    faultSieveMesh->updateClosure(*c_iter, forcesInitialVisitor);
  } // for
  // Close properties database
  _dbInitialTract->close();

  forcesInitial.complete(); // Assemble contributions

  // Rotate forces from fault coordinate system to global coordinate system
  const int orientationSize = spaceDim * spaceDim;
  const ALE::Obj<RealSection>& orientationSection =
    _fields->get("orientation").section();

  double_array forcesInitialVertexFault(spaceDim);
  double_array forcesInitialVertexGlobal(spaceDim);

  const int numVertices = _cohesiveVertices.size();
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_fault = _cohesiveVertices[iVertex].fault;

    assert(orientationSize == orientationSection->getFiberDimension(v_fault));
    assert(spaceDim == forcesInitialSection->getFiberDimension(v_fault));

    const double* orientationVertex = 
      orientationSection->restrictPoint(v_fault);

    forcesInitialSection->restrictPoint(v_fault, 
					&forcesInitialVertexFault[0], 
					forcesInitialVertexFault.size());

    forcesInitialVertexGlobal = 0.0;
    for (int iDim = 0; iDim < spaceDim; ++iDim)
      for (int kDim = 0; kDim < spaceDim; ++kDim)
        forcesInitialVertexGlobal[iDim] +=
          forcesInitialVertexFault[kDim] * 
	  orientationVertex[kDim*spaceDim+iDim];

    assert(forcesInitialVertexGlobal.size() == 
	   forcesInitialSection->getFiberDimension(v_fault));
    forcesInitialSection->updatePoint(v_fault, &forcesInitialVertexGlobal[0]);
  } // for

  //forcesInitial.view("INITIAL FORCES"); // DEBUGGING
} // _setupInitialTractions

// ----------------------------------------------------------------------
// Compute tractions on fault surface using solution.
void
pylith::faults::FaultCohesiveDyn::_calcTractions(
    topology::Field<topology::SubMesh>* tractions,
    const topology::Field<topology::Mesh>& dispT)
{ // _calcTractions
  assert(0 != tractions);
  assert(0 != _faultMesh);
  assert(0 != _fields);
  assert(0 != _normalizer);

  // Fiber dimension of tractions matches spatial dimension.
  const int fiberDim = _quadrature->spaceDim();
  double_array tractionsVertex(fiberDim);

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
  const double pressureScale = _normalizer->pressureScale();
  tractions->label("traction");
  tractions->scale(pressureScale);
  tractions->zero();

  const int numVertices = _cohesiveVertices.size();
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_lagrange = _cohesiveVertices[iVertex].lagrange;
    const int v_fault = _cohesiveVertices[iVertex].fault;

    assert(fiberDim == dispTSection->getFiberDimension(v_lagrange));
    assert(fiberDim == tractionsSection->getFiberDimension(v_fault));
    assert(1 == areaSection->getFiberDimension(v_fault));

    const double* dispTVertex = dispTSection->restrictPoint(v_lagrange);
    assert(0 != dispTVertex);
    const double* areaVertex = areaSection->restrictPoint(v_fault);
    assert(0 != areaVertex);

    for (int i=0; i < fiberDim; ++i)
      tractionsVertex[i] = dispTVertex[i] / areaVertex[0];

    assert(tractionsVertex.size() == 
	   tractionsSection->getFiberDimension(v_fault));
    tractionsSection->updatePoint(v_fault, &tractionsVertex[0]);
  } // for

  PetscLogFlops(numVertices * (1 + fiberDim) );

#if 0 // DEBUGGING
  tractions->view("TRACTIONS");
#endif

} // _calcTractions

// ----------------------------------------------------------------------
// Compute initial tractions on fault surface.
void
pylith::faults::FaultCohesiveDyn::_getInitialTractions(
    topology::Field<topology::SubMesh>* tractions)
{ // _getInitialTractions
  assert(0 != tractions);
  assert(0 != _faultMesh);
  assert(0 != _fields);
  assert(0 != _normalizer);

  // Fiber dimension of tractions matches spatial dimension.
  const int spaceDim = _quadrature->spaceDim();
  double_array tractionsVertexGlobal(spaceDim);
  double_array tractionsVertexFault(spaceDim);

  // Get sections.
  const ALE::Obj<RealSection>& areaSection = _fields->get("area").section();
  assert(!areaSection.isNull());
  const ALE::Obj<RealSection>& orientationSection = 
    _fields->get("orientation").section();
  assert(!orientationSection.isNull());
  const ALE::Obj<RealSection>& forcesInitialSection = 
    _fields->get("initial forces").section();
  assert(!forcesInitialSection.isNull());

  // Allocate buffer for tractions field (if necessary).
  const ALE::Obj<RealSection>& tractionsSection = tractions->section();
  if (tractionsSection.isNull()) {
    ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
    //logger.stagePush("Fault");

    const topology::Field<topology::SubMesh>& slip = _fields->get("slip");
    tractions->cloneSection(slip);

    //logger.stagePop();
  } // if
  const double pressureScale = _normalizer->pressureScale();
  tractions->label("initial_traction");
  tractions->scale(pressureScale);
  tractions->zero();

  const int numVertices = _cohesiveVertices.size();
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_fault = _cohesiveVertices[iVertex].fault;

    assert(spaceDim == forcesInitialSection->getFiberDimension(v_fault));
    assert(spaceDim == tractionsSection->getFiberDimension(v_fault));
    assert(1 == areaSection->getFiberDimension(v_fault));

    const double* forcesInitialVertex = 
      forcesInitialSection->restrictPoint(v_fault);
    assert(0 != forcesInitialVertex);
    const double* areaVertex = areaSection->restrictPoint(v_fault);
    assert(0 != areaVertex);
    const double* orientationVertex = 
      orientationSection->restrictPoint(v_fault);
    assert(0 != orientationVertex);

    for (int i = 0; i < spaceDim; ++i)
      tractionsVertexGlobal[i] = forcesInitialVertex[i] / areaVertex[0];

    // Rotate from global coordinate system to local coordinate system
    tractionsVertexFault = 0.0;
    for (int iDim = 0; iDim < spaceDim; ++iDim)
      for (int kDim = 0; kDim < spaceDim; ++kDim)
        tractionsVertexFault[iDim] +=
          tractionsVertexGlobal[kDim] * orientationVertex[iDim*spaceDim+kDim];
    
    assert(tractionsVertexFault.size() == 
	   tractionsSection->getFiberDimension(v_fault));
    tractionsSection->updatePoint(v_fault, &tractionsVertexFault[0]);
  } // for

  PetscLogFlops(numVertices * (1 + spaceDim) );

#if 0 // DEBUGGING
  tractions->view("TRACTIONS");
#endif

} // _getInitialTractions

// ----------------------------------------------------------------------
// Update slip rate associated with Lagrange vertex k corresponding
// to diffential velocity between conventional vertices i and j.
void
pylith::faults::FaultCohesiveDyn::_updateSlipRate(const topology::SolutionFields& fields)
{ // _updateSlipRate
  assert(0 != _fields);

  const int spaceDim = _quadrature->spaceDim();

  // Get section information
  double_array velocityVertexN(spaceDim);
  double_array velocityVertexP(spaceDim);
  const ALE::Obj<RealSection>& velocitySection =
      fields.get("velocity(t)").section();
  assert(!velocitySection.isNull());

  double_array slipRateVertex(spaceDim);
  const ALE::Obj<RealSection>& slipRateSection =
      _fields->get("slip rate").section();

  double_array orientationVertex(spaceDim*spaceDim);
  const ALE::Obj<RealSection>& orientationSection =
      _fields->get("orientation").section();
  assert(!orientationSection.isNull());

  const int numVertices = _cohesiveVertices.size();
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_lagrange = _cohesiveVertices[iVertex].lagrange;
    const int v_fault = _cohesiveVertices[iVertex].fault;
    const int v_negative = _cohesiveVertices[iVertex].negative;
    const int v_positive = _cohesiveVertices[iVertex].positive;

    // Get values
    const double* velocityVertexN = velocitySection->restrictPoint(v_negative);
    assert(0 != velocityVertexN);
    assert(spaceDim == velocitySection->getFiberDimension(v_negative));

    const double* velocityVertexP = velocitySection->restrictPoint(v_positive);
    assert(0 != velocityVertexP);
    assert(spaceDim == velocitySection->getFiberDimension(v_positive));

    const double* orientationVertex = orientationSection->restrictPoint(v_fault);
    assert(0 != orientationVertex);
    assert(spaceDim*spaceDim == orientationSection->getFiberDimension(v_fault));

    slipRateVertex = 0.0;
    // Velocity for negative vertex.
    for (int iDim = 0; iDim < spaceDim; ++iDim)
      for (int kDim = 0; kDim < spaceDim; ++kDim)
        slipRateVertex[iDim] +=
          velocityVertexN[kDim] * -orientationVertex[iDim*spaceDim+kDim];

    // Velocity for positive vertex.
    for (int iDim = 0; iDim < spaceDim; ++iDim)
      for (int kDim = 0; kDim < spaceDim; ++kDim)
        slipRateVertex[iDim] +=
          velocityVertexP[kDim] * +orientationVertex[iDim*spaceDim+kDim];

    // Update slip rate field.
    assert(slipRateVertex.size() == 
	   slipRateSection->getFiberDimension(v_fault));
    slipRateSection->updatePoint(v_fault, &slipRateVertex[0]);
  } // for

  PetscLogFlops(numVertices*spaceDim*spaceDim*4);
} // _updateSlipRate

// ----------------------------------------------------------------------
// Setup sensitivity problem to compute change in slip given change in Lagrange multipliers.
void
pylith::faults::FaultCohesiveDyn::_sensitivitySetup(const topology::Jacobian& jacobian)
{ // _sensitivitySetup
  assert(0 != _fields);
  assert(0 != _quadrature);

  const int spaceDim = _quadrature->spaceDim();

  // Setup fields involved in sensitivity solve.
  if (!_fields->hasField("sensitivity solution")) {
    _fields->add("sensitivity solution", "sensitivity_soln");
    topology::Field<topology::SubMesh>& solution =
        _fields->get("sensitivity solution");
    const topology::Field<topology::SubMesh>& slip =
        _fields->get("slip");
    solution.cloneSection(slip);
    solution.createVector();
    solution.createScatter();
  } // if
  const topology::Field<topology::SubMesh>& solution =
      _fields->get("sensitivity solution");

  if (!_fields->hasField("sensitivity residual")) {
    _fields->add("sensitivity residual", "sensitivity_residual");
    topology::Field<topology::SubMesh>& residual =
        _fields->get("sensitivity residual");
    residual.cloneSection(solution);
    residual.createVector();
    residual.createScatter();
  } // if

  if (!_fields->hasField("sensitivity dispRel")) {
    _fields->add("sensitivity dispRel", "sensitivity_disprel");
    topology::Field<topology::SubMesh>& dispRel =
        _fields->get("sensitivity dispRel");
    dispRel.cloneSection(solution);
  } // if
  topology::Field<topology::SubMesh>& dispRel =
    _fields->get("sensitivity dispRel");
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
  assert(0 != _jacobian);
  _jacobian->zero();

  // Setup PETSc KSP linear solver.
  if (0 == _ksp) {
    PetscErrorCode err = 0;
    err = KSPCreate(_faultMesh->comm(), &_ksp); CHECK_PETSC_ERROR(err);
    err = KSPSetInitialGuessNonzero(_ksp, PETSC_FALSE); CHECK_PETSC_ERROR(err);
    double rtol = 0.0;
    double atol = 0.0;
    double dtol = 0.0;
    int maxIters = 0;
    err = KSPGetTolerances(_ksp, &rtol, &atol, &dtol, &maxIters); 
    CHECK_PETSC_ERROR(err);
    rtol = 1.0e-15;
    atol = 1.0e-25;
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
pylith::faults::FaultCohesiveDyn::_sensitivityUpdateJacobian(const bool negativeSide,
                                                             const topology::Jacobian& jacobian,
                                                             const topology::SolutionFields& fields)
{ // _sensitivityUpdateJacobian
  assert(0 != _quadrature);
  assert(0 != _fields);

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
  const ALE::Obj<SieveMesh::label_sequence>& cellsCohesive =
    sieveMesh->getLabelStratum("material-id", id());
  assert(!cellsCohesive.isNull());
  const SieveMesh::label_sequence::iterator cellsCohesiveBegin =
    cellsCohesive->begin();
  const SieveMesh::label_sequence::iterator cellsCohesiveEnd =
    cellsCohesive->end();

  // Visitor for Jacobian matrix associated with domain.
  const PetscMat jacobianDomainMatrix = jacobian.matrix();
  assert(0 != jacobianDomainMatrix);
  const ALE::Obj<SieveMesh::order_type>& globalOrderDomain =
    sieveMesh->getFactory()->getGlobalOrder(sieveMesh, "default", solutionDomainSection);
  assert(!globalOrderDomain.isNull());
  const ALE::Obj<SieveMesh::sieve_type>& sieve = sieveMesh->getSieve();
  assert(!sieve.isNull());
  ALE::ISieveVisitor::NConeRetriever<SieveMesh::sieve_type> ncV(*sieve,
      (size_t) pow(sieve->getMaxConeSize(), std::max(0, sieveMesh->depth())));
  int_array indicesGlobal(subnrows);

  // Get fault Sieve mesh
  const ALE::Obj<SieveSubMesh>& faultSieveMesh = _faultMesh->sieveMesh();
  assert(!faultSieveMesh.isNull());

  // Get sensitivity solution field
  const ALE::Obj<RealSection>& solutionFaultSection =
    _fields->get("sensitivity solution").section();
  assert(!solutionFaultSection.isNull());

  // Visitor for Jacobian matrix associated with fault.
  double_array jacobianSubCell(submatrixSize);
  assert(0 != _jacobian);
  const PetscMat jacobianFaultMatrix = _jacobian->matrix();
  assert(0 != jacobianFaultMatrix);
  const ALE::Obj<SieveSubMesh::order_type>& globalOrderFault =
    faultSieveMesh->getFactory()->getGlobalOrder(faultSieveMesh, "default", solutionFaultSection);
  assert(!globalOrderFault.isNull());
  // We would need to request unique points here if we had an interpolated mesh
  topology::SubMesh::IndicesVisitor jacobianFaultVisitor(*solutionFaultSection,
                                                 *globalOrderFault,
                           (int) pow(faultSieveMesh->getSieve()->getMaxConeSize(),
                                     faultSieveMesh->depth())*spaceDim);

  const int iCone = (negativeSide) ? 0 : 1;
  for (SieveMesh::label_sequence::iterator c_iter=cellsCohesiveBegin;
       c_iter != cellsCohesiveEnd;
       ++c_iter) {
    // Get cone for cohesive cell
    ncV.clear();
    ALE::ISieveTraversal<SieveMesh::sieve_type>::orientedClosure(*sieve,
								 *c_iter, ncV);
    const int coneSize = ncV.getSize();
    assert(coneSize == 3*numBasis);
    const SieveMesh::point_type *cohesiveCone = ncV.getPoints();
    assert(0 != cohesiveCone);

    const SieveMesh::point_type c_fault = _cohesiveToFault[*c_iter];
    jacobianSubCell = 0.0;

    // Get indices
    for (int iBasis = 0; iBasis < numBasis; ++iBasis) {
      // negative side of the fault: iCone=0
      // positive side of the fault: iCone=1
      const int v_domain = cohesiveCone[iCone*numBasis+iBasis];
      
      for (int iDim=0, iB=iBasis*spaceDim; iDim < spaceDim; ++iDim) {
	if (globalOrderDomain->isLocal(v_domain))
	  indicesGlobal[iB+iDim] = globalOrderDomain->getIndex(v_domain) + iDim;
	else
	  indicesGlobal[iB+iDim] = -1;

	// Set matrix diagonal entries to 1.0 (used when vertex is not local).
	jacobianSubCell[(iB+iDim)*numBasis*spaceDim+iB+iDim] = 1.0;
      } // for
    } // for
    
    PetscErrorCode err = MatGetValues(jacobianDomainMatrix, 
				      indicesGlobal.size(), &indicesGlobal[0],
				      indicesGlobal.size(), &indicesGlobal[0],
				      &jacobianSubCell[0]);
    CHECK_PETSC_ERROR_MSG(err, "Restrict from PETSc Mat failed.");

    // Insert cell contribution into PETSc Matrix
    jacobianFaultVisitor.clear();
    err = updateOperator(jacobianFaultMatrix, *faultSieveMesh->getSieve(),
			 jacobianFaultVisitor, c_fault,
			 &jacobianSubCell[0], INSERT_VALUES);
    CHECK_PETSC_ERROR_MSG(err, "Update to PETSc Mat failed.");
  } // for

  _jacobian->assemble("final_assembly");

  //_jacobian->view(); // DEBUGGING
} // _sensitivityUpdateJacobian

// ----------------------------------------------------------------------
// Reform residual for sensitivity problem.
void
pylith::faults::FaultCohesiveDyn::_sensitivityReformResidual(const bool negativeSide)
{ // _sensitivityReformResidual
  assert(0 != _fields);
  assert(0 != _quadrature);

  const int spaceDim = _quadrature->spaceDim();

  // Compute residual -C^T dLagrange
  double_array residualVertex(spaceDim);
  topology::Field<topology::SubMesh>& residual =
      _fields->get("sensitivity residual");
  const ALE::Obj<RealSection>& residualSection = residual.section();
  residual.zero();

  const ALE::Obj<RealSection>& dLagrangeSection =
      _fields->get("sensitivity dLagrange").section();

  const ALE::Obj<RealSection>& orientationSection =
      _fields->get("orientation").section();
  assert(!orientationSection.isNull());

  const double sign = (negativeSide) ? -1.0 : 1.0;

  const int numVertices = _cohesiveVertices.size();
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_fault = _cohesiveVertices[iVertex].fault;

    const double* dLagrangeVertex = dLagrangeSection->restrictPoint(v_fault);
    assert(0 != dLagrangeVertex);
    assert(spaceDim == dLagrangeSection->getFiberDimension(v_fault));

    const double* orientationVertex = orientationSection->restrictPoint(v_fault);
    assert(0 != orientationVertex);
    assert(spaceDim*spaceDim == orientationSection->getFiberDimension(v_fault));

    residualVertex = 0.0;
    for (int iDim = 0; iDim < spaceDim; ++iDim)
      for (int kDim = 0; kDim < spaceDim; ++kDim)
        residualVertex[iDim] +=
          sign * dLagrangeVertex[kDim] * -orientationVertex[kDim*spaceDim+iDim];

    assert(residualVertex.size() == residualSection->getFiberDimension(v_fault));
    residualSection->updatePoint(v_fault, &residualVertex[0]);
  } // for

  PetscLogFlops(numVertices*spaceDim*spaceDim*4);
} // _sensitivityReformResidual

// ----------------------------------------------------------------------
// Solve sensitivity problem.
void
pylith::faults::FaultCohesiveDyn::_sensitivitySolve(void)
{ // _sensitivitySolve
  assert(0 != _fields);
  assert(0 != _jacobian);
  assert(0 != _ksp);

  const topology::Field<topology::SubMesh>& residual =
      _fields->get("sensitivity residual");
  const topology::Field<topology::SubMesh>& solution =
      _fields->get("sensitivity solution");

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

  //solution.view("SENSITIVITY SOLUTION"); // DEBUGGING
} // _sensitivitySolve

// ----------------------------------------------------------------------
// Update the relative displacement field values based on the
// sensitivity solve.
void
pylith::faults::FaultCohesiveDyn::_sensitivityUpdateSoln(const bool negativeSide)
{ // _sensitivityUpdateSoln
  assert(0 != _fields);
  assert(0 != _quadrature);

  const int spaceDim = _quadrature->spaceDim();

  double_array dispVertex(spaceDim);
  const ALE::Obj<RealSection>& solutionSection =
      _fields->get("sensitivity solution").section();
  const ALE::Obj<RealSection>& dispRelSection =
    _fields->get("sensitivity dispRel").section();

  const double sign = (negativeSide) ? -1.0 : 1.0;

  const int numVertices = _cohesiveVertices.size();
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_fault = _cohesiveVertices[iVertex].fault;

    solutionSection->restrictPoint(v_fault, &dispVertex[0], dispVertex.size());

    dispVertex *= sign;

    assert(dispVertex.size() == dispRelSection->getFiberDimension(v_fault));
    dispRelSection->updateAddPoint(v_fault, &dispVertex[0]);
  } // for
} // _sensitivityUpdateSoln

// ----------------------------------------------------------------------
// Solve slip/Lagrange multiplier sensitivity problem for case of lumped Jacobian in 1-D.
void
pylith::faults::FaultCohesiveDyn::_sensitivitySolveLumped1D(
                                     double_array* slip,
				     const double_array& dLagrangeTpdt,
				     const double_array& jacobianN,
				     const double_array& jacobianP)
{ // _sensitivitySolveLumped1D
  assert(0 != slip);

  // Sensitivity of slip to changes in the Lagrange multipliers
  const double Spp = 1.0 / jacobianN[0] + 1.0
    / jacobianP[0];

  const double dlp = dLagrangeTpdt[0];
  (*slip)[0] -= Spp * dlp;

  PetscLogFlops(2);
} // _sensitivitySolveLumped1D

// ----------------------------------------------------------------------
// Solve slip/Lagrange multiplier sensitivity problem for case of lumped Jacobian in 2-D.
void
pylith::faults::FaultCohesiveDyn::_sensitivitySolveLumped2D(
                                     double_array* slip,
				     const double_array& dLagrangeTpdt,
				     const double_array& jacobianN,
				     const double_array& jacobianP)
{ // _sensitivitySolveLumped2D
  assert(0 != slip);

  // Sensitivity of slip to changes in the Lagrange multipliers
  // Spp = Sqq = 1.0/Aix + 1.0/Ajx
  assert(jacobianN[0] > 0.0);
  assert(jacobianP[0] > 0.0);

  // Check to make sure Jacobian is same at all DOF for
  // vertices i and j (means S is diagonal with equal enties).
  assert(jacobianN[0] == jacobianN[1]);
  assert(jacobianP[0] == jacobianP[1]);
  
  const double Spp = 1.0 / jacobianN[0] + 1.0 / jacobianP[0];
  const double Sqq = Spp;

  const double dlp = dLagrangeTpdt[0];
  const double dlq = dLagrangeTpdt[1];
  (*slip)[0] -= Spp * dlp;
  (*slip)[1] -= Sqq * dlq;

  PetscLogFlops(7);
} // _sensitivitySolveLumped2D

// ----------------------------------------------------------------------
// Solve slip/Lagrange multiplier sensitivity problem for case of lumped Jacobian in 3-D.
void
pylith::faults::FaultCohesiveDyn::_sensitivitySolveLumped3D(
                                     double_array* slip,
				     const double_array& dLagrangeTpdt,
				     const double_array& jacobianN,
				     const double_array& jacobianP)
{ // _sensitivitySolveLumped3D
  assert(0 != slip);

  // Sensitivity of slip to changes in the Lagrange multipliers
  // Aixjx = 1.0/Aix + 1.0/Ajx
  assert(jacobianN[0] > 0.0);
  assert(jacobianP[0] > 0.0);

  // Check to make sure Jacobian is same at all DOF for
  // vertices i and j (means S is diagonal with equal enties).
  assert(jacobianN[0] == jacobianN[1] && 
	 jacobianN[0] == jacobianN[2]);
  assert(jacobianP[0] == jacobianP[1] && 
	 jacobianP[0] == jacobianP[2]);


  const double Spp = 1.0 / jacobianN[0] + 1.0 / jacobianP[0];
  const double Sqq = Spp;
  const double Srr = Spp;

  const double dlp = dLagrangeTpdt[0];
  const double dlq = dLagrangeTpdt[1];
  const double dlr = dLagrangeTpdt[2];
  (*slip)[0] -= Spp * dlp;
  (*slip)[1] -= Sqq * dlq;
  (*slip)[2] -= Srr * dlr;

  PetscLogFlops(9);
} // _sensitivitySolveLumped3D

// ----------------------------------------------------------------------
// Constrain solution space with lumped Jacobian in 1-D.
void
pylith::faults::FaultCohesiveDyn::_constrainSolnSpace1D(double_array* dLagrangeTpdt,
         const double_array& slip,
         const double_array& sliprate,
         const double_array& tractionTpdt,
         const double area)
{ // _constrainSolnSpace1D
  assert(0 != dLagrangeTpdt);

    if (tractionTpdt[0] < 0) {
      // if compression, then no changes to solution
    } else {
      // if tension, then traction is zero.

      const double dlp = -tractionTpdt[0] * area;
      (*dLagrangeTpdt)[0] = dlp;
    } // else

    PetscLogFlops(2);
} // _constrainSolnSpace1D

// ----------------------------------------------------------------------
// Constrain solution space with lumped Jacobian in 2-D.
void
pylith::faults::FaultCohesiveDyn::_constrainSolnSpace2D(double_array* dLagrangeTpdt,
         const double_array& slip,
         const double_array& slipRate,
         const double_array& tractionTpdt,
         const double area)
{ // _constrainSolnSpace2D
  assert(0 != dLagrangeTpdt);

  const double slipMag = fabs(slip[0]);
  const double slipRateMag = fabs(slipRate[0]);

  const double tractionNormal = tractionTpdt[1];
  const double tractionShearMag = fabs(tractionTpdt[0]);

  if (tractionNormal < 0 && 0.0 == slip[1]) {
    // if in compression and no opening
    const double frictionStress = _friction->calcFriction(slipMag, slipRateMag,
                tractionNormal);
    if (tractionShearMag > frictionStress) {
      // traction is limited by friction, so have sliding
      
      // Update slip based on value required to stick versus friction
      const double dlp = -(tractionShearMag - frictionStress) * area *
  tractionTpdt[0] / tractionShearMag;
      (*dLagrangeTpdt)[0] = dlp;
      (*dLagrangeTpdt)[1] = 0.0;
    } else {
      // else friction exceeds value necessary, so stick
      // no changes to solution
    } // if/else
  } else {
    // if in tension, then traction is zero.
    (*dLagrangeTpdt)[0] = -tractionTpdt[0] * area;
    (*dLagrangeTpdt)[1] = -tractionTpdt[1] * area;
  } // else

  PetscLogFlops(8);
} // _constrainSolnSpace2D

// ----------------------------------------------------------------------
// Constrain solution space with lumped Jacobian in 3-D.
void
pylith::faults::FaultCohesiveDyn::_constrainSolnSpace3D(double_array* dLagrangeTpdt,
         const double_array& slip,
         const double_array& slipRate,
         const double_array& tractionTpdt,
         const double area)
{ // _constrainSolnSpace3D
  assert(0 != dLagrangeTpdt);

  const double slipShearMag = sqrt(slip[0] * slip[0] +
             slip[1] * slip[1]);
  double slipRateMag = sqrt(slipRate[0]*slipRate[0] + 
            slipRate[1]*slipRate[1]);
  
  const double tractionNormal = tractionTpdt[2];
  const double tractionShearMag = 
    sqrt(tractionTpdt[0] * tractionTpdt[0] +
	 tractionTpdt[1] * tractionTpdt[1]);
  
  if (tractionNormal < 0.0 && 0.0 == slip[2]) {
    // if in compression and no opening
    const double frictionStress = 
      _friction->calcFriction(slipShearMag, slipRateMag, tractionNormal);
    if (tractionShearMag > frictionStress) {
      // traction is limited by friction, so have sliding
      // Update slip based on value required to stick versus friction
      const double dlp = -(tractionShearMag - frictionStress) * area *
  tractionTpdt[0] / tractionShearMag;
      const double dlq = -(tractionShearMag - frictionStress) * area *
  tractionTpdt[1] / tractionShearMag;

      (*dLagrangeTpdt)[0] = dlp;
      (*dLagrangeTpdt)[1] = dlq;
      (*dLagrangeTpdt)[2] = 0.0;
      
    } else {
      // else friction exceeds value necessary, so stick
      // no changes to solution
    } // if/else
  } else {
    // if in tension, then traction is zero.
    (*dLagrangeTpdt)[0] = -tractionTpdt[0] * area;
    (*dLagrangeTpdt)[1] = -tractionTpdt[1] * area;
    (*dLagrangeTpdt)[2] = -tractionTpdt[2] * area;
  } // else

  PetscLogFlops(22);
} // _constrainSolnSpace3D


// End of file 
