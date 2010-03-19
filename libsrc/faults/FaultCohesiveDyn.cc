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
  _friction(0)
{ // constructor
  _needJacobianDiag = true;
  _needVelocity = true;
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
					      const double upDir[3],
					      const double normalDir[3])
{ // initialize
  assert(0 != upDir);
  assert(0 != normalDir);
  assert(0 != _quadrature);
  assert(0 != _normalizer);

  FaultCohesiveLagrange::initialize(mesh, upDir, normalDir);

  // Get initial tractions using a spatial database.
  _setupInitialTractions();

  // Setup fault constitutive model.
  assert(0 != _friction);
  _friction->initialize(*_faultMesh, _quadrature, _fields->get("area"));

  const spatialdata::geocoords::CoordSys* cs = mesh.coordsys();
  assert(0 != cs);

  topology::Field<topology::SubMesh>& slip = _fields->get("slip");

  // Create field for diagonal entries of Jacobian at conventional
  // vertices i and j associated with Lagrange vertex k
  _fields->add("Jacobian diagonal", "jacobian_diagonal");
  topology::Field<topology::SubMesh>& jacobianDiag = _fields->get(
    "Jacobian diagonal");
  jacobianDiag.newSection(slip, 2 * cs->spaceDim());
  jacobianDiag.allocate();
  jacobianDiag.vectorFieldType(topology::FieldBase::OTHER);

  // Create field for slip rate associated with Lagrange vertex k
  _fields->add("slip rate", "slip_rate");
  topology::Field<topology::SubMesh>& slipRate = _fields->get("slip rate");
  slipRate.cloneSection(slip);
  slipRate.vectorFieldType(topology::FieldBase::OTHER);

  //logger.stagePop();
} // initialize

// ----------------------------------------------------------------------
// Integrate contributions to residual term (r) for operator.
void
pylith::faults::FaultCohesiveDyn::integrateResidualAssembled(
			     const topology::Field<topology::Mesh>& residual,
			     const double t,
			     topology::SolutionFields* const fields)
{ // integrateResidualAssembled
  assert(0 != fields);
  assert(0 != _fields);

  FaultCohesiveLagrange::integrateResidualAssembled(residual, t, fields);

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

  const int numVertices = _cohesiveVertices.size();
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_fault = _cohesiveVertices[iVertex].fault;
    const int v_lagrange = _cohesiveVertices[iVertex].lagrange;
    const int v_negative = _cohesiveVertices[iVertex].negative;
    const int v_positive = _cohesiveVertices[iVertex].positive;

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
} // integrateResidualAssembled

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
      _friction->updateStateVars(slipMag, slipRateMag, tractionNormal);
      break;
    } // case 1
    case 2: { // case 2
      const double slipMag = fabs(slipVertex[0]);
      const double slipRateMag = fabs(slipRateVertex[0]);
      const double tractionNormal = tractionTpdtVertex[1];
      _friction->updateStateVars(slipMag, slipRateMag, tractionNormal);
      break;
    } // case 2
    case 3: { // case 3
      const double slipMag = 
	sqrt(slipVertex[0]*slipVertex[0] + slipVertex[1]*slipVertex[1]);
      const double slipRateMag = 
	sqrt(slipRateVertex[0]*slipRateVertex[0] + 
	     slipRateVertex[1]*slipRateVertex[1]);
      const double tractionNormal = tractionTpdtVertex[2];
      _friction->updateStateVars(slipMag, slipRateMag, tractionNormal);
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
  /// Member prototype for _constrainSolnSpaceLumpedXD()
  typedef void (pylith::faults::FaultCohesiveDyn::*constrainSolnSpace_fn_type)
    (double_array*, double_array*, const double_array&,
     const double_array&, const double_array&, 
     const double_array&, const double_array&, 
     const double);

  assert(0 != fields);
  assert(0 != _quadrature);
  assert(0 != _fields);
  assert(0 != _friction);

  _updateSlipRate(*fields);

  const int spaceDim = _quadrature->spaceDim();

  // Allocate arrays for vertex values
  double_array tractionTVertex(spaceDim);
  double_array tractionTpdtVertex(spaceDim);
  double_array slipTpdtVertex(spaceDim);
  double_array lagrangeTpdtVertex(spaceDim);
  double_array dLagrangeTpdtVertex(spaceDim);

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

  double_array orientationVertex(spaceDim * spaceDim);
  const ALE::Obj<RealSection>& orientationSection =
      _fields->get("orientation").section();
  assert(!orientationSection.isNull());

  double_array lagrangeTVertex(spaceDim);
  const ALE::Obj<RealSection>& dispTSection = fields->get("disp(t)").section();
  assert(!dispTSection.isNull());

  double_array lagrangeTIncrVertex(spaceDim);
  const ALE::Obj<RealSection>& dispTIncrSection =
      fields->get("dispIncr(t->t+dt)").section();
  assert(!dispTIncrSection.isNull());

  double_array jacobianVertexN(spaceDim);
  double_array jacobianVertexP(spaceDim);
  const ALE::Obj<RealSection>& jacobianDiagSection =
        fields->get("Jacobian diagonal").section();
  assert(!jacobianDiagSection.isNull());

  constrainSolnSpace_fn_type constrainSolnSpaceFn;
  switch (spaceDim) { // switch
  case 1:
    constrainSolnSpaceFn = 
      &pylith::faults::FaultCohesiveDyn::_constrainSolnSpaceLumped1D;
    break;
  case 2: 
    constrainSolnSpaceFn = 
      &pylith::faults::FaultCohesiveDyn::_constrainSolnSpaceLumped2D;
    break;
  case 3:
    constrainSolnSpaceFn = 
      &pylith::faults::FaultCohesiveDyn::_constrainSolnSpaceLumped3D;
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

    // Get fault orientation
    orientationSection->restrictPoint(v_fault, &orientationVertex[0],
    orientationVertex.size());

    // Get diagonal of Jacobian at conventional vertices i and j
    // associated with Lagrange vertex k
    jacobianDiagSection->restrictPoint(v_negative, &jacobianVertexN[0],
      jacobianVertexN.size());
    jacobianDiagSection->restrictPoint(v_positive, &jacobianVertexP[0],
      jacobianVertexP.size());

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
					 &slipVertex, slipRateVertex, 
					 tractionTpdtVertex, orientationVertex,
					 jacobianVertexN, jacobianVertexP,
					 *areaVertex);

    // Update Lagrange multiplier values.
    lagrangeTIncrVertex += dLagrangeTpdtVertex;
    assert(lagrangeTIncrVertex.size() ==
        dispTIncrSection->getFiberDimension(v_lagrange));
    dispTIncrSection->updatePoint(v_lagrange, &lagrangeTIncrVertex[0]);

    // Update the slip estimate based on adjustment to the Lagrange
    // multiplier values.
    assert(slipVertex.size() ==
        slipSection->getFiberDimension(v_fault));
    slipSection->updatePoint(v_fault, &slipVertex[0]);
  } // if

#if 0 // DEBUGGING
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
  /// Member prototype for _constrainSolnSpaceLumpedXD()
  typedef void (pylith::faults::FaultCohesiveDyn::*constrainSolnSpace_fn_type)
    (double_array*, double_array*, const double_array&,
     const double_array&, const double_array&, 
     const double_array&, const double_array&, 
     const double);

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

  double_array jacobianVertexN(spaceDim);
  double_array jacobianVertexP(spaceDim);
  const ALE::Obj<RealSection>& jacobianSection = jacobian.section();
  assert(!jacobianSection.isNull());

  double_array residualVertexN(spaceDim);
  double_array residualVertexP(spaceDim);
  const ALE::Obj<RealSection>& residualSection =
      fields->get("residual").section();

  adjustSolnLumped_fn_type adjustSolnLumpedFn;
  constrainSolnSpace_fn_type constrainSolnSpaceFn;
  switch (spaceDim) { // switch
  case 1:
    adjustSolnLumpedFn = 
      &pylith::faults::FaultCohesiveDyn::_adjustSolnLumped1D;
    constrainSolnSpaceFn = 
      &pylith::faults::FaultCohesiveDyn::_constrainSolnSpaceLumped1D;
    break;
  case 2: 
    adjustSolnLumpedFn = 
      &pylith::faults::FaultCohesiveDyn::_adjustSolnLumped2D;
    constrainSolnSpaceFn = 
      &pylith::faults::FaultCohesiveDyn::_constrainSolnSpaceLumped2D;
    break;
  case 3:
    adjustSolnLumpedFn = 
      &pylith::faults::FaultCohesiveDyn::_adjustSolnLumped3D;
    constrainSolnSpaceFn = 
      &pylith::faults::FaultCohesiveDyn::_constrainSolnSpaceLumped3D;
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
					 &slipVertex, slipRateVertex,
					 tractionTpdtVertex, orientationVertex,
					 jacobianVertexN, jacobianVertexP,
					 *areaVertex);

    lagrangeTIncrVertex += dLagrangeTpdtVertex;

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

  } else if (0 == strncasecmp("initial_traction", name, slipStrLen)) {
    assert(0 != _dbInitialTract);
    topology::Field<topology::SubMesh>& buffer =
        _fields->get("buffer (vector)");
    _calcInitialTractions(&buffer);
    return buffer;

  } else if (0 == strcasecmp("traction", name)) {
    assert(0 != fields);
    const topology::Field<topology::Mesh>& dispT = fields->get("disp(t)");
    _allocateBufferVectorField();
    topology::Field<topology::SubMesh>& buffer =
        _fields->get("buffer (vector)");
    _calcTractions(&buffer, dispT);
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

    for (int i = 0; i < fiberDim; ++i)
      tractionsVertex[i] = dispTVertex[i] / areaVertex[0];

    assert(tractionsVertex.size() == tractionsSection->getFiberDimension(v_fault));
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
pylith::faults::FaultCohesiveDyn::_calcInitialTractions(
    topology::Field<topology::SubMesh>* tractions)
{ // _calcInitialTractions
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

} // _calcInitialTractions

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
          velocityVertexN[kDim] * -orientationVertex[kDim*spaceDim+iDim];

    // Velocity for positive vertex.
    for (int iDim = 0; iDim < spaceDim; ++iDim)
      for (int kDim = 0; kDim < spaceDim; ++kDim)
        slipRateVertex[iDim] +=
          velocityVertexP[kDim] * +orientationVertex[kDim*spaceDim+iDim];

    // Update slip rate field.
    assert(slipRateVertex.size() == slipRateSection->getFiberDimension(v_fault));
    slipRateSection->updatePoint(v_fault, &slipRateVertex[0]);
  } // for

  PetscLogFlops(numVertices*spaceDim*spaceDim*4);
} // _updateSlipRate

// ----------------------------------------------------------------------
// Constrain solution space in 1-D.
void
pylith::faults::FaultCohesiveDyn::_constrainSolnSpaceLumped1D(
				     double_array* dLagrangeTpdt,
				     double_array* slip,
				     const double_array& sliprate,
				     const double_array& tractionTpdt,
				     const double_array& orientation,
				     const double_array& jacobianN,
				     const double_array& jacobianP,
				     const double area)
{ // constrainSolnSpace1D
  assert(0 != dLagrangeTpdt);
  assert(0 != slip);

  // Sensitivity of slip to changes in the Lagrange multipliers
  // Aixjx = 1.0/Aix + 1.0/Ajx
  const double Aixjx = 1.0 / jacobianN[0] + 1.0
    / jacobianP[0];
  const double Spp = 1.0;
  
  if (tractionTpdt[0] < 0) {
    // if compression, then no changes to solution
  } else {
    // if tension, then traction is zero.
    
    // Update slip based on value required to stick versus
    // zero traction
    const double dlp = tractionTpdt[0] * area;
    (*dLagrangeTpdt)[0] = -dlp;
    (*slip)[0] += Spp * dlp;    
  } // else

  PetscLogFlops(6);
} // constrainSolnSpace1D

// ----------------------------------------------------------------------
// Constrain solution space in 2-D.
void
pylith::faults::FaultCohesiveDyn::_constrainSolnSpaceLumped2D(
				     double_array* dLagrangeTpdt,
				     double_array* slip,
				     const double_array& slipRate,
				     const double_array& tractionTpdt,
				     const double_array& orientation,
				     const double_array& jacobianN,
				     const double_array& jacobianP,
				     const double area)
{ // constrainSolnSpace2D
  assert(0 != dLagrangeTpdt);
  assert(0 != slip);

  std::cout << "Normal traction:" << tractionTpdt[1] << std::endl;

  // Sensitivity of slip to changes in the Lagrange multipliers
  // Aixjx = 1.0/Aix + 1.0/Ajx
  assert(jacobianN[0] > 0.0);
  assert(jacobianP[0] > 0.0);
  const double Aixjx = 1.0 / jacobianN[0] + 1.0
    / jacobianP[0];
  // Aiyjy = 1.0/Aiy + 1.0/Ajy
  assert(jacobianN[1] > 0.0);
  assert(jacobianP[1] > 0.0);
  const double Aiyjy = 1.0 / jacobianN[1] + 1.0
    / jacobianP[1];
  const double Cpx = orientation[0];
  const double Cpy = orientation[1];
  const double Cqx = orientation[2];
  const double Cqy = orientation[3];

  // Sensitivity matrix using only diagonal of A to approximate its inverse
   const double Spp = Cpx * Cpx * Aixjx + Cpy * Cpy * Aiyjy;
   const double Spq = Cpx * Cqx * Aixjx + Cpy * Cqy * Aiyjy;
   const double Sqq = Cqx * Cqx * Aixjx + Cqy * Cqy * Aiyjy;
  
  const double slipMag = fabs((*slip)[0]);
  const double slipRateMag = fabs(slipRate[0]);

  const double tractionNormal = tractionTpdt[1];
  const double tractionShearMag = fabs(tractionTpdt[0]);

  if (tractionNormal < 0 && 0.0 == (*slip)[1]) {
    // if in compression and no opening
    std::cout << "FAULT IN COMPRESSION" << std::endl;
    const double frictionStress = _friction->calcFriction(slipMag, slipRateMag,
							  tractionNormal);
    std::cout << "frictionStress: " << frictionStress << std::endl;
    if (tractionShearMag > frictionStress || 
	(tractionShearMag < frictionStress && slipMag > 0.0)) {
      // traction is limited by friction, so have sliding
      std::cout << "LIMIT TRACTION, HAVE SLIDING" << std::endl;
      
      // Update slip based on value required to stick versus friction
      const double dlp = (tractionShearMag - frictionStress) * area *
	tractionTpdt[0] / tractionShearMag;
      (*dLagrangeTpdt)[0] = -dlp;
      (*slip)[0] += Spp * dlp;
      std::cout << "Estimated slip: " << (*slip)[0] << std::endl;
    } else {
      // else friction exceeds value necessary, so stick
      std::cout << "STICK" << std::endl;
      // no changes to solution
    } // if/else
  } else {
    // if in tension, then traction is zero.
    std::cout << "FAULT IN TENSION" << std::endl;
    
    // Update slip based on value required to stick versus
    // zero traction
    const double dlp = tractionTpdt[0] * area;
    const double dlq = tractionTpdt[1] * area;

    (*dLagrangeTpdt)[0] = -dlp;
    (*dLagrangeTpdt)[1] = -dlq;
    (*slip)[0] += Spp * dlp + Spq * dlq;
    (*slip)[1] += Spq * dlp + Sqq * dlq;
  } // else

  PetscLogFlops(0); // :TODO: Fix this
} // constrainSolnSpace2D

// ----------------------------------------------------------------------
// Constrain solution space in 3-D.
void
pylith::faults::FaultCohesiveDyn::_constrainSolnSpaceLumped3D(
				     double_array* dLagrangeTpdt,
				     double_array* slip,
				     const double_array& slipRate,
				     const double_array& tractionTpdt,
				     const double_array& orientation,
				     const double_array& jacobianN,
				     const double_array& jacobianP,
				     const double area)
{ // constrainSolnSpace3D
  assert(0 != dLagrangeTpdt);
  assert(0 != slip);

  std::cout << "Normal traction:" << tractionTpdt[2] << std::endl;

  // Sensitivity of slip to changes in the Lagrange multipliers
  // Aixjx = 1.0/Aix + 1.0/Ajx
  assert(jacobianN[0] > 0.0);
  assert(jacobianP[0] > 0.0);
  const double Aixjx = 1.0 / jacobianN[0] + 1.0 / jacobianP[0];
  // Aiyjy = 1.0/Aiy + 1.0/Ajy
  assert(jacobianN[1] > 0.0);
  assert(jacobianP[1] > 0.0);
  const double Aiyjy = 1.0 / jacobianN[1] + 1.0 / jacobianP[1];
  // Aizjz = 1.0/Aiz + 1.0/Ajz
  assert(jacobianN[2] > 0.0);
  assert(jacobianP[2] > 0.0);
  const double Aizjz = 1.0 / jacobianN[2] + 1.0 / jacobianP[2];
  const double Cpx = orientation[0];
  const double Cpy = orientation[1];
  const double Cpz = orientation[2];
  const double Cqx = orientation[3];
  const double Cqy = orientation[4];
  const double Cqz = orientation[5];
  const double Crx = orientation[6];
  const double Cry = orientation[7];
  const double Crz = orientation[8];

  // Sensitivity matrix using only diagonal of A to approximate its inverse
    const double Spp = Cpx * Cpx * Aixjx + Cpy * Cpy * Aiyjy + Cpz * Cpz * Aizjz;
    const double Spq = Cpx * Cqx * Aixjx + Cpy * Cqy * Aiyjy + Cpz * Cqz * Aizjz;
    const double Spr = Cpx * Crx * Aixjx + Cpy * Cry * Aiyjy + Cpz * Crz * Aizjz;
    const double Sqq = Cqx * Cqx * Aixjx + Cqy * Cqy * Aiyjy + Cqz * Cqz * Aizjz;
    const double Sqr = Cqx * Crx * Aixjx + Cqy * Cry * Aiyjy + Cqz * Crz * Aizjz;
    const double Srr = Crx * Crx * Aixjx + Cry * Cry * Aiyjy + Crz * Crz * Aizjz;
  
  const double slipShearMag = sqrt((*slip)[0] * (*slip)[0] + 
				   (*slip)[1] * (*slip)[1]);
  double slipRateMag = sqrt(slipRate[0]*slipRate[0] + 
			    slipRate[1]*slipRate[1]);
  
  const double tractionNormal = tractionTpdt[2];
  const double tractionShearMag = 
    sqrt(tractionTpdt[0] * tractionTpdt[0] +
	 tractionTpdt[1] * tractionTpdt[1]);
  
  if (tractionNormal < 0.0 && 0.0 == (*slip)[2]) {
    // if in compression and no opening
    std::cout << "FAULT IN COMPRESSION" << std::endl;
    const double frictionStress = 
      _friction->calcFriction(slipShearMag, slipRateMag, tractionNormal);
    std::cout << "frictionStress: " << frictionStress << std::endl;
    if (tractionShearMag > frictionStress || 
	(tractionShearMag < frictionStress && slipShearMag > 0.0)) {
      // traction is limited by friction, so have sliding
      std::cout << "LIMIT TRACTION, HAVE SLIDING" << std::endl;
      
      // Update slip based on value required to stick versus friction
      const double dlp = (tractionShearMag - frictionStress) * area *
	tractionTpdt[0] / tractionShearMag;
      const double dlq = (tractionShearMag - frictionStress) * area *
	tractionTpdt[1] / tractionShearMag;

      (*dLagrangeTpdt)[0] = -dlp;
      (*dLagrangeTpdt)[1] = -dlq;
      (*slip)[0] += Spp * dlp + Spq * dlq;
      (*slip)[1] += Spq * dlp + Sqq * dlq;
      
      std::cout << "Estimated slip: " << "  " << (*slip)[0] << "  "
		<< (*slip)[1] << "  " << (*slip)[2] << std::endl;
    } else {
      // else friction exceeds value necessary, so stick
      std::cout << "STICK" << std::endl;
      // no changes to solution
    } // if/else
  } else {
    // if in tension, then traction is zero.
    std::cout << "FAULT IN TENSION" << std::endl;
    
    // Update slip based on value required to stick versus
    // zero traction
    const double dlp = tractionTpdt[0] * area;
    const double dlq = tractionTpdt[1] * area;
    const double dlr = tractionTpdt[2] * area;

    (*dLagrangeTpdt)[0] = -dlp;
    (*dLagrangeTpdt)[1] = -dlq;
    (*dLagrangeTpdt)[2] = -dlr;
    (*slip)[0] += Spp * dlp + Spq * dlq + Spr * dlr;
    (*slip)[1] += Spq * dlp + Sqq * dlq + Sqr * dlr;
    (*slip)[2] += Spr * dlp + Sqr * dlq + Srr * dlr;
    
    std::cout << "Estimated slip: " << "  " << (*slip)[0] << "  "
	      << (*slip)[1] << "  " << (*slip)[2] << std::endl;
    
  } // else

  PetscLogFlops(0); // :TODO: Fix this
} // constrainSolnSpace3D

// ----------------------------------------------------------------------
// Constrain solution space in 2-D.
void
pylith::faults::FaultCohesiveDyn::_constrainSolnSpaceLumped2x22D(
				     double_array* dLagrangeTpdt,
				     double_array* slip,
				     const double_array& slipRate,
				     const double_array& tractionTpdt,
				     const double_array& orientation,
				     const double_array& jacobianN,
				     const double_array& jacobianP,
				     const double area)
{ // constrainSolnSpace2x22D
 assert(0 != dLagrangeTpdt);
 assert(0 != slip);

 std::cout << "Normal traction:" << tractionTpdt[1] << std::endl;

 // Sensitivity of slip to changes in the Lagrange multipliers
 assert(jacobianN[0] > 0.0);
 assert(jacobianP[0] > 0.0);
 assert(jacobianN[1] > 0.0);
 assert(jacobianP[1] > 0.0);
 assert(jacobianN[2] > 0.0);
 assert(jacobianP[2] > 0.0);
 assert(jacobianN[3] > 0.0);
 assert(jacobianP[3] > 0.0);

 // Fault orientation matrix
 const double Cpx = orientation[0];
 const double Cpy = orientation[1];
 const double Cqx = orientation[2];
 const double Cqy = orientation[3];

 // 3 x 3 Jacobian block of i side of fault
 const double Aixx = jacobianN[0];
 const double Aixy = jacobianN[1];
 const double Aiyx = jacobianN[2];
 const double Aiyy = jacobianN[3];

 // 3 x 3 Jacobian block of j side of fault
 const double Ajxx = jacobianP[0];
 const double Ajxy = jacobianP[1];
 const double Ajyx = jacobianP[2];
 const double Ajyy = jacobianP[3];

 // Determinant of 2 X 2 block of Jacobian on each side of the fault
 const double Deti = Aixx * Aiyy - Aixy*Aiyx;
 const double Detj = Ajxx * Ajyy - Ajxy*Ajyx;
 assert(Deti > 0.0);
 assert(Detj > 0.0);

 // Co-factor matrix for i side of fault
 const double Ci11 = Aiyy;
 const double Ci12 = -Aiyx;
 const double Ci21 = -Aixy;
 const double Ci22 = Aixx;

 // Co-factor matrix for j side of fault
 const double Cj11 = Ajyy;
 const double Cj12 = -Ajyx;
 const double Cj21 = -Ajxy;
 const double Cj22 = Ajxx;

 // Contribution to sensitivity from i side of fault
 const double Sipp = Ci11 * Cpx * Cpx + Ci22 * Cpy * Cpy + (Ci12 + Ci21) * Cpx * Cpy;
 const double Sipq = (Ci11 * Cpx + Ci12 * Cpy) * Cqx + (Ci21 * Cpx + Ci22 * Cpy) * Cqy;
 const double Siqp = (Ci11 * Cqx + Ci12 * Cqy) * Cpx + (Ci21 * Cqx + Ci22 * Cqy) * Cpy;
 const double Siqq = Ci11 * Cqx * Cqx + Ci22 * Cqy * Cqy + (Ci12 + Ci21) * Cqx * Cqy;

 // Contribution to sensitivity from i side of fault
 const double Sjpp = Cj11 * Cpx * Cpx + Cj22 * Cpy * Cpy + (Cj12 + Cj21) * Cpx * Cpy;
 const double Sjpq = (Cj11 * Cpx + Cj12 * Cpy) * Cqx + (Cj21 * Cpx + Cj22 * Cpy) * Cqy;
 const double Sjqp = (Cj11 * Cqx + Cj12 * Cqy) * Cpx + (Cj21 * Cqx + Cj22 * Cqy) * Cpy;
 const double Sjqq = Cj11 * Cqx * Cqx + Cj22 * Cqy * Cqy + (Cj12 + Cj21) * Cqx * Cqy;

 // Sensitivity Matrix
 const double Spp = Sipp/Deti + Sjpp/Detj;
 const double Spq = Sipq/Deti + Sjpq/Detj;
 const double Sqp = Siqp/Deti + Sjqp/Detj;
 const double Sqq = Siqq/Deti + Sjqq/Detj;

 const double slipMag = fabs((*slip)[0]);
 const double slipRateMag = fabs(slipRate[0]);

 const double tractionNormal = tractionTpdt[1];
 const double tractionShearMag = fabs(tractionTpdt[0]);

 if (tractionNormal < 0 && 0.0 == (*slip)[1]) {
   // if in compression and no opening
   std::cout << "FAULT IN COMPRESSION" << std::endl;
   const double frictionStress = _friction->calcFriction(slipMag, slipRateMag,
							  tractionNormal);
   std::cout << "frictionStress: " << frictionStress << std::endl;
   if (tractionShearMag > frictionStress ||
	(tractionShearMag < frictionStress && slipMag > 0.0)) {
     // traction is limited by friction, so have sliding
     std::cout << "LIMIT TRACTION, HAVE SLIDING" << std::endl;

     // Update slip based on value required to stick versus friction
     const double dlp = (tractionShearMag - frictionStress) * area *
	tractionTpdt[0] / tractionShearMag;
     (*dLagrangeTpdt)[0] = -dlp;
     (*slip)[0] += Spp * dlp;
     std::cout << "Estimated slip: " << (*slip)[0] << std::endl;
   } else {
     // else friction exceeds value necessary, so stick
     std::cout << "STICK" << std::endl;
     // no changes to solution
   } // if/else
 } else {
   // if in tension, then traction is zero.
   std::cout << "FAULT IN TENSION" << std::endl;

   // Update slip based on value required to stick versus
   // zero traction
   const double dlp = tractionTpdt[0] * area;
   const double dlq = tractionTpdt[1] * area;

   (*dLagrangeTpdt)[0] = -dlp;
   (*dLagrangeTpdt)[1] = -dlq;
   (*slip)[0] += Spp * dlp + Spq * dlq;
   (*slip)[1] += Spq * dlp + Sqq * dlq;
 } // else

 PetscLogFlops(0); // :TODO: Fix this
} // constrainSolnSpace2x22D

// ----------------------------------------------------------------------
// Constrain solution space in 3-D.
void
pylith::faults::FaultCohesiveDyn::_constrainSolnSpaceLumped3x33D(
				     double_array* dLagrangeTpdt,
				     double_array* slip,
				     const double_array& slipRate,
				     const double_array& tractionTpdt,
				     const double_array& orientation,
				     const double_array& jacobianN,
				     const double_array& jacobianP,
				     const double area)
{ // constrainSolnSpace3x33D
 assert(0 != dLagrangeTpdt);
 assert(0 != slip);

 std::cout << "Normal traction:" << tractionTpdt[2] << std::endl;

 // Sensitivity of slip to changes in the Lagrange multipliers
 assert(jacobianN[0] > 0.0);
 assert(jacobianP[0] > 0.0);
 assert(jacobianN[1] > 0.0);
 assert(jacobianP[1] > 0.0);
 assert(jacobianN[2] > 0.0);
 assert(jacobianP[2] > 0.0);
 assert(jacobianN[3] > 0.0);
 assert(jacobianP[3] > 0.0);
 assert(jacobianN[4] > 0.0);
 assert(jacobianP[4] > 0.0);
 assert(jacobianN[5] > 0.0);
 assert(jacobianP[5] > 0.0);
 assert(jacobianN[6] > 0.0);
 assert(jacobianP[6] > 0.0);
 assert(jacobianN[7] > 0.0);
 assert(jacobianP[7] > 0.0);
 assert(jacobianN[8] > 0.0);
 assert(jacobianP[8] > 0.0);

 // 3 x 3 Jacobian block of i side of fault
 const double Aixx = jacobianN[0];
 const double Aixy = jacobianN[1];
 const double Aixz = jacobianN[2];
 const double Aiyx = jacobianN[3];
 const double Aiyy = jacobianN[4];
 const double Aiyz = jacobianN[5];
 const double Aizx = jacobianN[6];
 const double Aizy = jacobianN[7];
 const double Aizz = jacobianN[8];

 // 3 x 3 Jacobian block of j side of fault
 const double Ajxx = jacobianP[0];
 const double Ajxy = jacobianP[1];
 const double Ajxz = jacobianP[2];
 const double Ajyx = jacobianP[3];
 const double Ajyy = jacobianP[4];
 const double Ajyz = jacobianP[5];
 const double Ajzx = jacobianP[6];
 const double Ajzy = jacobianP[7];
 const double Ajzz = jacobianP[8];

 // Fault orientation matrix
 const double Cpx = orientation[0];
 const double Cpy = orientation[1];
 const double Cpz = orientation[2];
 const double Cqx = orientation[3];
 const double Cqy = orientation[4];
 const double Cqz = orientation[5];
 const double Crx = orientation[6];
 const double Cry = orientation[7];
 const double Crz = orientation[8];

 // Determinant of 3 X 3 block of Jacobian on each side of the fault
 const double Deti = Aixz * (-Aiyy * Aizx + Aiyx * Aizy) +
   Aixy * (Aiyz * Aizx - Aiyx * Aizz) + Aixx * (-Aiyz * Aizy + Aiyy * Aizz);
 const double Detj = Ajxz * (-Ajyy * Ajzx + Ajyx * Ajzy) +
   Ajxy * (Ajyz * Ajzx - Ajyx * Ajzz) + Ajxx * (-Ajyz * Ajzy + Ajyy * Ajzz);
 assert(Deti > 0.0);
 assert(Detj > 0.0);

 // Co-factor matrix for i side of fault
 const double Ci11 = Aiyz * Aizy + Aiyy * Aizz;
 const double Ci12 = Aiyz * Aizx - Aiyx * Aizz;
 const double Ci13 = -Aiyy * Aizx + Aiyx * Aizy;
 const double Ci21 = Aixz * Aizy - Aixy * Aizz;
 const double Ci22 = -Aixz * Aizx + Aixx * Aizz;
 const double Ci23 = Aixy * Aizx - Aixx * Aizy;
 const double Ci31 = Aixz * Aiyy + Aixy * Aiyz;
 const double Ci32 = Aixz * Aiyx - Aixx * Aiyz;
 const double Ci33 = -Aixy * Aiyx + Aixx * Aiyy;

 // Co-factor matrix for j side of fault
 const double Cj11 = Ajyz * Ajzy + Ajyy * Ajzz;
 const double Cj12 = Ajyz * Ajzx - Ajyx * Ajzz;
 const double Cj13 = -Ajyy * Ajzx + Ajyx * Ajzy;
 const double Cj21 = Ajxz * Ajzy - Ajxy * Ajzz;
 const double Cj22 = -Ajxz * Ajzx + Ajxx * Ajzz;
 const double Cj23 = Ajxy * Ajzx - Ajxx * Ajzy;
 const double Cj31 = Ajxz * Ajyy + Ajxy * Ajyz;
 const double Cj32 = Ajxz * Ajyx - Ajxx * Ajyz;
 const double Cj33 = -Ajxy * Ajyx + Ajxx * Ajyy;

 // Contribution to sensitivity from i side of fault
 const double Sipp = Ci11 * Cpx * Cpx + Ci22 * Cpy * Cpy + Ci33 * Cpz * Cpz +
   (Ci12 + Ci21) * Cpx * Cpy + (Ci13 + Ci31) * Cpx * Cpz + (Ci23 + Ci32) * Cpy * Cpz;
 const double Siqq = Ci11 * Cqx * Cqx + Ci22 * Cqy * Cqy + Ci33 * Cqz * Cqz +
   (Ci12 + Ci21) * Cqx * Cqy + (Ci13 + Ci31) * Cqx * Cqz + (Ci23 + Ci32) * Cqy * Cqz;
 const double Sirr = Ci11 * Crx * Crx + Ci22 * Cry * Cry + Ci33 * Crz * Crz +
   (Ci12 + Ci21) * Crx * Cry + (Ci13 + Ci31) * Crx * Crz + (Ci23 + Ci32) * Cry * Crz;
 const double Sipq = (Ci11 * Cpx + Ci12 * Cpy + Ci13 * Cpz) * Cqx +
   (Ci21 * Cpx + Ci22 * Cpy + Ci23 * Cpz) * Cqy +
   (Ci31 * Cpx + Ci32 * Cpy + Ci33 * Cpz) * Cqz;
 const double Sipr = (Ci11 * Cpx + Ci12 * Cpy + Ci13 * Cpz) * Crx +
   (Ci21 * Cpx + Ci22 * Cpy + Ci23 * Cpz) * Cry +
   (Ci31 * Cpx + Ci32 * Cpy + Ci33 * Cpz) * Crz;
 const double Siqp = (Ci11 * Cqx + Ci12 * Cqy + Ci13 * Cqz) * Cpx +
   (Ci21 * Cqx + Ci22 * Cqy + Ci23 * Cqz) * Cpy +
   (Ci31 * Cqx + Ci32 * Cqy + Ci33 * Cqz) * Cpz;
 const double Siqr = (Ci11 * Cqx + Ci12 * Cqy + Ci13 * Cqz) * Crx +
   (Ci21 * Cqx + Ci22 * Cqy + Ci23 * Cqz) * Cry +
   (Ci31 * Cqx + Ci32 * Cqy + Ci33 * Cqz) * Crz;
 const double Sirp = (Ci11 * Crx + Ci12 * Cry + Ci13 * Crz) * Cpx +
   (Ci21 * Crx + Ci22 * Cry + Ci23 * Crz) * Cpy +
   (Ci31 * Crx + Ci32 * Cry + Ci33 * Crz) * Cpz;
 const double Sirq = (Ci11 * Crx + Ci12 * Cry + Ci13 * Crz) * Cqx +
   (Ci21 * Crx + Ci22 * Cry + Ci23 * Crz) * Cqy +
   (Ci31 * Crx + Ci32 * Cry + Ci33 * Crz) * Cqz;

 // Contribution to sensitivity from i side of fault
 const double Sjpp = Cj11 * Cpx * Cpx + Cj22 * Cpy * Cpy + Cj33 * Cpz * Cpz +
   (Cj12 + Cj21) * Cpx * Cpy + (Cj13 + Cj31) * Cpx * Cpz + (Cj23 + Cj32) * Cpy * Cpz;
 const double Sjqq = Cj11 * Cqx * Cqx + Cj22 * Cqy * Cqy + Cj33 * Cqz * Cqz +
   (Cj12 + Cj21) * Cqx * Cqy + (Cj13 + Cj31) * Cqx * Cqz + (Cj23 + Cj32) * Cqy * Cqz;
 const double Sjrr = Cj11 * Crx * Crx + Cj22 * Cry * Cry + Cj33 * Crz * Crz +
   (Cj12 + Cj21) * Crx * Cry + (Cj13 + Cj31) * Crx * Crz + (Cj23 + Cj32) * Cry * Crz;
 const double Sjpq = (Cj11 * Cpx + Cj12 * Cpy + Cj13 * Cpz) * Cqx +
   (Cj21 * Cpx + Cj22 * Cpy + Cj23 * Cpz) * Cqy +
   (Cj31 * Cpx + Cj32 * Cpy + Cj33 * Cpz) * Cqz;
 const double Sjpr = (Cj11 * Cpx + Cj12 * Cpy + Cj13 * Cpz) * Crx +
   (Cj21 * Cpx + Cj22 * Cpy + Cj23 * Cpz) * Cry +
   (Cj31 * Cpx + Cj32 * Cpy + Cj33 * Cpz) * Crz;
 const double Sjqp = (Cj11 * Cqx + Cj12 * Cqy + Cj13 * Cqz) * Cpx +
   (Cj21 * Cqx + Cj22 * Cqy + Cj23 * Cqz) * Cpy +
   (Cj31 * Cqx + Cj32 * Cqy + Cj33 * Cqz) * Cpz;
 const double Sjqr = (Cj11 * Cqx + Cj12 * Cqy + Cj13 * Cqz) * Crx +
   (Cj21 * Cqx + Cj22 * Cqy + Cj23 * Cqz) * Cry +
   (Cj31 * Cqx + Cj32 * Cqy + Cj33 * Cqz) * Crz;
 const double Sjrp = (Cj11 * Crx + Cj12 * Cry + Cj13 * Crz) * Cpx +
   (Cj21 * Crx + Cj22 * Cry + Cj23 * Crz) * Cpy +
   (Cj31 * Crx + Cj32 * Cry + Cj33 * Crz) * Cpz;
 const double Sjrq = (Cj11 * Crx + Cj12 * Cry + Cj13 * Crz) * Cqx +
   (Cj21 * Crx + Cj22 * Cry + Cj23 * Crz) * Cqy +
   (Cj31 * Crx + Cj32 * Cry + Cj33 * Crz) * Cqz;

 // Sensitivity Matrix
 const double Spp = Sipp/Deti + Sjpp/Detj;
 const double Spq = Sipq/Deti + Sjpq/Detj;
 const double Spr = Sipr/Deti + Sjpr/Detj;
 const double Sqp = Siqp/Deti + Sjqp/Detj;
 const double Sqq = Siqq/Deti + Sjqq/Detj;
 const double Sqr = Siqr/Deti + Sjqr/Detj;
 const double Srp = Sirp/Deti + Sjrp/Detj;
 const double Srq = Sirq/Deti + Sjrq/Detj;
 const double Srr = Sirr/Deti + Sjrr/Detj;


 const double slipShearMag = sqrt((*slip)[0] * (*slip)[0] +
				   (*slip)[1] * (*slip)[1]);
 double slipRateMag = sqrt(slipRate[0]*slipRate[0] +
			    slipRate[1]*slipRate[1]);

 const double tractionNormal = tractionTpdt[2];
 const double tractionShearMag =
   sqrt(tractionTpdt[0] * tractionTpdt[0] +
	 tractionTpdt[1] * tractionTpdt[1]);

 if (tractionNormal < 0.0 && 0.0 == (*slip)[2]) {
   // if in compression and no opening
   std::cout << "FAULT IN COMPRESSION" << std::endl;
   const double frictionStress =
     _friction->calcFriction(slipShearMag, slipRateMag, tractionNormal);
   std::cout << "frictionStress: " << frictionStress << std::endl;
   if (tractionShearMag > frictionStress ||
	(tractionShearMag < frictionStress && slipShearMag > 0.0)) {
     // traction is limited by friction, so have sliding
     std::cout << "LIMIT TRACTION, HAVE SLIDING" << std::endl;

     // Update slip based on value required to stick versus friction
     const double dlp = (tractionShearMag - frictionStress) * area *
	tractionTpdt[0] / tractionShearMag;
     const double dlq = (tractionShearMag - frictionStress) * area *
	tractionTpdt[1] / tractionShearMag;

     (*dLagrangeTpdt)[0] = -dlp;
     (*dLagrangeTpdt)[1] = -dlq;
     (*slip)[0] += Spp * dlp + Spq * dlq;
     (*slip)[1] += Spq * dlp + Sqq * dlq;

     std::cout << "Estimated slip: " << "  " << (*slip)[0] << "  "
		<< (*slip)[1] << "  " << (*slip)[2] << std::endl;
   } else {
     // else friction exceeds value necessary, so stick
     std::cout << "STICK" << std::endl;
     // no changes to solution
   } // if/else
 } else {
   // if in tension, then traction is zero.
   std::cout << "FAULT IN TENSION" << std::endl;

   // Update slip based on value required to stick versus
   // zero traction
   const double dlp = tractionTpdt[0] * area;
   const double dlq = tractionTpdt[1] * area;
   const double dlr = tractionTpdt[2] * area;

   (*dLagrangeTpdt)[0] = -dlp;
   (*dLagrangeTpdt)[1] = -dlq;
   (*dLagrangeTpdt)[2] = -dlr;
   (*slip)[0] += Spp * dlp + Spq * dlq + Spr * dlr;
   (*slip)[1] += Spq * dlp + Sqq * dlq + Sqr * dlr;
   (*slip)[2] += Spr * dlp + Sqr * dlq + Srr * dlr;

   std::cout << "Estimated slip: " << "  " << (*slip)[0] << "  "
	      << (*slip)[1] << "  " << (*slip)[2] << std::endl;

 } // else

 PetscLogFlops(0); // :TODO: Fix this
} // constrainSolnSpace3x33D


// End of file 
