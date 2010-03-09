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
  _getInitialTractions();

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

  // No contribution if no initial tractions are specified.
  if (0 == _dbInitialTract)
    return;

  const int spaceDim = _quadrature->spaceDim();

  // Get sections
  double_array forcesInitialVertex(spaceDim);
  const ALE::Obj<RealSection>& forcesInitialSection = 
    _fields->get("initial force").section();
  assert(!forcesInitialSection.isNull());

  double_array residualVertex(spaceDim);
  const ALE::Obj<RealSection>& residualSection = residual.section();
  assert(!residualSection.isNull());

  const int numVertices = _cohesiveVertices.size();
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_lagrange = _cohesiveVertices[iVertex].lagrange;
    const int v_fault = _cohesiveVertices[iVertex].fault;
    const int v_negative = _cohesiveVertices[iVertex].negative;
    const int v_positive = _cohesiveVertices[iVertex].positive;

    // Get initial forces at fault vertex
    forcesInitialSection->restrictPoint(v_fault, 
					&forcesInitialVertex[0],
					forcesInitialVertex.size());

    residualVertex = forcesInitialVertex;
    assert(residualVertex.size() == 
	   residualSection->getFiberDimension(v_positive));
    residualSection->updatePoint(v_positive, &residualVertex[0]);

    residualVertex *= -1.0;
    assert(residualVertex.size() == 
	   residualSection->getFiberDimension(v_negative));
    residualSection->updatePoint(v_negative, &residualVertex[0]);
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
      const double slipMag = slipVertex[0];
      const double slipRateMag = slipRateVertex[0];
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
      throw std::logic_error("Unknown spatial dimension in updateStateVars().");
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

  slipSection->view("SLIP");
  areaSection->view("AREA");
  dispTSection->view("DISP (t)");
  dispTIncrSection->view("DISP INCR (t->t+dt)");

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
    switch (spaceDim) { // switch
    case 1: { // case 1
      // Sensitivity of slip to changes in the Lagrange multipliers
      // Aixjx = 1.0/Aix + 1.0/Ajx
      const double Aixjx = 1.0 / jacobianVertexN[0] + 1.0
          / jacobianVertexP[0];
      const double Spp = 1.0;

      if (tractionTpdtVertex[0] < 0) {
        // if compression, then no changes to solution
      } else {
        // if tension, then traction is zero.

        // Update slip based on value required to stick versus
        // zero traction
        dLagrangeTpdtVertex[0] = tractionTpdtVertex[0] * (*areaVertex);
        slipVertex[0] += Spp * dLagrangeTpdtVertex[0];

        // Set traction to zero.
        tractionTpdtVertex[0] = 0.0;
      } // else
      PetscLogFlops(0); // :TODO: Fix this
      break;
    } // case 1
    case 2: { // case 2
      std::cout << "Normal traction:" << tractionTpdtVertex[1] << std::endl;

      // Sensitivity of slip to changes in the Lagrange multipliers
      // Aixjx = 1.0/Aix + 1.0/Ajx
      assert(jacobianVertexN[0] > 0.0);
      assert(jacobianVertexP[0] > 0.0);
      const double Aixjx = 1.0 / jacobianVertexN[0] + 1.0
          / jacobianVertexP[0];
      // Aiyjy = 1.0/Aiy + 1.0/Ajy
      assert(jacobianVertexN[1] > 0.0);
      assert(jacobianVertexP[1] > 0.0);
      const double Aiyjy = 1.0 / jacobianVertexN[1] + 1.0
          / jacobianVertexP[1];
      const double Cpx = orientationVertex[0];
      const double Cpy = orientationVertex[1];
      const double Cqx = orientationVertex[2];
      const double Cqy = orientationVertex[3];
      const double Spp = Cpx * Cpx * Aixjx + Cpy * Cpy * Aiyjy;
      const double Spq = Cpx * Cqx * Aixjx + Cpy * Cqy * Aiyjy;
      const double Sqq = Cqx * Cqx * Aixjx + Cqy * Cqy * Aiyjy;

      const double tractionNormal = tractionTpdtVertex[1];
      const double slip = slipVertex[0];
      const double slipRate = slipRateVertex[0];

      if (tractionTpdtVertex[1] < 0 && 0.0 == slipVertex[1]) {
        // if in compression and no opening
        std::cout << "FAULT IN COMPRESSION" << std::endl;
        const double frictionStress = _friction->calcFriction(slip, slipRate,
          tractionNormal);
        std::cout << "frictionStress: " << frictionStress << std::endl;
        if (tractionTpdtVertex[0] > frictionStress || (tractionTpdtVertex[0]
            < frictionStress && slipVertex[0] > 0.0)) {
          // traction is limited by friction, so have sliding
          std::cout << "LIMIT TRACTION, HAVE SLIDING" << std::endl;

          // Update slip based on value required to stick versus friction
          dLagrangeTpdtVertex[0] = (tractionTpdtVertex[0] - frictionStress)
              * (*areaVertex);
          slipVertex[0] += Spp * dLagrangeTpdtVertex[0];
          std::cout << "Estimated slip: " << slipVertex[0] << std::endl;
          // Limit traction
          tractionTpdtVertex[0] = frictionStress;
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
        dLagrangeTpdtVertex[0] = tractionTpdtVertex[0] * (*areaVertex);
        dLagrangeTpdtVertex[1] = tractionTpdtVertex[1] * (*areaVertex);
        slipVertex[0] += Spp * dLagrangeTpdtVertex[0] + Spq
            * dLagrangeTpdtVertex[1];
        slipVertex[1] += Spq * dLagrangeTpdtVertex[0] + Sqq
            * dLagrangeTpdtVertex[1];

        // Set traction to zero
        tractionTpdtVertex = 0.0;
      } // else
      PetscLogFlops(0); // :TODO: Fix this
      break;
    } // case 2
    case 3: { // case 3
      std::cout << "Normal traction:" << tractionTpdtVertex[2] << std::endl;

      // Sensitivity of slip to changes in the Lagrange multipliers
      // Aixjx = 1.0/Aix + 1.0/Ajx
      assert(jacobianVertexN[0] > 0.0);
      assert(jacobianVertexP[0] > 0.0);
      const double Aixjx = 1.0 / jacobianVertexN[0] + 1.0
          / jacobianVertexP[0];
      // Aiyjy = 1.0/Aiy + 1.0/Ajy
      assert(jacobianVertexN[1] > 0.0);
      assert(jacobianVertexP[1] > 0.0);
      const double Aiyjy = 1.0 / jacobianVertexN[1] + 1.0
          / jacobianVertexP[1];
      // Aizjz = 1.0/Aiz + 1.0/Ajz
      assert(jacobianVertexN[2] > 0.0);
      assert(jacobianVertexP[2] > 0.0);
      const double Aizjz = 1.0 / jacobianVertexN[2] + 1.0
          / jacobianVertexP[2];
      const double Cpx = orientationVertex[0];
      const double Cpy = orientationVertex[1];
      const double Cpz = orientationVertex[2];
      const double Cqx = orientationVertex[3];
      const double Cqy = orientationVertex[4];
      const double Cqz = orientationVertex[5];
      const double Crx = orientationVertex[6];
      const double Cry = orientationVertex[7];
      const double Crz = orientationVertex[8];
      const double Spp = Cpx * Cpx * Aixjx + Cpy * Cpy * Aiyjy + Cpz * Cpz
          * Aizjz;
      const double Spq = Cpx * Cqx * Aixjx + Cpy * Cqy * Aiyjy + Cpz * Cqz
          * Aizjz;
      const double Spr = Cpx * Crx * Aixjx + Cpy * Cry * Aiyjy + Cpz * Crz
          * Aizjz;
      const double Sqq = Cqx * Cqx * Aixjx + Cqy * Cqy * Aiyjy + Cqz * Cqz
          * Aizjz;
      const double Sqr = Cqx * Crx * Aixjx + Cqy * Cry * Aiyjy + Cqz * Crz
          * Aizjz;
      const double Srr = Crx * Crx * Aixjx + Cry * Cry * Aiyjy + Crz * Crz
          * Aizjz;

      double slip = sqrt(pow(slipVertex[1], 2) + pow(slipVertex[0], 2));
      double slipRate = sqrt(pow(slipRateVertex[1], 2) + pow(slipRateVertex[0],
        2));

      const double tractionNormalVertex = tractionTpdtVertex[2];
      const double tractionShearVertex = sqrt(pow(tractionTpdtVertex[1], 2) + pow(
        tractionTpdtVertex[0], 2));
      const double slipShearVertex = sqrt(pow(slipVertex[1], 2) + pow(slipVertex[0], 2));

      if (tractionNormalVertex < 0 && 0 == slipVertex[2]) {
        // if in compression and no opening
        std::cout << "FAULT IN COMPRESSION" << std::endl;
        const double frictionStress = _friction->calcFriction(slip, slipRate,
          tractionNormalVertex);
        std::cout << "frictionStress: " << frictionStress << std::endl;
        if (tractionShearVertex > frictionStress || (tractionShearVertex
            < frictionStress && slipShearVertex > 0.0)) {
          // traction is limited by friction, so have sliding
          std::cout << "LIMIT TRACTION, HAVE SLIDING" << std::endl;

          // Update slip based on value required to stick versus friction
          dLagrangeTpdtVertex[0] = (tractionShearVertex - frictionStress)
              * tractionTpdtVertex[0] / tractionShearVertex * (*areaVertex);
          dLagrangeTpdtVertex[1] = (tractionShearVertex - frictionStress)
              * tractionTpdtVertex[1] / tractionShearVertex * (*areaVertex);
          slipVertex[0] += Spp * dLagrangeTpdtVertex[0] + Spq
              * dLagrangeTpdtVertex[1];

          slipVertex[1] += Spq * dLagrangeTpdtVertex[0] + Sqq
              * dLagrangeTpdtVertex[1];

          std::cout << "Estimated slip: " << "  " << slipVertex[0] << "  "
              << slipVertex[1] << "  " << slipVertex[2] << std::endl;

          // Limit traction
          tractionTpdtVertex[0] = frictionStress * tractionTpdtVertex[0]
              / tractionShearVertex;
          tractionTpdtVertex[1] = frictionStress * tractionTpdtVertex[1]
              / tractionShearVertex;
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
        dLagrangeTpdtVertex[0] = tractionTpdtVertex[0] * (*areaVertex);
        dLagrangeTpdtVertex[1] = tractionTpdtVertex[1] * (*areaVertex);
        dLagrangeTpdtVertex[2] = tractionTpdtVertex[2] * (*areaVertex);
        slipVertex[0] += Spp * dLagrangeTpdtVertex[0] + Spq
            * dLagrangeTpdtVertex[1] + Spr * dLagrangeTpdtVertex[2];
        slipVertex[1] += Spq * dLagrangeTpdtVertex[0] + Sqq
            * dLagrangeTpdtVertex[1] + Sqr * dLagrangeTpdtVertex[2];
        slipVertex[2] += Spr * dLagrangeTpdtVertex[0] + Sqr
            * dLagrangeTpdtVertex[1] + Srr * dLagrangeTpdtVertex[2];

        std::cout << "Estimated slip: " << "  " << slipVertex[0] << "  "
            << slipVertex[1] << "  " << slipVertex[2] << std::endl;

        // Set traction to zero
        tractionTpdtVertex = 0.0;
      } // else
      PetscLogFlops(0); // :TODO: Fix this
      break;
    } // case 3
    default:
      assert(0);
      throw std::logic_error("Unknown spatial dimension in updateStateVars().");
    } // switch

    // Update Lagrange multiplier values.
    // :KLUDGE: (TEMPORARY) Solution at Lagrange constraint vertices
    // is the Lagrange multiplier value, which is currently the
    // force.  Compute force by multipling traction by area
    lagrangeTIncrVertex =
      (tractionTpdtVertex - tractionTVertex) * (*areaVertex);
    assert(lagrangeTIncrVertex.size() ==
        dispTIncrSection->getFiberDimension(v_lagrange));
    dispTIncrSection->updatePoint(v_lagrange, &lagrangeTIncrVertex[0]);

    // Update the slip estimate based on adjustment to the Lagrange
    // multiplier values.
    assert(slipVertex.size() ==
        slipSection->getFiberDimension(v_fault));
    slipSection->updatePoint(v_fault, &slipVertex[0]);
  } // if

  dispTIncrSection->view("AFTER DISP INCR (t->t+dt)");
  slipSection->view("AFTER SLIP");
} // constrainSolnSpace

// ----------------------------------------------------------------------
// Adjust solution from solver with lumped Jacobian to match Lagrange
// multiplier constraints.
void
pylith::faults::FaultCohesiveDyn::adjustSolnLumped(topology::SolutionFields* const fields,
                                                        const topology::Field<
              topology::Mesh>& jacobian)
{ // adjustSolnLumped
  // :TODO: Update this to constrain solution space using friction
  assert(false);
#if 0
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

  double_array solutionVertexN(spaceDim);
  double_array solutionVertexP(spaceDim);
  double_array solutionVertexL(spaceDim);
  const ALE::Obj<RealSection>& solutionSection = fields->get(
    "dispIncr(t->t+dt)").section();

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

    // Get slip at fault cell's vertices.
    slipSection->restrictPoint(v_fault, &slipVertex[0], slipVertex.size());

    // Get residual at cohesive cell's vertices.
    residualSection->restrictPoint(v_negative, &residualVertexN[0], residualVertexN.size());
    residualSection->restrictPoint(v_positive, &residualVertexP[0], residualVertexP.size());

    // Get jacobian at cohesive cell's vertices.
    jacobianSection->restrictPoint(v_negative, &jacobianVertexN[0], jacobianVertexN.size());
    jacobianSection->restrictPoint(v_positive, &jacobianVertexP[0], jacobianVertexP.size());

    // Get disp(t) at cohesive cell's vertices.
    dispTSection->restrictPoint(v_negative, &dispTVertexN[0], dispTVertexN.size());
    dispTSection->restrictPoint(v_positive, &dispTVertexP[0], dispTVertexP.size());

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(restrictEvent);
    _logger->eventBegin(computeEvent);
#endif

      switch (spaceDim) { // switch
    case 1: {
      assert(jacobianVertexN[0] > 0.0);
      assert(jacobianVertexP[0] > 0.0);

      const double Sinv = jacobianVertexN[0] * jacobianVertexP[0]
          / (jacobianVertexN[0] + jacobianVertexP[0]);

      // Aru = A_i^{-1} r_i - A_j^{-1} r_j + u_i - u_j
      const double Aru = residualVertexN[0] / jacobianVertexN[0]
          - residualVertexP[0] / jacobianVertexP[0] + dispTVertexN[0]
          - dispTVertexP[0];

      // dl_k = D^{-1}( - C_{ki} Aru - d_k)
      const double Aruslip = -Aru - slipVertex[0];
      const double dlp = Sinv * Aruslip;

      // Update displacements at negative vertex
      solutionVertexN[0] = +1.0 / jacobianVertexN[0] * dlp;

      // Update displacements at positive vertex
      solutionVertexP[0] = -1.0 / jacobianVertexP[0] * dlp;

      // Update Lagrange multiplier.
      solutionVertexL[0] = dlp;

      break;
    } // case 1
    case 2: {
      assert(jacobianVertexN[0] > 0.0);
      assert(jacobianVertexN[1] > 0.0);
      assert(jacobianVertexP[0] > 0.0);
      assert(jacobianVertexP[1] > 0.0);

      const double Cpx = orientationVertex[0];
      const double Cpy = orientationVertex[1];
      const double Cqx = orientationVertex[2];
      const double Cqy = orientationVertex[3];

      // Check to make sure Jacobian is same at all DOF for
      // vertices i and j (means S is diagonal with equal enties).
      assert(jacobianVertexN[0] == jacobianVertexN[1]);
      assert(jacobianVertexP[0] == jacobianVertexP[1]);

      const double Sinv = jacobianVertexN[0] * jacobianVertexP[0]
          / (jacobianVertexN[0] + jacobianVertexP[0]);

      // Aru = A_i^{-1} r_i - A_j^{-1} r_j + u_i - u_j
      const double Arux = residualVertexN[0] / jacobianVertexN[0]
          - residualVertexP[0] / jacobianVertexP[0] + dispTVertexN[0]
          - dispTVertexP[0];
      const double Aruy = residualVertexN[1] / jacobianVertexN[1]
          - residualVertexP[1] / jacobianVertexP[1] + dispTVertexN[1]
          - dispTVertexP[1];

      // dl_k = S^{-1}(-C_{ki} Aru - d_k)
      const double Arup = Cpx * Arux + Cpy * Aruy;
      const double Aruq = Cqx * Arux + Cqy * Aruy;
      const double Arupslip = -Arup - slipVertex[0];
      const double Aruqslip = -Aruq - slipVertex[1];
      const double dlp = Sinv * Arupslip;
      const double dlq = Sinv * Aruqslip;

      const double dlx = Cpx * dlp + Cqx * dlq;
      const double dly = Cpy * dlp + Cqy * dlq;

      // Update displacements at negative vertex.
      solutionVertexN[0] = dlx / jacobianVertexN[0];
      solutionVertexN[1] = dly / jacobianVertexN[1];

      // Update displacements at positive vertex.
      solutionVertexP[0] = -dlx / jacobianVertexP[0];
      solutionVertexP[1] = -dly / jacobianVertexP[0];

      // Update Lagrange multiplier.
      solutionVertexL[0] = dlp;
      solutionVertexL[1] = dlq;

      break;
    } // case 2
    case 3: {
      assert(jacobianVertexN[0] > 0.0);
      assert(jacobianVertexN[1] > 0.0);
      assert(jacobianVertexN[2] > 0.0);
      assert(jacobianVertexP[0] > 0.0);
      assert(jacobianVertexP[1] > 0.0);
      assert(jacobianVertexP[2] > 0.0);

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
      assert(jacobianVertexN[0] == jacobianVertexN[1] && jacobianVertexN[0] == jacobianVertexN[2]);
      assert(jacobianVertexP[0] == jacobianVertexP[1] && jacobianVertexP[0] == jacobianVertexP[2]);

      const double Sinv = jacobianVertexN[0] * jacobianVertexP[0]
          / (jacobianVertexN[0] + jacobianVertexP[0]);

      // Aru = A_i^{-1} r_i - A_j^{-1} r_j + u_i - u_j
      const double Arux = residualVertexN[0] / jacobianVertexN[0]
          - residualVertexP[0] / jacobianVertexP[0] + dispTVertexN[0]
          - dispTVertexP[0];
      const double Aruy = residualVertexN[1] / jacobianVertexN[1]
          - residualVertexP[1] / jacobianVertexP[1] + dispTVertexN[1]
          - dispTVertexP[1];
      const double Aruz = residualVertexN[2] / jacobianVertexN[2]
          - residualVertexP[2] / jacobianVertexP[2] + dispTVertexN[2]
          - dispTVertexP[2];

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

      const double dlx = Cpx * dlp + Cqx * dlq + Crx * dlr;
      const double dly = Cpy * dlp + Cqy * dlq + Cry * dlr;
      const double dlz = Cpz * dlp + Cqz * dlq + Crz * dlr;

      // Update displacements at negative vertex.
      solutionVertexN[0] = dlx / jacobianVertexN[0];
      solutionVertexN[1] = dly / jacobianVertexN[1];
      solutionVertexN[2] = dlz / jacobianVertexN[2];

      // Update displacements at positive vertex.
      solutionVertexP[0] = -dlx / jacobianVertexP[0];
      solutionVertexP[1] = -dly / jacobianVertexP[1];
      solutionVertexP[2] = -dlz / jacobianVertexP[2];

      // Update Lagrange multiplier.
      solutionVertexL[0] = dlp;
      solutionVertexL[1] = dlq;
      solutionVertexL[2] = dlr;

      break;
    } // case 3
    default:
      assert(0);
      throw std::logic_error("Unknown spatial dimension.");
    } // switch

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(computeEvent);
    _logger->eventBegin(updateEvent);
#endif

    assert(solutionVertexN.size() == solutionSection->getFiberDimension(v_negative));
    solutionSection->updateAddPoint(v_negative, &solutionVertexN[0]);

    assert(solutionVertexP.size() == solutionSection->getFiberDimension(v_positive));
    solutionSection->updateAddPoint(v_positive, &solutionVertexP[0]);

    assert(solutionVertexL.size() == solutionSection->getFiberDimension(v_lagrange));
    solutionSection->updateAddPoint(v_lagrange, &solutionVertexL[0]);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(updateEvent);
#endif
  } // for

  switch(spaceDim) {
  case 1:
    PetscLogFlops(numVertices*17);
    break;
  case 2:
    PetscLogFlops(numVertices*41);
    break;
  case 3:
    PetscLogFlops(numVertices*72);
    break;
  default:
    assert(0);
    throw std::logic_error("Unknown spatial dimension.");
  } // switch

#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventEnd(computeEvent);
#endif
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
    // :TODO: Need to use scratch buffer to convert initial forces to
    // initial tractions.
    assert(0 != _dbInitialTract);
    const topology::Field<topology::SubMesh>& initialTraction = _fields->get(
      "initial traction");
    return initialTraction;

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
pylith::faults::FaultCohesiveDyn::_getInitialTractions(void)
{ // _getInitialTractions
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

      // Compute action for traction bc terms
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

  //intialForces.view("INITIAL FORCES"); // DEBUGGING
} // _getInitialTractions

// ----------------------------------------------------------------------
// Compute change in tractions on fault surface using solution.
// NOTE: We must convert vertex labels to fault vertex labels
void
pylith::faults::FaultCohesiveDyn::_calcTractions(
    topology::Field<topology::SubMesh>* tractions,
    const topology::Field<topology::Mesh>& dispT)
{ // _calcTractionsChange
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
    tractions->newSection(slip, fiberDim);
    tractions->allocate();

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
    tractionsSection->updateAddPoint(v_fault, &tractionsVertex[0]);
  } // for

  PetscLogFlops(numVertices * (1 + fiberDim) );

#if 0 // DEBUGGING
  tractions->view("TRACTIONS");
#endif

} // _calcTractions

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


// End of file 
