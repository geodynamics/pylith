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
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestFaultCohesiveDyn.hh" // Implementation of class methods

#include "pylith/faults/FaultCohesiveDyn.hh" // USES FaultCohesiveDyn

#include "data/CohesiveDynData.hh" // USES CohesiveDynData

#include "pylith/faults/EqKinSrc.hh" // USES EqKinSrc
#include "pylith/faults/BruneSlipFn.hh" // USES BruneSlipFn
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/SubMesh.hh" // USES SubMesh
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields
#include "pylith/topology/Jacobian.hh" // USES Jacobian
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/friction/StaticFriction.hh" // USES StaticFriction

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii

#include <stdexcept> // USES runtime_error

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesiveDyn );

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::RealSection RealSection;
typedef pylith::topology::SubMesh::SieveMesh SieveSubMesh;

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::faults::TestFaultCohesiveDyn::setUp(void)
{ // setUp
  _data = 0;
  _quadrature = new feassemble::Quadrature<topology::SubMesh>();
  CPPUNIT_ASSERT(0 != _quadrature);
  _dbInitialTract = 0;
  _friction = 0;
  _dbFriction = 0;
  _flipFault = false;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::faults::TestFaultCohesiveDyn::tearDown(void)
{ // tearDown
  delete _data; _data = 0;
  delete _quadrature; _quadrature = 0;
  delete _dbInitialTract; _dbInitialTract = 0;
  delete _friction; _friction = 0;
  delete _dbFriction; _dbFriction = 0;
} // tearDown

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::faults::TestFaultCohesiveDyn::testConstructor(void)
{ // testConstructor
  FaultCohesiveDyn fault;
} // testConstructor

// ----------------------------------------------------------------------
// Test dbInitialTract().
void
pylith::faults::TestFaultCohesiveDyn::testDBInitialTract(void)
{ // testDBInitialTract
  FaultCohesiveDyn fault;

  const std::string& label = "test database";
  spatialdata::spatialdb::SimpleDB db;
  db.label(label.c_str());
  fault.dbInitialTract(&db);
  CPPUNIT_ASSERT(0 != fault._dbInitialTract);
  CPPUNIT_ASSERT_EQUAL(label, std::string(fault._dbInitialTract->label()));
 } // testDBInitialTract

// ----------------------------------------------------------------------
// Test zeroTolerance().
void
pylith::faults::TestFaultCohesiveDyn::testZeroTolerance(void)
{ // testZeroTolerance
  FaultCohesiveDyn fault;

  CPPUNIT_ASSERT_EQUAL(1.0e-10, fault._zeroTolerance); // default

  const double value = 1.0e-20;
  fault.zeroTolerance(value);
  CPPUNIT_ASSERT_EQUAL(value, fault._zeroTolerance);
 } // zeroTolerance

// ----------------------------------------------------------------------
// Test initialize().
void
pylith::faults::TestFaultCohesiveDyn::testInitialize(void)
{ // testInitialize
  CPPUNIT_ASSERT(0 != _data);

  topology::Mesh mesh;
  FaultCohesiveDyn fault;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &fault, &fields);

  const ALE::Obj<SieveSubMesh>& faultSieveMesh = fault._faultMesh->sieveMesh();
  CPPUNIT_ASSERT(!faultSieveMesh.isNull());
  SieveSubMesh::renumbering_type& renumbering = 
    faultSieveMesh->getRenumbering();
  const ALE::Obj<SieveSubMesh::label_sequence>& vertices = 
    faultSieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());
  const SieveSubMesh::label_sequence::iterator verticesBegin = vertices->begin();
  const SieveSubMesh::label_sequence::iterator verticesEnd = vertices->end();
  int iVertex = 0;
  for (SieveSubMesh::label_sequence::iterator v_iter=verticesBegin;
       v_iter != verticesEnd;
       ++v_iter, ++iVertex) {
    CPPUNIT_ASSERT(renumbering.find(_data->constraintVertices[iVertex]) !=
		   renumbering.end());
    CPPUNIT_ASSERT_EQUAL(renumbering[_data->constraintVertices[iVertex]],
			 *v_iter);
  } // for
  CPPUNIT_ASSERT_EQUAL(_data->numConstraintVert, iVertex);

  // Check orientation
  //fault._fields->get("orientation").view("ORIENTATION"); // DEBUGGING

  const ALE::Obj<RealSection>& orientationSection = 
    fault._fields->get("orientation").section();
  CPPUNIT_ASSERT(!orientationSection.isNull());
  const int spaceDim = _data->spaceDim;
  const int orientationSize = spaceDim*spaceDim;
  iVertex = 0;
  for (SieveSubMesh::label_sequence::iterator v_iter=verticesBegin;
       v_iter != verticesEnd;
       ++v_iter, ++iVertex) {
    const int fiberDim = orientationSection->getFiberDimension(*v_iter);
    CPPUNIT_ASSERT_EQUAL(orientationSize, fiberDim);
    const double* orientationVertex =
      orientationSection->restrictPoint(*v_iter);
    CPPUNIT_ASSERT(0 != orientationVertex);

    const double tolerance = 1.0e-06;
    for (int i=0; i < orientationSize; ++i) {
      const int index = iVertex*orientationSize+i;
      CPPUNIT_ASSERT_DOUBLES_EQUAL(_data->orientation[index],
				   orientationVertex[i], tolerance);
    } // for
  } // for

  // Initial tractions
  if (0 != fault._dbInitialTract) {
    //fault._fields->get("initial traction").view("INITIAL TRACTIONS"); // DEBUGGING
    const ALE::Obj<RealSection>& initialTractionsSection = 
      fault._fields->get("initial traction").section();
    CPPUNIT_ASSERT(!initialTractionsSection.isNull());
    const int spaceDim = _data->spaceDim;
    iVertex = 0;
    for (SieveSubMesh::label_sequence::iterator v_iter = verticesBegin;
        v_iter != verticesEnd;
        ++v_iter, ++iVertex) {
      const int fiberDim = initialTractionsSection->getFiberDimension(*v_iter);
      CPPUNIT_ASSERT_EQUAL(spaceDim, fiberDim);
      const double* initialTractionsVertex = 
	initialTractionsSection->restrictPoint(*v_iter);
      CPPUNIT_ASSERT(initialTractionsVertex);

      const double tolerance = 1.0e-06;
      for (int i = 0; i < spaceDim; ++i) {
        const int index = iVertex * spaceDim + i;
        CPPUNIT_ASSERT_DOUBLES_EQUAL(_data->initialTractions[index], 
				     initialTractionsVertex[i], tolerance);
      } // for
    } // for
  } // if
} // testInitialize

// ----------------------------------------------------------------------
// Test constrainSolnSpace() for sticking case.
void
pylith::faults::TestFaultCohesiveDyn::testConstrainSolnSpaceStick(void)
{ // testConstrainSolnSpaceStick
  assert(0 != _data);

  topology::Mesh mesh;
  FaultCohesiveDyn fault;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &fault, &fields);
  topology::Jacobian jacobian(fields.solution(), "seqdense");
  _setFieldsJacobian(&mesh, &fault, &fields, &jacobian, _data->fieldIncrStick);

  const int spaceDim = _data->spaceDim;

  const double t = 2.134;
  const double dt = 0.01;
  fault.timeStep(dt);
  fault.constrainSolnSpace(&fields, t, jacobian);

  { // Check solution values
    // No change to Lagrange multipliers for stick case.
    const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
    CPPUNIT_ASSERT(!sieveMesh.isNull());
    const ALE::Obj<SieveMesh::label_sequence>& vertices =
      sieveMesh->depthStratum(0);
    CPPUNIT_ASSERT(!vertices.isNull());
    const SieveMesh::label_sequence::iterator verticesBegin = vertices->begin();
    const SieveMesh::label_sequence::iterator verticesEnd = vertices->end();

    // Get section containing solution (disp + Lagrange multipliers)
    const ALE::Obj<RealSection>& dispIncrSection =
      fields.get("dispIncr(t->t+dt)").section();
    CPPUNIT_ASSERT(!dispIncrSection.isNull());

    //dispIncrSection->view("DISP INCREMENT"); // DEBUGGING

    // Get expected values
    const double* valsE = _data->fieldIncrStick; // No change in dispIncr
    int iVertex = 0; // variable to use as index into valsE array
    const int fiberDimE = spaceDim; // number of values per point
    const double tolerance = 1.0e-06;
    for (SieveMesh::label_sequence::iterator v_iter = verticesBegin;
	 v_iter != verticesEnd;
	 ++v_iter, ++iVertex) { // loop over all vertices in mesh
      // Check fiber dimension (number of values at point)
      const int fiberDim = dispIncrSection->getFiberDimension(*v_iter);
      CPPUNIT_ASSERT_EQUAL(fiberDimE, fiberDim);
      const double* vals = dispIncrSection->restrictPoint(*v_iter);
      CPPUNIT_ASSERT(0 != vals);

      // Check values at point
      for (int i = 0; i < fiberDimE; ++i) {
        const int index = iVertex * spaceDim + i;
        const double valE = valsE[index];
        if (fabs(valE) > tolerance)
          CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, vals[i]/valE, tolerance);
        else
          CPPUNIT_ASSERT_DOUBLES_EQUAL(valE, vals[i], tolerance);
      } // for
    } // for
  } // Check solution values

  { // Check slip values
    // Slip should be zero for the stick case.

    // Get fault vertex info
    const ALE::Obj<SieveSubMesh>& faultSieveMesh = fault._faultMesh->sieveMesh();
    CPPUNIT_ASSERT(!faultSieveMesh.isNull());
    const ALE::Obj<SieveSubMesh::label_sequence>& vertices =
      faultSieveMesh->depthStratum(0);
    CPPUNIT_ASSERT(!vertices.isNull());
    const SieveSubMesh::label_sequence::iterator verticesBegin = vertices->begin();
    const SieveSubMesh::label_sequence::iterator verticesEnd = vertices->end();

    // Get section containing slip
    const ALE::Obj<RealSection>& slipSection =
      fault.vertexField("slip").section();
    CPPUNIT_ASSERT(!slipSection.isNull());

    const double valE = 0.0; // slip should be zero
    int iVertex = 0; // variable to use as index into valsE array
    const int fiberDimE = spaceDim; // number of values per point
    const double tolerance = 1.0e-06;
    for (SieveMesh::label_sequence::iterator v_iter = verticesBegin; v_iter
	   != verticesEnd;
	 ++v_iter, ++iVertex) { // loop over fault vertices
      // Check fiber dimension (number of values at point)
      const int fiberDim = slipSection->getFiberDimension(*v_iter);
      CPPUNIT_ASSERT_EQUAL(fiberDimE, fiberDim);
      const double* vals = slipSection->restrictPoint(*v_iter);
      CPPUNIT_ASSERT(0 != vals);

      // Check values at point
      for (int i = 0; i < fiberDimE; ++i) {
        const int index = iVertex * spaceDim + i;
	CPPUNIT_ASSERT_DOUBLES_EQUAL(valE, vals[i], tolerance);
      } // for
    } // for
  } // Check slip values

} // testConstrainSolnSpaceStick

// ----------------------------------------------------------------------
// Test constrainSolnSpace() for slipping case.
void
pylith::faults::TestFaultCohesiveDyn::testConstrainSolnSpaceSlip(void)
{ // testConstrainSolnSpaceSlip
  assert(0 != _data);

  topology::Mesh mesh;
  FaultCohesiveDyn fault;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &fault, &fields);
  topology::Jacobian jacobian(fields.solution(), "seqdense");
  _setFieldsJacobian(&mesh, &fault, &fields, &jacobian, _data->fieldIncrSlip);

  const int spaceDim = _data->spaceDim;

  const double t = 2.134;
  const double dt = 0.01;
  fault.timeStep(dt);
  fault.constrainSolnSpace(&fields, t, jacobian);

  { // Check slip values
    // Slip values should be adjusted based on the change in the
    // Lagrange multipliers as reflected in the slipSlipE data member.

    // Get fault vertex info
    const ALE::Obj<SieveSubMesh>& faultSieveMesh = fault._faultMesh->sieveMesh();
    CPPUNIT_ASSERT(!faultSieveMesh.isNull());
    const ALE::Obj<SieveSubMesh::label_sequence>& vertices =
      faultSieveMesh->depthStratum(0);
    CPPUNIT_ASSERT(!vertices.isNull());
    const SieveSubMesh::label_sequence::iterator verticesBegin = vertices->begin();
    const SieveSubMesh::label_sequence::iterator verticesEnd = vertices->end();

    // Get section containing slip
    const ALE::Obj<RealSection>& slipSection =
      fault.vertexField("slip").section();
    CPPUNIT_ASSERT(!slipSection.isNull());

    //slipSection->view("SLIP"); // DEBUGGING

    // Get expected values
    const double* valsE = _data->slipSlipE;
    int iVertex = 0; // variable to use as index into valsE array
    const int fiberDimE = spaceDim; // number of values per point
    const double tolerance = 1.0e-06;
    for (SieveMesh::label_sequence::iterator v_iter = verticesBegin; v_iter
	   != verticesEnd;
	 ++v_iter, ++iVertex) { // loop over fault vertices
      // Check fiber dimension (number of values at point)
      const int fiberDim = slipSection->getFiberDimension(*v_iter);
      CPPUNIT_ASSERT_EQUAL(fiberDimE, fiberDim);
      const double* vals = slipSection->restrictPoint(*v_iter);
      CPPUNIT_ASSERT(0 != vals);

      // Check values at point
      for (int i = 0; i < fiberDimE; ++i) {
        const int index = iVertex * spaceDim + i;
        const double valE = valsE[index];
#if 0 // DEBUGGING
	std::cout << "SLIP valE: " << valE
		  << ", val: " << vals[i]
		  << ", error: " << fabs(1.0-vals[i]/valE)
		  << std::endl;
#endif // DEBUGGING
        if (fabs(valE) > tolerance)
	  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, vals[i]/valE, tolerance);
        else
	  CPPUNIT_ASSERT_DOUBLES_EQUAL(valE, vals[i], tolerance);
      } // for
    } // for
  } // Check slip values

  { // Check solution values
    // Lagrange multipliers should be adjusted according to friction
    // as reflected in the fieldIncrSlipE data member.
    const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
    CPPUNIT_ASSERT(!sieveMesh.isNull());
    const ALE::Obj<SieveMesh::label_sequence>& vertices =
      sieveMesh->depthStratum(0);
    CPPUNIT_ASSERT(!vertices.isNull());
    const SieveMesh::label_sequence::iterator verticesBegin = vertices->begin();
    const SieveMesh::label_sequence::iterator verticesEnd = vertices->end();

    // Get section containing solution (disp + Lagrange multipliers)
    const ALE::Obj<RealSection>& dispIncrSection =
      fields.get("dispIncr(t->t+dt)").section();
    CPPUNIT_ASSERT(!dispIncrSection.isNull());

    // Get expected values
    const double* valsE = _data->fieldIncrSlipE; // Expected values for dispIncr
    int iVertex = 0; // variable to use as index into valsE array
    const int fiberDimE = spaceDim; // number of values per point
    const double tolerance = 1.0e-06;
    for (SieveMesh::label_sequence::iterator v_iter = verticesBegin;
	 v_iter != verticesEnd;
	 ++v_iter, ++iVertex) { // loop over all vertices in mesh
      // Check fiber dimension (number of values at point)
      const int fiberDim = dispIncrSection->getFiberDimension(*v_iter);
      CPPUNIT_ASSERT_EQUAL(fiberDimE, fiberDim);
      const double* vals = dispIncrSection->restrictPoint(*v_iter);
      CPPUNIT_ASSERT(0 != vals);

      // Check values at point
      for (int i = 0; i < fiberDimE; ++i) {
        const int index = iVertex * spaceDim + i;
        const double valE = valsE[index];
#if 0 // DEBUGGING
	std::cout << "SOLUTION valE: " << valE
		  << ", val: " << vals[i]
		  << ", error: " << fabs(1.0-vals[i]/valE)
		  << std::endl;
#endif // DEBUGGING
	if (fabs(valE) > tolerance)
	  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, vals[i]/valE, tolerance);
        else
	  CPPUNIT_ASSERT_DOUBLES_EQUAL(valE, vals[i], tolerance);
      } // for
    } // for
  } // Check solution values

} // testConstrainSolnSpaceSlip

// ----------------------------------------------------------------------
// Test constrainSolnSpace() for opening case.
void
pylith::faults::TestFaultCohesiveDyn::testConstrainSolnSpaceOpen(void)
{ // testConstrainSolnSpaceOpen
  assert(0 != _data);

  topology::Mesh mesh;
  FaultCohesiveDyn fault;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &fault, &fields);
  topology::Jacobian jacobian(fields.solution(), "seqdense");
  _setFieldsJacobian(&mesh, &fault, &fields, &jacobian, _data->fieldIncrOpen);

  const int spaceDim = _data->spaceDim;

  const double t = 2.134;
  const double dt = 0.01;
  fault.timeStep(dt);
  fault.constrainSolnSpace(&fields, t, jacobian);

  //residual.view("RESIDUAL"); // DEBUGGING

  { // Check solution values
    // Lagrange multipliers should be set to zero as reflected in the
    // fieldIncrOpenE data member.
    const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
    CPPUNIT_ASSERT(!sieveMesh.isNull());
    const ALE::Obj<SieveMesh::label_sequence>& vertices =
      sieveMesh->depthStratum(0);
    CPPUNIT_ASSERT(!vertices.isNull());
    const SieveMesh::label_sequence::iterator verticesBegin = vertices->begin();
    const SieveMesh::label_sequence::iterator verticesEnd = vertices->end();

    // Get section containing solution (disp + Lagrange multipliers)
    const ALE::Obj<RealSection>& dispIncrSection =
      fields.get("dispIncr(t->t+dt)").section();
    CPPUNIT_ASSERT(!dispIncrSection.isNull());

    // Get expected values
    const double* valsE = _data->fieldIncrOpenE; // Expected values for dispIncr
    int iVertex = 0; // variable to use as index into valsE array
    const int fiberDimE = spaceDim; // number of values per point
    const double tolerance = 1.0e-06;
    for (SieveMesh::label_sequence::iterator v_iter = verticesBegin;
	 v_iter != verticesEnd;
	 ++v_iter, ++iVertex) { // loop over all vertices in mesh
      // Check fiber dimension (number of values at point)
      const int fiberDim = dispIncrSection->getFiberDimension(*v_iter);
      CPPUNIT_ASSERT_EQUAL(fiberDimE, fiberDim);
      const double* vals = dispIncrSection->restrictPoint(*v_iter);
      CPPUNIT_ASSERT(0 != vals);

      // Check values at point
      for (int i = 0; i < fiberDimE; ++i) {
        const int index = iVertex * spaceDim + i;
        const double valE = valsE[index];
#if 0 // DEBUGGING
	std::cout << "valE: " << valE
		  << ", val: " << vals[i]
		  << std::endl;
#endif // DEBUGGING
        if (fabs(valE) > tolerance)
	  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, vals[i]/valE, tolerance);
        else
	  CPPUNIT_ASSERT_DOUBLES_EQUAL(valE, vals[i], tolerance);
      } // for
    } // for
  } // Check solution values

  { // Check slip values
    // Slip values should be adjusted based on the change in the
    // Lagrange multipliers as reflected in the slipOpenE data member.

    // Get fault vertex info
    const ALE::Obj<SieveSubMesh>& faultSieveMesh = fault._faultMesh->sieveMesh();
    CPPUNIT_ASSERT(!faultSieveMesh.isNull());
    const ALE::Obj<SieveSubMesh::label_sequence>& vertices =
      faultSieveMesh->depthStratum(0);
    CPPUNIT_ASSERT(!vertices.isNull());
    const SieveSubMesh::label_sequence::iterator verticesBegin = vertices->begin();
    const SieveSubMesh::label_sequence::iterator verticesEnd = vertices->end();

    // Get section containing slip
    const ALE::Obj<RealSection>& slipSection =
      fault.vertexField("slip").section();
    CPPUNIT_ASSERT(!slipSection.isNull());

    // Get expected values
    const double* valsE = _data->slipOpenE;
    int iVertex = 0; // variable to use as index into valsE array
    const int fiberDimE = spaceDim; // number of values per point
    const double tolerance = 1.0e-06;
    for (SieveMesh::label_sequence::iterator v_iter = verticesBegin; v_iter
	   != verticesEnd;
	 ++v_iter, ++iVertex) { // loop over fault vertices
      // Check fiber dimension (number of values at point)
      const int fiberDim = slipSection->getFiberDimension(*v_iter);
      CPPUNIT_ASSERT_EQUAL(fiberDimE, fiberDim);
      const double* vals = slipSection->restrictPoint(*v_iter);
      CPPUNIT_ASSERT(0 != vals);

      // Check values at point
      for (int i = 0; i < fiberDimE; ++i) {
        const int index = iVertex * spaceDim + i;
        const double valE = valsE[index];
#if 0 // DEBUGGING
	std::cout << "valE: " << valE
		  << ", val: " << vals[i]
		  << std::endl;
#endif // DEBUGGING
        if (fabs(valE) > tolerance)
	  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, vals[i]/valE, tolerance);
        else
	  CPPUNIT_ASSERT_DOUBLES_EQUAL(valE, vals[i], tolerance);
      } // for
    } // for
  } // Check slip values

} // testConstrainSolnSpaceOpen

// ----------------------------------------------------------------------
// Test updateStateVars().
void
pylith::faults::TestFaultCohesiveDyn::testUpdateStateVars(void)
{ // testUpdateStateVars
  assert(0 != _data);

  topology::Mesh mesh;
  FaultCohesiveDyn fault;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &fault, &fields);
  topology::Jacobian jacobian(fields.solution(), "seqdense");
  _setFieldsJacobian(&mesh, &fault, &fields, &jacobian, _data->fieldIncrSlip);

  const int spaceDim = _data->spaceDim;

  const double t = 2.134;
  const double dt = 0.01;
  fault.timeStep(dt);
  fault.updateStateVars(t, &fields);

  // :TODO: Need to verify that fault constitutive updateStateVars is called.
  // We don't have a way to verify state variables inside friction object.
} // testUpdateStateVars

// ----------------------------------------------------------------------
// Test calcTractions().
void
pylith::faults::TestFaultCohesiveDyn::testCalcTractions(void)
{ // testCalcTractions
  CPPUNIT_ASSERT(0 != _data);

  topology::Mesh mesh;
  FaultCohesiveDyn fault;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &fault, &fields);
  topology::Jacobian jacobian(fields.solution(), "seqdense");
  _setFieldsJacobian(&mesh, &fault, &fields, &jacobian, _data->fieldIncrStick);
  
  CPPUNIT_ASSERT(0 != fault._faultMesh);
  const int spaceDim = _data->spaceDim;
  topology::Field<topology::SubMesh> tractions(*fault._faultMesh);
  tractions.newSection(topology::FieldBase::VERTICES_FIELD, spaceDim);
  tractions.allocate();
  tractions.zero();
  const ALE::Obj<RealSection>& tractionsSection = tractions.section();
  CPPUNIT_ASSERT(!tractionsSection.isNull());

  const double t = 0;
  fault.updateStateVars(t, &fields);
  fault._calcTractions(&tractions, fields.get("disp(t)"));

  //tractions.view("TRACTIONS"); // DEBUGGING

  const ALE::Obj<SieveSubMesh>& faultSieveMesh = fault._faultMesh->sieveMesh();
  CPPUNIT_ASSERT(!faultSieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& vertices =
    faultSieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());
  const SieveSubMesh::label_sequence::iterator verticesBegin =
    vertices->begin();
  const SieveSubMesh::label_sequence::iterator verticesEnd = vertices->end();
  SieveSubMesh::renumbering_type& renumbering =
    faultSieveMesh->getRenumbering();
  const SieveMesh::renumbering_type::const_iterator renumberingBegin =
    renumbering.begin();
  const SieveMesh::renumbering_type::const_iterator renumberingEnd =
    renumbering.end();

  const ALE::Obj<RealSection>& dispSection = fields.get("disp(t)").section();
  CPPUNIT_ASSERT(!dispSection.isNull());

  int iVertex = 0;
  const double tolerance = 1.0e-06;
  for (SieveMesh::label_sequence::iterator v_iter=verticesBegin;
       v_iter != verticesEnd;
       ++v_iter, ++iVertex) {
    SieveMesh::point_type meshVertex = -1;
    bool found = false;

    for (SieveMesh::renumbering_type::const_iterator r_iter = renumberingBegin;
      r_iter != renumberingEnd;
      ++r_iter) {
      if (r_iter->second == *v_iter) {
        meshVertex = r_iter->first;
        found = true;
        break;
      } // if
    } // for
    CPPUNIT_ASSERT(found);
    int fiberDim = tractionsSection->getFiberDimension(*v_iter);
    CPPUNIT_ASSERT_EQUAL(spaceDim, fiberDim);
    const double* tractionsVertex = tractionsSection->restrictPoint(*v_iter);
    CPPUNIT_ASSERT(tractionsVertex);

    const double* tractionsVertexGlobalE = 
      dispSection->restrictPoint(meshVertex);
    CPPUNIT_ASSERT(tractionsVertexGlobalE);
    const double* orientationVertex = 
      &_data->orientation[iVertex*spaceDim*spaceDim];
    CPPUNIT_ASSERT(orientationVertex);

    for (int iDim=0; iDim < spaceDim; ++iDim) {
      double tractionE = 0.0;
      for (int jDim=0; jDim < spaceDim; ++jDim) {
	tractionE += 
	  orientationVertex[iDim*spaceDim+jDim]*tractionsVertexGlobalE[jDim];
      } // for
      if (tractionE != 0.0)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, tractionsVertex[iDim]/tractionE,
				     tolerance);
      else
        CPPUNIT_ASSERT_DOUBLES_EQUAL(tractionE, tractionsVertex[iDim],
				     tolerance);
    } // for
  } // for
} // testCalcTractions

// ----------------------------------------------------------------------
// Initialize FaultCohesiveDyn interface condition.
void
pylith::faults::TestFaultCohesiveDyn::_initialize(
					topology::Mesh* const mesh,
					FaultCohesiveDyn* const fault,
					topology::SolutionFields* const fields)
{ // _initialize
  CPPUNIT_ASSERT(0 != mesh);
  CPPUNIT_ASSERT(0 != fault);
  CPPUNIT_ASSERT(0 != fields);
  CPPUNIT_ASSERT(0 != _data);
  CPPUNIT_ASSERT(0 != _quadrature);

  meshio::MeshIOAscii iohandler;
  iohandler.filename(_data->meshFilename);
  iohandler.read(mesh);
  
  //mesh->debug(true); // DEBUGGING
  
  _quadrature->initialize(_data->basis, _data->numQuadPts, _data->numBasis,
			  _data->basisDeriv,
			  _data->numQuadPts, _data->numBasis, _data->cellDim,
			  _data->quadPts, _data->numQuadPts, _data->cellDim,
			  _data->quadWts, _data->numQuadPts,
			  _data->spaceDim);
  
  // Setup initial tractions
  spatialdata::spatialdb::SimpleDB* db =
      new spatialdata::spatialdb::SimpleDB("initial tractions");
  CPPUNIT_ASSERT(0 != db);
  spatialdata::spatialdb::SimpleIOAscii ioInitialTract;
  ioInitialTract.filename(_data->initialTractFilename);
  db->ioHandler(&ioInitialTract);
  delete _dbInitialTract; _dbInitialTract = db;
  fault->dbInitialTract(db);

  // Setup friction
  spatialdata::spatialdb::SimpleDB* dbFriction =
      new spatialdata::spatialdb::SimpleDB("static friction");
  CPPUNIT_ASSERT(0 != dbFriction);
  spatialdata::spatialdb::SimpleIOAscii ioFriction;
  if (2 == _data->spaceDim)
    ioFriction.filename("data/static_friction_2d.spatialdb");
  else if (3 == _data->spaceDim)
    ioFriction.filename("data/static_friction_3d.spatialdb");
  dbFriction->ioHandler(&ioFriction);
  delete _dbFriction; _dbFriction = dbFriction;
  friction::StaticFriction* friction = new pylith::friction::StaticFriction();
  CPPUNIT_ASSERT(0 != friction);
  friction->label("static friction");
  friction->dbProperties(dbFriction);
  _friction = friction;
  fault->frictionModel(friction);

  int firstFaultVertex    = 0;
  int firstLagrangeVertex = mesh->sieveMesh()->getIntSection(_data->label)->size();
  int firstFaultCell      = mesh->sieveMesh()->getIntSection(_data->label)->size();
  if (fault->useLagrangeConstraints())
    firstFaultCell += mesh->sieveMesh()->getIntSection(_data->label)->size();
  fault->id(_data->id);
  fault->label(_data->label);
  fault->quadrature(_quadrature);
  
  fault->adjustTopology(mesh, &firstFaultVertex, &firstLagrangeVertex,
			&firstFaultCell, _flipFault);
  
  spatialdata::geocoords::CSCart cs;
  spatialdata::units::Nondimensional normalizer;
  cs.setSpaceDim(mesh->dimension());
  cs.initialize();
  mesh->coordsys(&cs);
  mesh->nondimensionalize(normalizer);
  
  const double upDir[] = { 0.0, 0.0, 1.0 };
  
  fault->initialize(*mesh, upDir);
  
  // Setup fields
  fields->add("disp(t)", "displacement");
  fields->add("dispIncr(t->t+dt)", "displacement_increment");
  fields->add("velocity(t)", "velocity");
  fields->solutionName("dispIncr(t->t+dt)");
  
  const int spaceDim = _data->spaceDim;
  topology::Field<topology::Mesh>& disp = fields->get("disp(t)");
  disp.newSection(topology::FieldBase::VERTICES_FIELD, spaceDim);
  disp.allocate();
  fields->copyLayout("disp(t)");

  fault->verifyConfiguration(*mesh);
} // _initialize

// ----------------------------------------------------------------------
// Set values for fields and Jacobian.
void
pylith::faults::TestFaultCohesiveDyn::_setFieldsJacobian(
          topology::Mesh* const mesh,
          FaultCohesiveDyn* const fault,
          topology::SolutionFields* const fields,
          topology::Jacobian* const jacobian,
          const double* const fieldIncr)
{ // _initialize
  CPPUNIT_ASSERT(0 != mesh);
  CPPUNIT_ASSERT(0 != fault);
  CPPUNIT_ASSERT(0 != fields);
  CPPUNIT_ASSERT(0 != jacobian);
  CPPUNIT_ASSERT(0 != _data);
  CPPUNIT_ASSERT(0 != fieldIncr);

  const int spaceDim = _data->spaceDim;

  // Get vertices in mesh
  const ALE::Obj<SieveMesh>& sieveMesh = mesh->sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& vertices =
    sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());
  const SieveMesh::label_sequence::iterator verticesBegin = vertices->begin();
  const SieveMesh::label_sequence::iterator verticesEnd = vertices->end();

  // Set displacement values
  const ALE::Obj<RealSection>& dispSection =
    fields->get("disp(t)").section();
  CPPUNIT_ASSERT(!dispSection.isNull());
  int iVertex = 0;
  for (SieveMesh::label_sequence::iterator v_iter=verticesBegin;
       v_iter != verticesEnd;
       ++v_iter, ++iVertex)
    dispSection->updatePoint(*v_iter, &_data->fieldT[iVertex*spaceDim]);

  // Set increment values
  const ALE::Obj<RealSection>& dispIncrSection =
    fields->get("dispIncr(t->t+dt)").section();
  CPPUNIT_ASSERT(!dispIncrSection.isNull());
  iVertex = 0;
  for (SieveMesh::label_sequence::iterator v_iter=verticesBegin;
       v_iter != verticesEnd;
       ++v_iter, ++iVertex)
    dispIncrSection->updatePoint(*v_iter, &fieldIncr[iVertex*spaceDim]);

  // Setup Jacobian matrix
  const int nrows = dispIncrSection->sizeWithBC();
  const int ncols = nrows;
  int nrowsM = 0;
  int ncolsM = 0;
  PetscMat jacobianMat = jacobian->matrix();
  MatGetSize(jacobianMat, &nrowsM, &ncolsM);
  CPPUNIT_ASSERT_EQUAL(nrows, nrowsM);
  CPPUNIT_ASSERT_EQUAL(ncols, ncolsM);

  int_array rows(nrows);
  int_array cols(ncols);
  for (int iRow=0; iRow < nrows; ++iRow)
    rows[iRow] = iRow;
  for (int iCol=0; iCol < ncols; ++iCol)
    cols[iCol] = iCol;
  MatSetValues(jacobianMat, nrows, &rows[0], ncols, &cols[0], 
	       _data->jacobian, INSERT_VALUES);
  jacobian->assemble("final_assembly");

} // _setFieldsJacobian

// ----------------------------------------------------------------------
// Determine if vertex is a Lagrange multiplier constraint vertex.
bool
pylith::faults::TestFaultCohesiveDyn::_isConstraintVertex(const int vertex) const
{ // _isConstraintVertex
  assert(0 != _data);

  const int numConstraintVert = _data->numConstraintVert;
  bool isFound = false;
  for (int i=0; i < _data->numConstraintVert; ++i)
    if (_data->constraintVertices[i] == vertex) {
      isFound = true;
      break;
    } // if
  return isFound;
} // _isConstraintVertex


// End of file 
