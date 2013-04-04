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

#include "TestFaultCohesiveDyn.hh" // Implementation of class methods

#include "pylith/faults/FaultCohesiveDyn.hh" // USES FaultCohesiveDyn
#include "pylith/faults/TractPerturbation.hh" // USES TractPerturbation

#include "data/CohesiveDynData.hh" // USES CohesiveDynData

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
// Setup testing data.
void
pylith::faults::TestFaultCohesiveDyn::setUp(void)
{ // setUp
  _data = 0;
  _quadrature = new feassemble::Quadrature<topology::SubMesh>();
  CPPUNIT_ASSERT(0 != _quadrature);
  _tractPerturbation = 0;
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
  delete _tractPerturbation; _tractPerturbation = 0;
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
// Test tractPerturbation().
void
pylith::faults::TestFaultCohesiveDyn::testTractPerturbation(void)
{ // testTractPerturbation
  FaultCohesiveDyn fault;

  const std::string& label = "test database";
  TractPerturbation tract;
  tract.label(label.c_str());
  fault.tractPerturbation(&tract);
  CPPUNIT_ASSERT(fault._tractPerturbation);
 } // testTractPerturbation

// ----------------------------------------------------------------------
// Test zeroTolerance().
void
pylith::faults::TestFaultCohesiveDyn::testZeroTolerance(void)
{ // testZeroTolerance
  FaultCohesiveDyn fault;

  CPPUNIT_ASSERT_EQUAL(PylithScalar(1.0e-10), fault._zeroTolerance); // default

  const PylithScalar value = 1.0e-20;
  fault.zeroTolerance(value);
  CPPUNIT_ASSERT_EQUAL(value, fault._zeroTolerance);
 } // zeroTolerance

// ----------------------------------------------------------------------
// Test openFreeSurf().
void
pylith::faults::TestFaultCohesiveDyn::testOpenFreeSurf(void)
{ // testOpenFreeSurf
  FaultCohesiveDyn fault;

  CPPUNIT_ASSERT_EQUAL(true, fault._openFreeSurf); // default

  const bool value = false;
  fault.openFreeSurf(value);
  CPPUNIT_ASSERT_EQUAL(value, fault._openFreeSurf);
 } // testOpenFreeSurf

// ----------------------------------------------------------------------
// Test initialize().
void
pylith::faults::TestFaultCohesiveDyn::testInitialize(void)
{ // testInitialize
  CPPUNIT_ASSERT(_data);

  topology::Mesh mesh;
  FaultCohesiveDyn fault;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &fault, &fields);

  PetscDM dmMesh = fault._faultMesh->dmMesh();CPPUNIT_ASSERT(dmMesh);
  PetscIS subpointIS = NULL;
  const PetscInt *points = NULL;
  PetscInt vStart, vEnd, numPoints;
  PetscErrorCode  err;
  err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);
  err = DMPlexCreateSubpointIS(dmMesh, &subpointIS);CHECK_PETSC_ERROR(err);CPPUNIT_ASSERT(subpointIS);
  err = ISGetSize(subpointIS, &numPoints);CHECK_PETSC_ERROR(err);
  err = ISGetIndices(subpointIS, &points);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt faultPoint;
    err = PetscFindInt(_data->negativeVertices[v-vStart], numPoints, points, &faultPoint);CHECK_PETSC_ERROR(err);
    CPPUNIT_ASSERT(faultPoint >= 0);
    CPPUNIT_ASSERT_EQUAL(faultPoint, v);
  } // for
  err = ISRestoreIndices(subpointIS, &points);CHECK_PETSC_ERROR(err);
  err = ISDestroy(&subpointIS);CHECK_PETSC_ERROR(err);
  CPPUNIT_ASSERT_EQUAL(_data->numConstraintVert, vEnd-vStart);

  // Check orientation
  //fault._fields->get("orientation").view("ORIENTATION"); // DEBUGGING
  PetscSection orientationSection = fault._fields->get("orientation").petscSection();CPPUNIT_ASSERT(orientationSection);
  PetscVec orientationVec = fault._fields->get("orientation").localVector();CPPUNIT_ASSERT(orientationVec);
  PetscScalar *orientationArray = NULL;
  err = VecGetArray(orientationVec, &orientationArray);CHECK_PETSC_ERROR(err);

  const int spaceDim = _data->spaceDim;
  const int orientationSize = spaceDim*spaceDim;
  int iVertex = 0;
  for(PetscInt v = vStart; v < vEnd; ++v, ++iVertex) {
    PetscInt dof, off;
    err = PetscSectionGetDof(orientationSection, v, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(orientationSection, v, &off);CHECK_PETSC_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(orientationSize, dof);

    const PylithScalar tolerance = 1.0e-06;
    for(PetscInt d = 0; d < orientationSize; ++d) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(_data->orientation[iVertex*orientationSize+d], orientationArray[off+d], tolerance);
    } // for
  } // for
  err = VecRestoreArray(orientationVec, &orientationArray);CHECK_PETSC_ERROR(err);

  // Prescribed traction perturbation
  if (fault._tractPerturbation) {
    // :KLUDGE: Only check initial value
    PetscSection initialTractionsSection = fault.vertexField("traction_initial_value").petscSection();CPPUNIT_ASSERT(initialTractionsSection);
    PetscVec initialTractionsVec = fault.vertexField("traction_initial_value").localVector();CPPUNIT_ASSERT(initialTractionsVec);
    PetscScalar *initialTractionsArray = NULL;
    err = VecGetArray(initialTractionsVec, &initialTractionsArray);CHECK_PETSC_ERROR(err);

    const PylithScalar pressureScale = _data->pressureScale;

    iVertex = 0;
    for(PetscInt v = vStart; v < vEnd; ++v, ++iVertex) {
      PetscInt dof, off;
      err = PetscSectionGetDof(initialTractionsSection, v, &dof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(initialTractionsSection, v, &off);CHECK_PETSC_ERROR(err);
      CPPUNIT_ASSERT_EQUAL(spaceDim, dof);

      const PylithScalar tolerance = 1.0e-06;
      for(PetscInt d = 0; d < spaceDim; ++d) {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(_data->initialTractions[iVertex * spaceDim + d], initialTractionsArray[off+d]*pressureScale, tolerance);
      } // for
    } // for
    err = VecRestoreArray(initialTractionsVec, &initialTractionsArray);CHECK_PETSC_ERROR(err);
  } // if
} // testInitialize

// ----------------------------------------------------------------------
// Test constrainSolnSpace() for sticking case.
void
pylith::faults::TestFaultCohesiveDyn::testConstrainSolnSpaceStick(void)
{ // testConstrainSolnSpaceStick
  assert(_data);

  topology::Mesh mesh;
  FaultCohesiveDyn fault;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &fault, &fields);
  topology::Jacobian jacobian(fields.solution());
  _setFieldsJacobian(&mesh, &fault, &fields, &jacobian, _data->fieldIncrStick);

  const int spaceDim = _data->spaceDim;

  const PylithScalar t = 2.134 / _data->timeScale;
  const PylithScalar dt = 0.01;
  fault.timeStep(dt);
  fault.constrainSolnSpace(&fields, t, jacobian);
  
  topology::Field<topology::Mesh>& solution = fields.solution();
  const topology::Field<topology::Mesh>& dispIncrAdj = fields.get("dispIncr adjust");
  solution += dispIncrAdj;

  fault.updateStateVars(t, &fields);

  { // Check solution values
    // No change to Lagrange multipliers for stick case.
    PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
    PetscInt vStart, vEnd;
    PetscErrorCode err;
    err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);

    // Get section containing solution (disp + Lagrange multipliers)
    PetscSection dispTIncrSection = fields.get("dispIncr(t->t+dt)").petscSection();CPPUNIT_ASSERT(dispTIncrSection);
    PetscVec dispTIncrVec = fields.get("dispIncr(t->t+dt)").localVector();CPPUNIT_ASSERT(dispTIncrVec);
    PetscScalar *dispTIncrArray = NULL;
    err = VecGetArray(dispTIncrVec, &dispTIncrArray);CHECK_PETSC_ERROR(err);

    // Get expected values
    const PylithScalar* valsE = _data->fieldIncrStick; // No change in dispIncr
    int iVertex = 0; // variable to use as index into valsE array
    const int fiberDimE = spaceDim; // number of values per point
    const PylithScalar tolerance = 1.0e-06;
    for(PetscInt v = vStart; v < vEnd; ++v, ++iVertex) {
      PetscInt dof, off;
      err = PetscSectionGetDof(dispTIncrSection, v, &dof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(dispTIncrSection, v, &off);CHECK_PETSC_ERROR(err);
      assert(fiberDimE == dof);

      // Check values at point
      for(PetscInt d = 0; d < fiberDimE; ++d) {
        const PylithScalar valE = valsE[iVertex*spaceDim+d];

        if (fabs(valE) > tolerance)
          CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, dispTIncrArray[off+d]/valE, tolerance);
        else
          CPPUNIT_ASSERT_DOUBLES_EQUAL(valE, dispTIncrArray[off+d], tolerance);
      } // for
    } // for
    err = VecRestoreArray(dispTIncrVec, &dispTIncrArray);CHECK_PETSC_ERROR(err);
  } // Check solution values

  { // Check slip values
    // Slip should be zero for the stick case.

    // Get fault vertex info
    PetscDM faultDMMesh = fault._faultMesh->dmMesh();CPPUNIT_ASSERT(faultDMMesh);
    PetscInt vStart, vEnd;
    PetscErrorCode err;
    err = DMPlexGetDepthStratum(faultDMMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);

    // Get section containing slip
    PetscSection slipSection = fault.vertexField("slip").petscSection();CPPUNIT_ASSERT(slipSection);
    PetscVec slipVec = fault.vertexField("slip").localVector();CPPUNIT_ASSERT(slipVec);
    PetscScalar *slipArray = NULL;
    err = VecGetArray(slipVec, &slipArray);CHECK_PETSC_ERROR(err);

    // Get expected values
    CPPUNIT_ASSERT(_data->slipStickE);
    const PylithScalar* valsE = _data->slipStickE;
    int iVertex = 0; // variable to use as index into valsE array
    const int fiberDimE = spaceDim; // number of values per point
    const PylithScalar tolerance = 1.0e-06;
    for(PetscInt v = vStart; v < vEnd; ++v, ++iVertex) {
      PetscInt dof, off;
      err = PetscSectionGetDof(slipSection, v, &dof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(slipSection, v, &off);CHECK_PETSC_ERROR(err);
      CPPUNIT_ASSERT_EQUAL(fiberDimE, dof);

      // Check values at point
      for(PetscInt d = 0; d < dof; ++d) {
        const PylithScalar valE = valsE[iVertex*spaceDim+d];

        if (fabs(valE) > tolerance)
          CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, slipArray[off+d]/valE, tolerance);
        else
          CPPUNIT_ASSERT_DOUBLES_EQUAL(valE, slipArray[off+d], tolerance);
      }
    } // for
    err = VecRestoreArray(slipVec, &slipArray);CHECK_PETSC_ERROR(err);
  } // Check slip values

} // testConstrainSolnSpaceStick

// ----------------------------------------------------------------------
// Test constrainSolnSpace() for slipping case.
void
pylith::faults::TestFaultCohesiveDyn::testConstrainSolnSpaceSlip(void)
{ // testConstrainSolnSpaceSlip
  assert(_data);

  topology::Mesh mesh;
  FaultCohesiveDyn fault;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &fault, &fields);
  topology::Jacobian jacobian(fields.solution());
  _setFieldsJacobian(&mesh, &fault, &fields, &jacobian, _data->fieldIncrSlip);

  const int spaceDim = _data->spaceDim;

  const PylithScalar t = 2.134 / _data->timeScale;
  const PylithScalar dt = 0.01;
  fault.timeStep(dt);
  fault.constrainSolnSpace(&fields, t, jacobian);

  topology::Field<topology::Mesh>& solution = fields.solution();
  const topology::Field<topology::Mesh>& dispIncrAdj = fields.get("dispIncr adjust");
  solution += dispIncrAdj;

  fault.updateStateVars(t, &fields);

  { // Check solution values
    // Lagrange multipliers should be adjusted according to friction
    // as reflected in the fieldIncrSlipE data member.
    PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
    PetscInt vStart, vEnd;
    PetscErrorCode  err;
    err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);

    // Get section containing solution (disp + Lagrange multipliers)
    PetscSection dispIncrSection = fields.get("dispIncr(t->t+dt)").petscSection();CPPUNIT_ASSERT(dispIncrSection);
    PetscVec dispIncrVec = fields.get("dispIncr(t->t+dt)").localVector();CPPUNIT_ASSERT(dispIncrVec);
    PetscScalar *dispIncrArray = NULL;
    err = VecGetArray(dispIncrVec, &dispIncrArray);CHECK_PETSC_ERROR(err);

    // Get expected values
    const PylithScalar* valsE = _data->fieldIncrSlipE; // Expected values for dispIncr
    int iVertex = 0; // variable to use as index into valsE array
    const int fiberDimE = spaceDim; // number of values per point
    const PylithScalar tolerance = 1.0e-06;
    for(PetscInt v = vStart; v < vEnd; ++v, ++iVertex) {
      PetscInt dof, off;
      err = PetscSectionGetDof(dispIncrSection, v, &dof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(dispIncrSection, v, &off);CHECK_PETSC_ERROR(err);
      CPPUNIT_ASSERT_EQUAL(fiberDimE, dof);

      for(PetscInt d = 0; d < dof; ++d) {
        const PylithScalar valE = valsE[iVertex*spaceDim+d];
        if (fabs(valE) > tolerance)
          CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, dispIncrArray[off+d]/valE, tolerance);
        else
          CPPUNIT_ASSERT_DOUBLES_EQUAL(valE, dispIncrArray[off+d], tolerance);
      } // for
    } // for
    err = VecRestoreArray(dispIncrVec, &dispIncrArray);CHECK_PETSC_ERROR(err);
  } // Check solution values

  { // Check slip values
    // Slip values should be adjusted based on the change in the
    // Lagrange multipliers as reflected in the slipSlipE data member.

    // Get fault vertex info
    PetscDM faultDMMesh = fault._faultMesh->dmMesh();CPPUNIT_ASSERT(faultDMMesh);
    PetscInt vStart, vEnd;
    PetscErrorCode err;
    err = DMPlexGetDepthStratum(faultDMMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);

    // Get section containing slip
    PetscSection slipSection = fault.vertexField("slip").petscSection();CPPUNIT_ASSERT(slipSection);
    PetscVec slipVec = fault.vertexField("slip").localVector();CPPUNIT_ASSERT(slipVec);
    PetscScalar *slipArray = NULL;
    err = VecGetArray(slipVec, &slipArray);CHECK_PETSC_ERROR(err);

    // Get expected values
    const PylithScalar* valsE = _data->slipSlipE;
    int iVertex = 0; // variable to use as index into valsE array
    const int fiberDimE = spaceDim; // number of values per point
    const PylithScalar tolerance = (sizeof(double) == sizeof(PylithScalar)) ? 1.0e-06 : 1.0e-5;
    for(PetscInt v = vStart; v < vEnd; ++v, ++iVertex) {
      PetscInt dof, off;
      err = PetscSectionGetDof(slipSection, v, &dof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(slipSection, v, &off);CHECK_PETSC_ERROR(err);
      CPPUNIT_ASSERT_EQUAL(fiberDimE, dof);

      for(PetscInt d = 0; d < dof; ++d) {
        const PylithScalar valE = valsE[iVertex*spaceDim+d];
        if (fabs(valE) > tolerance)
          CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, slipArray[off+d]/valE, tolerance);
        else
          CPPUNIT_ASSERT_DOUBLES_EQUAL(valE, slipArray[off+d], tolerance);
      } // for
    } // for
    err = VecRestoreArray(slipVec, &slipArray);CHECK_PETSC_ERROR(err);
  } // Check slip values

} // testConstrainSolnSpaceSlip

// ----------------------------------------------------------------------
// Test constrainSolnSpace() for opening case.
void
pylith::faults::TestFaultCohesiveDyn::testConstrainSolnSpaceOpen(void)
{ // testConstrainSolnSpaceOpen
  assert(_data);

  topology::Mesh mesh;
  FaultCohesiveDyn fault;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &fault, &fields);
  topology::Jacobian jacobian(fields.solution());
  _setFieldsJacobian(&mesh, &fault, &fields, &jacobian, _data->fieldIncrOpen);

  const int spaceDim = _data->spaceDim;

  const PylithScalar t = 2.134 / _data->timeScale;
  const PylithScalar dt = 0.01;
  fault.timeStep(dt);
  fault.constrainSolnSpace(&fields, t, jacobian);

  topology::Field<topology::Mesh>& solution = fields.solution();
  const topology::Field<topology::Mesh>& dispIncrAdj = fields.get("dispIncr adjust");
  solution += dispIncrAdj;

  fault.updateStateVars(t, &fields);

  //residual.view("RESIDUAL"); // DEBUGGING

  { // Check solution values
    // Lagrange multipliers should be set to zero as reflected in the
    // fieldIncrOpenE data member.
    PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
    PetscInt vStart, vEnd;
    PetscErrorCode  err;
    err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);

    // Get section containing solution (disp + Lagrange multipliers)
    PetscSection dispIncrSection = fields.get("dispIncr(t->t+dt)").petscSection();CPPUNIT_ASSERT(dispIncrSection);
    PetscVec dispIncrVec = fields.get("dispIncr(t->t+dt)").localVector();CPPUNIT_ASSERT(dispIncrVec);
    PetscScalar *dispIncrArray = NULL;
    err = VecGetArray(dispIncrVec, &dispIncrArray);CHECK_PETSC_ERROR(err);

    // Get expected values
    const PylithScalar* valsE = _data->fieldIncrOpenE; // Expected values for dispIncr
    int iVertex = 0; // variable to use as index into valsE array
    const int fiberDimE = spaceDim; // number of values per point
    const PylithScalar tolerance = (sizeof(double) == sizeof(PylithScalar)) ? 1.0e-06 : 1.0e-05;
    for(PetscInt v = vStart; v < vEnd; ++v, ++iVertex) {
      PetscInt dof, off;
      err = PetscSectionGetDof(dispIncrSection, v, &dof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(dispIncrSection, v, &off);CHECK_PETSC_ERROR(err);
      CPPUNIT_ASSERT_EQUAL(fiberDimE, dof);

      for(PetscInt d = 0; d < dof; ++d) {
        const PylithScalar valE = valsE[iVertex*spaceDim+d];
        if (fabs(valE) > tolerance)
          CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, dispIncrArray[off+d]/valE, tolerance);
        else
          CPPUNIT_ASSERT_DOUBLES_EQUAL(valE, dispIncrArray[off+d], tolerance);
      } // for
    } // for
    err = VecRestoreArray(dispIncrVec, &dispIncrArray);CHECK_PETSC_ERROR(err);
  } // Check solution values

  { // Check slip values
    // Slip values should be adjusted based on the change in the
    // Lagrange multipliers as reflected in the slipOpenE data member.

    // Get fault vertex info
    PetscDM faultDMMesh = fault._faultMesh->dmMesh();CPPUNIT_ASSERT(faultDMMesh);
    PetscInt vStart, vEnd;
    PetscErrorCode err;
    err = DMPlexGetDepthStratum(faultDMMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);

    // Get section containing slip
    PetscSection slipSection = fault.vertexField("slip").petscSection();CPPUNIT_ASSERT(slipSection);
    PetscVec slipVec = fault.vertexField("slip").localVector();CPPUNIT_ASSERT(slipVec);
    PetscScalar *slipArray = NULL;
    err = VecGetArray(slipVec, &slipArray);CHECK_PETSC_ERROR(err);

    // Get expected values
    const PylithScalar* valsE = _data->slipOpenE;
    int iVertex = 0; // variable to use as index into valsE array
    const int fiberDimE = spaceDim; // number of values per point
    const PylithScalar tolerance = (sizeof(double) == sizeof(PylithScalar)) ? 1.0e-06 : 1.0e-05;
    for(PetscInt v = vStart; v < vEnd; ++v, ++iVertex) {
      PetscInt dof, off;
      err = PetscSectionGetDof(slipSection, v, &dof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(slipSection, v, &off);CHECK_PETSC_ERROR(err);
      CPPUNIT_ASSERT_EQUAL(fiberDimE, dof);

      for(PetscInt d = 0; d < dof; ++d) {
        const PylithScalar valE = valsE[iVertex*spaceDim+d];
        if (fabs(valE) > tolerance)
          CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, slipArray[off+d]/valE, tolerance);
        else
          CPPUNIT_ASSERT_DOUBLES_EQUAL(valE, slipArray[off+d], tolerance);
      } // for
    } // for
    err = VecRestoreArray(slipVec, &slipArray);CHECK_PETSC_ERROR(err);
  } // Check slip values

} // testConstrainSolnSpaceOpen

// ----------------------------------------------------------------------
// Test updateStateVars().
void
pylith::faults::TestFaultCohesiveDyn::testUpdateStateVars(void)
{ // testUpdateStateVars
  assert(_data);

  topology::Mesh mesh;
  FaultCohesiveDyn fault;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &fault, &fields);
  topology::Jacobian jacobian(fields.solution());
  _setFieldsJacobian(&mesh, &fault, &fields, &jacobian, _data->fieldIncrSlip);

  const int spaceDim = _data->spaceDim;

  const PylithScalar t = 2.134;
  const PylithScalar dt = 0.01;
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
  CPPUNIT_ASSERT(_data);

  topology::Mesh mesh;
  FaultCohesiveDyn fault;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &fault, &fields);
  topology::Jacobian jacobian(fields.solution());
  _setFieldsJacobian(&mesh, &fault, &fields, &jacobian, _data->fieldIncrStick);
  
  CPPUNIT_ASSERT(fault._faultMesh);
  const int spaceDim = _data->spaceDim;
  topology::Field<topology::SubMesh> tractions(*fault._faultMesh);
  tractions.newSection(topology::FieldBase::VERTICES_FIELD, spaceDim);
  tractions.allocate();
  tractions.zero();

  PetscSection tractionsSection = tractions.petscSection();CPPUNIT_ASSERT(tractionsSection);
  PetscVec tractionsVec = tractions.localVector();CPPUNIT_ASSERT(tractionsVec);
  PetscScalar *tractionsArray = NULL;

  const PylithScalar t = 0;
  fault.updateStateVars(t, &fields);
  fault._calcTractions(&tractions, fields.get("disp(t)"));

  PetscDM faultDMMesh = fault._faultMesh->dmMesh();CPPUNIT_ASSERT(faultDMMesh);
  PetscIS subpointIS = NULL;
  const PetscInt *points = NULL;
  PetscInt vStart, vEnd, numPoints;
  PetscErrorCode err;
  err = DMPlexGetDepthStratum(faultDMMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);
  err = DMPlexCreateSubpointIS(faultDMMesh, &subpointIS);CHECK_PETSC_ERROR(err);CPPUNIT_ASSERT(subpointIS);
  err = ISGetSize(subpointIS, &numPoints);CHECK_PETSC_ERROR(err);

  PetscSection dispSection = fields.get("disp(t)").petscSection();CPPUNIT_ASSERT(dispSection);
  PetscVec dispVec = fields.get("disp(t)").localVector();CPPUNIT_ASSERT(dispVec);
  PetscScalar *dispArray = NULL;

  err = VecGetArray(tractionsVec, &tractionsArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(dispVec, &dispArray);CHECK_PETSC_ERROR(err);
  err = ISGetIndices(subpointIS, &points);CHECK_PETSC_ERROR(err);

  int iVertex = 0;
  const PylithScalar tolerance = 1.0e-06;
  for(PetscInt v = vStart; v < vEnd; ++v, ++iVertex) {
    PetscInt tdof, toff, ddof, doff;
    err = PetscSectionGetDof(tractionsSection, v, &tdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(tractionsSection, v, &toff);CHECK_PETSC_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(spaceDim, tdof);
    err = PetscSectionGetDof(dispSection, points[v], &ddof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dispSection, points[v], &doff);CHECK_PETSC_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(spaceDim, ddof);

    const PylithScalar *orientationVertex = &_data->orientation[iVertex*spaceDim*spaceDim];
    CPPUNIT_ASSERT(orientationVertex);

    for(PetscInt d = 0; d < spaceDim; ++d) {
      PylithScalar tractionE = 0.0;
      for(PetscInt e = 0; e < spaceDim; ++e) {
        tractionE += orientationVertex[d*spaceDim+e]*dispArray[doff+e];
      } // for
      if (tractionE != 0.0)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, tractionsArray[toff+d]/tractionE, tolerance);
      else
        CPPUNIT_ASSERT_DOUBLES_EQUAL(tractionE, tractionsArray[toff+d], tolerance);
    } // for
  } // for
  err = VecRestoreArray(tractionsVec, &tractionsArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(dispVec, &dispArray);CHECK_PETSC_ERROR(err);
  err = ISRestoreIndices(subpointIS, &points);CHECK_PETSC_ERROR(err);
  err = ISDestroy(&subpointIS);CHECK_PETSC_ERROR(err);
} // testCalcTractions

// ----------------------------------------------------------------------
// Initialize FaultCohesiveDyn interface condition.
void
pylith::faults::TestFaultCohesiveDyn::_initialize(topology::Mesh* const mesh,
						  FaultCohesiveDyn* const fault,
						  topology::SolutionFields* const fields)
{ // _initialize
  CPPUNIT_ASSERT(mesh);
  CPPUNIT_ASSERT(fault);
  CPPUNIT_ASSERT(fields);
  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT(_quadrature);

  meshio::MeshIOAscii iohandler;
  iohandler.filename(_data->meshFilename);
  iohandler.read(mesh);

  // Set coordinate system
  spatialdata::geocoords::CSCart cs;
  cs.setSpaceDim(mesh->dimension());
  cs.initialize();
  mesh->coordsys(&cs);

  // Set scales
  // Most test data is insensitive to the scales because we set the fields directly.
  spatialdata::units::Nondimensional normalizer;
  normalizer.lengthScale(_data->lengthScale);
  normalizer.pressureScale(_data->pressureScale);
  normalizer.densityScale(_data->densityScale);
  normalizer.timeScale(_data->timeScale);
  mesh->nondimensionalize(normalizer);
  
  _quadrature->initialize(_data->basis, _data->numQuadPts, _data->numBasis,
			  _data->basisDeriv,
			  _data->numQuadPts, _data->numBasis, _data->cellDim,
			  _data->quadPts, _data->numQuadPts, _data->cellDim,
			  _data->quadWts, _data->numQuadPts,
			  _data->spaceDim);
  
  // Setup prescribed traction perturbation
  delete _tractPerturbation; _tractPerturbation = new TractPerturbation();
  _tractPerturbation->label("traction perturbation");
  spatialdata::spatialdb::SimpleDB* db = new spatialdata::spatialdb::SimpleDB("initial tractions");
  CPPUNIT_ASSERT(db);
  spatialdata::spatialdb::SimpleIOAscii ioInitialTract;
  ioInitialTract.filename(_data->initialTractFilename);
  db->ioHandler(&ioInitialTract);
  delete _dbInitialTract; _dbInitialTract = db;
  _tractPerturbation->dbInitial(db);
  fault->tractPerturbation(_tractPerturbation);

  // Setup friction
  spatialdata::spatialdb::SimpleDB* dbFriction =
      new spatialdata::spatialdb::SimpleDB("static friction");
  CPPUNIT_ASSERT(dbFriction);
  spatialdata::spatialdb::SimpleIOAscii ioFriction;
  if (2 == _data->spaceDim)
    ioFriction.filename("data/static_friction_2d.spatialdb");
  else if (3 == _data->spaceDim)
    ioFriction.filename("data/static_friction_3d.spatialdb");
  dbFriction->ioHandler(&ioFriction);
  delete _dbFriction; _dbFriction = dbFriction;
  friction::StaticFriction* friction = new pylith::friction::StaticFriction();
  CPPUNIT_ASSERT(friction);
  friction->label("static friction");
  friction->dbProperties(dbFriction);
  friction->normalizer(normalizer);
  _friction = friction;
  fault->frictionModel(friction);

  PetscInt labelSize;
  PetscErrorCode err;
  err = DMPlexGetStratumSize(mesh->dmMesh(), _data->label, 1, &labelSize);CHECK_PETSC_ERROR(err);

  PetscInt firstFaultVertex = 0;
  PetscInt firstLagrangeVertex = labelSize;
  PetscInt firstFaultCell = labelSize;
  if (fault->useLagrangeConstraints())
    firstFaultCell += labelSize;
  fault->id(_data->id);
  fault->label(_data->label);
  fault->quadrature(_quadrature);
  
  fault->adjustTopology(mesh, &firstFaultVertex, &firstLagrangeVertex,
			&firstFaultCell, _flipFault);
  
  const PylithScalar upDir[3] = { 0.0, 0.0, 1.0 };
  
  fault->normalizer(normalizer);
  fault->initialize(*mesh, upDir);
  
  // Setup fields
  fields->add("disp(t)", "displacement");
  fields->add("dispIncr(t->t+dt)", "displacement_increment");
  fields->add("velocity(t)", "velocity");
  fields->add("dispIncr adjust", "dispIncr_adjust");
  fields->solutionName("dispIncr(t->t+dt)");
  
  const int spaceDim = _data->spaceDim;
  topology::Field<topology::Mesh>& disp = fields->get("disp(t)");
  disp.newSection(topology::FieldBase::VERTICES_FIELD, spaceDim);
  disp.allocate();
  disp.scale(_data->lengthScale);
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
          const PylithScalar* const fieldIncr)
{ // _initialize
  CPPUNIT_ASSERT(mesh);
  CPPUNIT_ASSERT(fault);
  CPPUNIT_ASSERT(fields);
  CPPUNIT_ASSERT(jacobian);
  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT(fieldIncr);

  const int spaceDim = _data->spaceDim;

  // Get vertices in mesh
  PetscDM dmMesh = mesh->dmMesh();CPPUNIT_ASSERT(dmMesh);
  PetscInt vStart, vEnd;
  PetscErrorCode err;
  err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);

  // Set displacement values
  PetscSection dispSection = fields->get("disp(t)").petscSection();CPPUNIT_ASSERT(dispSection);
  PetscVec dispVec = fields->get("disp(t)").localVector();CPPUNIT_ASSERT(dispVec);
  PetscScalar *dispArray = NULL;
  err = VecGetArray(dispVec, &dispArray);CHECK_PETSC_ERROR(err);

  int iVertex = 0;
  for(PetscInt v = vStart; v < vEnd; ++v, ++iVertex) {
    PetscInt off;
    err = PetscSectionGetOffset(dispSection, v, &off);CHECK_PETSC_ERROR(err);
    for(PetscInt d = 0; d < spaceDim; ++d) dispArray[off+d] = _data->fieldT[iVertex*spaceDim+d];
  }
  err = VecRestoreArray(dispVec, &dispArray);CHECK_PETSC_ERROR(err);

  // Set increment values
  PetscSection dispIncrSection = fields->get("dispIncr(t->t+dt)").petscSection();CPPUNIT_ASSERT(dispIncrSection);
  PetscVec dispIncrVec = fields->get("dispIncr(t->t+dt)").localVector();CPPUNIT_ASSERT(dispIncrVec);
  PetscScalar *dispIncrArray = NULL;
  err = VecGetArray(dispIncrVec, &dispIncrArray);CHECK_PETSC_ERROR(err);

  iVertex = 0;
  for(PetscInt v = vStart; v < vEnd; ++v, ++iVertex) {
    PetscInt off;
    err = PetscSectionGetOffset(dispIncrSection, v, &off);CHECK_PETSC_ERROR(err);
    for(PetscInt d = 0; d < spaceDim; ++d) dispIncrArray[off+d] = fieldIncr[iVertex*spaceDim+d];
  }
  err = VecRestoreArray(dispIncrVec, &dispIncrArray);CHECK_PETSC_ERROR(err);

  // Setup Jacobian matrix
  PetscInt nrows;
  err = PetscSectionGetStorageSize(dispIncrSection, &nrows);CHECK_PETSC_ERROR(err);
  PetscInt ncols = nrows;
  int nrowsM = 0;
  int ncolsM = 0;
  PetscMat jacobianMat = jacobian->matrix();
  err = MatGetSize(jacobianMat, &nrowsM, &ncolsM);CHECK_PETSC_ERROR(err);
  CPPUNIT_ASSERT_EQUAL(nrows, nrowsM);
  CPPUNIT_ASSERT_EQUAL(ncols, ncolsM);
  // We ignore the sparsity patterns in our tests
  err = MatSetOption(jacobianMat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);CHECK_PETSC_ERROR(err);

  int_array rows(nrows);
  int_array cols(ncols);
  for (int iRow=0; iRow < nrows; ++iRow)
    rows[iRow] = iRow;
  for (int iCol=0; iCol < ncols; ++iCol)
    cols[iCol] = iCol;
  err = MatSetValues(jacobianMat, nrows, &rows[0], ncols, &cols[0], _data->jacobian, INSERT_VALUES);CHECK_PETSC_ERROR(err);
  jacobian->assemble("final_assembly");

} // _setFieldsJacobian

// ----------------------------------------------------------------------
// Determine if vertex is a Lagrange multiplier constraint vertex.
bool
pylith::faults::TestFaultCohesiveDyn::_isConstraintVertex(const int vertex) const
{ // _isConstraintVertex
  assert(_data);

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
