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

#include "TestFaultCohesiveImpulses.hh" // Implementation of class methods

#include "pylith/faults/FaultCohesiveImpulses.hh" // USES FaultCohesiveImpulses

#include "data/CohesiveImpulsesData.hh" // USES CohesiveImpulsesData

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps::nondimensionalize()
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <stdexcept> // USES runtime_error

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesiveImpulses );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::faults::TestFaultCohesiveImpulses::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  _data = 0;
  _quadrature = new feassemble::Quadrature();CPPUNIT_ASSERT(_quadrature);
  _dbImpulseAmp = 0;

  PYLITH_METHOD_END;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::faults::TestFaultCohesiveImpulses::tearDown(void)
{ // tearDown
  PYLITH_METHOD_BEGIN;

  delete _data; _data = 0;
  delete _quadrature; _quadrature = 0;
  delete _dbImpulseAmp; _dbImpulseAmp = 0;

  PYLITH_METHOD_END;
} // tearDown

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::faults::TestFaultCohesiveImpulses::testConstructor(void)
{ // testConstructor
  PYLITH_METHOD_BEGIN;

  FaultCohesiveImpulses fault;

  PYLITH_METHOD_END;
} // testConstructor

// ----------------------------------------------------------------------
// Test dbImpulseAmp().
void
pylith::faults::TestFaultCohesiveImpulses::testDBImpulseAmp(void)
{ // testDBImpulseAmp
  PYLITH_METHOD_BEGIN;

  FaultCohesiveImpulses fault;

  const std::string& label = "test database";
  spatialdata::spatialdb::SimpleDB db;
  db.label(label.c_str());
  fault.dbImpulseAmp(&db);
  CPPUNIT_ASSERT(fault._dbImpulseAmp);
  CPPUNIT_ASSERT_EQUAL(label, std::string(fault._dbImpulseAmp->label()));

  PYLITH_METHOD_END;
 } // testDBImpulseAmp

// ----------------------------------------------------------------------
// Test threshold().
void
pylith::faults::TestFaultCohesiveImpulses::testThreshold(void)
{ // testThreshold
  PYLITH_METHOD_BEGIN;

  FaultCohesiveImpulses fault;

  CPPUNIT_ASSERT_EQUAL(PylithScalar(1.0e-6), fault._threshold); // default

  const PylithScalar value = 1.0e-20;
  fault.threshold(value);
  CPPUNIT_ASSERT_EQUAL(value, fault._threshold);
 } // testThreshold

// ----------------------------------------------------------------------
// Test impulseDOF().
void
pylith::faults::TestFaultCohesiveImpulses::testImpulseDOF(void)
{ // testImpulseDOF
  PYLITH_METHOD_BEGIN;

  FaultCohesiveImpulses fault;

  const int ncomps = 2;
  const int dof[ncomps] = { 0, 2 };
  fault.impulseDOF(dof, ncomps);

  CPPUNIT_ASSERT_EQUAL(ncomps, fault.numComponents());
  CPPUNIT_ASSERT_EQUAL(ncomps, int(fault._impulseDOF.size()));
  for (int i=0; i < ncomps; ++i) {
    CPPUNIT_ASSERT_EQUAL(dof[i], fault._impulseDOF[i]);
  } // for
 } // testImpulseDOF

// ----------------------------------------------------------------------
// Test numImpulses().
void
pylith::faults::TestFaultCohesiveImpulses::testNumImpulses(void)
{ // testNumImpulses
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_data);

  topology::Mesh mesh;
  FaultCohesiveImpulses fault;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &fault, &fields);

  const int ncomps = 2;
  const int dof[ncomps] = { 0, 2 };
  fault.impulseDOF(dof, ncomps);

  CPPUNIT_ASSERT_EQUAL(_data->numImpulses, fault.numImpulses());
 } // testNumImpulses

// ----------------------------------------------------------------------
// Test initialize().
void
pylith::faults::TestFaultCohesiveImpulses::testInitialize(void)
{ // testInitialize
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_data);

  topology::Mesh mesh;
  FaultCohesiveImpulses fault;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &fault, &fields);

  PetscDM dmMesh = fault._faultMesh->dmMesh();CPPUNIT_ASSERT(dmMesh);
  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();

  PetscIS subpointIS = NULL;
  const PetscInt *points = NULL;
  PetscInt numPoints = 0;
  PetscErrorCode err;
  err = DMPlexCreateSubpointIS(dmMesh, &subpointIS);PYLITH_CHECK_ERROR(err);CPPUNIT_ASSERT(subpointIS);
  err = ISGetSize(subpointIS, &numPoints);PYLITH_CHECK_ERROR(err);
  err = ISGetIndices(subpointIS, &points);PYLITH_CHECK_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt faultPoint;

    err = PetscFindInt(_data->negativeVertices[v-vStart], numPoints, points, &faultPoint);PYLITH_CHECK_ERROR(err);
    CPPUNIT_ASSERT(faultPoint >= 0);
    CPPUNIT_ASSERT_EQUAL(faultPoint, v);
  } // for
  err = ISRestoreIndices(subpointIS, &points);PYLITH_CHECK_ERROR(err);
  err = ISDestroy(&subpointIS);PYLITH_CHECK_ERROR(err);
  CPPUNIT_ASSERT_EQUAL(_data->numConstraintEdges, vEnd-vStart);

  // Check orientation
  //fault._fields->get("orientation").view("ORIENTATION"); // DEBUGGING

  topology::VecVisitorMesh orientationVisitor(fault._fields->get("orientation"));
  const PetscScalar* orientationArray = orientationVisitor.localArray();

  const int spaceDim = _data->spaceDim;
  const int orientationSize = spaceDim*spaceDim;
  for(PetscInt v = vStart, iVertex = 0; v < vEnd; ++v, ++iVertex) {
    const PetscInt off = orientationVisitor.sectionOffset(v);
    CPPUNIT_ASSERT_EQUAL(orientationSize, orientationVisitor.sectionDof(v));

    const PylithScalar tolerance = 1.0e-06;
    for(PetscInt d = 0; d < orientationSize; ++d) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(_data->orientation[iVertex*orientationSize+d], orientationArray[off+d], tolerance);
    } // for
  } // for

  // Initial tractions
  if (fault._dbImpulseAmp) {
    topology::VecVisitorMesh amplitudeVisitor(fault._fields->get("impulse amplitude"));
    const PetscScalar* amplitudeArray = amplitudeVisitor.localArray();
    const PylithScalar amplitudeScale = _data->lengthScale;

    for(PetscInt v = vStart, iVertex = 0; v < vEnd; ++v, ++iVertex) {
      const PetscInt off = amplitudeVisitor.sectionOffset(v);
      CPPUNIT_ASSERT_EQUAL(1, amplitudeVisitor.sectionDof(v));

      const PylithScalar tolerance = 1.0e-06;
      CPPUNIT_ASSERT_DOUBLES_EQUAL(_data->amplitude[iVertex], amplitudeArray[off]*amplitudeScale, tolerance);
    } // for
  } // if

  PYLITH_METHOD_END;
} // testInitialize

// ----------------------------------------------------------------------
// Test integrateResidual().
void
pylith::faults::TestFaultCohesiveImpulses::testIntegrateResidual(void)
{ // testIntegrateResidual
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT(_data->fieldT);
  CPPUNIT_ASSERT(_data->residual);

  topology::Mesh mesh;
  FaultCohesiveImpulses fault;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &fault, &fields);

  const int spaceDim = _data->spaceDim;

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  PetscInt pStart, pEnd;

  // Set displacement values
  topology::Field& disp = fields.get("disp(t)");
  topology::VecVisitorMesh dispVisitor(disp);
  PetscErrorCode err = PetscSectionGetChart(disp.localSection(), &pStart, &pEnd);CPPUNIT_ASSERT(!err);
  PetscScalar* dispArray = dispVisitor.localArray();CPPUNIT_ASSERT(dispArray);
  for (PetscInt p = pStart, iVertex = 0; p < pEnd; ++p) {
    if (dispVisitor.sectionDof(p) > 0) {
      const PetscInt off = dispVisitor.sectionOffset(p);
      CPPUNIT_ASSERT_EQUAL(spaceDim, dispVisitor.sectionDof(p));
      for(PetscInt d = 0; d < spaceDim; ++d) {
	dispArray[off+d] = _data->fieldT[iVertex*spaceDim+d] / _data->lengthScale;
      } // for
      ++iVertex;
    } // if
  } // for

  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();

  const PylithScalar t = 1.0;
  const PylithScalar dt = 1.0;
  fault.timeStep(dt);

  topology::Field& residual = fields.get("residual");
  residual.zero();
  fault.integrateResidual(residual, t, &fields);
  //residual.view("RESIDUAL"); // DEBUGGING

  topology::VecVisitorMesh residualVisitor(residual);
  const PetscScalar* residualArray = residualVisitor.localArray();CPPUNIT_ASSERT(residualArray);

  // Check values
  const PylithScalar tolerance = (sizeof(double) == sizeof(PylithScalar)) ? 1.0e-06 : 1.0e-05;
  const PylithScalar residualScale = pow(_data->lengthScale, spaceDim);

  for (PetscInt p = pStart, iPoint = 0; p < pEnd; ++p) {
    if (residualVisitor.sectionDof(p) > 0) {
      const PetscInt off = residualVisitor.sectionOffset(p);
      CPPUNIT_ASSERT_EQUAL(spaceDim, residualVisitor.sectionDof(p));
      for(PetscInt d = 0; d < spaceDim; ++d) {
	const PylithScalar valE = _data->residual[iPoint*spaceDim+d];
	if (fabs(valE) > tolerance)
	  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, residualArray[off+d]/valE*residualScale, tolerance);
	else
	  CPPUNIT_ASSERT_DOUBLES_EQUAL(valE, residualArray[off+d]*residualScale, tolerance);
      } // for
      ++iPoint;
    } // if
  } // for

  PYLITH_METHOD_END;
} // testIntegrateResidual

// ----------------------------------------------------------------------
// Initialize FaultCohesiveImpulses interface condition.
void
pylith::faults::TestFaultCohesiveImpulses::_initialize(topology::Mesh* const mesh,
						       FaultCohesiveImpulses* const fault,
						       topology::SolutionFields* const fields)
{ // _initialize
  PYLITH_METHOD_BEGIN;

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
  spatialdata::units::Nondimensional normalizer;
  normalizer.lengthScale(_data->lengthScale);
  normalizer.pressureScale(_data->pressureScale);
  normalizer.densityScale(_data->densityScale);
  normalizer.timeScale(_data->timeScale);
  topology::MeshOps::nondimensionalize(mesh, normalizer);
  
  _quadrature->initialize(_data->basis, _data->numQuadPts, _data->numBasis,
			  _data->basisDeriv,
			  _data->numQuadPts, _data->numBasis, _data->cellDim,
			  _data->quadPts, _data->numQuadPts, _data->cellDim,
			  _data->quadWts, _data->numQuadPts,
			  _data->spaceDim);
  
  // Setup impulse amplitude
  spatialdata::spatialdb::SimpleDB* db = new spatialdata::spatialdb::SimpleDB("impulse amplitude");CPPUNIT_ASSERT(db);
  spatialdata::spatialdb::SimpleIOAscii ioImpulseAmp;
  ioImpulseAmp.filename(_data->impulseAmpFilename);
  db->ioHandler(&ioImpulseAmp);
  delete _dbImpulseAmp; _dbImpulseAmp = db;
  fault->dbImpulseAmp(db);

  CPPUNIT_ASSERT(_data->impulseDOF);
  fault->impulseDOF(_data->impulseDOF, _data->numComponents);

  PetscInt labelSize;
  PetscErrorCode err;
  err = DMGetStratumSize(mesh->dmMesh(), _data->label, 1, &labelSize);PYLITH_CHECK_ERROR(err);

  PetscInt firstFaultVertex = 0;
  PetscInt firstLagrangeVertex = labelSize;
  PetscInt firstFaultCell = labelSize;
  if (fault->useLagrangeConstraints())
    firstFaultCell += labelSize;
  fault->id(_data->id);
  fault->label(_data->label);
  fault->quadrature(_quadrature);
  
  fault->adjustTopology(mesh, &firstFaultVertex, &firstLagrangeVertex, &firstFaultCell);
  
  const PylithScalar upDir[3] = { 0.0, 0.0, 1.0 };
  
  fault->normalizer(normalizer);
  fault->initialize(*mesh, upDir);
  
  // Setup fields
  fields->add("residual", "residual");
  fields->add("disp(t)", "displacement");
  fields->add("dispIncr(t->t+dt)", "displacement_increment");
  fields->add("velocity(t)", "velocity");
  fields->solutionName("dispIncr(t->t+dt)");
  
  const int spaceDim = _data->spaceDim;
  topology::Field& residual = fields->get("residual");
  residual.subfieldAdd("displacement", spaceDim, topology::Field::VECTOR);
  residual.subfieldAdd("lagrange_multiplier", spaceDim, topology::Field::VECTOR);
  residual.subfieldsSetup();
  residual.setupSolnChart();
  residual.setupSolnDof(spaceDim);
  fault->setupSolnDof(&residual);
  residual.allocate();
  residual.zero();

  fields->copyLayout("residual");

  fault->verifyConfiguration(*mesh);

  PYLITH_METHOD_END;
} // _initialize

// ----------------------------------------------------------------------
// Determine if point is a Lagrange multiplier constraint point.
bool
pylith::faults::TestFaultCohesiveImpulses::_isConstraintEdge(const int point) const
{ // _isConstraintEdge
  PYLITH_METHOD_BEGIN;

  assert(_data);

  const int numConstraintEdges = _data->numConstraintEdges;
  bool isFound = false;
  for (int i=0; i < _data->numConstraintEdges; ++i)
    if (_data->constraintEdges[i] == point) {
      isFound = true;
      break;
    } // if
  PYLITH_METHOD_RETURN(isFound);
} // _isConstraintEdge


// End of file 
