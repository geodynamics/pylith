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

#include "TestPointForce.hh" // Implementation of class methods

#include "pylith/bc/PointForce.hh" // USES PointForce

#include "data/PointForceData.hh" // USES PointForceData

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii
#include "spatialdata/spatialdb/UniformDB.hh" // USES UniformDB

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestPointForce );

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::RealSection RealSection;
typedef pylith::topology::Mesh::RealUniformSection RealUniformSection;

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestPointForce::setUp(void)
{ // setUp
  _data = 0;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::bc::TestPointForce::tearDown(void)
{ // tearDown
  delete _data; _data = 0;
} // tearDown

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::bc::TestPointForce::testConstructor(void)
{ // testConstructor
  PointForce bc;
} // testConstructor

// ----------------------------------------------------------------------
// Test normalizer.
void
pylith::bc::TestPointForce::testNormalizer(void)
{ // testNormalizer
  PointForce bc;

  spatialdata::units::Nondimensional normalizer;
  const double scale = 4.0;
  normalizer.lengthScale(4.0);

  bc.normalizer(normalizer);
  CPPUNIT_ASSERT_EQUAL(scale, bc._getNormalizer().lengthScale());
} // testNormalizer

// ----------------------------------------------------------------------
// Test initialize().
void
pylith::bc::TestPointForce::testInitialize(void)
{ // testInitialize
  topology::Mesh mesh;
  PointForce bc;
  _initialize(&mesh, &bc);
  CPPUNIT_ASSERT(0 != _data);

  DM             dmMesh = mesh.dmMesh();
  PetscInt       cStart, cEnd;
  PetscErrorCode err;

  CPPUNIT_ASSERT(dmMesh);
  err = DMPlexGetHeightStratum(dmMesh, 0, &cStart, &cEnd);CHECK_PETSC_ERROR(err);
  const int numCells = cEnd-cStart;
  const int numForceDOF = _data->numForceDOF;
  const size_t numPoints = _data->numForcePts;

  // Check points
  const int offset = numCells;
  if (numForceDOF > 0) {
    CPPUNIT_ASSERT_EQUAL(numPoints, bc._points.size());
    for (int i=0; i < numPoints; ++i)
      CPPUNIT_ASSERT_EQUAL(_data->forcePoints[i]+offset, bc._points[i]);
  } // if

  CPPUNIT_ASSERT(0 != bc._parameters);
  PetscSection initialSection = bc._parameters->get("initial").petscSection();
  Vec          initialVec     = bc._parameters->get("initial").localVector();
  PetscScalar *initialArray;
  CPPUNIT_ASSERT(initialSection);CPPUNIT_ASSERT(initialVec);
  
  // Check values
  const PylithScalar tolerance = 1.0e-06;
  err = VecGetArray(initialVec, &initialArray);CHECK_PETSC_ERROR(err);
  for (int i=0; i < numPoints; ++i) {
    const int p_force = _data->forcePoints[i]+offset;
    PetscInt  dof, off;

    err = PetscSectionGetDof(initialSection, p_force, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(initialSection, p_force, &off);CHECK_PETSC_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(numForceDOF, dof);
    for (int iDOF=0; iDOF < numForceDOF; ++iDOF) 
      CPPUNIT_ASSERT_DOUBLES_EQUAL(_data->forceInitial[i*numForceDOF+iDOF], initialArray[off+iDOF], tolerance);
  } // for
  err = VecRestoreArray(initialVec, &initialArray);CHECK_PETSC_ERROR(err);

  // Check rate of change
  PetscSection rateSection = bc._parameters->get("rate").petscSection();
  Vec          rateVec     = bc._parameters->get("rate").localVector();
  PetscScalar *rateArray;
  CPPUNIT_ASSERT(rateSection);CPPUNIT_ASSERT(rateVec);

  err = VecGetArray(rateVec, &rateArray);CHECK_PETSC_ERROR(err);
  for (int i=0; i < numPoints; ++i) {
    const int p_force = _data->forcePoints[i]+offset;
    PetscInt  dof, off;

    err = PetscSectionGetDof(rateSection, p_force, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(rateSection, p_force, &off);CHECK_PETSC_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(numForceDOF, dof);
    for (int iDOF=0; iDOF < numForceDOF; ++iDOF) 
      CPPUNIT_ASSERT_DOUBLES_EQUAL(_data->forceRate, rateArray[off+iDOF], tolerance);
  } // for
  err = VecRestoreArray(rateVec, &rateArray);CHECK_PETSC_ERROR(err);
} // testInitialize

// ----------------------------------------------------------------------
// Test integrateResidual().
void
pylith::bc::TestPointForce::testIntegrateResidual(void)
{ // testIntegrateResidual
  topology::Mesh mesh;
  PointForce bc;
  _initialize(&mesh, &bc);

  topology::Field<topology::Mesh> residual(mesh);
  const spatialdata::geocoords::CoordSys* cs = mesh.coordsys();
  CPPUNIT_ASSERT(0 != cs);
  const int spaceDim = cs->spaceDim();
  residual.newSection(topology::FieldBase::VERTICES_FIELD, spaceDim);
  residual.allocate();
  residual.zero();

  topology::SolutionFields fields(mesh);

  const PylithScalar t = _data->tResidual;
  bc.integrateResidual(residual, t, &fields);

  DM             dmMesh = mesh.dmMesh();
  PetscInt       vStart, vEnd;
  PetscErrorCode err;

  CPPUNIT_ASSERT(dmMesh);
  err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);

  const PylithScalar* valsE = _data->residual;
  const int totalNumVertices = vEnd-vStart;
  const int sizeE = spaceDim * totalNumVertices;

  PetscSection residualSection = residual.petscSection();
  Vec          residualVec     = residual.localVector();
  PetscScalar *vals;
  PetscInt     size;

  CPPUNIT_ASSERT(residualSection);
  err = PetscSectionGetStorageSize(residualSection, &size);CHECK_PETSC_ERROR(err);
  CPPUNIT_ASSERT_EQUAL(sizeE, size);

  //residual.view("RESIDUAL");

  const PylithScalar tolerance = 1.0e-06;
  err = VecGetArray(residualVec, &vals);CHECK_PETSC_ERROR(err);
  for (int i=0; i < size; ++i)
    if (fabs(valsE[i]) > 1.0)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, vals[i]/valsE[i], tolerance);
    else
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valsE[i], vals[i], tolerance);
  err = VecRestoreArray(residualVec, &vals);CHECK_PETSC_ERROR(err);
} // testIntegrateResidual

// ----------------------------------------------------------------------
// Test verifyConfiguration().
void
pylith::bc::TestPointForce::testVerifyConfiguration(void)
{ // testVerifyConfiguration
  topology::Mesh mesh;
  PointForce bc;
  _initialize(&mesh, &bc);

  bc.verifyConfiguration(mesh);
} // testVerifyConfiguration

// ----------------------------------------------------------------------
void
pylith::bc::TestPointForce::_initialize(topology::Mesh* mesh,
					 PointForce* const bc) const
{ // _initialize
  CPPUNIT_ASSERT(0 != _data);
  CPPUNIT_ASSERT(0 != bc);

  meshio::MeshIOAscii iohandler;
  iohandler.filename(_data->meshFilename);
  iohandler.read(mesh);

  spatialdata::geocoords::CSCart cs;
  spatialdata::units::Nondimensional normalizer;
  cs.setSpaceDim(mesh->dimension());
  cs.initialize();
  mesh->coordsys(&cs);
  mesh->nondimensionalize(normalizer);

  spatialdata::spatialdb::SimpleDB dbInitial("TestPointForce initial");
  spatialdata::spatialdb::SimpleIOAscii dbInitialIO;
  dbInitialIO.filename(_data->dbFilename);
  dbInitial.ioHandler(&dbInitialIO);
  dbInitial.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

  spatialdata::spatialdb::UniformDB dbRate("TestPointForce rate");
  { // rate db
    const int numValues = 4;
    const char* names[numValues] = {
      "force-rate-x",
      "force-rate-y",
      "force-rate-z",
      "rate-start-time",
    };
    const char* units[numValues] = {
      "newton",
      "newton",
      "newton",
      "s",
    };
    const double values[numValues] = { _data->forceRate,
				       _data->forceRate,
				       _data->forceRate,
				       _data->tRef};
    dbRate.setData(names, units, values, numValues);
  } // rate db

  const PylithScalar upDir[] = { 0.0, 0.0, 1.0 };

  bc->label(_data->label);
  bc->dbInitial(&dbInitial);
  bc->dbRate(&dbRate);
  bc->bcDOF(_data->forceDOF, _data->numForceDOF);
  bc->initialize(*mesh, upDir);
} // _initialize


// End of file 
