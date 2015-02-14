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

#include "TestDirichletBCMulti.hh" // Implementation of class methods

#include "pylith/bc/DirichletBC.hh" // USES DirichletBC

#include "data/DirichletDataMulti.hh" // USES DirichletData

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps::nondimensionalize()
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorMesh.hh" // USES VisitorMesh
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestDirichletBCMulti::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  _data = 0;

  PYLITH_METHOD_END;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::bc::TestDirichletBCMulti::tearDown(void)
{ // tearDown
  PYLITH_METHOD_BEGIN;

  delete _data; _data = 0;

  PYLITH_METHOD_END;
} // tearDown

// ----------------------------------------------------------------------
// Test setConstraintSizes().
void
pylith::bc::TestDirichletBCMulti::testSetConstraintSizes(void)
{ // testSetConstraintSizes
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_data);

  topology::Mesh mesh;
  DirichletBC bcA;
  DirichletBC bcB;
  DirichletBC bcC;
  _initialize(&mesh, &bcA, &bcB, &bcC);

  const int fiberDim = _data->numDOF;
  topology::Field field(mesh);
  field.subfieldAdd("bc", fiberDim, topology::Field::VECTOR);
  field.subfieldsSetup();
  field.newSection(pylith::topology::FieldBase::VERTICES_FIELD, fiberDim);
  field.subfieldSetDof("bc", pylith::topology::FieldBase::VERTICES_FIELD, fiberDim);
  bcA.setConstraintSizes(field);
  bcB.setConstraintSizes(field);
  bcC.setConstraintSizes(field);
  field.allocate();

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);

  // Vertices
  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();

  // Cells
  topology::Stratum cellsStratum(dmMesh, topology::Stratum::HEIGHT, 0);
  const PetscInt numCells = cellsStratum.size();
  const PetscInt offset = numCells;

  topology::VecVisitorMesh fieldVisitor(field);

  int iConstraint = 0;
  for(PetscInt v = vStart; v < vEnd; ++v) {
    const PetscInt dof = fieldVisitor.sectionDof(v);
    const PetscInt cdof = fieldVisitor.sectionConstraintDof(v);
    CPPUNIT_ASSERT_EQUAL(_data->numDOF, dof);
    CPPUNIT_ASSERT_EQUAL(_data->constraintSizes[v-offset], cdof);
  } // for

  PYLITH_METHOD_END;
} // testSetConstraintSizes

// ----------------------------------------------------------------------
// Test setConstraints().
void
pylith::bc::TestDirichletBCMulti::testSetConstraints(void)
{ // testSetConstraints
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_data);

  topology::Mesh mesh;
  DirichletBC bcA;
  DirichletBC bcB;
  DirichletBC bcC;
  _initialize(&mesh, &bcA, &bcB, &bcC);

  const int fiberDim = _data->numDOF;
  topology::Field field(mesh);
  field.subfieldAdd("bc", fiberDim, topology::Field::VECTOR);
  field.subfieldsSetup();
  field.newSection(pylith::topology::FieldBase::VERTICES_FIELD, fiberDim);
  field.subfieldSetDof("bc", pylith::topology::FieldBase::VERTICES_FIELD, fiberDim);
  bcA.setConstraintSizes(field);
  bcB.setConstraintSizes(field);
  bcC.setConstraintSizes(field);
  field.allocate();
  bcA.setConstraints(field);
  bcB.setConstraints(field);
  bcC.setConstraints(field);

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);

  // Vertices
  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();

  // Cells
  topology::Stratum cellsStratum(dmMesh, topology::Stratum::HEIGHT, 0);
  const PetscInt numCells = cellsStratum.size();

  PetscSection fieldSection = field.localSection();CPPUNIT_ASSERT(fieldSection);
  const PetscInt offset = numCells;

  PetscErrorCode err = 0;
  int index = 0;
  int iConstraint = 0;
  for(PetscInt v = vStart; v < vEnd; ++v) {
    const int numConstrainedDOF = _data->constraintSizes[v-offset];
    const PetscInt *cInd;
    PetscInt dof, cdof;

    err = PetscSectionGetDof(fieldSection, v, &dof);PYLITH_CHECK_ERROR(err);
    err = PetscSectionGetConstraintDof(fieldSection, v, &cdof);PYLITH_CHECK_ERROR(err);
    err = PetscSectionGetConstraintIndices(fieldSection, v, &cInd);PYLITH_CHECK_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(numConstrainedDOF, cdof);
    if (numConstrainedDOF > 0) {
      for (int iDOF=0; iDOF < numConstrainedDOF; ++iDOF)
        CPPUNIT_ASSERT_EQUAL(_data->constrainedDOF[index++], cInd[iDOF]);
    } // if
  } // for

  PYLITH_METHOD_END;
} // testSetConstraints

// ----------------------------------------------------------------------
// Test setField().
void
pylith::bc::TestDirichletBCMulti::testSetField(void)
{ // testSetField
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_data);

  topology::Mesh mesh;
  DirichletBC bcA;
  DirichletBC bcB;
  DirichletBC bcC;
  _initialize(&mesh, &bcA, &bcB, &bcC);

  const int fiberDim = _data->numDOF;
  topology::Field field(mesh);
  field.subfieldAdd("bc", fiberDim, topology::Field::VECTOR);
  field.subfieldsSetup();
  field.newSection(pylith::topology::FieldBase::VERTICES_FIELD, fiberDim);
  field.subfieldSetDof("bc", pylith::topology::FieldBase::VERTICES_FIELD, fiberDim);
  bcA.setConstraintSizes(field);
  bcB.setConstraintSizes(field);
  bcC.setConstraintSizes(field);
  field.allocate();
  bcA.setConstraints(field);
  bcB.setConstraints(field);
  bcC.setConstraints(field);

  // Scales
  const PylithScalar tolerance = 1.0e-06;
  const PylithScalar valueScale = _data->lengthScale;

  // All values should be zero.
  field.zero();

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);

  // Vertices
  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();

  // Cells
  topology::Stratum cellsStratum(dmMesh, topology::Stratum::HEIGHT, 0);
  const PetscInt numCells = cellsStratum.size();
  const PetscInt offset = numCells;

  topology::VecVisitorMesh fieldVisitor(field);
  PetscScalar* fieldArray = fieldVisitor.localArray();
  for(PetscInt v = vStart; v < vEnd; ++v) {
    const PetscInt off = fieldVisitor.sectionOffset(v);
    CPPUNIT_ASSERT_EQUAL(fiberDim, fieldVisitor.sectionDof(v));

    for(int d = 0; d < fiberDim; ++d)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, fieldArray[off+d], tolerance);
  } // for
  fieldVisitor.clear();

  // Only unconstrained values should be zero.
  // Expected values set in _data->field
  const PylithScalar t = 10.0/_data->timeScale;
  bcA.setField(t, field);
  bcB.setField(t, field);
  bcC.setField(t, field);

  fieldVisitor.initialize(field);
  fieldArray = fieldVisitor.localArray();
  int i = 0;
  for(PetscInt v = vStart; v < vEnd; ++v) {
    const PetscInt off = fieldVisitor.sectionOffset(v);
    CPPUNIT_ASSERT_EQUAL(fiberDim, fieldVisitor.sectionDof(v));

    for(int iDOF = 0; iDOF < fiberDim; ++iDOF)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(_data->field[i++], fieldArray[off+iDOF]*valueScale, tolerance);
  } // for

  PYLITH_METHOD_END;
} // testSetField

// ----------------------------------------------------------------------
// Test setFieldIncr().
void
pylith::bc::TestDirichletBCMulti::testSetFieldIncr(void)
{ // testSetFieldIncr
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_data);

  topology::Mesh mesh;
  DirichletBC bcA;
  DirichletBC bcB;
  DirichletBC bcC;
  _initialize(&mesh, &bcA, &bcB, &bcC);

  const int fiberDim = _data->numDOF;
  topology::Field field(mesh);
  field.subfieldAdd("bc", fiberDim, topology::Field::VECTOR);
  field.subfieldsSetup();
  field.newSection(pylith::topology::FieldBase::VERTICES_FIELD, fiberDim);
  field.subfieldSetDof("bc", pylith::topology::FieldBase::VERTICES_FIELD, fiberDim);
  bcA.setConstraintSizes(field);
  bcB.setConstraintSizes(field);
  bcC.setConstraintSizes(field);
  field.allocate();
  bcA.setConstraints(field);
  bcB.setConstraints(field);
  bcC.setConstraints(field);

  // Scales
  const PylithScalar tolerance = 1.0e-06;
  const PylithScalar valueScale = _data->lengthScale;

  // All values should be zero.
  field.zero();

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);

  // Vertices
  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();

  // Cells
  topology::Stratum cellsStratum(dmMesh, topology::Stratum::HEIGHT, 0);
  const PetscInt numCells = cellsStratum.size();
  const PetscInt offset = numCells;

  topology::VecVisitorMesh fieldVisitor(field);
  PetscScalar* fieldArray = fieldVisitor.localArray();
  for(PetscInt v = vStart; v < vEnd; ++v) {
    const PetscInt off = fieldVisitor.sectionOffset(v);
    CPPUNIT_ASSERT_EQUAL(fiberDim, fieldVisitor.sectionDof(v));

    for(int d = 0; d < fiberDim; ++d)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, fieldArray[off+d], tolerance);
  } // for
  fieldVisitor.clear();

  // Only unconstrained values should be zero.
  // Expected values set in _data->field
  const PylithScalar t0 = 10.0/_data->timeScale;
  const PylithScalar t1 = 14.0/_data->timeScale;
  bcA.setFieldIncr(t0, t1, field);
  bcB.setFieldIncr(t0, t1, field);
  bcC.setFieldIncr(t0, t1, field);

  fieldVisitor.initialize(field);
  fieldArray = fieldVisitor.localArray();
  int i = 0;
  for(PetscInt v = vStart; v < vEnd; ++v) {
    const PetscInt off = fieldVisitor.sectionOffset(v);
    CPPUNIT_ASSERT_EQUAL(fiberDim, fieldVisitor.sectionDof(v));

    for(int iDOF = 0; iDOF < fiberDim; ++iDOF)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(_data->fieldIncr[i++], fieldArray[off+iDOF]*valueScale, tolerance);
  } // for

  PYLITH_METHOD_END;
} // testSetFieldIncr

// ----------------------------------------------------------------------
void
pylith::bc::TestDirichletBCMulti::_initialize(topology::Mesh* mesh,
					      DirichletBC* const bcA,
					      DirichletBC* const bcB,
					      DirichletBC* const bcC) const
{ // _initialize
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT(bcA);
  CPPUNIT_ASSERT(bcB);
  CPPUNIT_ASSERT(bcC);

  meshio::MeshIOAscii iohandler;
  iohandler.filename(_data->meshFilename);
  iohandler.read(mesh);

  spatialdata::geocoords::CSCart cs;
  cs.setSpaceDim(mesh->dimension());
  cs.initialize();
  mesh->coordsys(&cs);

  spatialdata::units::Nondimensional normalizer;
  normalizer.lengthScale(_data->lengthScale);
  normalizer.pressureScale(_data->pressureScale);
  normalizer.densityScale(_data->densityScale);
  normalizer.timeScale(_data->timeScale);
  topology::MeshOps::nondimensionalize(mesh, normalizer);

  // Setup boundary condition A
  spatialdata::spatialdb::SimpleDB db("TestDirichletBCMulti initial");
  spatialdata::spatialdb::SimpleIOAscii dbIO;
  dbIO.filename(_data->dbFilenameA);
  db.ioHandler(&dbIO);

  spatialdata::spatialdb::SimpleDB dbRate("TestDirichletBCMulti rate");
  spatialdata::spatialdb::SimpleIOAscii dbIORate;
  dbIORate.filename(_data->dbFilenameARate);
  dbRate.ioHandler(&dbIORate);

  const PylithScalar upDir[] = { 0.0, 0.0, 1.0 };

  bcA->label(_data->labelA);
  bcA->dbInitial(&db);
  bcA->dbRate(&dbRate);
  bcA->bcDOF(_data->fixedDOFA, _data->numFixedDOFA);
  bcA->normalizer(normalizer);
  bcA->initialize(*mesh, upDir);

  // Setup boundary condition B
  dbIO.filename(_data->dbFilenameB);
  db.ioHandler(&dbIO);

  dbIORate.filename(_data->dbFilenameBRate);
  dbRate.ioHandler(&dbIORate);

  bcB->label(_data->labelB);
  bcB->dbInitial(&db);
  bcB->dbRate(&dbRate);
  bcB->bcDOF(_data->fixedDOFB, _data->numFixedDOFB);
  bcB->normalizer(normalizer);
  bcB->initialize(*mesh, upDir);

  // Setup boundary condition C
  if (_data->numFixedDOFC > 0.0) {
    dbIO.filename(_data->dbFilenameC);
    db.ioHandler(&dbIO);
    
    dbIORate.filename(_data->dbFilenameCRate);
    dbRate.ioHandler(&dbIORate);
    
    bcC->label(_data->labelC);
    bcC->dbInitial(&db);
    bcC->dbRate(&dbRate);
    bcC->bcDOF(_data->fixedDOFC, _data->numFixedDOFC);
    bcC->normalizer(normalizer);
    bcC->initialize(*mesh, upDir);
  } // if

  PYLITH_METHOD_END;
} // _initialize


// End of file 
