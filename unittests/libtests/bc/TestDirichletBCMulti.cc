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

#include "TestDirichletBCMulti.hh" // Implementation of class methods

#include "pylith/bc/DirichletBC.hh" // USES DirichletBC

#include "data/DirichletDataMulti.hh" // USES DirichletData

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestDirichletBCMulti::setUp(void)
{ // setUp
  _data = 0;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::bc::TestDirichletBCMulti::tearDown(void)
{ // tearDown
  delete _data; _data = 0;
} // tearDown

// ----------------------------------------------------------------------
// Test setConstraintSizes().
void
pylith::bc::TestDirichletBCMulti::testSetConstraintSizes(void)
{ // testSetConstraintSizes
  topology::Mesh mesh;
  DirichletBC bcA;
  DirichletBC bcB;
  DirichletBC bcC;
  _initialize(&mesh, &bcA, &bcB, &bcC);
  CPPUNIT_ASSERT(0 != _data);

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  PetscInt cStart, cEnd, vStart, vEnd;
  PetscErrorCode err = 0;

  err = DMPlexGetHeightStratum(dmMesh, 0, &cStart, &cEnd);CHECK_PETSC_ERROR(err);
  err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);
  
  const int fiberDim = _data->numDOF;
  topology::Field<topology::Mesh> field(mesh);
  field.addField("bc", fiberDim);
  field.setupFields();
  field.newSection(pylith::topology::FieldBase::VERTICES_FIELD, fiberDim);
  field.updateDof("bc", pylith::topology::FieldBase::VERTICES_FIELD, fiberDim);
  bcA.setConstraintSizes(field);
  bcB.setConstraintSizes(field);
  bcC.setConstraintSizes(field);
  field.allocate();

  PetscSection fieldSection = field.petscSection();CPPUNIT_ASSERT(fieldSection);
  PetscVec fieldVec = field.localVector();CPPUNIT_ASSERT(fieldVec);

  const PetscInt numCells = cEnd - cStart;
  const PetscInt offset = numCells;
  int iConstraint = 0;
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt dof, cdof;

    err = PetscSectionGetDof(fieldSection, v, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetConstraintDof(fieldSection, v, &cdof);CHECK_PETSC_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(_data->numDOF, dof);
    CPPUNIT_ASSERT_EQUAL(_data->constraintSizes[v-offset], cdof);
  } // for
} // testSetConstraintSizes

// ----------------------------------------------------------------------
// Test setConstraints().
void
pylith::bc::TestDirichletBCMulti::testSetConstraints(void)
{ // testSetConstraints
  topology::Mesh mesh;
  DirichletBC bcA;
  DirichletBC bcB;
  DirichletBC bcC;
  _initialize(&mesh, &bcA, &bcB, &bcC);
  CPPUNIT_ASSERT(0 != _data);

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  PetscInt cStart, cEnd, vStart, vEnd;
  PetscErrorCode err = 0;

  err = DMPlexGetHeightStratum(dmMesh, 0, &cStart, &cEnd);CHECK_PETSC_ERROR(err);
  err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);
  
  const int fiberDim = _data->numDOF;
  topology::Field<topology::Mesh> field(mesh);
  field.addField("bc", fiberDim);
  field.setupFields();
  field.newSection(pylith::topology::FieldBase::VERTICES_FIELD, fiberDim);
  field.updateDof("bc", pylith::topology::FieldBase::VERTICES_FIELD, fiberDim);
  bcA.setConstraintSizes(field);
  bcB.setConstraintSizes(field);
  bcC.setConstraintSizes(field);
  field.allocate();
  bcA.setConstraints(field);
  bcB.setConstraints(field);
  bcC.setConstraints(field);

  PetscSection fieldSection = field.petscSection();CPPUNIT_ASSERT(fieldSection);
  PetscVec fieldVec     = field.localVector();CPPUNIT_ASSERT(fieldVec);

  const PetscInt numCells = cEnd - cStart;
  const PetscInt offset = numCells;
  int index = 0;
  int iConstraint = 0;
  for(PetscInt v = vStart; v < vEnd; ++v) {
    const int numConstrainedDOF = _data->constraintSizes[v-offset];
    const PetscInt *cInd;
    PetscInt dof, cdof;

    err = PetscSectionGetDof(fieldSection, v, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetConstraintDof(fieldSection, v, &cdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetConstraintIndices(fieldSection, v, &cInd);CHECK_PETSC_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(numConstrainedDOF, cdof);
    if (numConstrainedDOF > 0) {
      for (int iDOF=0; iDOF < numConstrainedDOF; ++iDOF)
        CPPUNIT_ASSERT_EQUAL(_data->constrainedDOF[index++], cInd[iDOF]);
    } // if
  } // for
} // testSetConstraints

// ----------------------------------------------------------------------
// Test setField().
void
pylith::bc::TestDirichletBCMulti::testSetField(void)
{ // testSetField
  topology::Mesh mesh;
  DirichletBC bcA;
  DirichletBC bcB;
  DirichletBC bcC;
  _initialize(&mesh, &bcA, &bcB, &bcC);
  CPPUNIT_ASSERT(0 != _data);

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  PetscInt cStart, cEnd, vStart, vEnd;
  PetscErrorCode err = 0;

  err = DMPlexGetHeightStratum(dmMesh, 0, &cStart, &cEnd);CHECK_PETSC_ERROR(err);
  err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);
  
  const int fiberDim = _data->numDOF;
  topology::Field<topology::Mesh> field(mesh);
  field.addField("bc", fiberDim);
  field.setupFields();
  field.newSection(pylith::topology::FieldBase::VERTICES_FIELD, fiberDim);
  field.updateDof("bc", pylith::topology::FieldBase::VERTICES_FIELD, fiberDim);
  bcA.setConstraintSizes(field);
  bcB.setConstraintSizes(field);
  bcC.setConstraintSizes(field);
  field.allocate();
  bcA.setConstraints(field);
  bcB.setConstraints(field);
  bcC.setConstraints(field);

  PetscSection fieldSection = field.petscSection();CPPUNIT_ASSERT(fieldSection);
  PetscVec fieldVec = field.localVector();CPPUNIT_ASSERT(fieldVec);

  const PylithScalar tolerance = 1.0e-06;
  const PylithScalar valueScale = _data->lengthScale;

  // All values should be zero.
  PetscScalar *values;
  field.zero();
  err = VecGetArray(fieldVec, &values);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt dof, off;

    err = PetscSectionGetDof(fieldSection, v, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(fieldSection, v, &off);CHECK_PETSC_ERROR(err);
    for(int d = 0; d < dof; ++d)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, values[off+d], tolerance);
  } // for
  err = VecRestoreArray(fieldVec, &values);CHECK_PETSC_ERROR(err);

  // Only unconstrained values should be zero.
  // Expected values set in _data->field
  const PylithScalar t = 10.0/_data->timeScale;
  bcA.setField(t, field);
  bcB.setField(t, field);
  bcC.setField(t, field);

  int i = 0;
  err = VecGetArray(fieldVec, &values);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt dof, cdof, off;

    err = PetscSectionGetDof(fieldSection, v, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(fieldSection, v, &off);CHECK_PETSC_ERROR(err);
    for(int iDOF = 0; iDOF < dof; ++iDOF)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(_data->field[i++], values[off+iDOF]*valueScale, tolerance);
  } // for
  err = VecRestoreArray(fieldVec, &values);CHECK_PETSC_ERROR(err);
} // testSetField

// ----------------------------------------------------------------------
// Test setFieldIncr().
void
pylith::bc::TestDirichletBCMulti::testSetFieldIncr(void)
{ // testSetFieldIncr
  topology::Mesh mesh;
  DirichletBC bcA;
  DirichletBC bcB;
  DirichletBC bcC;
  _initialize(&mesh, &bcA, &bcB, &bcC);
  CPPUNIT_ASSERT(0 != _data);

  PetscDM dmMesh = mesh.dmMesh();
  CPPUNIT_ASSERT(dmMesh);
  PetscInt       cStart, cEnd, vStart, vEnd;
  PetscErrorCode err;

  err = DMPlexGetHeightStratum(dmMesh, 0, &cStart, &cEnd);CHECK_PETSC_ERROR(err);
  err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);
  
  const int fiberDim = _data->numDOF;
  topology::Field<topology::Mesh> field(mesh);
  field.addField("bc", fiberDim);
  field.setupFields();
  field.newSection(pylith::topology::FieldBase::VERTICES_FIELD, fiberDim);
  field.updateDof("bc", pylith::topology::FieldBase::VERTICES_FIELD, fiberDim);
  bcA.setConstraintSizes(field);
  bcB.setConstraintSizes(field);
  bcC.setConstraintSizes(field);
  field.allocate();
  bcA.setConstraints(field);
  bcB.setConstraints(field);
  bcC.setConstraints(field);

  PetscSection fieldSection = field.petscSection();
  PetscVec          fieldVec     = field.localVector();
  CPPUNIT_ASSERT(fieldSection);CPPUNIT_ASSERT(fieldVec);

  const PylithScalar tolerance = 1.0e-06;
  const PylithScalar valueScale = _data->lengthScale;

  // All values should be zero.
  PetscScalar *values;
  field.zero();
  err = VecGetArray(fieldVec, &values);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt dof, off;

    err = PetscSectionGetDof(fieldSection, v, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(fieldSection, v, &off);CHECK_PETSC_ERROR(err);
    for(int d = 0; d < dof; ++d)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, values[off+d], tolerance);
  } // for
  err = VecRestoreArray(fieldVec, &values);CHECK_PETSC_ERROR(err);

  // Only unconstrained values should be zero.
  // Expected values set in _data->field
  const PylithScalar t0 = 10.0/_data->timeScale;
  const PylithScalar t1 = 14.0/_data->timeScale;
  bcA.setFieldIncr(t0, t1, field);
  bcB.setFieldIncr(t0, t1, field);
  bcC.setFieldIncr(t0, t1, field);

  int i = 0;
  err = VecGetArray(fieldVec, &values);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt dof, cdof, off;

    err = PetscSectionGetDof(fieldSection, v, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(fieldSection, v, &off);CHECK_PETSC_ERROR(err);
    for(int iDOF = 0; iDOF < dof; ++iDOF) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(_data->fieldIncr[i++], values[off+iDOF]*valueScale, tolerance);
    } // for
  } // for
  err = VecRestoreArray(fieldVec, &values);CHECK_PETSC_ERROR(err);
} // testSetFieldIncr

// ----------------------------------------------------------------------
void
pylith::bc::TestDirichletBCMulti::_initialize(topology::Mesh* mesh,
					      DirichletBC* const bcA,
					      DirichletBC* const bcB,
					      DirichletBC* const bcC) const
{ // _initialize
  CPPUNIT_ASSERT(0 != _data);
  CPPUNIT_ASSERT(0 != bcA);
  CPPUNIT_ASSERT(0 != bcB);
  CPPUNIT_ASSERT(0 != bcC);

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
  mesh->nondimensionalize(normalizer);

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
} // _initialize


// End of file 
