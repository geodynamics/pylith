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

#include "TestDirichletBC.hh" // Implementation of class methods

#include "pylith/bc/DirichletBC.hh" // USES DirichletBC

#include "data/DirichletData.hh" // USES DirichletData

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii
#include "spatialdata/spatialdb/UniformDB.hh" // USES UniformDB

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestDirichletBC );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestDirichletBC::setUp(void)
{ // setUp
  _data = 0;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::bc::TestDirichletBC::tearDown(void)
{ // tearDown
  delete _data; _data = 0;
} // tearDown

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::bc::TestDirichletBC::testConstructor(void)
{ // testConstructor
  DirichletBC bc;
} // testConstructor

// ----------------------------------------------------------------------
// Test initialize().
void
pylith::bc::TestDirichletBC::testInitialize(void)
{ // testInitialize
  CPPUNIT_ASSERT(_data);

  topology::Mesh mesh;
  DirichletBC bc;
  _initialize(&mesh, &bc);

  const PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  topology::Stratum heightStratum(dmMesh, topology::Stratum::HEIGHT, 0);
  const PetscInt numCells = heightStratum.size();

  const int numFixedDOF = _data->numFixedDOF;
  const size_t numPoints = _data->numConstrainedPts;

  // Check points
  const int offset = numCells;
  if (numFixedDOF > 0) {
    CPPUNIT_ASSERT_EQUAL(numPoints, bc._points.size());
    for (int i=0; i < numPoints; ++i)
      CPPUNIT_ASSERT_EQUAL(_data->constrainedPoints[i]+offset, bc._points[i]);
  } // if

  if (numFixedDOF > 0) {
    // Check values
    CPPUNIT_ASSERT(0 != bc._parameters);
    PetscSection initialSection = bc._parameters->get("initial").petscSection();
    PetscVec initialVec = bc._parameters->get("initial").localVector();
    PetscScalar* initialArray;
    PetscErrorCode err = 0;
    CPPUNIT_ASSERT(initialSection);CPPUNIT_ASSERT(initialVec);

    const PylithScalar tolerance = 1.0e-06;
    const PylithScalar dispScale = _data->lengthScale;
    const PylithScalar velocityScale = _data->lengthScale / _data->timeScale;
    err = VecGetArray(initialVec, &initialArray);CHECK_PETSC_ERROR(err);
    for (int i=0; i < numPoints; ++i) {
      const PetscInt p_value = _data->constrainedPoints[i]+offset;
      PetscInt dof, off;

      err = PetscSectionGetDof(initialSection, p_value, &dof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(initialSection, p_value, &off);CHECK_PETSC_ERROR(err);
      CPPUNIT_ASSERT_EQUAL(numFixedDOF, dof);
      for(int iDOF = 0; iDOF < numFixedDOF; ++iDOF) {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(_data->valuesInitial[i*numFixedDOF+iDOF]/dispScale, initialArray[off+iDOF], tolerance);
      } // for
    } // for
    err = VecRestoreArray(initialVec, &initialArray);CHECK_PETSC_ERROR(err);
    
    // Check rate of change
    PetscSection rateSection = bc._parameters->get("rate").petscSection();
    PetscVec rateVec = bc._parameters->get("rate").localVector();
    PetscScalar *rateArray;
    
    err = VecGetArray(rateVec, &rateArray);CHECK_PETSC_ERROR(err);
    for (int i=0; i < numPoints; ++i) {
      const PetscInt p_value = _data->constrainedPoints[i]+offset;
      PetscInt dof, off;

      err = PetscSectionGetDof(rateSection, p_value, &dof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(rateSection, p_value, &off);CHECK_PETSC_ERROR(err);
      CPPUNIT_ASSERT_EQUAL(numFixedDOF, dof);
      for(int iDOF = 0; iDOF < numFixedDOF; ++iDOF) 
        CPPUNIT_ASSERT_DOUBLES_EQUAL(_data->valueRate/velocityScale, rateArray[off+iDOF], tolerance);
    } // for
    err = VecRestoreArray(rateVec, &rateArray);CHECK_PETSC_ERROR(err);
  } // if
} // testInitialize

// ----------------------------------------------------------------------
// Test numDimConstrained().
void
pylith::bc::TestDirichletBC::testNumDimConstrained(void)
{ // testNumDimConstrained
  topology::Mesh mesh;
  DirichletBC bc;
  _initialize(&mesh, &bc);
  CPPUNIT_ASSERT(0 != _data);

  CPPUNIT_ASSERT_EQUAL(_data->numFixedDOF, bc.numDimConstrained());
} // testNumDimConstrained

// ----------------------------------------------------------------------
// Test setConstraintSizes().
void
pylith::bc::TestDirichletBC::testSetConstraintSizes(void)
{ // testSetConstraintSizes
  topology::Mesh mesh;
  DirichletBC bc;
  _initialize(&mesh, &bc);
  CPPUNIT_ASSERT(0 != _data);

  PetscDM dmMesh = mesh.dmMesh();
  CPPUNIT_ASSERT(dmMesh);
  PetscInt cStart, cEnd, vStart, vEnd;
  PetscErrorCode err = 0;

  err = DMPlexGetHeightStratum(dmMesh, 0, &cStart, &cEnd);CHECK_PETSC_ERROR(err);
  err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);
  
  const int fiberDim = _data->numDOF;
  const int spaceDim = mesh.dimension();
  topology::Field<topology::Mesh> field(mesh);
  field.addField("bc", spaceDim);
  field.setupFields();
  field.newSection(pylith::topology::FieldBase::VERTICES_FIELD, fiberDim);
  field.updateDof("bc", pylith::topology::FieldBase::VERTICES_FIELD, fiberDim);
  bc.setConstraintSizes(field); // Does not handle fields right now
  field.allocate();

  PetscSection fieldSection = field.petscSection();
  PetscVec fieldVec     = field.localVector();
  CPPUNIT_ASSERT(fieldSection);CPPUNIT_ASSERT(fieldVec);

  const PetscInt numCells = cEnd - cStart;
  const PetscInt offset   = numCells;
  int iConstraint = 0;
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt dof, cdof, fdof, fcdof;

    err = PetscSectionGetDof(fieldSection, v, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetConstraintDof(fieldSection, v, &cdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetFieldDof(fieldSection, v, 0, &fdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetFieldConstraintDof(fieldSection, v, 0, &fcdof);CHECK_PETSC_ERROR(err);
    if (v != _data->constrainedPoints[iConstraint] + offset) {
      CPPUNIT_ASSERT_EQUAL(_data->numDOF, dof);
      CPPUNIT_ASSERT_EQUAL(0, cdof);
      CPPUNIT_ASSERT_EQUAL(_data->numDOF, fdof);
      CPPUNIT_ASSERT_EQUAL(0, fcdof);
    } else {
      CPPUNIT_ASSERT_EQUAL(_data->numDOF, dof);
      CPPUNIT_ASSERT_EQUAL(_data->numFixedDOF, cdof);
      CPPUNIT_ASSERT_EQUAL(_data->numDOF, fdof);
      //CPPUNIT_ASSERT_EQUAL(_data->numFixedDOF, fcdof);
      ++iConstraint;
    } // if/else
  } // for
} // testSetConstraintSizes

// ----------------------------------------------------------------------
// Test setConstraints().
void
pylith::bc::TestDirichletBC::testSetConstraints(void)
{ // testSetConstraints
  topology::Mesh mesh;
  DirichletBC bc;
  _initialize(&mesh, &bc);
  CPPUNIT_ASSERT(0 != _data);

  PetscDM dmMesh = mesh.dmMesh();
  CPPUNIT_ASSERT(dmMesh);
  PetscInt cStart, cEnd, vStart, vEnd;
  PetscErrorCode err;

  err = DMPlexGetHeightStratum(dmMesh, 0, &cStart, &cEnd);CHECK_PETSC_ERROR(err);
  err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);
  
  const int spaceDim = mesh.dimension();
  const int fiberDim = _data->numDOF;
  topology::Field<topology::Mesh> field(mesh);
  field.addField("bc", spaceDim);
  field.setupFields();
  field.newSection(pylith::topology::FieldBase::VERTICES_FIELD, fiberDim);
  field.updateDof("bc", pylith::topology::FieldBase::VERTICES_FIELD, fiberDim);
  bc.setConstraintSizes(field);
  field.allocate();
  bc.setConstraints(field);

  PetscSection fieldSection = field.petscSection();
  PetscVec fieldVec = field.localVector();
  CPPUNIT_ASSERT(fieldSection);CPPUNIT_ASSERT(fieldVec);

  const PetscInt numCells = cEnd - cStart;
  const PetscInt offset = numCells;
  int iConstraint = 0;
  for(PetscInt v = vStart; v < vEnd; ++v) {
    const PetscInt *cInd, *fcInd;
    PetscInt dof, cdof, fdof, fcdof;

    err = PetscSectionGetDof(fieldSection, v, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetConstraintDof(fieldSection, v, &cdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetFieldDof(fieldSection, v, 0, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetFieldConstraintDof(fieldSection, v, 0, &fcdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetConstraintIndices(fieldSection, v, &cInd);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetFieldConstraintIndices(fieldSection, v, 0, &fcInd);CHECK_PETSC_ERROR(err);
    if (v != _data->constrainedPoints[iConstraint] + offset) {
      CPPUNIT_ASSERT_EQUAL(0, cdof);
      CPPUNIT_ASSERT_EQUAL(0, fcdof);
    } else {
      if (_data->numFixedDOF) {CPPUNIT_ASSERT(cInd);}
      CPPUNIT_ASSERT_EQUAL(_data->numFixedDOF, cdof);
      //CPPUNIT_ASSERT_EQUAL(_data->numFixedDOF, fcdof);
      for(PetscInt iDOF = 0; iDOF < _data->numFixedDOF; ++iDOF) {
        CPPUNIT_ASSERT_EQUAL(_data->fixedDOF[iDOF], cInd[iDOF]);
        //CPPUNIT_ASSERT_EQUAL(_data->fixedDOF[iDOF], fcInd[iDOF]);
      }
      ++iConstraint;
    } // if/else
  } // for
} // testSetConstraints

// ----------------------------------------------------------------------
// Test setField().
void
pylith::bc::TestDirichletBC::testSetField(void)
{ // testSetField
  topology::Mesh mesh;
  DirichletBC bc;
  _initialize(&mesh, &bc);
  CPPUNIT_ASSERT(0 != _data);

  PetscDM dmMesh = mesh.dmMesh();
  CPPUNIT_ASSERT(dmMesh);
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
  bc.setConstraintSizes(field);
  field.allocate();
  bc.setConstraints(field);

  PetscSection fieldSection = field.petscSection();
  PetscVec fieldVec = field.localVector();
  CPPUNIT_ASSERT(fieldSection);CPPUNIT_ASSERT(fieldVec);
  const PylithScalar tolerance = 1.0e-06;
  const PylithScalar dispScale = _data->lengthScale;
  const PylithScalar velocityScale = _data->lengthScale / _data->timeScale;
  const PylithScalar timeScale = _data->timeScale;

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
  const PylithScalar t = 1.0 / timeScale;
  bc.setField(t, field);

  // Create list of unconstrained DOF at constrained DOF
  const int numFreeDOF = _data->numDOF - _data->numFixedDOF;
  int_array freeDOF(numFreeDOF);
  int index = 0;
  for(int iDOF = 0; iDOF < _data->numDOF; ++iDOF) {
    bool free = true;
    for(int iFixed = 0; iFixed < _data->numFixedDOF; ++iFixed)
      if (iDOF == _data->fixedDOF[iFixed])
        free = false;
    if (free)
      freeDOF[index++] = iDOF;
  } // for
  assert(index == numFreeDOF);

  const PetscInt numCells = cEnd - cStart;
  const PetscInt offset = numCells;
  const PetscInt numFixedDOF = _data->numFixedDOF;
  int iConstraint = 0;

  const PylithScalar tRef = _data->tRef / timeScale;
  
  err = VecGetArray(fieldVec, &values);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt dof, cdof, off;

    err = PetscSectionGetDof(fieldSection, v, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(fieldSection, v, &off);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetConstraintDof(fieldSection, v, &cdof);CHECK_PETSC_ERROR(err);
    if (v != _data->constrainedPoints[iConstraint] + offset) {
      // unconstrained point
      for(PetscInt d = 0; d < dof; ++d)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, values[off+d], tolerance);
    } else {
      // constrained point

      // check unconstrained DOF
      for (int iDOF=0; iDOF < numFreeDOF; ++iDOF)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, values[off+freeDOF[iDOF]], tolerance);

      // check constrained DOF
      for (int iDOF=0; iDOF < numFixedDOF; ++iDOF) {
        const int index = iConstraint * numFixedDOF + iDOF;
        const PylithScalar valueE = (t > tRef) ?
          _data->valuesInitial[index]/dispScale + (t-tRef)*_data->valueRate/velocityScale :
          _data->valuesInitial[index]/dispScale;
        CPPUNIT_ASSERT_DOUBLES_EQUAL(valueE, values[off+_data->fixedDOF[iDOF]], tolerance);
      } // for
      ++iConstraint;
    } // if/else
  } // for
  err = VecRestoreArray(fieldVec, &values);CHECK_PETSC_ERROR(err);
} // testSetField

// ----------------------------------------------------------------------
// Test setFieldIncr().
void
pylith::bc::TestDirichletBC::testSetFieldIncr(void)
{ // testSetFieldIncr
  topology::Mesh mesh;
  DirichletBC bc;
  _initialize(&mesh, &bc);
  CPPUNIT_ASSERT(0 != _data);

  PetscDM dmMesh = mesh.dmMesh();
  CPPUNIT_ASSERT(dmMesh);
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
  bc.setConstraintSizes(field);
  field.allocate();
  bc.setConstraints(field);

  PetscSection fieldSection = field.petscSection();
  PetscVec fieldVec = field.localVector();
  CPPUNIT_ASSERT(fieldSection);CPPUNIT_ASSERT(fieldVec);
  const PylithScalar tolerance = 1.0e-06;
  const PylithScalar dispScale = _data->lengthScale;
  const PylithScalar velocityScale = _data->lengthScale / _data->timeScale;
  const PylithScalar timeScale = _data->timeScale;

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
  const PylithScalar t0 = 1.0 / timeScale;
  const PylithScalar t1 = 2.0 / timeScale;
  bc.setFieldIncr(t0, t1, field);

  // Create list of unconstrained DOF at constrained DOF
  const int numFreeDOF = _data->numDOF - _data->numFixedDOF;
  int_array freeDOF(numFreeDOF);
  int index = 0;
  for(int iDOF = 0; iDOF < _data->numDOF; ++iDOF) {
    bool free = true;
    for(int iFixed = 0; iFixed < _data->numFixedDOF; ++iFixed)
      if (iDOF == _data->fixedDOF[iFixed])
        free = false;
    if (free)
      freeDOF[index++] = iDOF;
  } // for
  assert(index == numFreeDOF);

  const PetscInt numCells = cEnd - cStart;
  const PetscInt offset = numCells;
  const PetscInt numFixedDOF = _data->numFixedDOF;
  int iConstraint = 0;

  const PylithScalar tRef = _data->tRef / timeScale;

  err = VecGetArray(fieldVec, &values);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt dof, cdof, off;

    err = PetscSectionGetDof(fieldSection, v, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(fieldSection, v, &off);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetConstraintDof(fieldSection, v, &cdof);CHECK_PETSC_ERROR(err);
    if (v != _data->constrainedPoints[iConstraint] + offset) {
      // unconstrained point
      for(PetscInt d = 0; d < dof; ++d)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, values[off+d], tolerance);
    } else {
      // constrained point

      // check unconstrained DOF
      for (int iDOF=0; iDOF < numFreeDOF; ++iDOF)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, values[off+freeDOF[iDOF]], tolerance);

      // check constrained DOF
      for (int iDOF=0; iDOF < numFixedDOF; ++iDOF) {
        const PylithScalar valueE = (t0 > tRef) ? (t1-t0)*_data->valueRate/velocityScale :
          (t1 > tRef) ? (t1-tRef)*_data->valueRate/velocityScale : 0.0;
        CPPUNIT_ASSERT_DOUBLES_EQUAL(valueE, values[off+_data->fixedDOF[iDOF]], tolerance);
      } // for
      ++iConstraint;
    } // if/else
  } // for
} // testSetFieldIncr

// ----------------------------------------------------------------------
void
pylith::bc::TestDirichletBC::_initialize(topology::Mesh* mesh,
					 DirichletBC* const bc) const
{ // _initialize
  CPPUNIT_ASSERT(0 != _data);
  CPPUNIT_ASSERT(0 != bc);

  meshio::MeshIOAscii iohandler;
  iohandler.filename(_data->meshFilename);
  iohandler.read(mesh);

  spatialdata::geocoords::CSCart cs;
  spatialdata::units::Nondimensional normalizer;
  normalizer.lengthScale(_data->lengthScale);
  normalizer.pressureScale(_data->pressureScale);
  normalizer.densityScale(_data->densityScale);
  normalizer.timeScale(_data->timeScale);
  cs.setSpaceDim(mesh->dimension());
  cs.initialize();
  mesh->coordsys(&cs);
  mesh->nondimensionalize(normalizer);

  spatialdata::spatialdb::SimpleDB db("TestDirichletBC initial");
  spatialdata::spatialdb::SimpleIOAscii dbIO;
  dbIO.filename(_data->dbFilename);
  db.ioHandler(&dbIO);
  db.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

  spatialdata::spatialdb::UniformDB dbRate("TestDirichletBC rate");
  const int numValues = 4;
  const char* names[] = { 
    "displacement-rate-x", 
    "displacement-rate-y", 
    "displacement-rate-z",
    "rate-start-time"};
  const char* units[] = { 
    "m", 
    "m", 
    "m",
    "s"};
  const double values[numValues] = { 
    _data->valueRate,
    _data->valueRate,
    _data->valueRate,
    _data->tRef,
  };
  dbRate.setData(names, units, values, numValues);

  const PylithScalar upDir[] = { 0.0, 0.0, 1.0 };

  bc->label(_data->label);
  bc->dbInitial(&db);
  bc->dbRate(&dbRate);
  bc->bcDOF(_data->fixedDOF, _data->numFixedDOF);
  bc->normalizer(normalizer);
  bc->initialize(*mesh, upDir);
} // _initialize


// End of file 
