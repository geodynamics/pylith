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

#include "TestDirichletBC.hh" // Implementation of class methods

#include "pylith/bc/DirichletBC.hh" // USES DirichletBC

#include "data/DirichletData.hh" // USES DirichletData

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps::nondimensionalize()
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorMesh.hh" // USES VisitorMesh
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii
#include "spatialdata/spatialdb/UniformDB.hh" // USES UniformDB
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestDirichletBC );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestDirichletBC::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  _data = 0;

  PYLITH_METHOD_END;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::bc::TestDirichletBC::tearDown(void)
{ // tearDown
  PYLITH_METHOD_BEGIN;

  delete _data; _data = 0;

  PYLITH_METHOD_END;
} // tearDown

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::bc::TestDirichletBC::testConstructor(void)
{ // testConstructor
  PYLITH_METHOD_BEGIN;

  DirichletBC bc;

  PYLITH_METHOD_END;
} // testConstructor

// ----------------------------------------------------------------------
// Test initialize().
void
pylith::bc::TestDirichletBC::testInitialize(void)
{ // testInitialize
  PYLITH_METHOD_BEGIN;

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
    CPPUNIT_ASSERT(bc._parameters);
    topology::VecVisitorMesh initialVisitor(bc._parameters->get("initial"));
    const PetscScalar* initialArray = initialVisitor.localArray();CPPUNIT_ASSERT(initialArray);

    const PylithScalar tolerance = 1.0e-06;
    const PylithScalar dispScale = _data->lengthScale;
    const PylithScalar velocityScale = _data->lengthScale / _data->timeScale;
    for (int i=0; i < numPoints; ++i) {
      const PetscInt p_value = _data->constrainedPoints[i]+offset;

      const PetscInt off = initialVisitor.sectionOffset(p_value);
      CPPUNIT_ASSERT_EQUAL(numFixedDOF, initialVisitor.sectionDof(p_value));
      for(int iDOF = 0; iDOF < numFixedDOF; ++iDOF) {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(_data->valuesInitial[i*numFixedDOF+iDOF]/dispScale, initialArray[off+iDOF], tolerance);
      } // for
    } // for
    
    // Check rate of change
    topology::VecVisitorMesh rateVisitor(bc._parameters->get("rate"));
    const PetscScalar* rateArray = rateVisitor.localArray();CPPUNIT_ASSERT(rateArray);

    for (int i=0; i < numPoints; ++i) {
      const PetscInt p_value = _data->constrainedPoints[i]+offset;

      const PetscInt off = rateVisitor.sectionOffset(p_value);
      CPPUNIT_ASSERT_EQUAL(numFixedDOF, rateVisitor.sectionDof(p_value));
      for(int iDOF = 0; iDOF < numFixedDOF; ++iDOF) 
        CPPUNIT_ASSERT_DOUBLES_EQUAL(_data->valueRate/velocityScale, rateArray[off+iDOF], tolerance);
    } // for
  } // if

  PYLITH_METHOD_END;
} // testInitialize

// ----------------------------------------------------------------------
// Test numDimConstrained().
void
pylith::bc::TestDirichletBC::testNumDimConstrained(void)
{ // testNumDimConstrained
  PYLITH_METHOD_BEGIN;

  topology::Mesh mesh;
  DirichletBC bc;
  _initialize(&mesh, &bc);
  CPPUNIT_ASSERT(_data);

  CPPUNIT_ASSERT_EQUAL(_data->numFixedDOF, bc.numDimConstrained());

  PYLITH_METHOD_END;
} // testNumDimConstrained

// ----------------------------------------------------------------------
// Test setConstraintSizes().
void
pylith::bc::TestDirichletBC::testSetConstraintSizes(void)
{ // testSetConstraintSizes
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_data);

  topology::Mesh mesh;
  DirichletBC bc;
  _initialize(&mesh, &bc);

  const int fiberDim = _data->numDOF;
  const int spaceDim = mesh.dimension();
  topology::Field field(mesh);
  field.subfieldAdd("bc", spaceDim, topology::Field::VECTOR);
  field.subfieldsSetup();
  field.newSection(pylith::topology::FieldBase::VERTICES_FIELD, fiberDim);
  field.subfieldSetDof("bc", pylith::topology::FieldBase::VERTICES_FIELD, fiberDim);
  bc.setConstraintSizes(field); // Does not handle fields right now
  field.allocate();

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

  int iConstraint = 0;
  PetscErrorCode err = 0;
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt dof, cdof, fdof, fcdof;

    err = PetscSectionGetDof(fieldSection, v, &dof);PYLITH_CHECK_ERROR(err);
    err = PetscSectionGetConstraintDof(fieldSection, v, &cdof);PYLITH_CHECK_ERROR(err);
    err = PetscSectionGetFieldDof(fieldSection, v, 0, &fdof);PYLITH_CHECK_ERROR(err);
    err = PetscSectionGetFieldConstraintDof(fieldSection, v, 0, &fcdof);PYLITH_CHECK_ERROR(err);
    if (v != _data->constrainedPoints[iConstraint] + offset) {
      CPPUNIT_ASSERT_EQUAL(_data->numDOF, dof);
      CPPUNIT_ASSERT_EQUAL(0, cdof);
      CPPUNIT_ASSERT_EQUAL(_data->numDOF, fdof);
      CPPUNIT_ASSERT_EQUAL(0, fcdof);
    } else {
      CPPUNIT_ASSERT_EQUAL(_data->numDOF, dof);
      CPPUNIT_ASSERT_EQUAL(_data->numFixedDOF, cdof);
      CPPUNIT_ASSERT_EQUAL(_data->numDOF, fdof);
      CPPUNIT_ASSERT_EQUAL(_data->numFixedDOF, fcdof);
      ++iConstraint;
    } // if/else
  } // for

  PYLITH_METHOD_END;
} // testSetConstraintSizes

// ----------------------------------------------------------------------
// Test setConstraints().
void
pylith::bc::TestDirichletBC::testSetConstraints(void)
{ // testSetConstraints
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_data);

  topology::Mesh mesh;
  DirichletBC bc;
  _initialize(&mesh, &bc);

  const int spaceDim = mesh.dimension();
  const int fiberDim = _data->numDOF;
  topology::Field field(mesh);
  field.subfieldAdd("bc", spaceDim, topology::Field::VECTOR);
  field.subfieldsSetup();
  field.newSection(pylith::topology::FieldBase::VERTICES_FIELD, fiberDim);
  field.subfieldSetDof("bc", pylith::topology::FieldBase::VERTICES_FIELD, fiberDim);
  bc.setConstraintSizes(field);
  field.allocate();
  bc.setConstraints(field);

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

  int iConstraint = 0;
  PetscErrorCode err = 0;
  for(PetscInt v = vStart; v < vEnd; ++v) {
    const PetscInt *cInd, *fcInd;
    PetscInt dof, cdof, fdof, fcdof;

    err = PetscSectionGetDof(fieldSection, v, &dof);PYLITH_CHECK_ERROR(err);
    err = PetscSectionGetConstraintDof(fieldSection, v, &cdof);PYLITH_CHECK_ERROR(err);
    err = PetscSectionGetFieldDof(fieldSection, v, 0, &dof);PYLITH_CHECK_ERROR(err);
    err = PetscSectionGetFieldConstraintDof(fieldSection, v, 0, &fcdof);PYLITH_CHECK_ERROR(err);
    err = PetscSectionGetConstraintIndices(fieldSection, v, &cInd);PYLITH_CHECK_ERROR(err);
    err = PetscSectionGetFieldConstraintIndices(fieldSection, v, 0, &fcInd);PYLITH_CHECK_ERROR(err);
    if (v != _data->constrainedPoints[iConstraint] + offset) {
      CPPUNIT_ASSERT_EQUAL(0, cdof);
      CPPUNIT_ASSERT_EQUAL(0, fcdof);
    } else {
      if (_data->numFixedDOF) {CPPUNIT_ASSERT(cInd);}
      CPPUNIT_ASSERT_EQUAL(_data->numFixedDOF, cdof);
      CPPUNIT_ASSERT_EQUAL(_data->numFixedDOF, fcdof);
      for(PetscInt iDOF = 0; iDOF < _data->numFixedDOF; ++iDOF) {
        CPPUNIT_ASSERT_EQUAL(_data->fixedDOF[iDOF], cInd[iDOF]);
        CPPUNIT_ASSERT_EQUAL(_data->fixedDOF[iDOF], fcInd[iDOF]);
      }
      ++iConstraint;
    } // if/else
  } // for

  PYLITH_METHOD_END;
} // testSetConstraints

// ----------------------------------------------------------------------
// Test setField().
void
pylith::bc::TestDirichletBC::testSetField(void)
{ // testSetField
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_data);

  topology::Mesh mesh;
  DirichletBC bc;
  _initialize(&mesh, &bc);

  const int fiberDim = _data->numDOF;
  topology::Field field(mesh);
  field.subfieldAdd("bc", fiberDim, topology::Field::VECTOR);
  field.subfieldsSetup();
  field.newSection(pylith::topology::FieldBase::VERTICES_FIELD, fiberDim);
  field.subfieldSetDof("bc", pylith::topology::FieldBase::VERTICES_FIELD, fiberDim);
  bc.setConstraintSizes(field);
  field.allocate();
  bc.setConstraints(field);

  // Scales
  const PylithScalar tolerance = 1.0e-06;
  const PylithScalar dispScale = _data->lengthScale;
  const PylithScalar velocityScale = _data->lengthScale / _data->timeScale;
  const PylithScalar timeScale = _data->timeScale;

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

  const PylithScalar tRef = _data->tRef / timeScale;
  const PetscInt numFixedDOF = _data->numFixedDOF;
  int iConstraint = 0;

  fieldVisitor.initialize(field);
  fieldArray = fieldVisitor.localArray();
  for(PetscInt v = vStart; v < vEnd; ++v) {
    const PetscInt off = fieldVisitor.sectionOffset(v);
    const PetscInt dof = fieldVisitor.sectionDof(v);
    const PetscInt cdof = fieldVisitor.sectionConstraintDof(v);
    if (v != _data->constrainedPoints[iConstraint] + offset) {
      // unconstrained point
      for(PetscInt d = 0; d < dof; ++d)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, fieldArray[off+d], tolerance);
    } else {
      // constrained point

      // check unconstrained DOF
      for (int iDOF=0; iDOF < numFreeDOF; ++iDOF)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, fieldArray[off+freeDOF[iDOF]], tolerance);

      // check constrained DOF
      for (int iDOF=0; iDOF < numFixedDOF; ++iDOF) {
        const int index = iConstraint * numFixedDOF + iDOF;
        const PylithScalar valueE = (t > tRef) ?
          _data->valuesInitial[index]/dispScale + (t-tRef)*_data->valueRate/velocityScale :
          _data->valuesInitial[index]/dispScale;
        CPPUNIT_ASSERT_DOUBLES_EQUAL(valueE, fieldArray[off+_data->fixedDOF[iDOF]], tolerance);
      } // for
      ++iConstraint;
    } // if/else
  } // for

  PYLITH_METHOD_END;
} // testSetField

// ----------------------------------------------------------------------
// Test setFieldIncr().
void
pylith::bc::TestDirichletBC::testSetFieldIncr(void)
{ // testSetFieldIncr
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_data);

  topology::Mesh mesh;
  DirichletBC bc;
  _initialize(&mesh, &bc);

  const int fiberDim = _data->numDOF;
  topology::Field field(mesh);
  field.subfieldAdd("bc", fiberDim, topology::Field::VECTOR);
  field.subfieldsSetup();
  field.newSection(pylith::topology::FieldBase::VERTICES_FIELD, fiberDim);
  field.subfieldSetDof("bc", pylith::topology::FieldBase::VERTICES_FIELD, fiberDim);
  bc.setConstraintSizes(field);
  field.allocate();
  bc.setConstraints(field);

  // Scales
  const PylithScalar tolerance = 1.0e-06;
  const PylithScalar dispScale = _data->lengthScale;
  const PylithScalar velocityScale = _data->lengthScale / _data->timeScale;
  const PylithScalar timeScale = _data->timeScale;

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

  const PylithScalar tRef = _data->tRef / timeScale;
  const PetscInt numFixedDOF = _data->numFixedDOF;
  int iConstraint = 0;

  fieldVisitor.initialize(field);
  fieldArray = fieldVisitor.localArray();
  for(PetscInt v = vStart; v < vEnd; ++v) {
    const PetscInt off = fieldVisitor.sectionOffset(v);
    const PetscInt dof = fieldVisitor.sectionDof(v);
    const PetscInt cdof = fieldVisitor.sectionConstraintDof(v);
    if (v != _data->constrainedPoints[iConstraint] + offset) {
      // unconstrained point
      for(PetscInt d = 0; d < dof; ++d)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, fieldArray[off+d], tolerance);
    } else {
      // constrained point

      // check unconstrained DOF
      for (int iDOF=0; iDOF < numFreeDOF; ++iDOF)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, fieldArray[off+freeDOF[iDOF]], tolerance);

      // check constrained DOF
      for (int iDOF=0; iDOF < numFixedDOF; ++iDOF) {
        const PylithScalar valueE = (t0 > tRef) ? (t1-t0)*_data->valueRate/velocityScale :
          (t1 > tRef) ? (t1-tRef)*_data->valueRate/velocityScale : 0.0;
        CPPUNIT_ASSERT_DOUBLES_EQUAL(valueE, fieldArray[off+_data->fixedDOF[iDOF]], tolerance);
      } // for
      ++iConstraint;
    } // if/else
  } // for

  PYLITH_METHOD_END;
} // testSetFieldIncr

// ----------------------------------------------------------------------
void
pylith::bc::TestDirichletBC::_initialize(topology::Mesh* mesh,
					 DirichletBC* const bc) const
{ // _initialize
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT(bc);

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
  topology::MeshOps::nondimensionalize(mesh, normalizer);

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

  PYLITH_METHOD_END;
} // _initialize


// End of file 
