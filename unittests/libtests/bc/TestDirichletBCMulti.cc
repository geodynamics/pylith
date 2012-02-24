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
#include "pylith/utils/sievetypes.hh" // USES PETSc Mesh

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::RealSection RealSection;

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

  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& vertices =
		 sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());
  
  const int fiberDim = _data->numDOF;
  topology::Field<topology::Mesh> field(mesh);
  field.newSection(vertices, fiberDim);
  const ALE::Obj<RealSection>& fieldSection = field.section();
  CPPUNIT_ASSERT(!fieldSection.isNull());

  bcA.setConstraintSizes(field);
  bcB.setConstraintSizes(field);
  bcC.setConstraintSizes(field);

  const int numCells = sieveMesh->heightStratum(0)->size();
  const int offset = numCells;
  for (SieveMesh::label_sequence::iterator v_iter = vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    CPPUNIT_ASSERT_EQUAL(_data->numDOF,
			 fieldSection->getFiberDimension(*v_iter));
    
    CPPUNIT_ASSERT_EQUAL(_data->constraintSizes[*v_iter-offset],
			 fieldSection->getConstraintDimension(*v_iter));
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

  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& vertices =
		 sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());
  
  const int fiberDim = _data->numDOF;
  topology::Field<topology::Mesh> field(mesh);
  field.newSection(vertices, fiberDim);
  const ALE::Obj<RealSection>& fieldSection = field.section();
  CPPUNIT_ASSERT(!fieldSection.isNull());

  bcA.setConstraintSizes(field);
  bcB.setConstraintSizes(field);
  bcC.setConstraintSizes(field);
  field.allocate();
  bcA.setConstraints(field);
  bcB.setConstraints(field);
  bcC.setConstraints(field);

  const int numCells = sieveMesh->heightStratum(0)->size();
  const int offset = numCells;
  int index = 0;
  for (SieveMesh::label_sequence::iterator v_iter = vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    const int numConstrainedDOF = _data->constraintSizes[*v_iter-offset];
    if (numConstrainedDOF > 0) {
      const int* fixedDOF = fieldSection->getConstraintDof(*v_iter);
      for (int iDOF=0; iDOF < numConstrainedDOF; ++iDOF)
	CPPUNIT_ASSERT_EQUAL(_data->constrainedDOF[index++], fixedDOF[iDOF]);
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

  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& vertices =
		 sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());
  
  const int fiberDim = _data->numDOF;
  topology::Field<topology::Mesh> field(mesh);
  field.newSection(vertices, fiberDim);
  const ALE::Obj<RealSection>& fieldSection = field.section();
  CPPUNIT_ASSERT(!fieldSection.isNull());

  bcA.setConstraintSizes(field);
  bcB.setConstraintSizes(field);
  bcC.setConstraintSizes(field);
  field.allocate();
  bcA.setConstraints(field);
  bcB.setConstraints(field);
  bcC.setConstraints(field);

  const double tolerance = 1.0e-06;

  // All values should be zero.
  field.zero();
  for (SieveMesh::label_sequence::iterator v_iter = vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    const int fiberDim = fieldSection->getFiberDimension(*v_iter);
    const RealSection::value_type* values = fieldSection->restrictPoint(*v_iter);
    for (int i=0; i < fiberDim; ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, values[i], tolerance);
  } // for

  // Only unconstrained values should be zero.
  // Expected values set in _data->field
  const double t = 10.0;
  bcA.setField(t, field);
  bcB.setField(t, field);
  bcC.setField(t, field);

  int i = 0;
  for (SieveMesh::label_sequence::iterator v_iter = vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    const int fiberDim = fieldSection->getFiberDimension(*v_iter);
    const RealSection::value_type* values = fieldSection->restrictPoint(*v_iter);
    for (int iDOF=0; iDOF < fiberDim; ++iDOF)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(_data->field[i++], values[iDOF], tolerance);
  } // for
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
  bcA.useSolnIncr(true);
  bcB.useSolnIncr(true);
  bcC.useSolnIncr(true);

  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& vertices =
		 sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());
  
  const int fiberDim = _data->numDOF;
  topology::Field<topology::Mesh> field(mesh);
  field.newSection(vertices, fiberDim);
  const ALE::Obj<RealSection>& fieldSection = field.section();
  CPPUNIT_ASSERT(!fieldSection.isNull());

  bcA.setConstraintSizes(field);
  bcB.setConstraintSizes(field);
  bcC.setConstraintSizes(field);
  field.allocate();
  bcA.setConstraints(field);
  bcB.setConstraints(field);
  bcC.setConstraints(field);

  const double tolerance = 1.0e-06;

  // All values should be zero.
  field.zero();
  for (SieveMesh::label_sequence::iterator v_iter = vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    const int fiberDim = fieldSection->getFiberDimension(*v_iter);
    const RealSection::value_type* values = fieldSection->restrictPoint(*v_iter);
    for (int i=0; i < fiberDim; ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, values[i], tolerance);
  } // for

  // Only unconstrained values should be zero.
  // Expected values set in _data->field
  const double t0 = 10.0;
  const double t1 = 14.0;
  bcA.setFieldIncr(t0, t1, field);
  bcB.setFieldIncr(t0, t1, field);
  bcC.setFieldIncr(t0, t1, field);

  int i = 0;
  for (SieveMesh::label_sequence::iterator v_iter = vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    const int fiberDim = fieldSection->getFiberDimension(*v_iter);
    const RealSection::value_type* values = fieldSection->restrictPoint(*v_iter);
    for (int iDOF=0; iDOF < fiberDim; ++iDOF)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(_data->fieldIncr[i++], values[iDOF], tolerance);
  } // for
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
  spatialdata::units::Nondimensional normalizer;
  cs.setSpaceDim(mesh->dimension());
  cs.initialize();
  mesh->coordsys(&cs);
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

  const double upDir[] = { 0.0, 0.0, 1.0 };

  bcA->label(_data->labelA);
  bcA->dbInitial(&db);
  bcA->dbRate(&dbRate);
  bcA->bcDOF(_data->fixedDOFA, _data->numFixedDOFA);
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
    bcC->initialize(*mesh, upDir);
  } // if
} // _initialize


// End of file 
