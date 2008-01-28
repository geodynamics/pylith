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

#include "TestDirichletPointsMulti.hh" // Implementation of class methods

#include "pylith/bc/DirichletPoints.hh" // USES DirichletPoints

#include "data/DirichletPointsDataMulti.hh" // USES DirichletPointsData
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/utils/sievetypes.hh" // USES PETSc Mesh

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestDirichletPointsMulti::setUp(void)
{ // setUp
  _data = 0;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::bc::TestDirichletPointsMulti::tearDown(void)
{ // tearDown
  delete _data; _data = 0;
} // tearDown

// ----------------------------------------------------------------------
// Test setConstraintSizes().
void
pylith::bc::TestDirichletPointsMulti::testSetConstraintSizes(void)
{ // testSetConstraintSizes
  ALE::Obj<Mesh> mesh;
  DirichletPoints bcA;
  DirichletPoints bcB;
  _initialize(&mesh, &bcA, &bcB);

  const ALE::Obj<real_section_type>& field = mesh->getRealSection("field");
  const ALE::Obj<Mesh::label_sequence>& vertices = mesh->depthStratum(0);
  field->setFiberDimension(vertices, _data->numDOF);
  bcA.setConstraintSizes(field, mesh);
  bcB.setConstraintSizes(field, mesh);

  CPPUNIT_ASSERT(0 != _data);

  const int numCells = mesh->heightStratum(0)->size();
  const int offset = numCells;
  for (Mesh::label_sequence::iterator v_iter = vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    CPPUNIT_ASSERT_EQUAL(_data->numDOF, field->getFiberDimension(*v_iter));
    
    CPPUNIT_ASSERT_EQUAL(_data->constraintSizes[*v_iter-offset],
			 field->getConstraintDimension(*v_iter));
  } // for
} // testSetConstraintSizes

// ----------------------------------------------------------------------
// Test setConstraints().
void
pylith::bc::TestDirichletPointsMulti::testSetConstraints(void)
{ // testSetConstraints
  ALE::Obj<Mesh> mesh;
  DirichletPoints bcA;
  DirichletPoints bcB;
  _initialize(&mesh, &bcA, &bcB);

  const ALE::Obj<real_section_type>& field = mesh->getRealSection("field");
  const ALE::Obj<Mesh::label_sequence>& vertices = mesh->depthStratum(0);
  field->setFiberDimension(vertices, _data->numDOF);
  bcA.setConstraintSizes(field, mesh);
  bcB.setConstraintSizes(field, mesh);
  mesh->allocate(field);
  bcA.setConstraints(field, mesh);
  bcB.setConstraints(field, mesh);

  CPPUNIT_ASSERT(0 != _data);

  const int numCells = mesh->heightStratum(0)->size();
  const int offset = numCells;
  int index = 0;
  for (Mesh::label_sequence::iterator v_iter = vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    const int numConstrainedDOF = _data->constraintSizes[*v_iter-offset];
    if (numConstrainedDOF > 0) {
      const int* fixedDOF = field->getConstraintDof(*v_iter);
      for (int iDOF=0; iDOF < numConstrainedDOF; ++iDOF)
	CPPUNIT_ASSERT_EQUAL(_data->constrainedDOF[index++], fixedDOF[iDOF]);
    } // if
  } // for
} // testSetConstraints

// ----------------------------------------------------------------------
// Test setField().
void
pylith::bc::TestDirichletPointsMulti::testSetField(void)
{ // testSetField
  ALE::Obj<Mesh> mesh;
  DirichletPoints bcA;
  DirichletPoints bcB;
  _initialize(&mesh, &bcA, &bcB);

  const ALE::Obj<real_section_type>& field = mesh->getRealSection("field");
  const ALE::Obj<Mesh::label_sequence>& vertices = mesh->depthStratum(0);
  field->setFiberDimension(vertices, _data->numDOF);
  bcA.setConstraintSizes(field, mesh);
  bcB.setConstraintSizes(field, mesh);
  mesh->allocate(field);
  bcA.setConstraints(field, mesh);
  bcB.setConstraints(field, mesh);

  CPPUNIT_ASSERT(0 != _data);
  const double tolerance = 1.0e-06;

  // All values should be zero.
  field->zero();
  for (Mesh::label_sequence::iterator v_iter = vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    const int fiberDim = field->getFiberDimension(*v_iter);
    const real_section_type::value_type* values = 
      mesh->restrict(field, *v_iter);
    for (int i=0; i < fiberDim; ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, values[i], tolerance);
  } // for

  // Only unconstrained values should be zero.
  // Expected values set in _data->field
  const double t = 1.0;
  bcA.setField(t, field, mesh);
  bcB.setField(t, field, mesh);

  int i = 0;
  for (Mesh::label_sequence::iterator v_iter = vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    const int fiberDim = field->getFiberDimension(*v_iter);
    const real_section_type::value_type* values = 
      mesh->restrict(field, *v_iter);
    for (int iDOF=0; iDOF < fiberDim; ++iDOF)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(_data->field[i++], values[iDOF], tolerance);
  } // for
} // testSetField

// ----------------------------------------------------------------------
void
pylith::bc::TestDirichletPointsMulti::_initialize(ALE::Obj<Mesh>* mesh,
					    DirichletPoints* const bcA,
					    DirichletPoints* const bcB) const
{ // _initialize
  CPPUNIT_ASSERT(0 != _data);
  CPPUNIT_ASSERT(0 != bcA);
  CPPUNIT_ASSERT(0 != bcB);

  meshio::MeshIOAscii iohandler;
  iohandler.filename(_data->meshFilename);
  iohandler.read(mesh);
  CPPUNIT_ASSERT(!mesh->isNull());

  spatialdata::geocoords::CSCart cs;
  cs.setSpaceDim((*mesh)->getDimension());
  cs.initialize();

  // Setup boundary condition A
  spatialdata::spatialdb::SimpleDB db("TestDirichletPointsMulti");
  spatialdata::spatialdb::SimpleIOAscii dbIO;
  dbIO.filename(_data->dbFilenameA);
  db.ioHandler(&dbIO);

  int_array fixedDOFA(_data->fixedDOFA, _data->numFixedDOFA);
  const double upDirVals[] = { 0.0, 0.0, 1.0 };
  double_array upDir(upDirVals, 3);

  bcA->id(_data->idA);
  bcA->label(_data->labelA);
  bcA->db(&db);
  bcA->fixedDOF(fixedDOFA);
  bcA->initialize(*mesh, &cs, upDir);

  // Setup boundary condition B
  dbIO.filename(_data->dbFilenameB);
  db.ioHandler(&dbIO);

  int_array fixedDOFB(_data->fixedDOFB, _data->numFixedDOFB);

  bcB->id(_data->idB);
  bcB->label(_data->labelB);
  bcB->db(&db);
  bcB->fixedDOF(fixedDOFB);
  bcB->initialize(*mesh, &cs, upDir);
} // _initialize


// End of file 
