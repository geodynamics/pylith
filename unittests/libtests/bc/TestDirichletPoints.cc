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

#include "TestDirichletPoints.hh" // Implementation of class methods

#include "pylith/bc/DirichletPoints.hh" // USES DirichletPoints

#include "data/DirichletData.hh" // USES DirichletData
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/utils/sievetypes.hh" // USES PETSc Mesh

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii
#include "spatialdata/spatialdb/UniformDB.hh" // USES UniformDB

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestDirichletPoints );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestDirichletPoints::setUp(void)
{ // setUp
  _data = 0;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::bc::TestDirichletPoints::tearDown(void)
{ // tearDown
  delete _data; _data = 0;
} // tearDown

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::bc::TestDirichletPoints::testConstructor(void)
{ // testConstructor
  DirichletPoints bc;
} // testConstructor

// ----------------------------------------------------------------------
// Test fixedDOF()
void
pylith::bc::TestDirichletPoints::testFixedDOF(void)
{ // testfixedDOF
  DirichletPoints bc;
  
  const size_t numDOF = 4;
  const int dof[] = { 0, 2, 3, 5 };
  int_array fixedDOF(dof, numDOF);
  bc.fixedDOF(fixedDOF);

  CPPUNIT_ASSERT_EQUAL(numDOF, bc._fixedDOF.size());
  for (int i=0; i < numDOF; ++i)
    CPPUNIT_ASSERT_EQUAL(fixedDOF[i], bc._fixedDOF[i]);
} // testFixedDOF

// ----------------------------------------------------------------------
// Test initialize().
void
pylith::bc::TestDirichletPoints::testInitialize(void)
{ // testInitialize
  ALE::Obj<Mesh> mesh;
  DirichletPoints bc;
  _initialize(&mesh, &bc);

  CPPUNIT_ASSERT(0 != _data);

  const int numCells = mesh->heightStratum(0)->size();

  const int numFixedDOF = _data->numFixedDOF;
  const size_t numPoints = _data->numConstrainedPts;

  // Check points
  const int offset = numCells;
  if (numFixedDOF > 0) {
    CPPUNIT_ASSERT_EQUAL(numPoints, bc._points.size());
    for (int i=0; i < numPoints; ++i)
      CPPUNIT_ASSERT_EQUAL(_data->constrainedPoints[i]+offset, bc._points[i]);
  } // if

  // Check values
  const size_t size = numPoints * numFixedDOF;
  CPPUNIT_ASSERT_EQUAL(size, bc._valuesInitial.size());
  const double tolerance = 1.0e-06;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(_data->valuesInitial[i], bc._valuesInitial[i], 
				 tolerance);

  CPPUNIT_ASSERT_EQUAL(size, bc._valuesRate.size());
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(_data->valueRate, bc._valuesRate[i], 
				 tolerance);
} // testInitialize

// ----------------------------------------------------------------------
// Test setConstraintSizes().
void
pylith::bc::TestDirichletPoints::testSetConstraintSizes(void)
{ // testSetConstraintSizes
  ALE::Obj<Mesh> mesh;
  DirichletPoints bc;
  _initialize(&mesh, &bc);

  const ALE::Obj<real_section_type>& field = mesh->getRealSection("field");
  const ALE::Obj<Mesh::label_sequence>& vertices = mesh->depthStratum(0);
  field->setFiberDimension(vertices, _data->numDOF);
  bc.setConstraintSizes(field, mesh);

  CPPUNIT_ASSERT(0 != _data);

  const int numCells = mesh->heightStratum(0)->size();
  const int offset = numCells;
  int iConstraint = 0;
  for (Mesh::label_sequence::iterator v_iter = vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    if (*v_iter != _data->constrainedPoints[iConstraint] + offset) {
      CPPUNIT_ASSERT_EQUAL(_data->numDOF, field->getFiberDimension(*v_iter));
      CPPUNIT_ASSERT_EQUAL(0, field->getConstraintDimension(*v_iter));
    } else {
      CPPUNIT_ASSERT_EQUAL(_data->numDOF, field->getFiberDimension(*v_iter));
      CPPUNIT_ASSERT_EQUAL(_data->numFixedDOF, 
			   field->getConstraintDimension(*v_iter));
      ++iConstraint;
    } // if/else
  } // for
} // testSetConstraintSizes

// ----------------------------------------------------------------------
// Test setConstraints().
void
pylith::bc::TestDirichletPoints::testSetConstraints(void)
{ // testSetConstraints
  ALE::Obj<Mesh> mesh;
  DirichletPoints bc;
  _initialize(&mesh, &bc);

  const ALE::Obj<real_section_type>& field = mesh->getRealSection("field");
  const ALE::Obj<Mesh::label_sequence>& vertices = mesh->depthStratum(0);
  field->setFiberDimension(vertices, _data->numDOF);
  bc.setConstraintSizes(field, mesh);
  mesh->allocate(field);
  bc.setConstraints(field, mesh);

  CPPUNIT_ASSERT(0 != _data);

  const int numCells = mesh->heightStratum(0)->size();
  const int offset = numCells;
  int iConstraint = 0;
  for (Mesh::label_sequence::iterator v_iter = vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    const int* fixedDOF = field->getConstraintDof(*v_iter);
    if (*v_iter != _data->constrainedPoints[iConstraint] + offset) {
      CPPUNIT_ASSERT_EQUAL(0, field->getConstraintDimension(*v_iter));
      //CPPUNIT_ASSERT(0 == fixedDOF);
    } else {
      CPPUNIT_ASSERT(0 != fixedDOF);
      CPPUNIT_ASSERT_EQUAL(_data->numFixedDOF, 
			   field->getConstraintDimension(*v_iter));
      for (int iDOF=0; iDOF < _data->numFixedDOF; ++iDOF)
	CPPUNIT_ASSERT_EQUAL(_data->fixedDOF[iDOF], fixedDOF[iDOF]);
      ++iConstraint;
    } // if/else
  } // for
} // testSetConstraints

// ----------------------------------------------------------------------
// Test setField().
void
pylith::bc::TestDirichletPoints::testSetField(void)
{ // testSetField
  ALE::Obj<Mesh> mesh;
  DirichletPoints bc;
  _initialize(&mesh, &bc);

  const ALE::Obj<real_section_type>& field = mesh->getRealSection("field");
  const ALE::Obj<Mesh::label_sequence>& vertices = mesh->depthStratum(0);
  field->setFiberDimension(vertices, _data->numDOF);
  bc.setConstraintSizes(field, mesh);
  mesh->allocate(field);
  bc.setConstraints(field, mesh);

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
  const double t = 1.0;
  bc.setField(t, field, mesh);

  // Create list of unconstrained DOF at constrained DOF
  const int numFreeDOF = _data->numDOF - _data->numFixedDOF;
  int_array freeDOF(numFreeDOF);
  int index = 0;
  for (int iDOF=0; iDOF < _data->numDOF; ++iDOF) {
    bool free = true;
    for (int iFixed=0; iFixed < _data->numFixedDOF; ++iFixed)
      if (iDOF == _data->fixedDOF[iFixed])
	free = false;
    if (free)
      freeDOF[index] = iDOF;
  } // for

  const int numCells = mesh->heightStratum(0)->size();
  const int offset = numCells;
  const int numFixedDOF = _data->numFixedDOF;
  int iConstraint = 0;
  for (Mesh::label_sequence::iterator v_iter = vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    const int fiberDim = field->getFiberDimension(*v_iter);
    const real_section_type::value_type* values = 
      mesh->restrict(field, *v_iter);

    if (*v_iter != _data->constrainedPoints[iConstraint] + offset) {
      // unconstrained point
      for (int i=0; i < fiberDim; ++i)
	CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, values[i], tolerance);
    } else {
      // constrained point

      // check unconstrained DOF
      for (int iDOF=0; iDOF < numFreeDOF; ++iDOF)
	CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, values[freeDOF[iDOF]], tolerance);

      // check constrained DOF
      for (int iDOF=0; iDOF < numFixedDOF; ++iDOF) {
	const int index = iConstraint * numFixedDOF + iDOF;
	const double valueE = (t > _data->tRef) ?
	  _data->valuesInitial[index] + (t-_data->tRef)*_data->valueRate :
	  _data->valuesInitial[index];
	CPPUNIT_ASSERT_DOUBLES_EQUAL(valueE, values[_data->fixedDOF[iDOF]],
				     tolerance);
      } // for
      ++iConstraint;
    } // if/else
  } // for
} // testSetField

// ----------------------------------------------------------------------
void
pylith::bc::TestDirichletPoints::_initialize(ALE::Obj<Mesh>* mesh,
				       DirichletPoints* const bc) const
{ // _initialize
  CPPUNIT_ASSERT(0 != _data);
  CPPUNIT_ASSERT(0 != bc);

  meshio::MeshIOAscii iohandler;
  iohandler.filename(_data->meshFilename);
  iohandler.read(mesh);
  CPPUNIT_ASSERT(!mesh->isNull());
  (*mesh)->getFactory()->clear();

  spatialdata::geocoords::CSCart cs;
  cs.setSpaceDim((*mesh)->getDimension());
  cs.initialize();

  spatialdata::spatialdb::SimpleDB db("TestDirichletPoints initial");
  spatialdata::spatialdb::SimpleIOAscii dbIO;
  dbIO.filename(_data->dbFilename);
  db.ioHandler(&dbIO);
  db.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

  spatialdata::spatialdb::UniformDB dbRate("TestDirichletPoints rate");
  const char* names[] = { "dof-0", "dof-1", "dof-2" };
  const double values[] = { _data->valueRate,
			    _data->valueRate,
			    _data->valueRate };
  const int numValues = 3;
  dbRate.setData(names, values, numValues);

  int_array fixedDOF(_data->fixedDOF, _data->numFixedDOF);
  const double upDirVals[] = { 0.0, 0.0, 1.0 };
  double_array upDir(upDirVals, 3);

  bc->label(_data->label);
  bc->db(&db);
  bc->dbRate(&dbRate);
  bc->referenceTime(_data->tRef);
  bc->fixedDOF(fixedDOF);
  bc->initialize(*mesh, &cs, upDir);
} // _initialize


// End of file 
