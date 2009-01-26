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

#include "TestDirichletBC.hh" // Implementation of class methods

#include "pylith/bc/DirichletBC.hh" // USES DirichletBC

#include "data/DirichletData.hh" // USES DirichletData

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/FieldUniform.hh" // USES Field
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/utils/sievetypes.hh" // USES PETSc Mesh

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
// Test fixedDOF()
void
pylith::bc::TestDirichletBC::testFixedDOF(void)
{ // testfixedDOF
  DirichletBC bc;
  
  const size_t numDOF = 4;
  const int fixedDOF[] = { 0, 2, 3, 5 };
  bc.fixedDOF(fixedDOF, numDOF);

  CPPUNIT_ASSERT_EQUAL(numDOF, bc._fixedDOF.size());
  for (int i=0; i < numDOF; ++i)
    CPPUNIT_ASSERT_EQUAL(fixedDOF[i], bc._fixedDOF[i]);
} // testFixedDOF

// ----------------------------------------------------------------------
// Test initialize().
void
pylith::bc::TestDirichletBC::testInitialize(void)
{ // testInitialize
  topology::Mesh mesh;
  DirichletBC bc;
  _initialize(&mesh, &bc);
  CPPUNIT_ASSERT(0 != _data);

  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  
  const int numCells = sieveMesh->heightStratum(0)->size();
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
pylith::bc::TestDirichletBC::testSetConstraintSizes(void)
{ // testSetConstraintSizes
  topology::Mesh mesh;
  DirichletBC bc;
  _initialize(&mesh, &bc);
  CPPUNIT_ASSERT(0 != _data);

  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& vertices =
		 sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());
  
  const int fiberDim = _data->numDOF;
  topology::Field field(sieveMesh);
  field.newSection(vertices, fiberDim);
  const ALE::Obj<SieveRealSection>& fieldSection = field.section();
  CPPUNIT_ASSERT(!fieldSection.isNull());

  bc.setConstraintSizes(field, mesh);

  const int numCells = sieveMesh->heightStratum(0)->size();
  const int offset = numCells;
  int iConstraint = 0;
  for (SieveMesh::label_sequence::iterator v_iter = vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    if (*v_iter != _data->constrainedPoints[iConstraint] + offset) {
      CPPUNIT_ASSERT_EQUAL(_data->numDOF,
			   fieldSection->getFiberDimension(*v_iter));
      CPPUNIT_ASSERT_EQUAL(0,
			   fieldSection->getConstraintDimension(*v_iter));
    } else {
      CPPUNIT_ASSERT_EQUAL(_data->numDOF,
			   fieldSection->getFiberDimension(*v_iter));
      CPPUNIT_ASSERT_EQUAL(_data->numFixedDOF, 
			   fieldSection->getConstraintDimension(*v_iter));
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

  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& vertices =
		 sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());
  
  const int fiberDim = _data->numDOF;
  topology::Field field(sieveMesh);
  field.newSection(vertices, fiberDim);
  const ALE::Obj<SieveRealSection>& fieldSection = field.section();
  CPPUNIT_ASSERT(!fieldSection.isNull());

  bc.setConstraintSizes(field, mesh);
  sieveMesh->allocate(fieldSection);
  bc.setConstraints(field, mesh);

  const int numCells = sieveMesh->heightStratum(0)->size();
  const int offset = numCells;
  int iConstraint = 0;
  for (SieveMesh::label_sequence::iterator v_iter = vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    const int* fixedDOF = fieldSection->getConstraintDof(*v_iter);
    if (*v_iter != _data->constrainedPoints[iConstraint] + offset) {
      CPPUNIT_ASSERT_EQUAL(0, fieldSection->getConstraintDimension(*v_iter));
      //CPPUNIT_ASSERT(0 == fixedDOF);
    } else {
      CPPUNIT_ASSERT(0 != fixedDOF);
      CPPUNIT_ASSERT_EQUAL(_data->numFixedDOF, 
			   fieldSection->getConstraintDimension(*v_iter));
      for (int iDOF=0; iDOF < _data->numFixedDOF; ++iDOF)
	CPPUNIT_ASSERT_EQUAL(_data->fixedDOF[iDOF], fixedDOF[iDOF]);
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

  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& vertices =
		 sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());
  
  const int fiberDim = _data->numDOF;
  topology::Field field(sieveMesh);
  field.newSection(vertices, fiberDim);
  const ALE::Obj<SieveRealSection>& fieldSection = field.section();
  CPPUNIT_ASSERT(!fieldSection.isNull());

  bc.setConstraintSizes(field, mesh);
  sieveMesh->allocate(fieldSection);
  bc.setConstraints(field, mesh);

  const double tolerance = 1.0e-06;

  // All values should be zero.
  field.zero();
  for (SieveMesh::label_sequence::iterator v_iter = vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    const int fiberDim = fieldSection->getFiberDimension(*v_iter);
    const SieveRealSection::value_type* values = 
      sieveMesh->restrictClosure(fieldSection, *v_iter);
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

  const int numCells = sieveMesh->heightStratum(0)->size();
  const int offset = numCells;
  const int numFixedDOF = _data->numFixedDOF;
  int iConstraint = 0;
  for (SieveMesh::label_sequence::iterator v_iter = vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    const int fiberDim = fieldSection->getFiberDimension(*v_iter);
    const SieveRealSection::value_type* values = 
      sieveMesh->restrictClosure(fieldSection, *v_iter);

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
pylith::bc::TestDirichletBC::_initialize(topology::Mesh* mesh,
					 DirichletBC* const bc) const
{ // _initialize
  CPPUNIT_ASSERT(0 != _data);
  CPPUNIT_ASSERT(0 != bc);

  meshio::MeshIOAscii iohandler;
  iohandler.filename(_data->meshFilename);
  iohandler.read(mesh);

  spatialdata::geocoords::CSCart cs;
  cs.setSpaceDim(mesh->dimension());
  cs.initialize();
  mesh->coordsys(&cs);

  spatialdata::spatialdb::SimpleDB db("TestDirichletBC initial");
  spatialdata::spatialdb::SimpleIOAscii dbIO;
  dbIO.filename(_data->dbFilename);
  db.ioHandler(&dbIO);
  db.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

  spatialdata::spatialdb::UniformDB dbRate("TestDirichletBC rate");
  const char* names[] = { "dof-0", "dof-1", "dof-2" };
  const double values[] = { _data->valueRate,
			    _data->valueRate,
			    _data->valueRate };
  const int numValues = 3;
  dbRate.setData(names, values, numValues);

  const double upDir[] = { 0.0, 0.0, 1.0 };

  bc->label(_data->label);
  bc->db(&db);
  bc->dbRate(&dbRate);
  bc->referenceTime(_data->tRef);
  bc->fixedDOF(_data->fixedDOF, _data->numFixedDOF);
  bc->initialize(*mesh, upDir);
} // _initialize


// End of file 
