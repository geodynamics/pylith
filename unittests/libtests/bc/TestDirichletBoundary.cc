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

#include "TestDirichletBoundary.hh" // Implementation of class methods

#include "pylith/bc/DirichletBoundary.hh" // USES DirichletBoundary

#include "data/DirichletData.hh" // USES DirichletData

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/SubMesh.hh" // USES SubMesh
#include "pylith/topology/FieldSubMesh.hh" // USES FieldSubMesh
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii
#include "spatialdata/spatialdb/UniformDB.hh" // USES UniformDB

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestDirichletBoundary );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestDirichletBoundary::setUp(void)
{ // setUp
  _data = 0;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::bc::TestDirichletBoundary::tearDown(void)
{ // tearDown
  delete _data; _data = 0;
} // tearDown

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::bc::TestDirichletBoundary::testConstructor(void)
{ // testConstructor
  DirichletBoundary bc;
} // testConstructor

// ----------------------------------------------------------------------
// Test initialize().
void
pylith::bc::TestDirichletBoundary::testInitialize(void)
{ // testInitialize
  topology::Mesh mesh;
  DirichletBoundary bc;
  _initialize(&mesh, &bc);

  CPPUNIT_ASSERT(0 != _data);

  const int numCells = mesh.sieveMesh()->heightStratum(0)->size();

  const int numFixedDOF = _data->numFixedDOF;
  const size_t numBoundary = _data->numConstrainedPts;
  // Check vertices in boundary mesh
  const ALE::Obj<SieveSubMesh>& sieveMesh = bc._boundaryMesh->sieveMesh();
  const ALE::Obj<SieveSubMesh::label_sequence>& vertices =
    sieveMesh->depthStratum(0);
  const SieveSubMesh::label_sequence::iterator verticesEnd = vertices->end();

  const int offset = numCells;
  if (numFixedDOF > 0) {
    int i = 0;
    for (SieveSubMesh::label_sequence::iterator v_iter=vertices->begin();
	 v_iter != verticesEnd;
	 ++v_iter, ++i) {
      CPPUNIT_ASSERT_EQUAL(_data->constrainedPoints[i]+offset, *v_iter);
    } // for
    CPPUNIT_ASSERT_EQUAL(int(numBoundary), i);
  } // if

  // Check values
  const size_t size = numBoundary * numFixedDOF;
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
void
pylith::bc::TestDirichletBoundary::_initialize(topology::Mesh* mesh,
					       DirichletBoundary* const bc) const
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

  spatialdata::spatialdb::SimpleDB db("TestDirichletBoundary initial");
  spatialdata::spatialdb::SimpleIOAscii dbIO;
  dbIO.filename(_data->dbFilename);
  db.ioHandler(&dbIO);
  db.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

  spatialdata::spatialdb::UniformDB dbRate("TestDirichletBoundary rate");
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
