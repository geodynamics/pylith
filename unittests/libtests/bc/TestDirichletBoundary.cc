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
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestDirichletBoundary.hh" // Implementation of class methods

#include "pylith/bc/DirichletBoundary.hh" // USES DirichletBoundary

#include "data/DirichletData.hh" // USES DirichletData

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/SubMesh.hh" // USES SubMesh
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional
#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii
#include "spatialdata/spatialdb/UniformDB.hh" // USES UniformDB

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestDirichletBoundary );

// ----------------------------------------------------------------------
typedef pylith::topology::SubMesh::SieveMesh SieveMesh;

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
  const ALE::Obj<SieveMesh>& sieveMesh = bc._boundaryMesh->sieveMesh();
  const ALE::Obj<SieveMesh::label_sequence>& vertices =
    sieveMesh->depthStratum(0);
  const SieveMesh::label_sequence::iterator verticesEnd = vertices->end();

  // Check cells
  const int offset = numCells;
  if (numFixedDOF > 0) {
    int i = 0;
    for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
	 v_iter != verticesEnd;
	 ++v_iter, ++i) {
      CPPUNIT_ASSERT_EQUAL(_data->constrainedPoints[i]+offset, *v_iter);
    } // for
    CPPUNIT_ASSERT_EQUAL(int(numBoundary), i);
  } // if
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
  spatialdata::units::Nondimensional normalizer;
  cs.setSpaceDim(mesh->dimension());
  cs.initialize();
  mesh->coordsys(&cs);
  mesh->nondimensionalize(normalizer);

  spatialdata::spatialdb::SimpleDB db("TestDirichletBoundary initial");
  spatialdata::spatialdb::SimpleIOAscii dbIO;
  dbIO.filename(_data->dbFilename);
  db.ioHandler(&dbIO);
  db.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

  spatialdata::spatialdb::UniformDB dbRate("TestDirichletBoundary rate");
  const int numValues = 4;
  const char* names[numValues] = { 
    "displacement-rate-x", 
    "displacement-rate-y", 
    "displacement-rate-z", 
    "rate-start-time",
  };
  const char* units[numValues] = { 
    "m", 
    "m", 
    "m", 
    "s",
  };
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
  bc->initialize(*mesh, upDir);
} // _initialize


// End of file 
