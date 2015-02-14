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

#include "TestDirichletBoundary.hh" // Implementation of class methods

#include "pylith/bc/DirichletBoundary.hh" // USES DirichletBoundary

#include "data/DirichletData.hh" // USES DirichletData

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps::nondimensionalize()
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional
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
  PYLITH_METHOD_BEGIN;

  _data = 0;

  PYLITH_METHOD_END;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::bc::TestDirichletBoundary::tearDown(void)
{ // tearDown
  PYLITH_METHOD_BEGIN;

  delete _data; _data = 0;

  PYLITH_METHOD_END;
} // tearDown

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::bc::TestDirichletBoundary::testConstructor(void)
{ // testConstructor
  PYLITH_METHOD_BEGIN;

  DirichletBoundary bc;

  PYLITH_METHOD_END;
} // testConstructor

// ----------------------------------------------------------------------
// Test initialize().
void
pylith::bc::TestDirichletBoundary::testInitialize(void)
{ // testInitialize
  PYLITH_METHOD_BEGIN;

  topology::Mesh mesh;
  DirichletBoundary bc;
  _initialize(&mesh, &bc);

  CPPUNIT_ASSERT(_data);

  const PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);

  const int numFixedDOF = _data->numFixedDOF;
  const size_t numBoundary = _data->numConstrainedPts;

  // Check vertices in boundary mesh.
  const PetscDM dmSubMesh = bc._boundaryMesh->dmMesh();CPPUNIT_ASSERT(dmSubMesh);
  topology::Stratum depthStratum(dmSubMesh, topology::Stratum::DEPTH, 0);
  CPPUNIT_ASSERT_EQUAL(PetscInt(numBoundary), depthStratum.size());

  // :TODO: Check cells in boundary mesh.

  PYLITH_METHOD_END;
} // testInitialize

// ----------------------------------------------------------------------
void
pylith::bc::TestDirichletBoundary::_initialize(topology::Mesh* mesh,
					       DirichletBoundary* const bc) const
{ // _initialize
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT(bc);

  meshio::MeshIOAscii iohandler;
  iohandler.filename(_data->meshFilename);
  iohandler.read(mesh);

  spatialdata::geocoords::CSCart cs;
  spatialdata::units::Nondimensional normalizer;
  cs.setSpaceDim(mesh->dimension());
  cs.initialize();
  mesh->coordsys(&cs);
  topology::MeshOps::nondimensionalize(mesh, normalizer);

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

  PYLITH_METHOD_END;
} // _initialize


// End of file 
