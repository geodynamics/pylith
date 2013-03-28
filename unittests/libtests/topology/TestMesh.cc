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

#include "TestMesh.hh" // Implementation of class methods
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor

#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::topology::TestMesh );

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::topology::TestMesh::testConstructor(void)
{ // testConstructor
  PYLITH_METHOD_BEGIN;

  int result = 0;

  Mesh mesh;
  CPPUNIT_ASSERT(mesh._mesh.isNull());
  CPPUNIT_ASSERT_EQUAL(0, mesh.dimension());
  CPPUNIT_ASSERT_EQUAL(false, mesh.debug());
  MPI_Comm_compare(PETSC_COMM_WORLD, mesh.comm(), &result);
  CPPUNIT_ASSERT_EQUAL(int(MPI_IDENT), result);
  
  Mesh mesh2(2);
  CPPUNIT_ASSERT(!mesh2._mesh.isNull());
  CPPUNIT_ASSERT_EQUAL(2, mesh2.dimension());
  MPI_Comm_compare(PETSC_COMM_WORLD, mesh2.comm(), &result);
  CPPUNIT_ASSERT_EQUAL(int(MPI_IDENT), result);

  Mesh mesh3(1, PETSC_COMM_SELF);
  CPPUNIT_ASSERT(!mesh3._mesh.isNull());
  CPPUNIT_ASSERT_EQUAL(1, mesh3.dimension());
  MPI_Comm_compare(PETSC_COMM_WORLD, mesh3.comm(), &result);
  CPPUNIT_ASSERT_EQUAL(int(MPI_CONGRUENT), result);

  PYLITH_METHOD_END;
} // testConstructor

// ----------------------------------------------------------------------
// Test createDMMesh().
void
pylith::topology::TestMesh::testCreateDMMesh(void)
{ // testCreateDMMesh
  PYLITH_METHOD_BEGIN;

  Mesh mesh;
  CPPUNIT_ASSERT(mesh._mesh.isNull());

  int dim = 2;
  mesh.createDMMesh(dim);
  CPPUNIT_ASSERT(mesh._newMesh);
  CPPUNIT_ASSERT_EQUAL(dim, mesh.dimension());

  dim = 1;
  mesh.createDMMesh(dim);
  CPPUNIT_ASSERT(mesh._newMesh);
  CPPUNIT_ASSERT_EQUAL(dim, mesh.dimension());

  PYLITH_METHOD_END;
} // testCreateDMMesh

// ----------------------------------------------------------------------
// Test dmMesh().
void
pylith::topology::TestMesh::testDMMesh(void)
{ // testDMMesh
  PYLITH_METHOD_BEGIN;

  const int dim = 2;
  PetscInt dmDim;
  Mesh mesh(dim);
  
  DM dmMesh = mesh.dmMesh();
  CPPUNIT_ASSERT(dmMesh);
  PetscErrorCode err = DMPlexGetDimension(dmMesh, &dmDim);CHECK_PETSC_ERROR(err);
  CPPUNIT_ASSERT_EQUAL(dim, dmDim);

  PYLITH_METHOD_END;
} // testDMMesh

// ----------------------------------------------------------------------
// Test coordsys().
void
pylith::topology::TestMesh::testCoordsys(void)
{ // testCoordsys
  PYLITH_METHOD_BEGIN;

  Mesh mesh;

  spatialdata::geocoords::CSCart cs;
  cs.setSpaceDim(2);

  mesh.coordsys(&cs);

  CPPUNIT_ASSERT_EQUAL(cs.spaceDim(), mesh.coordsys()->spaceDim());

  PYLITH_METHOD_END;
} // testCoordsys

// ----------------------------------------------------------------------
// Test debug().
void
pylith::topology::TestMesh::testDebug(void)
{ // testDebug
  PYLITH_METHOD_BEGIN;

  Mesh mesh;
  CPPUNIT_ASSERT_EQUAL(false, mesh.debug());

  mesh.debug(true);
  CPPUNIT_ASSERT_EQUAL(true, mesh.debug());

  PYLITH_METHOD_END;
} // testDebug

// ----------------------------------------------------------------------
// Test dimension().
void
pylith::topology::TestMesh::testDimension(void)
{ // testDimension
  PYLITH_METHOD_BEGIN;

  Mesh mesh;
  CPPUNIT_ASSERT_EQUAL(0, mesh.dimension());

  const int dim = 2;
  Mesh mesh2(dim);
  CPPUNIT_ASSERT_EQUAL(dim, mesh2.dimension());

  PYLITH_METHOD_END;
} // testDimension

// ----------------------------------------------------------------------
// Test comm().
void
pylith::topology::TestMesh::testComm(void)
{ // testComm
  PYLITH_METHOD_BEGIN;

  Mesh mesh;
  CPPUNIT_ASSERT_EQUAL(PETSC_COMM_WORLD, mesh.comm());

  mesh.comm(PETSC_COMM_SELF);
  CPPUNIT_ASSERT_EQUAL(PETSC_COMM_SELF, mesh.comm());

  PYLITH_METHOD_END;
} // testComm

// ----------------------------------------------------------------------
// Test nondimensionalize().
void
pylith::topology::TestMesh::testNondimensionalize(void)
{ // testNondimensionalizer
  PYLITH_METHOD_BEGIN;

  const PylithScalar lengthScale = 2.0;
  const int spaceDim = 2;
  const int numVertices = 4;
  const int coordinates[] = { 
    -1.0, 0.0,
    0.0, -1.0,
    0.0, 1.0,
    1.0, 0.0,
  };

  Mesh mesh;
  meshio::MeshIOAscii iohandler;
  iohandler.filename("data/tri3.mesh");
  iohandler.read(&mesh);

  spatialdata::geocoords::CSCart cs;
  cs.setSpaceDim(2);
  mesh.coordsys(&cs);
  spatialdata::units::Nondimensional normalizer;
  normalizer.lengthScale(lengthScale);
  mesh.nondimensionalize(normalizer);

  // Get vertices
  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  // Check nondimensional coordinates
  CoordsVisitor coordsVisitor(dmMesh);
  const PetscScalar* coordsArray = coordsVisitor.localArray();

  const PylithScalar tolerance = 1.0e-06;
  for(PetscInt v = vStart, i = 0; v < vEnd; ++v) {
    CPPUNIT_ASSERT_EQUAL(spaceDim, coordsVisitor.sectionDof(v));
    const PetscInt off = coordsVisitor.sectionOffset(v);
    for (int iDim=0; iDim < spaceDim; ++iDim, ++i) {
      const PylithScalar coordE = coordinates[i] / lengthScale;
      if (fabs(coordE) < 1.0) {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(coordE, coordsArray[off+iDim], tolerance);
      } else {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, coordsArray[off+iDim]/coordE, tolerance);
      } // if/else
    } // for
  } // for
  
  PYLITH_METHOD_END;
} // testNondimensionalize


// End of file 
