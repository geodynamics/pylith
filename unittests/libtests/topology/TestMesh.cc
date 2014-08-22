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
// Copyright (c) 2010-2014 University of California, Davis
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
  CPPUNIT_ASSERT(!mesh._dmMesh);
  CPPUNIT_ASSERT_EQUAL(0, mesh.dimension());
  CPPUNIT_ASSERT_EQUAL(false, mesh.debug());
  MPI_Comm_compare(PETSC_COMM_WORLD, mesh.comm(), &result);
  CPPUNIT_ASSERT_EQUAL(int(MPI_IDENT), result);

  int dim = 2;
  Mesh mesh2(dim);
  CPPUNIT_ASSERT(mesh2._dmMesh);
  CPPUNIT_ASSERT_EQUAL(dim, mesh2.dimension());
  MPI_Comm_compare(PETSC_COMM_WORLD, mesh2.comm(), &result);
  CPPUNIT_ASSERT_EQUAL(int(MPI_CONGRUENT), result);

  dim = 1;
  Mesh mesh3(dim, PETSC_COMM_SELF);
  CPPUNIT_ASSERT(mesh3._dmMesh);
  CPPUNIT_ASSERT_EQUAL(dim, mesh3.dimension());
  MPI_Comm_compare(PETSC_COMM_WORLD, mesh3.comm(), &result);
  CPPUNIT_ASSERT_EQUAL(int(MPI_CONGRUENT), result);

  PYLITH_METHOD_END;
} // testConstructor

// ----------------------------------------------------------------------
// Test dmMesh().
void
pylith::topology::TestMesh::testDMMesh(void)
{ // testDMMesh
  PYLITH_METHOD_BEGIN;

  const int dim = 2;
  PetscInt dmDim;
  Mesh mesh(dim);
  
  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  PetscErrorCode err = DMGetDimension(dmMesh, &dmDim);PYLITH_CHECK_ERROR(err);
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
  int result = 0;
  MPI_Comm_compare(PETSC_COMM_WORLD, mesh.comm(), &result);
  CPPUNIT_ASSERT_EQUAL(int(MPI_IDENT), result);


  Mesh mesh2(2, PETSC_COMM_SELF);
  result = 0;
  MPI_Comm_compare(PETSC_COMM_SELF, mesh2.comm(), &result);
  CPPUNIT_ASSERT_EQUAL(int(MPI_CONGRUENT), result);

  PYLITH_METHOD_END;
} // testComm


// ----------------------------------------------------------------------
// Test view().
void
pylith::topology::TestMesh::testView(void)
{ // testView
  PYLITH_METHOD_BEGIN;

  const char* filename = "data/tri3.mesh";

  Mesh mesh;
  meshio::MeshIOAscii iohandler;
  iohandler.filename(filename);
  iohandler.read(&mesh);

  mesh.view();
  mesh.view(":mesh.view:ascii_info_detail");
  mesh.view("vtk:mesh.vtk:ascii_vtk");
  
  PYLITH_METHOD_END;
} // testView


// End of file 
