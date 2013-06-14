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

#include "TestMesh.hh" // Implementation of class methods
#include "pylith/topology/Mesh.hh" // USES Mesh

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
  Mesh mesh;
  CPPUNIT_ASSERT(mesh._mesh.isNull());
  CPPUNIT_ASSERT_EQUAL(0, mesh.dimension());
  CPPUNIT_ASSERT_EQUAL(false, mesh.debug());
  CPPUNIT_ASSERT_EQUAL(PETSC_COMM_WORLD, mesh.comm());
  
  Mesh mesh2(2);
  CPPUNIT_ASSERT(!mesh2._mesh.isNull());
  CPPUNIT_ASSERT_EQUAL(2, mesh2.dimension());
  CPPUNIT_ASSERT_EQUAL(PETSC_COMM_WORLD, mesh2.comm());

  Mesh mesh3(1, PETSC_COMM_SELF);
  CPPUNIT_ASSERT(!mesh3._mesh.isNull());
  CPPUNIT_ASSERT_EQUAL(1, mesh3.dimension());
  CPPUNIT_ASSERT_EQUAL(PETSC_COMM_SELF, mesh3.comm());
} // testConstructor

// ----------------------------------------------------------------------
// Test createSieveMesh().
void
pylith::topology::TestMesh::testCreateSieveMesh(void)
{ // testCreateSieveMesh
  Mesh mesh;
  CPPUNIT_ASSERT(mesh._mesh.isNull());

  int dim = 2;
  mesh.createSieveMesh(dim);
  CPPUNIT_ASSERT(!mesh._mesh.isNull());
  CPPUNIT_ASSERT_EQUAL(dim, mesh.dimension());

  dim = 1;
  mesh.createSieveMesh(dim);
  CPPUNIT_ASSERT(!mesh._mesh.isNull());
  CPPUNIT_ASSERT_EQUAL(dim, mesh.dimension());
} // testCreateMeshSieve

// ----------------------------------------------------------------------
// Test sieveMesh().
void
pylith::topology::TestMesh::testSieveMesh(void)
{ // testSieveMesh
  const int dim = 2;
  Mesh mesh(dim);
  
  const ALE::Obj<Mesh::SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  CPPUNIT_ASSERT_EQUAL(dim, sieveMesh->getDimension());
} // testSieveMesh

// ----------------------------------------------------------------------
// Test coordsys().
void
pylith::topology::TestMesh::testCoordsys(void)
{ // testCoordsys
  Mesh mesh;

  spatialdata::geocoords::CSCart cs;
  cs.setSpaceDim(2);

  mesh.coordsys(&cs);

  CPPUNIT_ASSERT_EQUAL(cs.spaceDim(), mesh.coordsys()->spaceDim());
} // testCoordsys

// ----------------------------------------------------------------------
// Test debug().
void
pylith::topology::TestMesh::testDebug(void)
{ // testDebug
  Mesh mesh;
  CPPUNIT_ASSERT_EQUAL(false, mesh.debug());

  mesh.debug(true);
  CPPUNIT_ASSERT_EQUAL(true, mesh.debug());
} // testDebug

// ----------------------------------------------------------------------
// Test dimension().
void
pylith::topology::TestMesh::testDimension(void)
{ // testDimension
  Mesh mesh;
  CPPUNIT_ASSERT_EQUAL(0, mesh.dimension());

  const int dim = 2;
  Mesh mesh2(dim);
  CPPUNIT_ASSERT_EQUAL(dim, mesh2.dimension());
} // testDimension

// ----------------------------------------------------------------------
// Test comm().
void
pylith::topology::TestMesh::testComm(void)
{ // testComm
  Mesh mesh;
  CPPUNIT_ASSERT_EQUAL(PETSC_COMM_WORLD, mesh.comm());

  mesh.comm(PETSC_COMM_SELF);
  CPPUNIT_ASSERT_EQUAL(PETSC_COMM_SELF, mesh.comm());
} // testComm

// ----------------------------------------------------------------------
// Test nondimensionalize().
void
pylith::topology::TestMesh::testNondimensionalize(void)
{ // testNondimensionalizer
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

  const ALE::Obj<Mesh::SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  const ALE::Obj<Mesh::SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());
  const Mesh::SieveMesh::label_sequence::iterator verticesBegin =
    vertices->begin();
  const Mesh::SieveMesh::label_sequence::iterator verticesEnd =
    vertices->end();
  CPPUNIT_ASSERT_EQUAL(numVertices, int(vertices->size()));

  // Check nondimensional coordinates
  const ALE::Obj<Mesh::RealSection>& coordsField =
    sieveMesh->getRealSection("coordinates");
  CPPUNIT_ASSERT(!coordsField.isNull());
  CPPUNIT_ASSERT_EQUAL(spaceDim, 
		       coordsField->getFiberDimension(*verticesBegin));
  int i = 0;
  for(Mesh::SieveMesh::label_sequence::iterator v_iter = verticesBegin;
      v_iter != verticesEnd;
      ++v_iter) {
    const PylithScalar* coordsVertex = coordsField->restrictPoint(*v_iter);
    CPPUNIT_ASSERT(0 != coordsVertex);
    const PylithScalar tolerance = 1.0e-06;
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      const PylithScalar coordE = coordinates[i++] / lengthScale;
      if (coordE < 1.0)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(coordE, coordsVertex[iDim],
				     tolerance);
      else
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, coordsVertex[iDim]/coordE,
				     tolerance);
    } // for
  } // for
  
  // Check dimensioned coordinates
  const ALE::Obj<Mesh::RealSection>& coordsDimField =
    sieveMesh->getRealSection("coordinates_dimensioned");
  CPPUNIT_ASSERT(!coordsDimField.isNull());
  CPPUNIT_ASSERT_EQUAL(spaceDim, 
		       coordsDimField->getFiberDimension(*verticesBegin));
  i = 0;
  for(Mesh::SieveMesh::label_sequence::iterator v_iter = verticesBegin;
      v_iter != verticesEnd;
      ++v_iter) {
    const PylithScalar* coordsVertex = coordsDimField->restrictPoint(*v_iter);
    CPPUNIT_ASSERT(0 != coordsVertex);
    const PylithScalar tolerance = 1.0e-06;
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      const PylithScalar coordE = coordinates[i++];
      if (coordE < 1.0)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(coordE, coordsVertex[iDim],
				     tolerance);
      else
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, coordsVertex[iDim]/coordE,
				     tolerance);
    } // for
  } // for
} // testNondimensionalize


// End of file 
