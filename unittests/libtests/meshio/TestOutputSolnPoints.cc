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
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestOutputSolnPoints.hh" // Implementation of class methods

#include "pylith/meshio/OutputSolnPoints.hh"

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart

#include <string.h> // USES strcmp()

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestOutputSolnPoints );

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::RealSection RealSection;

// ----------------------------------------------------------------------
// Test constructor
void
pylith::meshio::TestOutputSolnPoints::testConstructor(void)
{ // testConstructor
  OutputSolnPoints output;
} // testConstructor


// ----------------------------------------------------------------------
// Test setupInterpolator() for 2D points.
void
pylith::meshio::TestOutputSolnPoints::testSetupInterpolator2D(void)
{ // testSetupInterpolator2D
  const char* meshFilename = "data/quad4.mesh";
  const int numPoints = 7;
  const PylithScalar points[14] = { 
    0.0,  0.1, // interior points
    0.3,  0.4,
   -0.6, -0.7,
   -1.0,  0.9,
   -0.3,  0.8,
    0.3,  0.99999, // edge point
   -0.999999,  0.99999, // vertex point
  };
  const int nvertices = numPoints;
  const int verticesE[7] = { 7, 8, 9, 10, 11, 12, 13 };
  const int ncells = numPoints;
  const int ncorners = 1;
  const int cellsE[7] = { 7, 8, 9, 10, 11, 12, 13 };
  const int spaceDim = 2;

  topology::Mesh mesh;
  spatialdata::geocoords::CSCart cs;
  cs.setSpaceDim(spaceDim);
  cs.initialize();
  mesh.coordsys(&cs);
  MeshIOAscii iohandler;
  iohandler.filename(meshFilename);
  iohandler.read(&mesh);

  OutputSolnPoints output;
  output.setupInterpolator(&mesh, points, numPoints, spaceDim);

  const topology::Mesh& pointsMesh = output.pointsMesh();
  const ALE::Obj<SieveMesh>& sievePointsMesh = pointsMesh.sieveMesh();
  CPPUNIT_ASSERT(!sievePointsMesh.isNull());

  pointsMesh.view("POINTS MESH");

  // Check vertices
  const ALE::Obj<SieveMesh::label_sequence>& vertices = 
    sievePointsMesh->depthStratum(0);
  const SieveMesh::label_sequence::iterator verticesEnd = vertices->end();
  CPPUNIT_ASSERT_EQUAL(nvertices, int(vertices->size()));
  int ipt = 0;
  for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != verticesEnd;
       ++v_iter, ++ipt)
    CPPUNIT_ASSERT_EQUAL(verticesE[ipt], *v_iter);

  // Check cells
  const ALE::Obj<SieveMesh::label_sequence>& cells = 
    sievePointsMesh->heightStratum(0);
  const SieveMesh::label_sequence::iterator cellsEnd = cells->end();
  const ALE::Obj<SieveMesh::sieve_type>& sieve = sievePointsMesh->getSieve();
  assert(!sieve.isNull());

  CPPUNIT_ASSERT_EQUAL(ncells, int(cells->size()));

  ALE::ISieveVisitor::PointRetriever<SieveMesh::sieve_type> pV(sieve->getMaxConeSize());
  int i = 0;
  for (SieveMesh::label_sequence::iterator c_iter=cells->begin(); c_iter != cellsEnd; ++c_iter) {
    sieve->cone(*c_iter, pV);
    const SieveMesh::point_type* cone = pV.getPoints();
    CPPUNIT_ASSERT_EQUAL(ncorners, (int) pV.getSize());
    for(int p = 0; p < pV.getSize(); ++p, ++i) {
      CPPUNIT_ASSERT_EQUAL(cellsE[i], cone[p]);
    }
    pV.clear();
  } // for
} // testSetupInterpolator2D


// ----------------------------------------------------------------------
// Test setupInterpolator() for 3D points.
void
pylith::meshio::TestOutputSolnPoints::testSetupInterpolator3D(void)
{ // testSetupInterpolator3D
  const char* meshFilename = "data/hex8.mesh";
  const int numPoints = 7;
  const PylithScalar points[21] = { 
    0.0001,  0.1,  0.2, // interior points
    0.3, -0.4,  0.5,
    0.6,  0.7, -0.8,
   -0.9999,  0.1, -0.3,
    0.3,  0.8,  0.5,
    0.3,  0.9999,  0.9999, // edge point
    0.000001,  0.9999,  -0.999999, // vertex point
  };
  const int nvertices = numPoints;
  const int verticesE[7] = { 7, 8, 9, 10, 11, 12, 13 };
  const int ncells = numPoints;
  const int ncorners = 1;
  const int cellsE[7] = { 7, 8, 9, 10, 11, 12, 13 };
  const int spaceDim = 3;

  topology::Mesh mesh;
  spatialdata::geocoords::CSCart cs;
  cs.setSpaceDim(spaceDim);
  cs.initialize();
  mesh.coordsys(&cs);
  MeshIOAscii iohandler;
  iohandler.filename(meshFilename);
  iohandler.read(&mesh);

  OutputSolnPoints output;
  output.setupInterpolator(&mesh, points, numPoints, spaceDim);

  const topology::Mesh& pointsMesh = output.pointsMesh();
  const ALE::Obj<SieveMesh>& sievePointsMesh = pointsMesh.sieveMesh();
  CPPUNIT_ASSERT(!sievePointsMesh.isNull());

  pointsMesh.view("POINTS MESH");

  // Check vertices
  const ALE::Obj<SieveMesh::label_sequence>& vertices = 
    sievePointsMesh->depthStratum(0);
  const SieveMesh::label_sequence::iterator verticesEnd = vertices->end();
  CPPUNIT_ASSERT_EQUAL(nvertices, int(vertices->size()));
  int ipt = 0;
  for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != verticesEnd;
       ++v_iter, ++ipt)
    CPPUNIT_ASSERT_EQUAL(verticesE[ipt], *v_iter);

  // Check cells
  const ALE::Obj<SieveMesh::label_sequence>& cells = 
    sievePointsMesh->heightStratum(0);
  const SieveMesh::label_sequence::iterator cellsEnd = cells->end();
  const ALE::Obj<SieveMesh::sieve_type>& sieve = sievePointsMesh->getSieve();
  assert(!sieve.isNull());

  CPPUNIT_ASSERT_EQUAL(ncells, int(cells->size()));

  ALE::ISieveVisitor::PointRetriever<SieveMesh::sieve_type> pV(sieve->getMaxConeSize());
  int i = 0;
  for (SieveMesh::label_sequence::iterator c_iter=cells->begin(); c_iter != cellsEnd; ++c_iter) {
    sieve->cone(*c_iter, pV);
    const SieveMesh::point_type *cone = pV.getPoints();
    CPPUNIT_ASSERT_EQUAL(ncorners, (int) pV.getSize());
    for(int p = 0; p < pV.getSize(); ++p, ++i) {
      CPPUNIT_ASSERT_EQUAL(cellsE[i], cone[p]);
    }
    pV.clear();
  } // for
} // testSetupInterpolator3D


// End of file 
