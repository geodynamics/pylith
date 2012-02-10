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

#include "data/OutputSolnPointsDataTri3.hh"
#include "data/OutputSolnPointsDataQuad4.hh"
#include "data/OutputSolnPointsDataTet4.hh"
#include "data/OutputSolnPointsDataHex8.hh"

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
// Test setupInterpolator for tri3 mesh.
void
pylith::meshio::TestOutputSolnPoints::testSetupInterpolatorTri3(void)
{ // testSetupInterpolatorTri3
  OutputSolnPoints output;
  OutputSolnPointsDataTri3 data;

  _testSetupInterpolator(data);
} // testSetupInterpolatorTri3


// ----------------------------------------------------------------------
// Test setupInterpolator for quad4 mesh.
void
pylith::meshio::TestOutputSolnPoints::testSetupInterpolatorQuad4(void)
{ // testSetupInterpolatorQuad4
  OutputSolnPoints output;
  OutputSolnPointsDataQuad4 data;

  _testSetupInterpolator(data);
} // testSetupInterpolatorQuad4


// ----------------------------------------------------------------------
// Test setupInterpolator for tet4 mesh.
void
pylith::meshio::TestOutputSolnPoints::testSetupInterpolatorTet4(void)
{ // testSetupInterpolatorTet4
  OutputSolnPoints output;
  OutputSolnPointsDataTet4 data;

  _testSetupInterpolator(data);
} // testSetupInterpolatorTet4


// ----------------------------------------------------------------------
// Test setupInterpolator for hex8 mesh.
void
pylith::meshio::TestOutputSolnPoints::testSetupInterpolatorHex8(void)
{ // testSetupInterpolatorHex8
  OutputSolnPoints output;
  OutputSolnPointsDataHex8 data;

  _testSetupInterpolator(data);
} // testSetupInterpolatorHex8


// ----------------------------------------------------------------------
// Test setupInterpolator().
void
pylith::meshio::TestOutputSolnPoints::_testSetupInterpolator(const OutputSolnPointsData& data)
{ // _testSetupInterpolator
  const int numPoints = data.numPoints;
  const int spaceDim = data.spaceDim;

  const int numVerticesE = numPoints;
  const int numCellsE = numPoints;
  const int numCornersE = 1;

  topology::Mesh mesh;
  spatialdata::geocoords::CSCart cs;
  cs.setSpaceDim(spaceDim);
  cs.initialize();
  mesh.coordsys(&cs);
  MeshIOAscii iohandler;
  iohandler.filename(data.meshFilename);
  iohandler.read(&mesh);

  OutputSolnPoints output;
  CPPUNIT_ASSERT(data.points);
  output.setupInterpolator(&mesh, data.points, numPoints, spaceDim);

  const topology::Mesh& pointsMesh = output.pointsMesh();
  const ALE::Obj<SieveMesh>& sievePointsMesh = pointsMesh.sieveMesh();
  CPPUNIT_ASSERT(!sievePointsMesh.isNull());

  //pointsMesh.view("POINTS MESH"); // DEBUGGING

  // Check vertices
  const ALE::Obj<SieveMesh::label_sequence>& vertices = 
    sievePointsMesh->depthStratum(0);
  const SieveMesh::label_sequence::iterator verticesEnd = vertices->end();
  CPPUNIT_ASSERT_EQUAL(numVerticesE, int(vertices->size()));
  int ipt = 0;
  for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != verticesEnd;
       ++v_iter, ++ipt) {
    const int vertexE = numCellsE + ipt;
    CPPUNIT_ASSERT_EQUAL(vertexE, *v_iter);
  } // for

  // Check cells
  const ALE::Obj<SieveMesh::label_sequence>& cells = 
    sievePointsMesh->heightStratum(0);
  const SieveMesh::label_sequence::iterator cellsEnd = cells->end();
  const ALE::Obj<SieveMesh::sieve_type>& sieve = sievePointsMesh->getSieve();
  assert(!sieve.isNull());

  CPPUNIT_ASSERT_EQUAL(numCellsE, int(cells->size()));

  ALE::ISieveVisitor::PointRetriever<SieveMesh::sieve_type> pV(sieve->getMaxConeSize());
  int i = 0;
  for (SieveMesh::label_sequence::iterator c_iter=cells->begin(); c_iter != cellsEnd; ++c_iter) {
    sieve->cone(*c_iter, pV);
    const SieveMesh::point_type* cone = pV.getPoints();
    CPPUNIT_ASSERT_EQUAL(numCornersE, (int) pV.getSize());
    for(int p = 0; p < pV.getSize(); ++p, ++i) {
      const int coneE = numCellsE+i;
      CPPUNIT_ASSERT_EQUAL(coneE, cone[p]);
    } // for
    pV.clear();
  } // for
} // _testSetupInterpolator


// End of file 
