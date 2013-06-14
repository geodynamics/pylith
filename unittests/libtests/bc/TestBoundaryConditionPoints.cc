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

#include "TestBoundaryConditionPoints.hh" // Implementation of class methods

#include "pylith/bc/PointForce.hh" // USES PointForce

#include "data/PointForceDataTri3.hh" // USES PointForceDataTri3

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestBoundaryConditionPoints );

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;

// ----------------------------------------------------------------------
// Test _getPoints().
void
pylith::bc::TestBoundaryConditionPoints::testGetPoints(void)
{ // testGetPoints
  topology::Mesh mesh;
  PointForce bc;
  PointForceDataTri3 data;

  meshio::MeshIOAscii iohandler;
  iohandler.filename(data.meshFilename);
  iohandler.read(&mesh);

  spatialdata::geocoords::CSCart cs;
  spatialdata::units::Nondimensional normalizer;
  cs.setSpaceDim(mesh.dimension());
  cs.initialize();
  mesh.coordsys(&cs);
  mesh.nondimensionalize(normalizer);

  bc.label(data.label);
  bc.BoundaryConditionPoints::_getPoints(mesh);

  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  
  const int numCells = sieveMesh->heightStratum(0)->size();
  const size_t numPoints = data.numForcePts;

  // Check points
  const int offset = numCells;
  CPPUNIT_ASSERT_EQUAL(numPoints, bc._points.size());
  for (int i=0; i < numPoints; ++i)
    CPPUNIT_ASSERT_EQUAL(data.forcePoints[i]+offset, bc._points[i]);
} // testGetPoints


// End of file 
