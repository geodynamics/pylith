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

#include "TestBoundaryConditionPoints.hh" // Implementation of class methods

#include "pylith/bc/PointForce.hh" // USES PointForce

#include "data/PointForceDataTri3.hh" // USES PointForceDataTri3

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps::nondimensionalize()
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestBoundaryConditionPoints );

// ----------------------------------------------------------------------
// Test _getPoints().
void
pylith::bc::TestBoundaryConditionPoints::testGetPoints(void)
{ // testGetPoints
  PYLITH_METHOD_BEGIN;

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
  topology::MeshOps::nondimensionalize(&mesh, normalizer);

  bc.label(data.label);
  bc.BoundaryConditionPoints::_getPoints(mesh);

  const PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  topology::Stratum heightStratum(dmMesh, topology::Stratum::HEIGHT, 0);
  const PetscInt numCells = heightStratum.size();
  const size_t numPoints = data.numForcePts;

  // Check points
  const int offset = numCells;
  CPPUNIT_ASSERT_EQUAL(numPoints, bc._points.size());
  for (int i=0; i < numPoints; ++i)
    CPPUNIT_ASSERT_EQUAL(data.forcePoints[i]+offset, bc._points[i]);

  PYLITH_METHOD_END;
} // testGetPoints


// End of file 
