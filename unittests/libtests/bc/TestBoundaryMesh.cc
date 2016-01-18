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

#include "TestBoundaryMesh.hh" // Implementation of class methods

#include "data/BoundaryMeshData.hh" // USES BoundaryMeshData

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps::nondimensionalize()
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorSubMesh.hh" // USES SubMeshIS
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/faults/FaultCohesiveKin.hh" // USES FaultsCohesiveKin

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestBoundaryMesh::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  _data = 0;

  PYLITH_METHOD_END;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::bc::TestBoundaryMesh::tearDown(void)
{ // tearDown
  PYLITH_METHOD_BEGIN;

  delete _data; _data = 0;

  PYLITH_METHOD_END;
} // tearDown

// ----------------------------------------------------------------------
// Test submesh() without fault.
void
pylith::bc::TestBoundaryMesh::testSubmesh(void)
{ // testSubmesh
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_data);
  PetscErrorCode err;

  topology::Mesh mesh;
  meshio::MeshIOAscii iohandler;
  iohandler.filename(_data->filename);
  iohandler.read(&mesh);

  // Set up coordinates
  spatialdata::geocoords::CSCart cs;
  spatialdata::units::Nondimensional normalizer;
  cs.setSpaceDim(mesh.dimension());
  cs.initialize();
  mesh.coordsys(&cs);
  topology::MeshOps::nondimensionalize(&mesh, normalizer);

  // Create submesh
  topology::Mesh submesh(mesh, _data->bcLabel);
  PetscDM dmMesh = submesh.dmMesh();CPPUNIT_ASSERT(dmMesh);

  // Check vertices
  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();
  CPPUNIT_ASSERT_EQUAL(_data->numVerticesNoFault, verticesStratum.size());

  // Check cells
  topology::Stratum cellsStratum(dmMesh, topology::Stratum::HEIGHT, 1);
  const PetscInt cStart = cellsStratum.begin();
  const PetscInt cEnd = cellsStratum.end();
  CPPUNIT_ASSERT_EQUAL(_data->numCells, cellsStratum.size());
  for (PetscInt c = cStart, index = 0; c < cEnd; ++c) {
    PetscInt  vertices[32];
    PetscInt *closure = NULL;
    PetscInt  closureSize, numVertices = 0;

    err = DMPlexGetTransitiveClosure(dmMesh, c, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
    for (PetscInt p = 0; p < closureSize*2; p += 2) {
      if ((closure[p] >= vStart) && (closure[p] < vEnd)) vertices[numVertices++] = closure[p];
    } // for
    err = DMPlexRestoreTransitiveClosure(dmMesh, c, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(_data->numCorners, numVertices);
  } // for

  PYLITH_METHOD_END;
} // testSubmesh

// ----------------------------------------------------------------------
// Test submesh() with fault.
void
pylith::bc::TestBoundaryMesh::testSubmeshFault(void)
{ // testSubmeshFault
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_data);
  PetscErrorCode err;

  topology::Mesh mesh;
  meshio::MeshIOAscii iohandler;
  iohandler.filename(_data->filename);
  iohandler.read(&mesh);

  // Set up coordinates
  spatialdata::geocoords::CSCart cs;
  spatialdata::units::Nondimensional normalizer;
  cs.setSpaceDim(mesh.dimension());
  cs.initialize();
  mesh.coordsys(&cs);
  topology::MeshOps::nondimensionalize(&mesh, normalizer);

  // Adjust topology
  faults::FaultCohesiveKin fault;
  PetscInt firstFaultVertex = 0;
  PetscInt firstLagrangeVertex, firstFaultCell;

  err = DMGetStratumSize(mesh.dmMesh(), _data->faultLabel, 1, &firstLagrangeVertex);PYLITH_CHECK_ERROR(err);
  firstFaultCell = firstLagrangeVertex;
  fault.label(_data->faultLabel);
  fault.id(_data->faultId);
  fault.adjustTopology(&mesh, &firstFaultVertex, &firstLagrangeVertex, &firstFaultCell);

  // Create submesh
  topology::Mesh submesh(mesh, _data->bcLabel);
  PetscDM dmMesh = submesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
#if 0 // DEBUGGING
  PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_INFO_DETAIL);
  DMView(mesh.dmMesh(), PETSC_VIEWER_STDOUT_WORLD);
  PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD);
#endif

  // Check vertices
  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();
  CPPUNIT_ASSERT_EQUAL(_data->numVerticesFault, verticesStratum.size());

  // Check cells
  topology::Stratum cellsStratum(dmMesh, topology::Stratum::HEIGHT, 1);
  const PetscInt cStart = cellsStratum.begin();
  const PetscInt cEnd = cellsStratum.end();
  CPPUNIT_ASSERT_EQUAL(_data->numCells, cellsStratum.size());

  for (PetscInt c = cStart, index = 0; c < cEnd; ++c) {
    PetscInt  vertices[32];
    PetscInt *closure = NULL;
    PetscInt  closureSize, numVertices = 0;

    err = DMPlexGetTransitiveClosure(dmMesh, c, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
    for (PetscInt p = 0; p < closureSize*2; p += 2) {
      if ((closure[p] >= vStart) && (closure[p] < vEnd)) vertices[numVertices++] = closure[p];
    } // for
    err = DMPlexRestoreTransitiveClosure(dmMesh, c, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(_data->numCorners, numVertices);
  } // for

  PYLITH_METHOD_END;
} // testSubmeshFault


// End of file 
