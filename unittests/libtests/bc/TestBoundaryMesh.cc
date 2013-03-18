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

#include "TestBoundaryMesh.hh" // Implementation of class methods

#include "data/BoundaryMeshData.hh" // USES BoundaryMeshData

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/SubMesh.hh" // USES SubMesh
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/faults/FaultCohesiveKin.hh" // USES FaultsCohesiveKin

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestBoundaryMesh::setUp(void)
{ // setUp
  _data = 0;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::bc::TestBoundaryMesh::tearDown(void)
{ // tearDown
  delete _data; _data = 0;
} // tearDown

// ----------------------------------------------------------------------
// Test submesh() without fault.
void
pylith::bc::TestBoundaryMesh::testSubmesh(void)
{ // testSubmesh
  CPPUNIT_ASSERT(0 != _data);

  topology::Mesh mesh;
  PetscErrorCode err;

  meshio::MeshIOAscii iohandler;
  iohandler.filename(_data->filename);
  iohandler.read(&mesh);

  // Set up coordinates
  spatialdata::geocoords::CSCart cs;
  spatialdata::units::Nondimensional normalizer;
  cs.setSpaceDim(mesh.dimension());
  cs.initialize();
  mesh.coordsys(&cs);
  mesh.nondimensionalize(normalizer);

  // Create submesh
  topology::SubMesh submesh(mesh, _data->bcLabel);
  DM                dmMesh = submesh.dmMesh();
  CPPUNIT_ASSERT(dmMesh);

  // Check vertices
  IS              subpointIS;
  const PetscInt *subpointMap;
  PetscInt        vStart, vEnd;

  err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);
  CPPUNIT_ASSERT_EQUAL(_data->numVerticesNoFault, vEnd-vStart);
  err = DMPlexCreateSubpointIS(dmMesh, &subpointIS);CHECK_PETSC_ERROR(err);
  err = ISGetIndices(subpointIS, &subpointMap);CHECK_PETSC_ERROR(err);
  for (PetscInt v = vStart; v < vEnd; ++v)
    CPPUNIT_ASSERT_EQUAL(_data->verticesNoFault[v-vStart], subpointMap[v]);

  // Check cells
  PetscInt cStart, cEnd;

  err = DMPlexGetHeightStratum(dmMesh, 1, &cStart, &cEnd);CHECK_PETSC_ERROR(err);
  CPPUNIT_ASSERT_EQUAL(_data->numCells, cEnd-cStart);

  for (PetscInt c = cStart, index = 0; c < cEnd; ++c) {
    PetscInt *closure = NULL;
    PetscInt  closureSize, numVertices = 0;

    err = DMPlexGetTransitiveClosure(dmMesh, c, PETSC_TRUE, &closureSize, &closure);CHECK_PETSC_ERROR(err);
    for (PetscInt p = 0; p < closureSize*2; p += 2) {
      if ((closure[p] >= vStart) && (closure[p] < vEnd)) closure[numVertices++] = closure[p];
    }
    CPPUNIT_ASSERT_EQUAL(_data->numCorners, numVertices);
    for (PetscInt v = 0; v < numVertices; ++v, ++index)
      CPPUNIT_ASSERT_EQUAL(_data->cellsNoFault[index], subpointMap[closure[v]]);
    err = DMPlexRestoreTransitiveClosure(dmMesh, c, PETSC_TRUE, &closureSize, &closure);CHECK_PETSC_ERROR(err);
  } // for
  err = ISRestoreIndices(subpointIS, &subpointMap);CHECK_PETSC_ERROR(err);
  err = ISDestroy(&subpointIS);CHECK_PETSC_ERROR(err);
} // testSubmesh

// ----------------------------------------------------------------------
// Test submesh() with fault.
void
pylith::bc::TestBoundaryMesh::testSubmeshFault(void)
{ // testSubmeshFault
  CPPUNIT_ASSERT(0 != _data);
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
  mesh.nondimensionalize(normalizer);

  // Adjust topology
  faults::FaultCohesiveKin fault;
  PetscInt firstFaultVertex = 0;
  PetscInt firstLagrangeVertex, firstFaultCell;

  err = DMPlexGetStratumSize(mesh.dmMesh(), _data->faultLabel, 1, &firstLagrangeVertex);CHECK_PETSC_ERROR(err);
  firstFaultCell = firstLagrangeVertex;
  fault.label(_data->faultLabel);
  fault.id(_data->faultId);
  fault.adjustTopology(&mesh, &firstFaultVertex, &firstLagrangeVertex, &firstFaultCell, _flipFault);

  // Create submesh
  topology::SubMesh submesh(mesh, _data->bcLabel);
  DM                dmMesh = submesh.dmMesh();
  CPPUNIT_ASSERT(dmMesh);

  // Check vertices
  IS              subpointIS;
  const PetscInt *subpointMap;
  PetscInt        vStart, vEnd;

  err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);
  CPPUNIT_ASSERT_EQUAL(_data->numVerticesFault, vEnd-vStart);
  err = DMPlexCreateSubpointIS(dmMesh, &subpointIS);CHECK_PETSC_ERROR(err);
  err = ISGetIndices(subpointIS, &subpointMap);CHECK_PETSC_ERROR(err);
  for (PetscInt v = vStart; v < vEnd; ++v)
    CPPUNIT_ASSERT_EQUAL(_data->verticesFault[v-vStart], subpointMap[v]);
  // Check cells
  PetscInt cStart, cEnd;

  err = DMPlexGetHeightStratum(dmMesh, 1, &cStart, &cEnd);CHECK_PETSC_ERROR(err);
  CPPUNIT_ASSERT_EQUAL(_data->numCells, cEnd-cStart);

  for (PetscInt c = cStart, index = 0; c < cEnd; ++c) {
    PetscInt *closure = NULL;
    PetscInt  closureSize, numVertices = 0;

    err = DMPlexGetTransitiveClosure(dmMesh, c, PETSC_TRUE, &closureSize, &closure);CHECK_PETSC_ERROR(err);
    for (PetscInt p = 0; p < closureSize*2; p += 2) {
      if ((closure[p] >= vStart) && (closure[p] < vEnd)) closure[numVertices++] = closure[p];
    }
    CPPUNIT_ASSERT_EQUAL(_data->numCorners, numVertices);
    for (PetscInt v = 0; v < numVertices; ++v, ++index)
      CPPUNIT_ASSERT_EQUAL(_data->cellsFault[index], subpointMap[closure[v]]);
    err = DMPlexRestoreTransitiveClosure(dmMesh, c, PETSC_TRUE, &closureSize, &closure);CHECK_PETSC_ERROR(err);
  } // for
  err = ISRestoreIndices(subpointIS, &subpointMap);CHECK_PETSC_ERROR(err);
  err = ISDestroy(&subpointIS);CHECK_PETSC_ERROR(err);
} // testSubmeshFault


// End of file 
