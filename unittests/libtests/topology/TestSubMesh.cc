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

#include "TestSubMesh.hh" // Implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/SubMesh.hh" // USES SubMesh
#include "pylith/topology/Stratum.hh" // USES Stratum

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::topology::TestSubMesh );

// ----------------------------------------------------------------------
namespace pylith {
  namespace topology {
    namespace _TestSubMesh {
      const int cellDim = 2;
      const int nvertices = 4;
      const int ncells = 2;
      const int ncorners = 3;
      const int cells[] = {
	0, 1, 3,
	0, 3, 2,
      };
      const PylithScalar coordinates[] = {
	0.0, 0.0,
	1.0, 0.0,
	0.0, 1.0,
	1.0, 1.0,
      };
      const char* label = "bc";
      const int groupSize = 3;
      const int groupVertices[groupSize] = {
	1, 2, 3,
      };
      const int submeshNumVertices = groupSize;
      const int submeshVertices[submeshNumVertices] = {
	2, 3, 4,
      };
      const int submeshNumCells = 2;
      const int submeshCells[submeshNumCells] = {
	0, 1,
      };
    } // _TestSubMesh
  } // topology
} // pylith

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::topology::TestSubMesh::testConstructor(void)
{ // testConstructor
  SubMesh mesh;
  CPPUNIT_ASSERT_EQUAL(0, mesh.dimension());
  CPPUNIT_ASSERT_EQUAL(false, mesh.debug());
} // testConstructor

// ----------------------------------------------------------------------
// Test constructor w/mesh.
void
pylith::topology::TestSubMesh::testConstructorMesh(void)
{ // testConstructorMesh
  Mesh mesh2D;
  _buildMesh(&mesh2D);
  
  SubMesh mesh(mesh2D, _TestSubMesh::label);
  CPPUNIT_ASSERT_EQUAL(_TestSubMesh::cellDim-1, mesh.dimension());
  MPI_Comm commA, commB;
  PetscObjectGetComm((PetscObject)mesh2D.dmMesh(), &commA);
  PetscObjectGetComm((PetscObject)mesh.dmMesh(), &commB);
  std::cout << "mesh2D comm: " << mesh2D.comm()
	    << ", dm comm: " << commA
	    << ", PETSC_COMM_WORLD: " << PETSC_COMM_WORLD
	    << ", PETSC_COMM_SELF: " << PETSC_COMM_SELF
	    << ", submesh comm: " << mesh.comm()
	    << ", dm comm: " << commB
	    << std::endl;
  CPPUNIT_ASSERT_EQUAL(PETSC_COMM_WORLD, mesh.comm());

  // Check vertices
  const PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();
  
  const PetscInt nvertices = _TestSubMesh::submeshNumVertices;
  CPPUNIT_ASSERT_EQUAL(nvertices, depthStratum.size());
  for (PetscInt v = vStart, iV=0; v < vEnd; ++v, ++iV) {
    CPPUNIT_ASSERT_EQUAL(_TestSubMesh::submeshVertices[iV], v);
  } // for

  // Check cells
  Stratum heightStratum(dmMesh, Stratum::HEIGHT, 0);
  const PetscInt cStart = heightStratum.begin();
  const PetscInt cEnd = heightStratum.end();
  
  const PetscInt ncells = _TestSubMesh::submeshNumCells;
  CPPUNIT_ASSERT_EQUAL(ncells, heightStratum.size());
  for (PetscInt c = cStart, iC=0; c < cEnd; ++c, ++iC) {
    CPPUNIT_ASSERT_EQUAL(_TestSubMesh::submeshCells[iC], c);
  } // for
} // testConstructorMesh

// ----------------------------------------------------------------------
// Test createSubMesh().
void
pylith::topology::TestSubMesh::testCreateSubMesh(void)
{ // testCreateSubMesh
  Mesh mesh2D;
  _buildMesh(&mesh2D);
  
  SubMesh mesh;
  mesh.createSubMesh(mesh2D, _TestSubMesh::label);
  CPPUNIT_ASSERT_EQUAL(_TestSubMesh::cellDim-1, mesh.dimension());
  CPPUNIT_ASSERT_EQUAL(PETSC_COMM_WORLD, mesh.comm());

  // Check vertices
  const PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();
  
  const PetscInt nvertices = _TestSubMesh::submeshNumVertices;
  CPPUNIT_ASSERT_EQUAL(nvertices, depthStratum.size());
  for (PetscInt v = vStart, iV=0; v < vEnd; ++v, ++iV) {
    CPPUNIT_ASSERT_EQUAL(_TestSubMesh::submeshVertices[iV], v);
  } // for

  // Check cells
  Stratum heightStratum(dmMesh, Stratum::HEIGHT, 0);
  const PetscInt cStart = heightStratum.begin();
  const PetscInt cEnd = heightStratum.end();
  
  const PetscInt ncells = _TestSubMesh::submeshNumCells;
  CPPUNIT_ASSERT_EQUAL(ncells, heightStratum.size());
  for (PetscInt c = cStart, iC=0; c < cEnd; ++c, ++iC) {
    CPPUNIT_ASSERT_EQUAL(_TestSubMesh::submeshCells[iC], c);
  } // for
} // testCreateSubMesh

// ----------------------------------------------------------------------
// Test coordsys().
void
pylith::topology::TestSubMesh::testCoordsys(void)
{ // testCoordsys
  Mesh mesh2D;
  _buildMesh(&mesh2D);

  SubMesh mesh(mesh2D, _TestSubMesh::label);

  CPPUNIT_ASSERT_EQUAL(_TestSubMesh::cellDim, mesh.coordsys()->spaceDim());
} // testCoordsys

// ----------------------------------------------------------------------
// Test debug().
void
pylith::topology::TestSubMesh::testDebug(void)
{ // testDebug
  SubMesh mesh;
  CPPUNIT_ASSERT_EQUAL(false, mesh.debug());

  mesh.debug(true);
  CPPUNIT_ASSERT_EQUAL(true, mesh.debug());
} // testDebug

// ----------------------------------------------------------------------
// Test dimension().
void
pylith::topology::TestSubMesh::testDimension(void)
{ // testDimension
  SubMesh mesh;
  CPPUNIT_ASSERT_EQUAL(0, mesh.dimension());

  Mesh mesh2D;
  _buildMesh(&mesh2D);
  SubMesh mesh2(mesh2D, _TestSubMesh::label);
  CPPUNIT_ASSERT_EQUAL(_TestSubMesh::cellDim-1, mesh2.dimension());
} // testDimension

// ----------------------------------------------------------------------
// Test comm().
void
pylith::topology::TestSubMesh::testComm(void)
{ // testComm
  SubMesh mesh;

  Mesh mesh2D;
  _buildMesh(&mesh2D);
  mesh.createSubMesh(mesh2D, _TestSubMesh::label);

  CPPUNIT_ASSERT_EQUAL(PETSC_COMM_WORLD, mesh.comm());
} // testComm

// ----------------------------------------------------------------------
// Test initialize().
void
pylith::topology::TestSubMesh::testInitialize(void)
{ // testInitialize
  Mesh mesh2D;
  _buildMesh(&mesh2D);
  SubMesh mesh(mesh2D, _TestSubMesh::label);

  mesh.initialize();
} // testInitialize

// ----------------------------------------------------------------------
void
pylith::topology::TestSubMesh::_buildMesh(Mesh* mesh)
{ // _buildMesh
  assert(mesh);

  const int cellDim = _TestSubMesh::cellDim;
  const int ncells = _TestSubMesh::ncells;
  const int* cells = _TestSubMesh::cells;
  const int nvertices = _TestSubMesh::nvertices;
  const int ncorners = _TestSubMesh::ncorners;
  const int spaceDim = _TestSubMesh::cellDim;
  const PylithScalar* coordinates = _TestSubMesh::coordinates;
  const bool interpolate = false;

  mesh->createDMMesh(_TestSubMesh::cellDim);
  PetscDM dmMesh = mesh->dmMesh();CPPUNIT_ASSERT(dmMesh);
  PetscErrorCode err;
  
  err = DMPlexSetChart(dmMesh, 0, ncells+nvertices);CHECK_PETSC_ERROR(err);
  for(PetscInt c = 0; c < ncells; ++c) {
    err = DMPlexSetConeSize(dmMesh, c, ncorners);CHECK_PETSC_ERROR(err);
  } // for
  err = DMSetUp(dmMesh);CHECK_PETSC_ERROR(err);
  PetscInt *cone = new PetscInt[ncorners];
  for(PetscInt c = 0; c < ncells; ++c) {
    for(PetscInt v = 0; v < ncorners; ++v) {
      cone[v] = cells[c*ncorners+v]+ncells;
    } // for
    err = DMPlexSetCone(dmMesh, c, cone);CHECK_PETSC_ERROR(err);
  } // for
  delete[] cone; cone = 0;
  err = DMPlexSymmetrize(dmMesh);CHECK_PETSC_ERROR(err);
  err = DMPlexStratify(dmMesh);CHECK_PETSC_ERROR(err);
  PetscSection coordSection;
  PetscVec coordVec;
  PetscScalar *coords = NULL;
  PetscInt coordSize;

  err = DMPlexGetCoordinateSection(dmMesh, &coordSection);CHECK_PETSC_ERROR(err);
  err = PetscSectionSetChart(coordSection, ncells, ncells+nvertices);CHECK_PETSC_ERROR(err);
  for(PetscInt v = ncells; v < ncells+nvertices; ++v) {
    err = PetscSectionSetDof(coordSection, v, spaceDim);CHECK_PETSC_ERROR(err);
  } // for
  err = PetscSectionSetUp(coordSection);CHECK_PETSC_ERROR(err);
  err = PetscSectionGetStorageSize(coordSection, &coordSize);CHECK_PETSC_ERROR(err);
  err = VecCreate(mesh->comm(), &coordVec);CHECK_PETSC_ERROR(err);
  err = VecSetSizes(coordVec, coordSize, PETSC_DETERMINE);CHECK_PETSC_ERROR(err);
  err = VecSetFromOptions(coordVec);CHECK_PETSC_ERROR(err);
  err = VecGetArray(coordVec, &coords);CHECK_PETSC_ERROR(err);
  for(PetscInt v = 0; v < nvertices; ++v) {
    PetscInt off;

    err = PetscSectionGetOffset(coordSection, v+ncells, &off);CHECK_PETSC_ERROR(err);
    for(PetscInt d = 0; d < spaceDim; ++d) {
      coords[off+d] = coordinates[v*spaceDim+d];
    } // for
  } // for
  err = VecRestoreArray(coordVec, &coords);CHECK_PETSC_ERROR(err);
  err = DMSetCoordinatesLocal(dmMesh, coordVec);CHECK_PETSC_ERROR(err);

  spatialdata::geocoords::CSCart cs;
  cs.setSpaceDim(spaceDim);
  cs.initialize();
  mesh->coordsys(&cs);

  const int numPoints = _TestSubMesh::groupSize;
  for(PetscInt i = 0; i < numPoints; ++i) {
    err = DMPlexSetLabelValue(dmMesh, _TestSubMesh::label, ncells+_TestSubMesh::groupVertices[i], 1);CHECK_PETSC_ERROR(err);
  } // for
} // _buildMesh


// End of file 
