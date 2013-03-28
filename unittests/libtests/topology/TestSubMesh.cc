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
  PYLITH_METHOD_BEGIN;

  SubMesh mesh;
  CPPUNIT_ASSERT_EQUAL(0, mesh.dimension());
  CPPUNIT_ASSERT_EQUAL(false, mesh.debug());

  PYLITH_METHOD_END;
} // testConstructor

// ----------------------------------------------------------------------
// Test constructor w/mesh.
void
pylith::topology::TestSubMesh::testConstructorMesh(void)
{ // testConstructorMesh
  PYLITH_METHOD_BEGIN;

 Mesh mesh2D;
  _buildMesh(&mesh2D);
  
  SubMesh mesh(mesh2D, _TestSubMesh::label);
  CPPUNIT_ASSERT_EQUAL(_TestSubMesh::cellDim-1, mesh.dimension());

  int result = 0;
  MPI_Comm_compare(PETSC_COMM_WORLD, mesh.comm(), &result);
  CPPUNIT_ASSERT_EQUAL(int(MPI_CONGRUENT), result);

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

  PYLITH_METHOD_END;
} // testConstructorMesh

// ----------------------------------------------------------------------
// Test createSubMesh().
void
pylith::topology::TestSubMesh::testCreateSubMesh(void)
{ // testCreateSubMesh
  PYLITH_METHOD_BEGIN;

  Mesh mesh2D;
  _buildMesh(&mesh2D);
  
  SubMesh mesh;
  mesh.createSubMesh(mesh2D, _TestSubMesh::label);
  CPPUNIT_ASSERT_EQUAL(_TestSubMesh::cellDim-1, mesh.dimension());

  int result = 0;
  MPI_Comm_compare(PETSC_COMM_WORLD, mesh.comm(), &result);
  CPPUNIT_ASSERT_EQUAL(int(MPI_CONGRUENT), result);

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

  PYLITH_METHOD_END;
} // testCreateSubMesh

// ----------------------------------------------------------------------
// Test coordsys().
void
pylith::topology::TestSubMesh::testCoordsys(void)
{ // testCoordsys
  PYLITH_METHOD_BEGIN;

  Mesh mesh2D;
  _buildMesh(&mesh2D);

  SubMesh mesh(mesh2D, _TestSubMesh::label);

  CPPUNIT_ASSERT_EQUAL(_TestSubMesh::cellDim, mesh.coordsys()->spaceDim());

  PYLITH_METHOD_END;
} // testCoordsys

// ----------------------------------------------------------------------
// Test debug().
void
pylith::topology::TestSubMesh::testDebug(void)
{ // testDebug
  PYLITH_METHOD_BEGIN;

  SubMesh mesh;
  CPPUNIT_ASSERT_EQUAL(false, mesh.debug());

  mesh.debug(true);
  CPPUNIT_ASSERT_EQUAL(true, mesh.debug());

  PYLITH_METHOD_END;
} // testDebug

// ----------------------------------------------------------------------
// Test dimension().
void
pylith::topology::TestSubMesh::testDimension(void)
{ // testDimension
  PYLITH_METHOD_BEGIN;

  SubMesh mesh;
  CPPUNIT_ASSERT_EQUAL(0, mesh.dimension());

  Mesh mesh2D;
  _buildMesh(&mesh2D);
  SubMesh mesh2(mesh2D, _TestSubMesh::label);
  CPPUNIT_ASSERT_EQUAL(_TestSubMesh::cellDim-1, mesh2.dimension());

  PYLITH_METHOD_END;
} // testDimension

// ----------------------------------------------------------------------
// Test comm().
void
pylith::topology::TestSubMesh::testComm(void)
{ // testComm
  PYLITH_METHOD_BEGIN;

  SubMesh mesh;

  Mesh mesh2D;
  _buildMesh(&mesh2D);
  mesh.createSubMesh(mesh2D, _TestSubMesh::label);

  int result = 0;
  MPI_Comm_compare(PETSC_COMM_WORLD, mesh.comm(), &result);
  CPPUNIT_ASSERT_EQUAL(int(MPI_CONGRUENT), result);

  PYLITH_METHOD_END;
} // testComm

// ----------------------------------------------------------------------
// Test initialize().
void
pylith::topology::TestSubMesh::testInitialize(void)
{ // testInitialize
  PYLITH_METHOD_BEGIN;

  Mesh mesh2D;
  _buildMesh(&mesh2D);
  SubMesh mesh(mesh2D, _TestSubMesh::label);

  mesh.initialize();

  PYLITH_METHOD_END;
} // testInitialize

// ----------------------------------------------------------------------
void
pylith::topology::TestSubMesh::_buildMesh(Mesh* mesh)
{ // _buildMesh
  PYLITH_METHOD_BEGIN;

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

  PYLITH_METHOD_END;
} // _buildMesh


// End of file 
