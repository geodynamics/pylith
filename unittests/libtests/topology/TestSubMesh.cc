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

#include "TestSubMesh.hh" // Implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/MeshOps.hh" // USES MeshOps::createDMMesh()

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
      const int cells[ncells*ncorners] = {
	0, 1, 3,
	0, 3, 2,
      };
      const PylithScalar coordinates[nvertices*cellDim] = {
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
      const int submeshNumCorners = 2;
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

  Mesh submesh;
  CPPUNIT_ASSERT_EQUAL(0, submesh.dimension());
  CPPUNIT_ASSERT_EQUAL(false, submesh.debug());

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
  
  Mesh submesh(mesh2D, _TestSubMesh::label);
  CPPUNIT_ASSERT_EQUAL(_TestSubMesh::cellDim-1, submesh.dimension());

  int result = 0;
  MPI_Comm_compare(PETSC_COMM_WORLD, submesh.comm(), &result);
  CPPUNIT_ASSERT_EQUAL(int(MPI_CONGRUENT), result);

  // Check vertices
  const PetscDM dmMesh = submesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
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
// Test coordsys().
void
pylith::topology::TestSubMesh::testCoordsys(void)
{ // testCoordsys
  PYLITH_METHOD_BEGIN;

  Mesh mesh2D;
  _buildMesh(&mesh2D);

  Mesh submesh(mesh2D, _TestSubMesh::label);

  CPPUNIT_ASSERT_EQUAL(_TestSubMesh::cellDim, submesh.coordsys()->spaceDim());

  PYLITH_METHOD_END;
} // testCoordsys

// ----------------------------------------------------------------------
// Test debug().
void
pylith::topology::TestSubMesh::testDebug(void)
{ // testDebug
  PYLITH_METHOD_BEGIN;

  Mesh mesh2D;
  _buildMesh(&mesh2D);

  Mesh submesh(mesh2D, _TestSubMesh::label);
  CPPUNIT_ASSERT_EQUAL(false, submesh.debug());

  submesh.debug(true);
  CPPUNIT_ASSERT_EQUAL(true, submesh.debug());

  PYLITH_METHOD_END;
} // testDebug

// ----------------------------------------------------------------------
// Test dimension().
void
pylith::topology::TestSubMesh::testDimension(void)
{ // testDimension
  PYLITH_METHOD_BEGIN;

  Mesh submesh;
  CPPUNIT_ASSERT_EQUAL(0, submesh.dimension());

  Mesh mesh2D;
  _buildMesh(&mesh2D);
  Mesh submesh2(mesh2D, _TestSubMesh::label);
  CPPUNIT_ASSERT_EQUAL(_TestSubMesh::cellDim-1, submesh2.dimension());

  PYLITH_METHOD_END;
} // testDimension

// ----------------------------------------------------------------------
// Test numCorners().
void
pylith::topology::TestSubMesh::testNumCorners(void)
{ // testNumCorners
  PYLITH_METHOD_BEGIN;

  Mesh mesh;
  _buildMesh(&mesh);

  Mesh submesh(mesh, _TestSubMesh::label);
  CPPUNIT_ASSERT_EQUAL(_TestSubMesh::submeshNumCorners, submesh.numCorners());

  PYLITH_METHOD_END;
} // testNumCorners

// ----------------------------------------------------------------------
// Test numVertices().
void
pylith::topology::TestSubMesh::testNumVertices(void)
{ // testNumVertices
  PYLITH_METHOD_BEGIN;

  Mesh mesh;
  _buildMesh(&mesh);
  Mesh submesh(mesh, _TestSubMesh::label);
  CPPUNIT_ASSERT_EQUAL(_TestSubMesh::submeshNumVertices, submesh.numVertices());

  PYLITH_METHOD_END;
} // testNumVertices

// ----------------------------------------------------------------------
// Test numCells().
void
pylith::topology::TestSubMesh::testNumCells(void)
{ // testNumCells
  PYLITH_METHOD_BEGIN;

  Mesh mesh;
  _buildMesh(&mesh);
  Mesh submesh(mesh, _TestSubMesh::label);
  CPPUNIT_ASSERT_EQUAL(_TestSubMesh::submeshNumCells, submesh.numCells());

  PYLITH_METHOD_END;
} // testNumCells

// ----------------------------------------------------------------------
// Test comm().
void
pylith::topology::TestSubMesh::testComm(void)
{ // testComm
  PYLITH_METHOD_BEGIN;

  Mesh mesh2D;
  _buildMesh(&mesh2D);

  Mesh submesh(mesh2D, _TestSubMesh::label);

  int result = 0;
  MPI_Comm_compare(PETSC_COMM_WORLD, submesh.comm(), &result);
  CPPUNIT_ASSERT_EQUAL(int(MPI_CONGRUENT), result);

  PYLITH_METHOD_END;
} // testComm

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

  MeshOps::createDMMesh(mesh, _TestSubMesh::cellDim);
  PetscDM dmMesh = mesh->dmMesh();CPPUNIT_ASSERT(dmMesh);
  PetscErrorCode err;
  
  err = DMPlexSetChart(dmMesh, 0, ncells+nvertices);PYLITH_CHECK_ERROR(err);
  for(PetscInt c = 0; c < ncells; ++c) {
    err = DMPlexSetConeSize(dmMesh, c, ncorners);PYLITH_CHECK_ERROR(err);
  } // for
  err = DMSetUp(dmMesh);PYLITH_CHECK_ERROR(err);
  PetscInt *cone = new PetscInt[ncorners];
  for(PetscInt c = 0; c < ncells; ++c) {
    for(PetscInt v = 0; v < ncorners; ++v) {
      cone[v] = cells[c*ncorners+v]+ncells;
    } // for
    err = DMPlexSetCone(dmMesh, c, cone);PYLITH_CHECK_ERROR(err);
  } // for
  delete[] cone; cone = 0;
  err = DMPlexSymmetrize(dmMesh);PYLITH_CHECK_ERROR(err);
  err = DMPlexStratify(dmMesh);PYLITH_CHECK_ERROR(err);
  PetscSection coordSection = NULL;
  PetscVec coordVec = NULL;
  PetscScalar *coords = NULL;
  PetscInt coordSize = 0;

  err = DMGetCoordinateSection(dmMesh, &coordSection);PYLITH_CHECK_ERROR(err);
  err = PetscSectionSetNumFields(coordSection, 1);PYLITH_CHECK_ERROR(err);
  err = PetscSectionSetFieldComponents(coordSection, 0, spaceDim);PYLITH_CHECK_ERROR(err);
  err = PetscSectionSetChart(coordSection, ncells, ncells+nvertices);PYLITH_CHECK_ERROR(err);
  for(PetscInt v = ncells; v < ncells+nvertices; ++v) {
    err = PetscSectionSetDof(coordSection, v, spaceDim);PYLITH_CHECK_ERROR(err);
  } // for
  err = PetscSectionSetUp(coordSection);PYLITH_CHECK_ERROR(err);
  err = PetscSectionGetStorageSize(coordSection, &coordSize);PYLITH_CHECK_ERROR(err);
  err = VecCreate(mesh->comm(), &coordVec);PYLITH_CHECK_ERROR(err);
  err = VecSetSizes(coordVec, coordSize, PETSC_DETERMINE);PYLITH_CHECK_ERROR(err);
  err = VecSetFromOptions(coordVec);PYLITH_CHECK_ERROR(err);
  err = VecGetArray(coordVec, &coords);PYLITH_CHECK_ERROR(err);
  for(PetscInt v = 0; v < nvertices; ++v) {
    PetscInt off;

    err = PetscSectionGetOffset(coordSection, v+ncells, &off);PYLITH_CHECK_ERROR(err);
    for(PetscInt d = 0; d < spaceDim; ++d) {
      coords[off+d] = coordinates[v*spaceDim+d];
    } // for
  } // for
  err = VecRestoreArray(coordVec, &coords);PYLITH_CHECK_ERROR(err);
  err = DMSetCoordinatesLocal(dmMesh, coordVec);PYLITH_CHECK_ERROR(err);
  err = VecDestroy(&coordVec);PYLITH_CHECK_ERROR(err);

  spatialdata::geocoords::CSCart cs;
  cs.setSpaceDim(spaceDim);
  cs.initialize();
  mesh->coordsys(&cs);

  const int numPoints = _TestSubMesh::groupSize;
  for(PetscInt i = 0; i < numPoints; ++i) {
    err = DMSetLabelValue(dmMesh, _TestSubMesh::label, ncells+_TestSubMesh::groupVertices[i], 1);PYLITH_CHECK_ERROR(err);
  } // for

  PYLITH_METHOD_END;
} // _buildMesh


// End of file 
