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
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestSubMesh.hh"   // Implementation of class methods

#include "pylith/topology/Mesh.hh"  // USES Mesh
#include "pylith/topology/Stratum.hh"   // USES Stratum
#include "pylith/meshio/MeshBuilder.hh" // USES MeshBuilder::buildMesh()

#include "spatialdata/geocoords/CSCart.hh"  // USES CSCart

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::topology::TestSubMesh::setUp(void)
{ // setUp
    PYLITH_METHOD_BEGIN;

    _data = new TestSubMesh_Data; CPPUNIT_ASSERT(_data);

    PYLITH_METHOD_END;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::topology::TestSubMesh::tearDown(void)
{ // tearDown
    PYLITH_METHOD_BEGIN;

    delete _data; _data = NULL;

    PYLITH_METHOD_END;
} // tearDown

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
{   // testConstructorMesh
    PYLITH_METHOD_BEGIN;
    CPPUNIT_ASSERT(_data);

    Mesh mesh2D;
    _buildMesh(&mesh2D);

    Mesh submesh(mesh2D, _data->label);
    CPPUNIT_ASSERT_EQUAL(_data->cellDim-1, submesh.dimension());

    int result = 0;
    MPI_Comm_compare(PETSC_COMM_WORLD, submesh.comm(), &result);
    CPPUNIT_ASSERT_EQUAL(int(MPI_CONGRUENT), result);

    // Check vertices
    const PetscDM dmMesh = submesh.dmMesh(); CPPUNIT_ASSERT(dmMesh);
    Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
    const PetscInt vStart = depthStratum.begin();
    const PetscInt vEnd = depthStratum.end();

    const PetscInt nvertices = _data->submeshNumVertices;
    CPPUNIT_ASSERT_EQUAL(nvertices, depthStratum.size());
    for (PetscInt v = vStart, iV=0; v < vEnd; ++v, ++iV) {
        CPPUNIT_ASSERT_EQUAL(_data->submeshVertices[iV], v);
    }   // for

    // Check cells
    Stratum heightStratum(dmMesh, Stratum::HEIGHT, 0);
    const PetscInt cStart = heightStratum.begin();
    const PetscInt cEnd = heightStratum.end();

    const PetscInt ncells = _data->submeshNumCells;
    CPPUNIT_ASSERT_EQUAL(ncells, heightStratum.size());
    for (PetscInt c = cStart, iC=0; c < cEnd; ++c, ++iC) {
        CPPUNIT_ASSERT_EQUAL(_data->submeshCells[iC], c);
    }   // for

    CPPUNIT_ASSERT_THROW(Mesh submesh(mesh2D, "zzyyxx"), std::runtime_error);

    PYLITH_METHOD_END;
}   // testConstructorMesh

// ----------------------------------------------------------------------
// Test coordsys().
void
pylith::topology::TestSubMesh::testCoordsys(void)
{   // testCoordsys
    PYLITH_METHOD_BEGIN;
    CPPUNIT_ASSERT(_data);

    Mesh mesh2D;
    _buildMesh(&mesh2D);

    Mesh submesh(mesh2D, _data->label);

    CPPUNIT_ASSERT_EQUAL(_data->cellDim, submesh.coordsys()->spaceDim());

    PYLITH_METHOD_END;
}   // testCoordsys

// ----------------------------------------------------------------------
// Test debug().
void
pylith::topology::TestSubMesh::testDebug(void)
{   // testDebug
    PYLITH_METHOD_BEGIN;
    CPPUNIT_ASSERT(_data);

    Mesh mesh2D;
    _buildMesh(&mesh2D);

    Mesh submesh(mesh2D, _data->label);
    CPPUNIT_ASSERT_EQUAL(false, submesh.debug());

    submesh.debug(true);
    CPPUNIT_ASSERT_EQUAL(true, submesh.debug());

    PYLITH_METHOD_END;
}   // testDebug

// ----------------------------------------------------------------------
// Test dimension().
void
pylith::topology::TestSubMesh::testDimension(void)
{   // testDimension
    PYLITH_METHOD_BEGIN;
    CPPUNIT_ASSERT(_data);

    Mesh submesh;
    CPPUNIT_ASSERT_EQUAL(0, submesh.dimension());

    Mesh mesh2D;
    _buildMesh(&mesh2D);
    Mesh submesh2(mesh2D, _data->label);
    CPPUNIT_ASSERT_EQUAL(_data->cellDim-1, submesh2.dimension());

    PYLITH_METHOD_END;
}   // testDimension

// ----------------------------------------------------------------------
// Test numCorners().
void
pylith::topology::TestSubMesh::testNumCorners(void)
{   // testNumCorners
    PYLITH_METHOD_BEGIN;
    CPPUNIT_ASSERT(_data);

    Mesh mesh;
    _buildMesh(&mesh);

    Mesh submesh(mesh, _data->label);
    CPPUNIT_ASSERT_EQUAL(_data->submeshNumCorners, submesh.numCorners());

    PYLITH_METHOD_END;
}   // testNumCorners

// ----------------------------------------------------------------------
// Test numVertices().
void
pylith::topology::TestSubMesh::testNumVertices(void)
{   // testNumVertices
    PYLITH_METHOD_BEGIN;
    CPPUNIT_ASSERT(_data);

    Mesh mesh;
    _buildMesh(&mesh);
    Mesh submesh(mesh, _data->label);
    CPPUNIT_ASSERT_EQUAL(_data->submeshNumVertices, submesh.numVertices());

    PYLITH_METHOD_END;
}   // testNumVertices

// ----------------------------------------------------------------------
// Test numCells().
void
pylith::topology::TestSubMesh::testNumCells(void)
{   // testNumCells
    PYLITH_METHOD_BEGIN;
    CPPUNIT_ASSERT(_data);

    Mesh mesh;
    _buildMesh(&mesh);
    Mesh submesh(mesh, _data->label);
    CPPUNIT_ASSERT_EQUAL(_data->submeshNumCells, submesh.numCells());

    PYLITH_METHOD_END;
}   // testNumCells

// ----------------------------------------------------------------------
// Test comm().
void
pylith::topology::TestSubMesh::testComm(void)
{   // testComm
    PYLITH_METHOD_BEGIN;
    CPPUNIT_ASSERT(_data);

    Mesh mesh2D;
    _buildMesh(&mesh2D);

    Mesh submesh(mesh2D, _data->label);

    int result = 0;
    MPI_Comm_compare(PETSC_COMM_WORLD, submesh.comm(), &result);
    CPPUNIT_ASSERT_EQUAL(int(MPI_CONGRUENT), result);

    PYLITH_METHOD_END;
}   // testComm

// ----------------------------------------------------------------------
void
pylith::topology::TestSubMesh::_buildMesh(Mesh* mesh)
{   // _buildMesh
    PYLITH_METHOD_BEGIN;
    CPPUNIT_ASSERT(_data);
    CPPUNIT_ASSERT(mesh);

    const int cellDim = _data->cellDim;
    const int numCells = _data->numCells;
    const int numVertices = _data->numVertices;
    const int numCorners = _data->numCorners;
    const int spaceDim = _data->cellDim;

    PetscInt size = numVertices * spaceDim;
    scalar_array coordinates(size);
    for (PetscInt i=0; i < size; ++i) {
        coordinates[i] = _data->coordinates[i];
    }   // for

    size = numCells * numCorners;
    int_array cells(size);
    for (PetscInt i=0; i < size; ++i) {
        cells[i] = _data->cells[i];
    }   // for

    pylith::meshio::MeshBuilder::buildMesh(mesh, &coordinates, numVertices, spaceDim, cells, numCells, numCorners,
                                           cellDim);

    spatialdata::geocoords::CSCart cs;
    cs.setSpaceDim(spaceDim);
    cs.initialize();
    mesh->coordsys(&cs);

    PetscErrorCode err;
    const int numPoints = _data->groupSize;
    for(PetscInt i = 0; i < numPoints; ++i) {
        err = DMSetLabelValue(mesh->dmMesh(), _data->label, numCells+_data->groupVertices[i], 1); CPPUNIT_ASSERT(!err);
    }   // for

    PYLITH_METHOD_END;
}   // _buildMesh


// ----------------------------------------------------------------------
// Constructor
pylith::topology::TestSubMesh_Data::TestSubMesh_Data(void) :
    cellDim(0),
    numVertices(0),
    numCells(0),
    numCorners(0),
    cells(NULL),
    coordinates(NULL),
    label(NULL),
    groupSize(0),
    groupVertices(NULL),
    submeshNumCorners(0),
    submeshNumVertices(0),
    submeshVertices(NULL),
    submeshNumCells(0),
    submeshCells(NULL)
{   // constructor
}   // constructor


// ----------------------------------------------------------------------
// Destructor
pylith::topology::TestSubMesh_Data::~TestSubMesh_Data(void)
{   // destructor
}   // destructor


// End of file
