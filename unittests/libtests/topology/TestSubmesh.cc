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
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestSubmesh.hh" // Implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES createLowerDimMesh()
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/meshio/MeshBuilder.hh" // USES MeshBuilder::buildMesh()

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::topology::TestSubmesh::setUp(void) {
    PYLITH_METHOD_BEGIN;

    _mesh = NULL;
    _submesh = NULL;
    _data = new TestSubmesh_Data;CPPUNIT_ASSERT(_data);

    PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::topology::TestSubmesh::tearDown(void) {
    PYLITH_METHOD_BEGIN;

    delete _mesh;_mesh = NULL;
    delete _submesh;_submesh = NULL;
    delete _data;_data = NULL;

    PYLITH_METHOD_END;
} // tearDown


// ----------------------------------------------------------------------
// Test constructor.
void
pylith::topology::TestSubmesh::testConstructor(void) {
    PYLITH_METHOD_BEGIN;

    Mesh submesh;
    CPPUNIT_ASSERT_EQUAL(0, submesh.dimension());
    CPPUNIT_ASSERT_EQUAL(false, submesh.debug());

    PYLITH_METHOD_END;
} // testConstructor


// ----------------------------------------------------------------------
// Test constructor w/mesh.
void
pylith::topology::TestSubmesh::testConstructorMesh(void) {
    PYLITH_METHOD_BEGIN;
    CPPUNIT_ASSERT(_data);

    _buildMesh();

    CPPUNIT_ASSERT_EQUAL(_data->cellDim-1, _submesh->dimension());

    int result = 0;
    MPI_Comm_compare(PETSC_COMM_WORLD, _submesh->comm(), &result);
    CPPUNIT_ASSERT_EQUAL(int(MPI_CONGRUENT), result);

    // Check vertices
    const PetscDM dmMesh = _submesh->dmMesh();CPPUNIT_ASSERT(dmMesh);
    Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
    const PetscInt vStart = depthStratum.begin();
    const PetscInt vEnd = depthStratum.end();

    const PetscInt nvertices = _data->submeshNumVertices;
    CPPUNIT_ASSERT_EQUAL(nvertices, depthStratum.size());
    for (PetscInt v = vStart, iV = 0; v < vEnd; ++v, ++iV) {
        CPPUNIT_ASSERT_EQUAL(_data->submeshVertices[iV], v);
    } // for

    // Check cells
    Stratum heightStratum(dmMesh, Stratum::HEIGHT, 0);
    const PetscInt cStart = heightStratum.begin();
    const PetscInt cEnd = heightStratum.end();

    const PetscInt ncells = _data->submeshNumCells;
    CPPUNIT_ASSERT_EQUAL(ncells, heightStratum.size());
    for (PetscInt c = cStart, iC = 0; c < cEnd; ++c, ++iC) {
        CPPUNIT_ASSERT_EQUAL(_data->submeshCells[iC], c);
    } // for

    CPPUNIT_ASSERT_THROW(MeshOps::createLowerDimMesh(*_mesh, "zzyyxx"), std::runtime_error);

    PYLITH_METHOD_END;
} // testConstructorMesh


// ----------------------------------------------------------------------
// Test coordsys(), debug(), comm().
void
pylith::topology::TestSubmesh::testAccessors(void) {
    PYLITH_METHOD_BEGIN;
    CPPUNIT_ASSERT(_data);

    _buildMesh();

    CPPUNIT_ASSERT_EQUAL(_data->cellDim, _submesh->coordsys()->spaceDim());

    CPPUNIT_ASSERT_EQUAL(false, _submesh->debug());
    _submesh->debug(true);
    CPPUNIT_ASSERT_EQUAL(true, _submesh->debug());

    int result = 0;
    MPI_Comm_compare(PETSC_COMM_WORLD, _submesh->comm(), &result);
    CPPUNIT_ASSERT_EQUAL(int(MPI_CONGRUENT), result);

    PYLITH_METHOD_END;
} // testAccessors


// ----------------------------------------------------------------------
// Test dimension(), numCorners(), numVertices(), numCells().
void
pylith::topology::TestSubmesh::testSizes(void) {
    PYLITH_METHOD_BEGIN;
    CPPUNIT_ASSERT(_data);

    Mesh submesh;
    CPPUNIT_ASSERT_EQUAL(0, submesh.dimension());

    _buildMesh();
    CPPUNIT_ASSERT_EQUAL(_data->cellDim-1, _submesh->dimension());
    CPPUNIT_ASSERT_EQUAL(_data->submeshNumCorners, _submesh->numCorners());
    CPPUNIT_ASSERT_EQUAL(_data->submeshNumVertices, _submesh->numVertices());
    CPPUNIT_ASSERT_EQUAL(_data->submeshNumCells, _submesh->numCells());

    PYLITH_METHOD_END;
} // testSizes


// ----------------------------------------------------------------------
void
pylith::topology::TestSubmesh::_buildMesh(void) {
    PYLITH_METHOD_BEGIN;
    CPPUNIT_ASSERT(_data);
    CPPUNIT_ASSERT(!_mesh);
    CPPUNIT_ASSERT(!_submesh);

    const int cellDim = _data->cellDim;
    const int numCells = _data->numCells;
    const int numVertices = _data->numVertices;
    const int numCorners = _data->numCorners;
    const int spaceDim = _data->cellDim;

    PetscInt size = numVertices * spaceDim;
    scalar_array coordinates(size);
    for (PetscInt i = 0; i < size; ++i) {
        coordinates[i] = _data->coordinates[i];
    } // for

    size = numCells * numCorners;
    int_array cells(size);
    for (PetscInt i = 0; i < size; ++i) {
        cells[i] = _data->cells[i];
    } // for

    delete _mesh;_mesh = new pylith::topology::Mesh();CPPUNIT_ASSERT(_mesh);
    pylith::meshio::MeshBuilder::buildMesh(_mesh, &coordinates, numVertices, spaceDim, cells, numCells, numCorners,
                                           cellDim);

    spatialdata::geocoords::CSCart cs;
    cs.setSpaceDim(spaceDim);
    cs.initialize();
    _mesh->coordsys(&cs);

    PetscErrorCode err;
    const int numPoints = _data->groupSize;
    for (PetscInt i = 0; i < numPoints; ++i) {
        err = DMSetLabelValue(_mesh->dmMesh(), _data->label, numCells+_data->groupVertices[i], 1);CPPUNIT_ASSERT(!err);
    } // for

    delete _submesh;_submesh = MeshOps::createLowerDimMesh(*_mesh, _data->label);CPPUNIT_ASSERT(_submesh);

    PYLITH_METHOD_END;
} // _buildMesh


// ----------------------------------------------------------------------
// Constructor
pylith::topology::TestSubmesh_Data::TestSubmesh_Data(void) :
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
    submeshCells(NULL) { // constructor
} // constructor


// ----------------------------------------------------------------------
// Destructor
pylith::topology::TestSubmesh_Data::~TestSubmesh_Data(void) {}


// End of file
