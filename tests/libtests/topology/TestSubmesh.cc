// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------

#include <portinfo>

#include "TestSubmesh.hh" // Implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES createLowerDimMesh()
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/meshio/MeshBuilder.hh" // USES MeshBuilder::buildMesh()

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart

// ---------------------------------------------------------------------------------------------------------------------
// Setup testing data.
void
pylith::topology::TestSubmesh::setUp(void) {
    PYLITH_METHOD_BEGIN;

    _domainMesh = NULL;
    _testMesh = NULL;
    _data = new TestSubmesh_Data;CPPUNIT_ASSERT(_data);

    PYLITH_METHOD_END;
} // setUp


// ---------------------------------------------------------------------------------------------------------------------
// Tear down testing data.
void
pylith::topology::TestSubmesh::tearDown(void) {
    PYLITH_METHOD_BEGIN;

    delete _domainMesh;_domainMesh = NULL;
    delete _testMesh;_testMesh = NULL;
    delete _data;_data = NULL;

    PYLITH_METHOD_END;
} // tearDown


// ---------------------------------------------------------------------------------------------------------------------
// Test constructor w/mesh.
void
pylith::topology::TestSubmesh::testCreateLowerDimMesh(void) {
    PYLITH_METHOD_BEGIN;
    CPPUNIT_ASSERT(_data);

    _buildMesh();
    delete _testMesh;_testMesh = MeshOps::createLowerDimMesh(*_domainMesh, _data->groupLabel);CPPUNIT_ASSERT(_testMesh);

    CPPUNIT_ASSERT_EQUAL(_data->cellDim-1, _testMesh->getDimension());

    int result = 0;
    MPI_Comm_compare(PETSC_COMM_WORLD, _testMesh->getComm(), &result);
    CPPUNIT_ASSERT_EQUAL(int(MPI_CONGRUENT), result);

    // Check vertices
    const PetscDM dmMesh = _testMesh->getDM();CPPUNIT_ASSERT(dmMesh);
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

    CPPUNIT_ASSERT_THROW(MeshOps::createLowerDimMesh(*_domainMesh, "zzyyxx"), std::runtime_error);

    PYLITH_METHOD_END;
} // testCreateLowerDimMesh


// ---------------------------------------------------------------------------------------------------------------------
// Test constructor w/mesh.
void
pylith::topology::TestSubmesh::testCreateSubdomainMesh(void) {
    PYLITH_METHOD_BEGIN;
    CPPUNIT_ASSERT(_data);

    _buildMesh();
    delete _testMesh;_testMesh = MeshOps::createSubdomainMesh(*_domainMesh, _data->subdomainLabel,
                                                              _data->subdomainLabelValue, "Test subdomain");
    CPPUNIT_ASSERT(_testMesh);

    CPPUNIT_ASSERT_EQUAL(_data->cellDim, _testMesh->getDimension());

    int result = 0;
    MPI_Comm_compare(PETSC_COMM_WORLD, _testMesh->getComm(), &result);
    CPPUNIT_ASSERT_EQUAL(int(MPI_CONGRUENT), result);

    // Check vertices
    const PetscDM dmMesh = _testMesh->getDM();CPPUNIT_ASSERT(dmMesh);
    Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
    const PetscInt vStart = depthStratum.begin();
    const PetscInt vEnd = depthStratum.end();

    const PetscInt nvertices = _data->subdomainNumVertices;
    CPPUNIT_ASSERT_EQUAL(nvertices, depthStratum.size());
    for (PetscInt v = vStart, iV = 0; v < vEnd; ++v, ++iV) {
        CPPUNIT_ASSERT_EQUAL(_data->subdomainVertices[iV], v);
    } // for

    // Check cells
    Stratum heightStratum(dmMesh, Stratum::HEIGHT, 0);
    const PetscInt cStart = heightStratum.begin();
    const PetscInt cEnd = heightStratum.end();

    const PetscInt ncells = _data->subdomainNumCells;
    CPPUNIT_ASSERT_EQUAL(ncells, heightStratum.size());
    for (PetscInt c = cStart, iC = 0; c < cEnd; ++c, ++iC) {
        CPPUNIT_ASSERT_EQUAL(_data->subdomainCells[iC], c);
    } // for

    CPPUNIT_ASSERT_THROW(MeshOps::createSubdomainMesh(*_domainMesh, "material-id", -9, "Test subdomain"), std::runtime_error);
    CPPUNIT_ASSERT_THROW(MeshOps::createSubdomainMesh(*_domainMesh, "zzyyxx", -9, "Test subdomain"), std::runtime_error);

    PYLITH_METHOD_END;
} // testCreateLowerDimMesh


// ---------------------------------------------------------------------------------------------------------------------
// Test getCoordSys(), debug(), comm().
void
pylith::topology::TestSubmesh::testAccessors(void) {
    PYLITH_METHOD_BEGIN;
    CPPUNIT_ASSERT(_data);

    _buildMesh();
    delete _testMesh;_testMesh = MeshOps::createLowerDimMesh(*_domainMesh, _data->groupLabel);CPPUNIT_ASSERT(_testMesh);
    CPPUNIT_ASSERT(_testMesh->getCoordSys());

    CPPUNIT_ASSERT_EQUAL(size_t(_data->cellDim), _testMesh->getCoordSys()->getSpaceDim());

    int result = 0;
    MPI_Comm_compare(PETSC_COMM_WORLD, _testMesh->getComm(), &result);
    CPPUNIT_ASSERT_EQUAL(int(MPI_CONGRUENT), result);

    PYLITH_METHOD_END;
} // testAccessors


// ---------------------------------------------------------------------------------------------------------------------
// Test dimension(), numCorners(), numVertices(), numCells().
void
pylith::topology::TestSubmesh::testSizes(void) {
    PYLITH_METHOD_BEGIN;
    CPPUNIT_ASSERT(_data);

    Mesh submesh;
    CPPUNIT_ASSERT_EQUAL(0, submesh.getDimension());

    _buildMesh();
    delete _testMesh;_testMesh = MeshOps::createLowerDimMesh(*_domainMesh, _data->groupLabel);CPPUNIT_ASSERT(_testMesh);
    CPPUNIT_ASSERT(_testMesh);

    CPPUNIT_ASSERT_EQUAL(_data->cellDim-1, _testMesh->getDimension());

    PYLITH_METHOD_END;
} // testSizes


// ---------------------------------------------------------------------------------------------------------------------
void
pylith::topology::TestSubmesh::_buildMesh(void) {
    PYLITH_METHOD_BEGIN;
    CPPUNIT_ASSERT(_data);
    CPPUNIT_ASSERT(!_domainMesh);
    CPPUNIT_ASSERT(!_testMesh);

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

    delete _domainMesh;_domainMesh = new pylith::topology::Mesh();CPPUNIT_ASSERT(_domainMesh);
    pylith::meshio::MeshBuilder::buildMesh(_domainMesh, &coordinates, numVertices, spaceDim, cells, numCells, numCorners,
                                           cellDim);

    spatialdata::geocoords::CSCart cs;
    cs.setSpaceDim(spaceDim);
    _domainMesh->setCoordSys(&cs);

    PetscErrorCode err;
    const int numPoints = _data->groupSize;
    for (PetscInt i = 0; i < numPoints; ++i) {
        const PylithInt groupLabelValue = 1;
        err = DMSetLabelValue(_domainMesh->getDM(), _data->groupLabel, numCells+_data->groupVertices[i],
                              groupLabelValue);CPPUNIT_ASSERT(!err);
    } // for

    for (PetscInt c = 0; c < numCells; ++c) {
        err = DMSetLabelValue(_domainMesh->getDM(), _data->subdomainLabel, c,
                              _data->subdomainLabelValues[c]);PYLITH_CHECK_ERROR(err);
    } // for

    PYLITH_METHOD_END;
} // _buildMesh


// ---------------------------------------------------------------------------------------------------------------------
// Constructor
pylith::topology::TestSubmesh_Data::TestSubmesh_Data(void) :
    cellDim(0),
    numVertices(0),
    numCells(0),
    numCorners(0),
    cells(NULL),
    coordinates(NULL),
    groupLabel(NULL),
    groupSize(0),
    groupVertices(NULL),
    submeshNumCorners(0),
    submeshNumVertices(0),
    submeshVertices(NULL),
    submeshNumCells(0),
    submeshCells(NULL),
    subdomainLabel(NULL),
    subdomainLabelValues(NULL),
    subdomainLabelValue(0),
    subdomainNumCorners(0),
    subdomainNumVertices(0),
    subdomainVertices(NULL),
    subdomainNumCells(0),
    subdomainCells(NULL) {}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor
pylith::topology::TestSubmesh_Data::~TestSubmesh_Data(void) {}


// End of file
