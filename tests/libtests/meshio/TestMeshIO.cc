// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "TestMeshIO.hh" // Implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps::nondimensionalize()
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor
#include "pylith/meshio/MeshIO.hh" // USES MeshIO

#include "pylith/utils/array.hh" // USES int_array
#include "pylith/utils/journals.hh" // USES journal::debug_t

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers_floating_point.hpp"

#include <strings.h> // USES strcasecmp()
#include <stdexcept> // USES std::logic_error
#include <cassert> // USES assert()

static PetscErrorCode
DMPlexInvertCell_Private(PetscInt dim,
                         PetscInt numCorners,
                         PetscInt cone[]) {
#define SWAPCONE(cone,i,j)  \
    do {                      \
        int _cone_tmp;          \
        _cone_tmp = (cone)[i];  \
        (cone)[i] = (cone)[j];  \
        (cone)[j] = _cone_tmp;  \
    } while (0)

    PetscFunctionBegin;
    if (dim != 3) { PetscFunctionReturn(0);}
    switch (numCorners) {
    case 4: SWAPCONE(cone,0,1);break;
    case 6: SWAPCONE(cone,0,1);break;
    case 8: SWAPCONE(cone,1,3);break;
    default: break;
    }
    PetscFunctionReturn(0);
#undef SWAPCONE
}


// ----------------------------------------------------------------------
// Constructor.
pylith::meshio::TestMeshIO::TestMeshIO(TestMeshIO_Data* data) :
    _data(data) {
    assert(_data);
    _mesh = NULL;
} // constructor


// ----------------------------------------------------------------------
// Tear down testing data.
pylith::meshio::TestMeshIO::~TestMeshIO(void) {
    PYLITH_METHOD_BEGIN;

    delete _mesh;_mesh = NULL;

    PYLITH_METHOD_END;
} // destructor


// ----------------------------------------------------------------------
// Get simple mesh for testing I/O.
void
pylith::meshio::TestMeshIO::_createMesh(void) {
    PYLITH_METHOD_BEGIN;
    assert(_data);

    // buildTopology() requires zero based index
    assert(_data->useIndexZero);

    assert(_data->topology);
    assert(_data->geometry);
    assert(_data->materialIds);
    if (_data->numVertexGroups > 0) {
        assert(_data->vertexGroups);
        assert(_data->vertexGroupSizes);
        assert(_data->vertexGroupNames);
    } // if
    if (_data->numFaceGroups > 0) {
        assert(_data->faceGroups);
        assert(_data->faceGroupSizes);
        assert(_data->faceGroupNames);
    } // if

    delete _mesh;_mesh = new topology::Mesh();assert(_mesh);

    int_array cellsCopy(_data->topology->cells); // Create copy because building mesh may change cells (invert)
    pylith::meshio::MeshBuilder::buildMesh(_mesh, *_data->topology, *_data->geometry);
    _data->topology->cells = cellsCopy;
    const size_t numCells = _data->topology->numCells;

    { // material ids
        int_array materialIds(numCells);
        for (size_t i = 0; i < numCells; ++i) {
            materialIds[i] = _data->materialIds[i];
        } // for
        pylith::meshio::MeshBuilder::setMaterials(_mesh, materialIds);
    } // material ids

    // Vertex groups
    for (size_t iGroup = 0, index = 0; iGroup < _data->numVertexGroups; ++iGroup) {
        const size_t groupSize = _data->vertexGroupSizes[iGroup];assert(groupSize > 0);
        int_array points(groupSize);
        for (size_t i = 0; i < groupSize; ++i, ++index) {
            points[i] = _data->vertexGroups[index];
        } // for
        pylith::meshio::MeshBuilder::setVertexGroup(_mesh, _data->vertexGroupNames[iGroup], points);
    } // for

    // Face groups
    for (size_t iGroup = 0, index = 0; iGroup < _data->numFaceGroups; ++iGroup) {
        const size_t numFaceVertices = _data->numFaceVertices;assert(numFaceVertices > 0);
        const size_t numFaces = _data->faceGroupSizes[iGroup];assert(numFaces > 0);
        const size_t totalSize = numFaces * (1 + numFaceVertices); // cell + vertices
        int_array faceValues(totalSize);
        for (size_t i = 0; i < totalSize; ++i, ++index) {
            faceValues[i] = _data->faceGroups[index];
        } // for
        pylith::meshio::MeshBuilder::setFaceGroupFromCellVertices(_mesh, _data->faceGroupNames[iGroup], faceValues, numFaceVertices);
    } // for

    // Set coordinate system
    spatialdata::geocoords::CSCart cs;
    cs.setSpaceDim(_data->geometry->spaceDim);
    _mesh->setCoordSys(&cs);

    spatialdata::units::Nondimensional normalizer;
    normalizer.setLengthScale(10.0);
    topology::MeshOps::nondimensionalize(_mesh, normalizer);

    pythia::journal::debug_t debug("TestMeshIO");
    if (debug.state()) {
        _mesh->view();
        _mesh->view(":mesh.tex:ascii_latex");
    } // if

    PYLITH_METHOD_END;
} // _createMesh


// ----------------------------------------------------------------------
// Check values in mesh against data.
void
pylith::meshio::TestMeshIO::_checkVals(void) {
    PYLITH_METHOD_BEGIN;
    assert(_mesh);
    assert(_data);

    // Check mesh dimension
    CHECK(_data->topology->dimension == size_t(_mesh->getDimension()));
    const int spaceDim = _data->geometry->spaceDim;

    PetscDM dmMesh = _mesh->getDM();assert(dmMesh);

    // Check vertices
    topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
    const PylithInt vStart = verticesStratum.begin();
    const PylithInt vEnd = verticesStratum.end();

    REQUIRE(_data->geometry->numVertices == size_t(verticesStratum.size()));

    topology::CoordsVisitor coordsVisitor(dmMesh);
    const PetscScalar* coordsArray = coordsVisitor.localArray();
    const PylithScalar tolerance = 1.0e-06;
    for (PylithInt v = vStart, index = 0; v < vEnd; ++v) {
        const PylithInt off = coordsVisitor.sectionOffset(v);
        REQUIRE(spaceDim == coordsVisitor.sectionDof(v));

        for (int iDim = 0; iDim < spaceDim; ++iDim, ++index) {
            const double vtolerance = std::max(tolerance, fabs(_data->geometry->vertices[index])*tolerance);
            CHECK_THAT(coordsArray[off+iDim], Catch::Matchers::WithinAbs(_data->geometry->vertices[index], vtolerance));
        } // for
    } // for

    // Check cells
    topology::Stratum cellsStratum(dmMesh, topology::Stratum::HEIGHT, 0);
    const PylithInt cStart = cellsStratum.begin();
    const PylithInt cEnd = cellsStratum.end();
    const size_t numCells = cellsStratum.size();

    REQUIRE(_data->topology->numCells == numCells);
    const int offset = numCells;
    PetscErrorCode err = 0;
    for (PylithInt c = cStart, index = 0; c < cEnd; ++c) {
        PylithInt *closure = NULL;
        PylithInt closureSize;
        size_t numCorners = 0;

        err = DMPlexGetTransitiveClosure(dmMesh, c, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
        for (PylithInt p = 0; p < closureSize*2; p += 2) {
            const PylithInt point = closure[p];
            if ((point >= vStart) && (point < vEnd)) {
                closure[numCorners++] = point;
            } // if
        } // for
        err = DMPlexInvertCell_Private(_data->topology->dimension, numCorners, closure);PYLITH_CHECK_ERROR(err);
        REQUIRE(_data->topology->numCorners == numCorners);
        for (size_t p = 0; p < numCorners; ++p) {
            CHECK(_data->topology->cells[index++] == closure[p]-offset);
        } // for
        err = DMPlexRestoreTransitiveClosure(dmMesh, c, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
    } // for

    // check materials
    PylithInt matId = 0;
    for (PylithInt c = cStart; c < cEnd; ++c) {
        err = DMGetLabelValue(dmMesh, pylith::topology::Mesh::cells_label_name, c, &matId);PYLITH_CHECK_ERROR(err);
        CHECK(_data->materialIds[c-cStart] == matId);
    } // for

    // Check vertex groups
    pylith::string_vector vertexGroupNames;
    pylith::meshio::MeshBuilder::getVertexGroupNames(&vertexGroupNames, *_mesh);
    const size_t numVertexGroups = vertexGroupNames.size();
    REQUIRE(_data->numVertexGroups == numVertexGroups);
    for (size_t index = 0, iGroup = 0; iGroup < numVertexGroups; ++iGroup) {
        const char* groupName = _data->vertexGroupNames[iGroup];
        INFO("Checking " << groupName);

        PetscBool hasLabel = PETSC_TRUE;
        err = DMHasLabel(dmMesh, groupName, &hasLabel);assert(!err);
        REQUIRE(hasLabel);
        PetscDMLabel label = NULL;
        err = DMGetLabel(dmMesh, groupName, &label);assert(!err);

        const PylithInt labelValue = _data->vertexGroupTags ? _data->vertexGroupTags[iGroup] : 1;
        int_array points;
        pylith::meshio::MeshBuilder::getVertexGroup(&points, *_mesh, groupName, labelValue);
        const size_t numPoints = _data->vertexGroupSizes[iGroup];
        REQUIRE(numPoints == points.size());
        for (size_t p = 0; p < numPoints; ++p) {
            CHECK(_data->vertexGroups[index++] == points[p]);
        } // for
    } // for

    // Check face groups
    pylith::string_vector faceGroupNames;
    pylith::meshio::MeshBuilder::getFaceGroupNames(&faceGroupNames, *_mesh);
    const size_t numFaceGroups = faceGroupNames.size();
    REQUIRE(_data->numFaceGroups == numFaceGroups);
    if (_data->numFaceGroups > 0) {
        REQUIRE(_data->numFaceVertices > 0);
    } // if
    for (size_t index = 0, iGroup = 0; iGroup < numFaceGroups; ++iGroup) {
        const char* groupName = _data->faceGroupNames[iGroup];
        INFO("Checking " << groupName);

        PetscBool hasLabel = PETSC_TRUE;
        err = DMHasLabel(dmMesh, groupName, &hasLabel);assert(!err);
        REQUIRE(hasLabel);
        PetscDMLabel label = NULL;
        err = DMGetLabel(dmMesh, groupName, &label);assert(!err);

        const PylithInt labelValue = _data->faceGroupTags ? _data->faceGroupTags[iGroup] : 1;
        int_array faceValues;
        pylith::meshio::MeshBuilder::getFaceGroup(&faceValues, *_mesh, groupName, labelValue);
        const size_t numFaces = _data->faceGroupSizes[iGroup];
        const size_t totalSize = numFaces*(1+_data->numFaceVertices);
        REQUIRE(totalSize == faceValues.size() );
        for (size_t i = 0; i < totalSize; ++i) {
            CHECK(_data->faceGroups[index++] == faceValues[i]);
        } // for
    } // for

    PYLITH_METHOD_END;
} // _checkVals


// ----------------------------------------------------------------------
// Constructor
pylith::meshio::TestMeshIO_Data::TestMeshIO_Data(void) :
    topology(NULL),
    geometry(NULL),
    materialIds(NULL),

    vertexGroups(NULL),
    vertexGroupSizes(NULL),
    vertexGroupTags(NULL),
    vertexGroupNames(NULL),
    numVertexGroups(0),

    faceGroups(NULL),
    faceGroupSizes(NULL),
    faceGroupTags(NULL),
    faceGroupNames(NULL),
    numFaceGroups(0),

    useIndexZero(true),
    filename("") {}


// ----------------------------------------------------------------------
// Destructor
pylith::meshio::TestMeshIO_Data::~TestMeshIO_Data(void) {
    delete topology;
    delete geometry;
}


// End of file
