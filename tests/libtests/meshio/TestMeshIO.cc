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
// Copyright (c) 2010-2022 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestMeshIO.hh" // Implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps::nondimensionalize()
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor
#include "pylith/meshio/MeshIO.hh" // USES MeshIO

#include "pylith/utils/array.hh" // USES int_array

#include "data/MeshData.hh"

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

    assert(_data->vertices);
    assert(_data->cells);
    assert(_data->materialIds);
    if (_data->numGroups > 0) {
        assert(_data->groups);
        assert(_data->groupSizes);
        assert(_data->groupNames);
        assert(_data->groupTypes);
    } // if

    delete _mesh;_mesh = new topology::Mesh(_data->cellDim);assert(_mesh);

    // Cells and vertices
    PetscDM dmMesh = NULL;
    const PetscBool interpolateMesh = PETSC_TRUE;
    PylithInt bound = _data->numCells * _data->numCorners;
    PetscErrorCode err = PETSC_SUCCESS;

    PylithInt *cells = new PylithInt[bound];
    for (PylithInt coff = 0; coff < bound; ++coff) {
        cells[coff] = _data->cells[coff];
    }
    for (PylithInt coff = 0; coff < bound; coff += _data->numCorners) {
        err = DMPlexInvertCell_Private(_data->cellDim, _data->numCorners, &cells[coff]);PYLITH_CHECK_ERROR(err);
    } // for
    err = DMPlexCreateFromCellListPetsc(_mesh->getComm(), _data->cellDim, _data->numCells, _data->numVertices, _data->numCorners, interpolateMesh, cells, _data->spaceDim, _data->vertices, &dmMesh);PYLITH_CHECK_ERROR(err);
    delete [] cells;
    _mesh->setDM(dmMesh);

    // Material ids
    PylithInt cStart, cEnd;
    err = DMPlexGetHeightStratum(dmMesh, 0, &cStart, &cEnd);PYLITH_CHECK_ERROR(err);
    for (PylithInt c = cStart; c < cEnd; ++c) {
        err = DMSetLabelValue(dmMesh, pylith::topology::Mesh::cells_label_name, c, _data->materialIds[c-cStart]);PYLITH_CHECK_ERROR(err);
    } // for

    // Groups
    for (int iGroup = 0, index = 0; iGroup < _data->numGroups; ++iGroup) {
        err = DMCreateLabel(dmMesh, _data->groupNames[iGroup]);PYLITH_CHECK_ERROR(err);

        const int numPoints = _data->groupSizes[iGroup];
        if (0 == strcasecmp("cell", _data->groupTypes[iGroup])) {
            for (int i = 0; i < numPoints; ++i, ++index) {
                err = DMSetLabelValue(dmMesh, _data->groupNames[iGroup], _data->groups[index], 1);PYLITH_CHECK_ERROR(err);
            } // for
        } else if (0 == strcasecmp("vertex", _data->groupTypes[iGroup])) {
            PylithInt numCells;
            err = DMPlexGetHeightStratum(dmMesh, 0, NULL, &numCells);PYLITH_CHECK_ERROR(err);
            for (int i = 0; i < numPoints; ++i, ++index) {
                err = DMSetLabelValue(dmMesh, _data->groupNames[iGroup], _data->groups[index]+numCells, 1);PYLITH_CHECK_ERROR(err);
            } // for
        } else {
            throw std::logic_error("Could not parse group type.");
        } // else
    } // for
      // Set coordinate system
    spatialdata::geocoords::CSCart cs;
    cs.setSpaceDim(_data->spaceDim);
    _mesh->setCoordSys(&cs);

    spatialdata::units::Nondimensional normalizer;
    normalizer.setLengthScale(10.0);
    topology::MeshOps::nondimensionalize(_mesh, normalizer);

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
    CHECK(_data->cellDim == _mesh->getDimension());
    const int spaceDim = _data->spaceDim;

    PetscDM dmMesh = _mesh->getDM();assert(dmMesh);

    // Check vertices
    topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
    const PylithInt vStart = verticesStratum.begin();
    const PylithInt vEnd = verticesStratum.end();

    REQUIRE(_data->numVertices == verticesStratum.size());

    topology::CoordsVisitor coordsVisitor(dmMesh);
    const PetscScalar* coordsArray = coordsVisitor.localArray();
    const PylithScalar tolerance = 1.0e-06;
    for (PylithInt v = vStart, index = 0; v < vEnd; ++v) {
        const PylithInt off = coordsVisitor.sectionOffset(v);
        REQUIRE(spaceDim == coordsVisitor.sectionDof(v));

        for (int iDim = 0; iDim < spaceDim; ++iDim, ++index) {
            const double vtolerance = std::max(tolerance, fabs(_data->vertices[index])*tolerance);
            CHECK_THAT(coordsArray[off+iDim], Catch::Matchers::WithinAbs(_data->vertices[index], vtolerance));
        } // for
    } // for

    // Check cells
    topology::Stratum cellsStratum(dmMesh, topology::Stratum::HEIGHT, 0);
    const PylithInt cStart = cellsStratum.begin();
    const PylithInt cEnd = cellsStratum.end();
    const PylithInt numCells = cellsStratum.size();

    REQUIRE(_data->numCells == numCells);
    const int offset = numCells;
    PetscErrorCode err = 0;
    for (PylithInt c = cStart, index = 0; c < cEnd; ++c) {
        PylithInt *closure = NULL;
        PylithInt closureSize, numCorners = 0;

        err = DMPlexGetTransitiveClosure(dmMesh, c, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
        for (PylithInt p = 0; p < closureSize*2; p += 2) {
            const PylithInt point = closure[p];
            if ((point >= vStart) && (point < vEnd)) {
                closure[numCorners++] = point;
            } // if
        } // for
        err = DMPlexInvertCell_Private(_data->cellDim, numCorners, closure);PYLITH_CHECK_ERROR(err);
        REQUIRE(_data->numCorners == numCorners);
        for (PylithInt p = 0; p < numCorners; ++p, ++index) {
            CHECK(_data->cells[index] == closure[p]-offset);
        } // for
        err = DMPlexRestoreTransitiveClosure(dmMesh, c, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
    } // for

    // check materials
    PylithInt matId = 0;
    for (PylithInt c = cStart; c < cEnd; ++c) {
        err = DMGetLabelValue(dmMesh, pylith::topology::Mesh::cells_label_name, c, &matId);PYLITH_CHECK_ERROR(err);
        CHECK(_data->materialIds[c-cStart] == matId);
    } // for

    // Check groups
    PylithInt numGroups, pStart, pEnd;
    err = DMPlexGetChart(dmMesh, &pStart, &pEnd);PYLITH_CHECK_ERROR(err);
    err = DMGetNumLabels(dmMesh, &numGroups);PYLITH_CHECK_ERROR(err);
    numGroups -= 3; // Remove depth, celltype and material labels.
    REQUIRE(_data->numGroups == numGroups);
    PylithInt index = 0;
    for (PylithInt iGroup = 0; iGroup < numGroups; ++iGroup) {
        const char* groupName = _data->groupNames[iGroup];
        INFO("Checking " << groupName);

        PetscBool hasLabel = PETSC_TRUE;
        err = DMHasLabel(dmMesh, groupName, &hasLabel);assert(!err);
        assert(hasLabel);
        PetscDMLabel label = NULL;
        err = DMGetLabel(dmMesh, groupName, &label);assert(!err);

        PylithInt firstPoint = 0;

        for (PylithInt p = pStart; p < pEnd; ++p) {
            PylithInt val;

            err = DMGetLabelValue(dmMesh, groupName, p, &val);PYLITH_CHECK_ERROR(err);
            if (val >= 0) {
                firstPoint = p;
                break;
            } // if
        } // for
        std::string groupType = (firstPoint >= cStart && firstPoint < cEnd) ? "cell" : "vertex";
        const PylithInt labelValue = _data->groupTags ? _data->groupTags[iGroup] : 1;
        REQUIRE(std::string(_data->groupTypes[iGroup]) == groupType);
        PylithInt numPoints, numVertices = 0;
        err = DMGetStratumSize(dmMesh, groupName, labelValue, &numPoints);PYLITH_CHECK_ERROR(err);
        assert(numPoints > 0);
        PetscIS pointIS = NULL;
        const PylithInt *points = NULL;
        const PylithInt offset = ("vertex" == groupType) ? numCells : 0;
        err = DMGetStratumIS(dmMesh, groupName, labelValue, &pointIS);PYLITH_CHECK_ERROR(err);
        err = ISGetIndices(pointIS, &points);PYLITH_CHECK_ERROR(err);
        for (PylithInt p = 0; p < numPoints; ++p) {
            const PylithInt pStart = ("vertex" == groupType) ? vStart : cStart;
            const PylithInt pEnd = ("vertex" == groupType) ? vEnd   : cEnd;
            if ((points[p] >= pStart) && (points[p] < pEnd)) {
                CHECK(_data->groups[index++] == points[p]-offset);
                ++numVertices;
            }
        } // for
        CHECK(_data->groupSizes[iGroup] == numVertices);
        err = ISRestoreIndices(pointIS, &points);PYLITH_CHECK_ERROR(err);
        err = ISDestroy(&pointIS);PYLITH_CHECK_ERROR(err);
    } // for

    PYLITH_METHOD_END;
} // _checkVals


// ----------------------------------------------------------------------
// Constructor
pylith::meshio::TestMeshIO_Data::TestMeshIO_Data(void) :
    numVertices(0),
    spaceDim(0),
    numCells(0),
    cellDim(0),
    numCorners(0),
    vertices(NULL),
    cells(NULL),
    materialIds(NULL),
    groups(NULL),
    groupSizes(NULL),
    groupTags(NULL),
    groupNames(NULL),
    groupTypes(NULL),
    numGroups(0),
    useIndexZero(true),
    filename("") {}


// ----------------------------------------------------------------------
// Destructor
pylith::meshio::TestMeshIO_Data::~TestMeshIO_Data(void) {}


// End of file
