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

#include <strings.h> // USES strcasecmp()
#include <stdexcept> // USES std::logic_error

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
// Setup testing data.
void
pylith::meshio::TestMeshIO::setUp(void) { // setUp
    PYLITH_METHOD_BEGIN;

    _mesh = NULL;

    PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::meshio::TestMeshIO::tearDown(void) { // tearDown
    PYLITH_METHOD_BEGIN;

    delete _mesh;_mesh = NULL;

    PYLITH_METHOD_END;
} // tearDown


// ----------------------------------------------------------------------
// Get simple mesh for testing I/O.
void
pylith::meshio::TestMeshIO::_createMesh(void) { // _createMesh
    PYLITH_METHOD_BEGIN;

    const TestMeshIO_Data* data = _getData();CPPUNIT_ASSERT(data);

    // buildTopology() requires zero based index
    CPPUNIT_ASSERT_EQUAL(true, data->useIndexZero);

    CPPUNIT_ASSERT(data->vertices);
    CPPUNIT_ASSERT(data->cells);
    CPPUNIT_ASSERT(data->materialIds);
    if (data->numGroups > 0) {
        CPPUNIT_ASSERT(data->groups);
        CPPUNIT_ASSERT(data->groupSizes);
        CPPUNIT_ASSERT(data->groupNames);
        CPPUNIT_ASSERT(data->groupTypes);
    } // if

    delete _mesh;_mesh = new topology::Mesh(data->cellDim);CPPUNIT_ASSERT(_mesh);

    // Cells and vertices
    const bool interpolate = false;
    PetscDM dmMesh = NULL;
    PetscBool interpolateMesh = interpolate ? PETSC_TRUE : PETSC_FALSE;
    PylithInt bound = data->numCells*data->numCorners;
    PetscErrorCode err;

    PylithInt *cells = new int[bound];
    for (PylithInt coff = 0; coff < bound; ++coff) {
        cells[coff] = data->cells[coff];
    }
    for (PylithInt coff = 0; coff < bound; coff += data->numCorners) {
        err = DMPlexInvertCell_Private(data->cellDim, data->numCorners, &cells[coff]);PYLITH_CHECK_ERROR(err);
    } // for
    err = DMPlexCreateFromCellListPetsc(_mesh->getComm(), data->cellDim, data->numCells, data->numVertices, data->numCorners, interpolateMesh, cells, data->spaceDim, data->vertices, &dmMesh);PYLITH_CHECK_ERROR(err);
    delete [] cells;
    _mesh->setDM(dmMesh);

    // Material ids
    PylithInt cStart, cEnd;
    err = DMPlexGetHeightStratum(dmMesh, 0, &cStart, &cEnd);PYLITH_CHECK_ERROR(err);
    for (PylithInt c = cStart; c < cEnd; ++c) {
        err = DMSetLabelValue(dmMesh, "material-id", c, data->materialIds[c-cStart]);PYLITH_CHECK_ERROR(err);
    } // for

    // Groups
    for (int iGroup = 0, index = 0; iGroup < data->numGroups; ++iGroup) {
        err = DMCreateLabel(dmMesh, data->groupNames[iGroup]);PYLITH_CHECK_ERROR(err);

        const int numPoints = data->groupSizes[iGroup];
        if (0 == strcasecmp("cell", data->groupTypes[iGroup])) {
            for (int i = 0; i < numPoints; ++i, ++index) {
                err = DMSetLabelValue(dmMesh, data->groupNames[iGroup], data->groups[index], 1);PYLITH_CHECK_ERROR(err);
            } // for
        } else if (0 == strcasecmp("vertex", data->groupTypes[iGroup])) {
            PylithInt numCells;
            err = DMPlexGetHeightStratum(dmMesh, 0, NULL, &numCells);PYLITH_CHECK_ERROR(err);
            for (int i = 0; i < numPoints; ++i, ++index) {
                err = DMSetLabelValue(dmMesh, data->groupNames[iGroup], data->groups[index]+numCells, 1);PYLITH_CHECK_ERROR(err);
            } // for
        } else {
            throw std::logic_error("Could not parse group type.");
        } // else
    } // for
      // Set coordinate system
    spatialdata::geocoords::CSCart cs;
    cs.setSpaceDim(data->spaceDim);
    _mesh->setCoordSys(&cs);

    spatialdata::units::Nondimensional normalizer;
    normalizer.setLengthScale(10.0);
    topology::MeshOps::nondimensionalize(_mesh, normalizer);

    PYLITH_METHOD_END;
} // _createMesh


// ----------------------------------------------------------------------
// Check values in mesh against data->
void
pylith::meshio::TestMeshIO::_checkVals(void) { // _checkVals
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_mesh);

    const TestMeshIO_Data* data = _getData();CPPUNIT_ASSERT(data);

    // Check mesh dimension
    CPPUNIT_ASSERT_EQUAL(data->cellDim, _mesh->getDimension());
    const int spaceDim = data->spaceDim;

    PetscDM dmMesh = _mesh->getDM();CPPUNIT_ASSERT(dmMesh);

    // Check vertices
    topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
    const PylithInt vStart = verticesStratum.begin();
    const PylithInt vEnd = verticesStratum.end();

    CPPUNIT_ASSERT_EQUAL(data->numVertices, verticesStratum.size());

    topology::CoordsVisitor coordsVisitor(dmMesh);
    const PetscScalar* coordsArray = coordsVisitor.localArray();
    const PylithScalar tolerance = 1.0e-06;
    for (PylithInt v = vStart, index = 0; v < vEnd; ++v) {
        const PylithInt off = coordsVisitor.sectionOffset(v);
        CPPUNIT_ASSERT_EQUAL(spaceDim, coordsVisitor.sectionDof(v));

        for (int iDim = 0; iDim < spaceDim; ++iDim, ++index) {
            if (data->vertices[index] < 1.0) {
                CPPUNIT_ASSERT_DOUBLES_EQUAL(data->vertices[index], coordsArray[off+iDim], tolerance);
            } else {
                CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, coordsArray[off+iDim]/data->vertices[index], tolerance);
            } // if/else
        } // for
    } // for

    // Check cells
    topology::Stratum cellsStratum(dmMesh, topology::Stratum::HEIGHT, 0);
    const PylithInt cStart = cellsStratum.begin();
    const PylithInt cEnd = cellsStratum.end();
    const PylithInt numCells = cellsStratum.size();

    CPPUNIT_ASSERT_EQUAL(data->numCells, numCells);
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
        err = DMPlexInvertCell_Private(data->cellDim, numCorners, closure);PYLITH_CHECK_ERROR(err);
        CPPUNIT_ASSERT_EQUAL(data->numCorners, numCorners);
        for (PylithInt p = 0; p < numCorners; ++p, ++index) {
            CPPUNIT_ASSERT_EQUAL(data->cells[index], closure[p]-offset);
        } // for
        err = DMPlexRestoreTransitiveClosure(dmMesh, c, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
    } // for

    // check materials
    PylithInt matId = 0;
    for (PylithInt c = cStart; c < cEnd; ++c) {
        err = DMGetLabelValue(dmMesh, "material-id", c, &matId);PYLITH_CHECK_ERROR(err);
        CPPUNIT_ASSERT_EQUAL(data->materialIds[c-cStart], matId);
    } // for

    // Check groups
    PylithInt numGroups, pStart, pEnd;
    err = DMPlexGetChart(dmMesh, &pStart, &pEnd);PYLITH_CHECK_ERROR(err);
    err = DMGetNumLabels(dmMesh, &numGroups);PYLITH_CHECK_ERROR(err);
    numGroups -= 3; // Remove depth, celltype and material labels.
    CPPUNIT_ASSERT_EQUAL(data->numGroups, numGroups);
    PylithInt index = 0;
    for (PylithInt iGroup = 0; iGroup < numGroups; ++iGroup) {
        const char* groupName = data->groupNames[iGroup];

        PetscBool hasLabel = PETSC_TRUE;
        err = DMHasLabel(dmMesh, groupName, &hasLabel);CPPUNIT_ASSERT(!err);
        if (!hasLabel) {
            std::ostringstream msg;
            msg << "Mesh missing label '" << groupName << "'.";
            CPPUNIT_ASSERT_MESSAGE(msg.str().c_str(), hasLabel);
        } // if
        PetscDMLabel label = NULL;
        err = DMGetLabel(dmMesh, groupName, &label);CPPUNIT_ASSERT(!err);

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
        CPPUNIT_ASSERT_EQUAL(std::string(data->groupTypes[iGroup]), groupType);
        PylithInt numPoints, numVertices = 0;
        err = DMGetStratumSize(dmMesh, groupName, 1, &numPoints);PYLITH_CHECK_ERROR(err);
        PetscIS pointIS = NULL;
        const PylithInt *points = NULL;
        const PylithInt offset = ("vertex" == groupType) ? numCells : 0;
        err = DMGetStratumIS(dmMesh, groupName, 1, &pointIS);PYLITH_CHECK_ERROR(err);
        err = ISGetIndices(pointIS, &points);PYLITH_CHECK_ERROR(err);
        for (PylithInt p = 0; p < numPoints; ++p) {
            const PylithInt pStart = ("vertex" == groupType) ? vStart : cStart;
            const PylithInt pEnd = ("vertex" == groupType) ? vEnd   : cEnd;
            if ((points[p] >= pStart) && (points[p] < pEnd)) {
                CPPUNIT_ASSERT_EQUAL(data->groups[index++], points[p]-offset);
                ++numVertices;
            }
        } // for
        CPPUNIT_ASSERT_EQUAL(data->groupSizes[iGroup], numVertices);
        err = ISRestoreIndices(pointIS, &points);PYLITH_CHECK_ERROR(err);
        err = ISDestroy(&pointIS);PYLITH_CHECK_ERROR(err);
    } // for

    PYLITH_METHOD_END;
} // _checkVals


// ----------------------------------------------------------------------
// Test debug()
void
pylith::meshio::TestMeshIO::_testDebug(MeshIO& iohandler) { // _testDebug
    PYLITH_METHOD_BEGIN;

    bool debug = false;
    iohandler.debug(debug);
    CPPUNIT_ASSERT_EQUAL(debug, iohandler.debug());

    debug = true;
    iohandler.debug(debug);
    CPPUNIT_ASSERT_EQUAL(debug, iohandler.debug());

    PYLITH_METHOD_END;
} // _testDebug


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
    groupNames(NULL),
    groupTypes(NULL),
    numGroups(0),
    useIndexZero(true) { // constructor
} // constructor


// ----------------------------------------------------------------------
// Destructor
pylith::meshio::TestMeshIO_Data::~TestMeshIO_Data(void) { // destructor
} // destructor


// End of file
