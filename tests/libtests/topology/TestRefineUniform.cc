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

#include "TestRefineUniform.hh" // Implementation of class methods

#include "pylith/topology/RefineUniform.hh" // USES RefineUniform

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/testing/FaultCohesiveStub.hh" // USES FaultCohesiveStub

#include "pylith/utils/array.hh" // USES int_array

#include <strings.h> // USES strcasecmp()
#include <stdexcept> // USES std::logic_error

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::topology::TestRefineUniform::setUp(void) {
    _data = new TestRefineUniform_Data;CPPUNIT_ASSERT(_data);
} // setUp


// ----------------------------------------------------------------------
// Deallocate testing data.
void
pylith::topology::TestRefineUniform::tearDown(void) {
    delete _data;_data = NULL;
} // tearDown


// ----------------------------------------------------------------------
// Test refine().
void
pylith::topology::TestRefineUniform::testRefine(void) {
    PYLITH_METHOD_BEGIN;
    CPPUNIT_ASSERT(_data);

    Mesh mesh(_data->cellDim);
    _initializeMesh(&mesh);

    RefineUniform refiner;
    Mesh newMesh(_data->cellDim);
    refiner.refine(&newMesh, mesh, _data->refineLevel);

    // Check mesh dimension
    CPPUNIT_ASSERT_EQUAL(_data->cellDim, newMesh.getDimension());

    const PetscDM& dmMesh = newMesh.getDM();CPPUNIT_ASSERT(dmMesh);

    // Check vertices
    pylith::topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
    CPPUNIT_ASSERT_EQUAL(_data->numVertices, verticesStratum.size());
    const PetscInt vStart = verticesStratum.begin();
    const PetscInt vEnd = verticesStratum.end();

    pylith::topology::CoordsVisitor coordsVisitor(dmMesh);
    const int spaceDim = _data->spaceDim;
    for (PetscInt v = vStart; v < vEnd; ++v) {
        CPPUNIT_ASSERT_EQUAL(spaceDim, coordsVisitor.sectionDof(v));
    } // for

    // Check cells
    pylith::topology::Stratum cellsStratum(dmMesh, topology::Stratum::HEIGHT, 0);
    const PetscInt cStart = cellsStratum.begin();
    const PetscInt cEnd = cellsStratum.end();
    const PetscInt numCells = cellsStratum.size();

    CPPUNIT_ASSERT_EQUAL(_data->numCells+_data->numCellsCohesive, numCells);
    PetscErrorCode err;
    // Normal cells
    for (PetscInt c = cStart; c < _data->numCells; ++c) {
        DMPolytopeType ct;
        PetscInt *closure = NULL;
        PetscInt closureSize, numCorners = 0;

        err = DMPlexGetCellType(dmMesh, c, &ct);CPPUNIT_ASSERT(!err);
        err = DMPlexGetTransitiveClosure(dmMesh, c, PETSC_TRUE, &closureSize, &closure);CPPUNIT_ASSERT(!err);
        for (PetscInt p = 0; p < closureSize*2; p += 2) {
            const PetscInt point = closure[p];
            if ((point >= vStart) && (point < vEnd)) {
                closure[numCorners++] = point;
            } // if
        } // for
        err = DMPlexInvertCell(ct, closure);CPPUNIT_ASSERT(!err);
        CPPUNIT_ASSERT_EQUAL(_data->numCorners, numCorners);
        err = DMPlexRestoreTransitiveClosure(dmMesh, c, PETSC_TRUE, &closureSize, &closure);CPPUNIT_ASSERT(!err);
    } // for

    // Cohesive cells
    for (PetscInt c = _data->numCells; c < cEnd; ++c) {
        DMPolytopeType ct;
        PetscInt *closure = NULL;
        PetscInt closureSize, numCorners = 0;

        err = DMPlexGetCellType(dmMesh, c, &ct);CPPUNIT_ASSERT(!err);
        err = DMPlexGetTransitiveClosure(dmMesh, c, PETSC_TRUE, &closureSize, &closure);CPPUNIT_ASSERT(!err);
        for (PetscInt p = 0; p < closureSize*2; p += 2) {
            const PetscInt point = closure[p];
            if ((point >= vStart) && (point < vEnd)) {
                closure[numCorners++] = point;
            } // if
        } // for
        err = DMPlexInvertCell(ct, closure);CPPUNIT_ASSERT(!err);
        CPPUNIT_ASSERT_EQUAL(_data->numCornersCohesive, numCorners);
        err = DMPlexRestoreTransitiveClosure(dmMesh, c, PETSC_TRUE, &closureSize, &closure);CPPUNIT_ASSERT(!err);
    } // for

    // check materials
    PetscInt matId = 0;
    PetscInt matIdSum = 0; // Use sum of material ids as simple checksum.
    for (PetscInt c = cStart; c < cEnd; ++c) {
        err = DMGetLabelValue(dmMesh, "material-id", c, &matId);CPPUNIT_ASSERT(!err);
        matIdSum += matId;
    } // for
    CPPUNIT_ASSERT_EQUAL(_data->matIdSum, matIdSum);

    // Check groups
    PetscInt numGroups, pStart, pEnd;
    err = DMPlexGetChart(dmMesh, &pStart, &pEnd);CPPUNIT_ASSERT(!err);
    err = DMGetNumLabels(dmMesh, &numGroups);CPPUNIT_ASSERT(!err);
    for (PetscInt iGroup = 0; iGroup < _data->numGroups; ++iGroup) {
        // Omit depth, vtk, ghost and material-id labels
        // Don't know order of labels, so do brute force linear search
        bool foundLabel = false;
        int iLabel = 0;
        const char *name = NULL;
        PetscInt firstPoint = 0;

        while (iLabel < numGroups) {
            err = DMGetLabelName(dmMesh, iLabel, &name);CPPUNIT_ASSERT(!err);
            if (0 == strcmp(_data->groupNames[iGroup], name)) {
                foundLabel = true;
                break;
            } else {
                ++iLabel;
            } // if/else
        } // while
        CPPUNIT_ASSERT(foundLabel);

        for (PetscInt p = pStart; p < pEnd; ++p) {
            PetscInt val;
            err = DMGetLabelValue(dmMesh, name, p, &val);CPPUNIT_ASSERT(!err);
            if (val >= 0) {
                firstPoint = p;
                break;
            } // if
        } // for
        std::string groupType = (firstPoint >= cStart && firstPoint < cEnd) ? "cell" : "vertex";
        CPPUNIT_ASSERT_EQUAL(std::string(_data->groupTypes[iGroup]), groupType);
        PetscInt numPoints;
        err = DMGetStratumSize(dmMesh, name, 1, &numPoints);CPPUNIT_ASSERT(!err);
        CPPUNIT_ASSERT_EQUAL(_data->groupSizes[iGroup], numPoints);
        PetscIS pointIS = NULL;
        const PetscInt *points = NULL;
        err = DMGetStratumIS(dmMesh, name, 1, &pointIS);CPPUNIT_ASSERT(!err);
        err = ISGetIndices(pointIS, &points);CPPUNIT_ASSERT(!err);
        if (groupType == std::string("vertex")) {
            for (PetscInt p = 0; p < numPoints; ++p) {
                CPPUNIT_ASSERT((points[p] >= 0 && points[p] < cStart) || (points[p] >= cEnd));
            } // for
        } else {
            for (PetscInt p = 0; p < numPoints; ++p) {
                CPPUNIT_ASSERT(points[p] >= cStart && points[p] < cEnd);
            } // for
        } // if/else
        err = ISRestoreIndices(pointIS, &points);CPPUNIT_ASSERT(!err);
        err = ISDestroy(&pointIS);CPPUNIT_ASSERT(!err);
    } // for

    PYLITH_METHOD_END;
} // testRefine


// ----------------------------------------------------------------------
void
pylith::topology::TestRefineUniform::_initializeMesh(Mesh* const mesh) {
    PYLITH_METHOD_BEGIN;
    CPPUNIT_ASSERT(_data);
    CPPUNIT_ASSERT(mesh);

    pylith::meshio::MeshIOAscii iohandler;
    iohandler.filename(_data->filename);
    iohandler.read(mesh);

    // Adjust topology if necessary.
    if (_data->faultA) {
        faults::FaultCohesiveStub faultA;
        faultA.setInterfaceId(100);
        faultA.setSurfaceMarkerLabel(_data->faultA);
        faultA.adjustTopology(mesh);
    } // if

    if (_data->faultB) {
        faults::FaultCohesiveStub faultB;
        faultB.setInterfaceId(101);
        faultB.setSurfaceMarkerLabel(_data->faultB);
        faultB.adjustTopology(mesh);
    } // if

    PYLITH_METHOD_END;
} // _initializeMesh


// ----------------------------------------------------------------------
// Constructor
pylith::topology::TestRefineUniform_Data::TestRefineUniform_Data(void) :
    filename(NULL),
    refineLevel(0),
    faultA(NULL),
    faultB(NULL),
    isSimplexMesh(true),
    numVertices(0),
    spaceDim(0),
    cellDim(0),
    numCells(0),
    numCorners(0),
    numCellsCohesive(0),
    numCornersCohesive(0),
    matIdSum(0),
    groupSizes(NULL),
    groupNames(NULL),
    groupTypes(NULL),
    numGroups(0) { // constructor
} // constructor


// ----------------------------------------------------------------------
// Destructor
pylith::topology::TestRefineUniform_Data::~TestRefineUniform_Data(void) {} // destructor


// End of file
