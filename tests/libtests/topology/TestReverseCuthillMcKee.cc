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

#include "TestReverseCuthillMcKee.hh" // Implementation of class methods

#include "pylith/topology/ReverseCuthillMcKee.hh" // USES ReverseCuthillMcKee

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "tests/src/FaultCohesiveStub.hh" // USES FaultCohesiveStub
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor

#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers_floating_point.hpp"
#include "catch2/matchers/catch_matchers_exception.hpp"

// ------------------------------------------------------------------------------------------------
// Setup testing data.
pylith::topology::TestReverseCuthillMcKee::TestReverseCuthillMcKee(TestReverseCuthillMcKee_Data* data) :
    _data(data) {
    PYLITH_METHOD_BEGIN;
    assert(_data);

    _mesh = NULL;

    PYLITH_METHOD_END;
} // setUp


// ------------------------------------------------------------------------------------------------
// Tear down testing data.
pylith::topology::TestReverseCuthillMcKee::~TestReverseCuthillMcKee(void) {
    PYLITH_METHOD_BEGIN;

    delete _data;_data = NULL;
    delete _mesh;_mesh = NULL;

    PYLITH_METHOD_END;
} // tearDown


// ------------------------------------------------------------------------------------------------
// Test reorder().
void
pylith::topology::TestReverseCuthillMcKee::testReorder(void) {
    PYLITH_METHOD_BEGIN;

    _initialize();
    assert(_mesh);

    // Get original DM and create Mesh for it
    const PetscDM dmOrig = _mesh->getDM();
    PetscObjectReference((PetscObject) dmOrig);
    Mesh meshOrig;
    meshOrig.setDM(dmOrig);

    ReverseCuthillMcKee::reorder(_mesh);

    const PetscDM& dmMesh = _mesh->getDM();assert(dmMesh);

    // Check vertices (size only)
    topology::Stratum verticesStratumE(dmOrig, topology::Stratum::DEPTH, 0);
    topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
    CHECK(verticesStratumE.size() == verticesStratum.size());

    // Check cells (size only)
    topology::Stratum cellsStratumE(dmOrig, topology::Stratum::HEIGHT, 0);
    topology::Stratum cellsStratum(dmMesh, topology::Stratum::HEIGHT, 0);
    CHECK(cellsStratumE.size() == cellsStratum.size());

    // Check groups
    PetscInt numGroupsE, numGroups;
    PetscErrorCode err;
    err = DMGetNumLabels(dmOrig, &numGroupsE);REQUIRE(!err);
    err = DMGetNumLabels(dmMesh, &numGroups);REQUIRE(!err);
    REQUIRE(numGroupsE == numGroups);

    for (PetscInt iGroup = 0; iGroup < numGroups; ++iGroup) {
        const char *name = NULL;
        err = DMGetLabelName(dmMesh, iGroup, &name);REQUIRE(!err);

        PetscInt numPointsE, numPoints;
        err = DMGetStratumSize(dmOrig, name, 1, &numPointsE);REQUIRE(!err);
        err = DMGetStratumSize(dmMesh, name, 1, &numPoints);REQUIRE(!err);
        CHECK(numPointsE == numPoints);
    } // for

    // Check element centroids
    PylithScalar coordsCheckOrig = 0.0;
    PylithInt numCellsOrig = 0;
    PylithInt totalClosureSizeOrig = 0;
    { // original
        Stratum cellsStratum(dmOrig, Stratum::HEIGHT, 0);
        const PetscInt cStart = cellsStratum.begin();
        const PetscInt cEnd = cellsStratum.end();
        numCellsOrig = cEnd - cStart;
        pylith::topology::CoordsVisitor coordsVisitor(dmOrig);
        for (PetscInt cell = cStart; cell < cEnd; ++cell) {
            PetscScalar* coordsCell = NULL;
            PetscInt coordsSize = 0;
            PylithScalar value = 0.0;
            coordsVisitor.getClosure(&coordsCell, &coordsSize, cell);
            totalClosureSizeOrig += coordsSize;
            for (int i = 0; i < coordsSize; ++i) {
                value += coordsCell[i];
            } // for
            coordsCheckOrig += value*value;
            coordsVisitor.restoreClosure(&coordsCell, &coordsSize, cell);
        } // for
    } // original
    PylithScalar coordsCheckReorder = 0.0;
    PylithInt numCellsReorder = 0;
    PylithInt totalClosureSizeReorder = 0;
    { // reordered
        Stratum cellsStratum(dmMesh, Stratum::HEIGHT, 0);
        const PetscInt cStart = cellsStratum.begin();
        const PetscInt cEnd = cellsStratum.end();
        numCellsReorder = cEnd - cStart;
        pylith::topology::CoordsVisitor coordsVisitor(dmMesh);
        for (PetscInt cell = cStart; cell < cEnd; ++cell) {
            PetscScalar* coordsCell = NULL;
            PetscInt coordsSize = 0;
            PylithScalar value = 0.0;
            coordsVisitor.getClosure(&coordsCell, &coordsSize, cell);
            totalClosureSizeReorder += coordsSize;
            for (int i = 0; i < coordsSize; ++i) {
                value += coordsCell[i];
            } // for
            coordsCheckReorder += value*value;
            coordsVisitor.restoreClosure(&coordsCell, &coordsSize, cell);
        } // for
    } // reordered
    CHECK(numCellsOrig == numCellsReorder);
    CHECK(totalClosureSizeOrig == totalClosureSizeReorder);
    const PylithScalar tolerance = 1.0e-6;
    CHECK_THAT(coordsCheckReorder, Catch::Matchers::WithinAbs(coordsCheckOrig, tolerance*coordsCheckOrig));

    // Verify reduction in Jacobian bandwidth
    Field fieldOrig(meshOrig);
    Field::Description description;
    description.label = "solution";
    description.vectorFieldType = FieldBase::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = "field";
    description.scale = 1.0;
    description.validator = NULL;

    Field::Discretization discretization;
    discretization.basisOrder = 1;
    discretization.quadOrder = 1;
    fieldOrig.subfieldAdd(description, discretization);
    fieldOrig.subfieldsSetup();
    fieldOrig.createDiscretization();
    fieldOrig.allocate();
    PetscMat matrix = NULL;
    PetscInt bandwidthOrig = 0;
    err = DMCreateMatrix(fieldOrig.getDM(), &matrix);REQUIRE(!err);
    err = MatComputeBandwidth(matrix, 0.0, &bandwidthOrig);REQUIRE(!err);
    err = MatDestroy(&matrix);REQUIRE(!err);

    Field field(*_mesh);
    field.subfieldAdd(description, discretization);
    field.subfieldsSetup();
    field.createDiscretization();
    field.allocate();
    PetscInt bandwidth = 0;
    err = DMCreateMatrix(field.getDM(), &matrix);REQUIRE(!err);
    err = MatComputeBandwidth(matrix, 0.0, &bandwidth);REQUIRE(!err);
    err = MatDestroy(&matrix);REQUIRE(!err);

    REQUIRE(bandwidthOrig > 0);
    REQUIRE(bandwidth > 0);
    REQUIRE(bandwidth <= bandwidthOrig);

    PYLITH_METHOD_END;
} // testReorder


// ------------------------------------------------------------------------------------------------
void
pylith::topology::TestReverseCuthillMcKee::_initialize() {
    PYLITH_METHOD_BEGIN;
    assert(_data);

    delete _mesh;_mesh = new Mesh;assert(_mesh);

    meshio::MeshIOAscii iohandler;
    iohandler.setFilename(_data->filename);
    iohandler.read(_mesh);
    assert(pylith::topology::MeshOps::getNumCells(*_mesh) > 0);
    assert(pylith::topology::MeshOps::getNumVertices(*_mesh) > 0);

    // Adjust topology if necessary.
    if (_data->faultLabel) {
        pylith::faults::FaultCohesiveStub fault;
        fault.setCohesiveLabelValue(100);
        fault.setSurfaceLabelName(_data->faultLabel);
        fault.adjustTopology(_mesh);
    } // if

    PYLITH_METHOD_END;
} // _initialize


// ------------------------------------------------------------------------------------------------
// Constructor
pylith::topology::TestReverseCuthillMcKee_Data::TestReverseCuthillMcKee_Data(void) :
    filename(NULL),
    faultLabel(NULL) {}


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::topology::TestReverseCuthillMcKee_Data::~TestReverseCuthillMcKee_Data(void) {}


// End of file
