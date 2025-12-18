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
    REQUIRE(_data);

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
    REQUIRE(_mesh);

    pylith::topology::Mesh* meshNew = ReverseCuthillMcKee::reorder(*_mesh);assert(meshNew);

    const PetscDM dmOrig = _mesh->getDM();
    const PetscDM dmNew = meshNew->getDM();

    // Check vertices (size only)
    topology::Stratum verticesStratumE(dmOrig, topology::Stratum::DEPTH, 0);
    topology::Stratum verticesStratum(dmNew, topology::Stratum::DEPTH, 0);
    CHECK(verticesStratumE.size() == verticesStratum.size());

    // Check cells (size only)
    topology::Stratum cellsStratumE(dmOrig, topology::Stratum::HEIGHT, 0);
    topology::Stratum cellsStratum(dmNew, topology::Stratum::HEIGHT, 0);
    CHECK(cellsStratumE.size() == cellsStratum.size());

    // Check groups
    PetscInt numGroupsE, numGroups;
    PylithCallPetscRequire(DMGetNumLabels(dmOrig, &numGroupsE));
    PylithCallPetscRequire(DMGetNumLabels(dmNew, &numGroups));
    REQUIRE(numGroupsE == numGroups);

    for (PetscInt iGroup = 0; iGroup < numGroups; ++iGroup) {
        const char *name = NULL;
        PylithCallPetscRequire(DMGetLabelName(dmNew, iGroup, &name));

        PetscInt numPointsE, numPoints;
        PylithCallPetscRequire(DMGetStratumSize(dmOrig, name, 1, &numPointsE));
        PylithCallPetscRequire(DMGetStratumSize(dmNew, name, 1, &numPoints));
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
        Stratum cellsStratum(dmNew, Stratum::HEIGHT, 0);
        const PetscInt cStart = cellsStratum.begin();
        const PetscInt cEnd = cellsStratum.end();
        numCellsReorder = cEnd - cStart;
        pylith::topology::CoordsVisitor coordsVisitor(dmNew);
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

    Field fieldOrig(*_mesh);
    fieldOrig.subfieldAdd(description, discretization);
    fieldOrig.subfieldsSetup();
    fieldOrig.createDiscretization();
    fieldOrig.allocate();
    PetscMat matrix = NULL;
    PetscInt bandwidthOrig = 0;
    PylithCallPetscRequire(DMCreateMatrix(fieldOrig.getDM(), &matrix));
    PylithCallPetscRequire(MatComputeBandwidth(matrix, 0.0, &bandwidthOrig));
    PylithCallPetscRequire(MatDestroy(&matrix));

    Field field(*meshNew);
    field.subfieldAdd(description, discretization);
    field.subfieldsSetup();
    field.createDiscretization();
    field.allocate();
    PetscInt bandwidth = 0;
    PylithCallPetscRequire(DMCreateMatrix(field.getDM(), &matrix));
    PylithCallPetscRequire(MatComputeBandwidth(matrix, 0.0, &bandwidth));
    PylithCallPetscRequire(MatDestroy(&matrix));

    CHECK(bandwidthOrig > 0);
    CHECK(bandwidth > 0);
    CHECK(bandwidth <= bandwidthOrig);

    delete meshNew;meshNew = nullptr;
    PYLITH_METHOD_END;
} // testReorder


// ------------------------------------------------------------------------------------------------
void
pylith::topology::TestReverseCuthillMcKee::_initialize(void) {
    PYLITH_METHOD_BEGIN;
    REQUIRE(_data);

    delete _mesh;_mesh = new Mesh;REQUIRE(_mesh);

    meshio::MeshIOAscii iohandler;
    iohandler.setFilename(_data->filename);
    iohandler.read(_mesh);
    REQUIRE(pylith::topology::MeshOps::getNumCells(*_mesh) > 0);
    REQUIRE(pylith::topology::MeshOps::getNumVertices(*_mesh) > 0);

    // Adjust topology if necessary.
    if (_data->faultLabel) {
        pylith::faults::FaultCohesiveStub fault;
        fault.setCohesiveLabelValue(100);
        fault.setSurfaceLabelName(_data->faultLabel);
        pylith::topology::Mesh* meshNew = fault.transformTopology(_mesh);
        delete _mesh;_mesh = meshNew;
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
