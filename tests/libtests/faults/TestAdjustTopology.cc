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

#include "TestAdjustTopology.hh"

#include "tests/src/FaultCohesiveStub.hh" // USES FaultsCohesiveStub

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps::nondimensionalize()
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor
#include "pylith/faults/TopologyOps.hh" // USES TopologyOps

#include "pylith/utils/array.hh" // USES int_array, scalar_array
#include "pylith/utils/journals.hh" // USES journals
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

#include "pylith/scales/Scales.hh" // USES Scales
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart

#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers_exception.hpp"

#include <set> // USES std::set
#include <string> // USES std:;string

// ------------------------------------------------------------------------------------------------
// Constructor
pylith::faults::TestAdjustTopology::TestAdjustTopology(TestAdjustTopology_Data* data) :
    _data(data),
    _mesh(nullptr) {
    REQUIRE(_data);
    GenericComponent::setName("TestAdjustTopology");
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::faults::TestAdjustTopology::~TestAdjustTopology(void) {
    delete _data;_data = nullptr;
    delete _mesh;_mesh = nullptr;
} // tearDown


// ------------------------------------------------------------------------------------------------
// Run adjustTopology().
void
pylith::faults::TestAdjustTopology::run(void) {
    _initialize();
    REQUIRE(_mesh);
    REQUIRE(_data);

    for (size_t i = 0; i < _data->numFaults; ++i) {
        FaultCohesiveStub fault;

        REQUIRE(_data->interfaceIds);
        REQUIRE(_data->faultSurfaceLabels);
        REQUIRE(_data->faultEdgeLabels);

        fault.setCohesiveLabelName(pylith::topology::Mesh::cells_label_name);
        fault.setCohesiveLabelValue(_data->interfaceIds[i]);
        fault.setSurfaceLabelName(_data->faultSurfaceLabels[i]);
        fault.setSurfaceLabelValue(1);
        if (_data->faultEdgeLabels[i]) {
            fault.setBuriedEdgesLabelName(_data->faultEdgeLabels[i]);
            fault.setBuriedEdgesLabelValue(1);
        } // if
        if (!_data->failureExpected) {
            fault.adjustTopology(_mesh);
        } else {
            REQUIRE_THROWS_AS(fault.adjustTopology(_mesh), std::runtime_error);
            return;
        } // if/else
    } // for

    pythia::journal::debug_t debug(GenericComponent::getName());
    if (debug.state()) {
        _mesh->view(":mesh_transformed.txt:ascii_info_detail");
        _mesh->view(":mesh_transformed.tex:ascii_latex");
        _mesh->view("vtk:mesh_transformed.vtu");
    } // if

    REQUIRE(_data->cellDim == size_t(_mesh->getDimension()));
    PetscDM dmMesh = _mesh->getDM();REQUIRE(dmMesh);

    // Check vertices
    topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
    const PetscInt vStart = verticesStratum.begin();
    const PetscInt vEnd = verticesStratum.end();
    REQUIRE(_data->numVertices == size_t(verticesStratum.size()));

    topology::CoordsVisitor coordsVisitor(dmMesh);
    const PetscInt spaceDim = _data->spaceDim;

    for (PetscInt v = vStart; v < vEnd; ++v) {
        REQUIRE(spaceDim == coordsVisitor.sectionDof(v));
    } // for

    // check cells
    REQUIRE(_data->numCorners);
    pylith::topology::Stratum cellsStratum(dmMesh, pylith::topology::Stratum::HEIGHT, 0);
    const PetscInt cStart = cellsStratum.begin();
    const PetscInt cEnd = cellsStratum.end();
    REQUIRE(_data->numCells == size_t(cellsStratum.size()));
    for (PetscInt c = cStart, cell = 0; c < cEnd; ++c, ++cell) {
        PetscInt coneSize = 0;
        PylithCallPetsc(DMPlexGetConeSize(dmMesh, c, &coneSize));
        REQUIRE(_data->numCorners[cell] == coneSize);
    } // for

    // check materials
    REQUIRE(_data->materialIds);
    PetscDMLabel labelMaterials = NULL;
    const char* const cellsLabelName = pylith::topology::Mesh::cells_label_name;
    PylithCallPetsc(DMGetLabel(dmMesh, cellsLabelName, &labelMaterials));
    REQUIRE(labelMaterials);
    const PetscInt idDefault = -999;
    for (PetscInt c = cStart, cell = 0; c < cEnd; ++c, ++cell) {
        PetscInt labelValue;

        PylithCallPetsc(DMLabelGetValue(labelMaterials, c, &labelValue));
        if (labelValue == -1) {
            labelValue = idDefault;
        } // if
        double centroid[3] = {0.0, 0.0, 0.0};
        err = DMPlexComputeCellGeometryFVM(dmMesh, c, NULL, centroid, NULL);PYLITH_CHECK_ERROR(err);
        const PetscInt labelValueE = _data->getMatId(c, _data->numNoncohesiveCells, centroid);
        INFO("cell="<<c<<" with centroid ("<<centroid[0]<<","<<centroid[1]<<","<<centroid[2]<<")");
        CHECK(labelValueE == labelValue);
    } // for

    // Check groups
    REQUIRE(_data->groupSizes);
    REQUIRE(_data->groupNames);
    REQUIRE(_data->groupTypes);
    PetscInt numLabels;
    PylithCallPetsc(DMGetNumLabels(dmMesh, &numLabels));
    for (PetscInt iLabel = 0; iLabel < numLabels; ++iLabel) {
        PetscDMLabel label = NULL;
        PetscIS pointIS = NULL;
        const PetscInt *points = NULL;
        PetscInt numPoints = 0, depth = 0;
        const char *labelName = NULL;
        std::set<std::string> ignoreLabels;
        ignoreLabels.insert("depth");
        ignoreLabels.insert(pylith::topology::Mesh::cells_label_name);
        ignoreLabels.insert("vtk");
        ignoreLabels.insert("ghost");
        ignoreLabels.insert("dim");
        ignoreLabels.insert("celltype");

        PylithCallPetsc(DMGetLabelName(dmMesh, iLabel, &labelName));
        if (ignoreLabels.count(labelName) > 0) { continue; }
        PylithCallPetsc(DMGetLabel(dmMesh, labelName, &label));REQUIRE(label);
        PylithCallPetsc(DMLabelGetStratumIS(label, 1, &pointIS));
        PylithCallPetsc(ISGetLocalSize(pointIS, &numPoints));
        PylithCallPetsc(ISGetIndices(pointIS, &points));
        PylithCallPetsc(DMGetLabelValue(dmMesh, "depth", points[0], &depth));
        PylithCallPetsc(ISRestoreIndices(pointIS, &points));
        PylithCallPetsc(ISDestroy(&pointIS));
        std::string groupType = depth ? "face" : "vertex";

        bool foundGroup = false;
        for (size_t i = 0; i < _data->numGroups; ++i) {
            if (std::string(_data->groupNames[i]) == std::string(labelName)) {
                std::ostringstream msg;
                msg << "Checking group '" <<labelName<<"'.";
                INFO(labelName);
                CHECK(std::string(_data->groupTypes[i]) == groupType);
                CHECK(_data->groupSizes[i] == numPoints);
                foundGroup = true;
                break;
            } // if
        } // for
        INFO("Could not find mesh group '" << labelName << "' in test data.");
        CHECK(foundGroup);
    } // for

} // run


// ------------------------------------------------------------------------------------------------
// Test creation of cohesive cells using TransformApply().
void
pylith::faults::TestAdjustTopology::run_transform(void) {
    _initialize();
    assert(_mesh);
    assert(_data);

    for (size_t i = 0; i < _data->numFaults; ++i) {
        FaultCohesiveStub fault;

        assert(_data->interfaceIds);
        assert(_data->faultSurfaceLabels);
        assert(_data->faultEdgeLabels);

        fault.setCohesiveLabelName(pylith::topology::Mesh::cells_label_name);
        fault.setCohesiveLabelValue(_data->interfaceIds[i]);
        fault.setSurfaceLabelName(_data->faultSurfaceLabels[i]);
        fault.setSurfaceLabelValue(1);
        if (_data->faultEdgeLabels[i]) {
            fault.setBuriedEdgesLabelName(_data->faultEdgeLabels[i]);
            fault.setBuriedEdgesLabelValue(1);
        } // if
        if (!_data->failureExpected) {
            fault.transformTopology(_mesh);
        } else {
            REQUIRE_THROWS_AS(fault.transformTopology(_mesh), std::runtime_error);
            return;
        } // if/else
    } // for
    pythia::journal::debug_t debug(GenericComponent::getName());
    if (debug.state()) {
        _mesh->view(":mesh_transformed.txt:ascii_info_detail");
        _mesh->view(":mesh_transformed.tex:ascii_latex");
        _mesh->view("vtk:mesh_transformed.vtu");
    } // if

    REQUIRE(_data->cellDim == size_t(_mesh->getDimension()));
    PetscDM dmMesh = _mesh->getDM();assert(dmMesh);

    // Check vertices
    topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
    const PetscInt vStart = verticesStratum.begin();
    const PetscInt vEnd = verticesStratum.end();
    REQUIRE(_data->numVertices == size_t(verticesStratum.size()));

    topology::CoordsVisitor coordsVisitor(dmMesh);
    const PetscInt spaceDim = _data->spaceDim;

    for (PetscInt v = vStart; v < vEnd; ++v) {
        CHECK(spaceDim == coordsVisitor.sectionDof(v));
    } // for

    // check cells
    assert(_data->numCorners);
    PetscErrorCode err = 0;
    pylith::topology::Stratum cellsStratum(dmMesh, pylith::topology::Stratum::HEIGHT, 0);
    const PetscInt cStart = cellsStratum.begin();
    const PetscInt cEnd = cellsStratum.end();
    REQUIRE(_data->numCells == size_t(cellsStratum.size()));
    for (PetscInt c = cStart, cell = 0; c < cEnd; ++c, ++cell) {
        PetscInt coneSize = 0;
        err = DMPlexGetConeSize(dmMesh, c, &coneSize);PYLITH_CHECK_ERROR(err);
        INFO("cell="<<cell);
        CHECK(_data->numCorners[cell] == coneSize);
    } // for

    // check materials
    assert(_data->getMatId);
    PetscDMLabel labelMaterials = NULL;
    const char* const cellsLabelName = pylith::topology::Mesh::cells_label_name;
    err = DMGetLabel(dmMesh, cellsLabelName, &labelMaterials);PYLITH_CHECK_ERROR(err);
    assert(labelMaterials);
    const PetscInt idDefault = -999;
    for (PetscInt c = cStart, cell = 0; c < cEnd; ++c, ++cell) {
        PetscInt labelValue;

        err = DMLabelGetValue(labelMaterials, c, &labelValue);PYLITH_CHECK_ERROR(err);
        if (labelValue == -1) {
            labelValue = idDefault;
        } // if
        double centroid[3];
        err = DMPlexComputeCellGeometryFVM(dmMesh, c, NULL, centroid, NULL);PYLITH_CHECK_ERROR(err);
        const PetscInt labelValueE = _data->getMatId(c, _data->numNoncohesiveCells, centroid);
        INFO("cell="<<c<<" with centroid ("<<centroid[0]<<","<<centroid[1]<<","<<centroid[2]<<")");
        CHECK(labelValueE == labelValue);
    } // for

    // Check groups
    assert(_data->groupSizes);
    assert(_data->groupNames);
    assert(_data->groupTypes);
    PetscInt numLabels;
    err = DMGetNumLabels(dmMesh, &numLabels);PYLITH_CHECK_ERROR(err);
    for (PetscInt iLabel = 0; iLabel < numLabels; ++iLabel) {
        PetscDMLabel label = NULL;
        PetscIS pointIS = NULL;
        const PetscInt *points = NULL;
        PetscInt numPoints = 0, depth = 0;
        const char *labelName = NULL;
        std::set<std::string> ignoreLabels;
        ignoreLabels.insert("depth");
        ignoreLabels.insert(pylith::topology::Mesh::cells_label_name);
        ignoreLabels.insert("vtk");
        ignoreLabels.insert("ghost");
        ignoreLabels.insert("dim");
        ignoreLabels.insert("celltype");

        err = DMGetLabelName(dmMesh, iLabel, &labelName);PYLITH_CHECK_ERROR(err);
        if (ignoreLabels.count(labelName) > 0) { continue; }
        INFO("Checking label '" <<labelName << "'");
        err = DMGetLabel(dmMesh, labelName, &label);PYLITH_CHECK_ERROR(err);assert(label);
        err = DMLabelGetStratumIS(label, 1, &pointIS);PYLITH_CHECK_ERROR(err);
        REQUIRE(pointIS);
        err = ISGetLocalSize(pointIS, &numPoints);PYLITH_CHECK_ERROR(err);
        err = ISGetIndices(pointIS, &points);PYLITH_CHECK_ERROR(err);
        err = DMGetLabelValue(dmMesh, "depth", points[0], &depth);PYLITH_CHECK_ERROR(err);
        err = ISRestoreIndices(pointIS, &points);PYLITH_CHECK_ERROR(err);
        err = ISDestroy(&pointIS);PYLITH_CHECK_ERROR(err);
        std::string groupType = depth ? "face" : "vertex";

        bool foundGroup = false;
        for (size_t i = 0; i < _data->numGroups; ++i) {
            if (std::string(_data->groupNames[i]) == std::string(labelName)) {
                INFO("group name="<<labelName);
                CHECK(std::string(_data->groupTypes[i]) == groupType);
                CHECK(_data->groupSizes[i] == numPoints);
                foundGroup = true;
                break;
            } // if
        } // for
        INFO("Could not find group '" << labelName << "' in test data.");
        CHECK(foundGroup);
    } // for

} // run_transform


// ------------------------------------------------------------------------------------------------
void
pylith::faults::TestAdjustTopology::_initialize(void) {
    PYLITH_METHOD_BEGIN;
    REQUIRE(_data);

    delete _mesh;_mesh = new pylith::topology::Mesh;REQUIRE(_mesh);

    pylith::meshio::MeshIOAscii iohandler;
    iohandler.setFilename(_data->filename);
    iohandler.read(_mesh);
    REQUIRE(pylith::topology::MeshOps::getNumCells(*_mesh) > 0);
    REQUIRE(pylith::topology::MeshOps::getNumVertices(*_mesh) > 0);

    pylith::topology::Stratum cellsStratum(_mesh->getDM(), pylith::topology::Stratum::HEIGHT, 0);
    _data->numNoncohesiveCells = cellsStratum.size(); // Number of cells in original mesh (no cohesive cells).

    PYLITH_METHOD_END;
} // _initialize


// ------------------------------------------------------------------------------------------------
// Constructor
pylith::faults::TestAdjustTopology_Data::TestAdjustTopology_Data(void) :
    filename(NULL),
    numFaults(0),
    faultSurfaceLabels(NULL),
    faultEdgeLabels(NULL),
    interfaceIds(NULL),
    cellDim(0),
    spaceDim(0),
    numVertices(0),
    numCells(0),
    numNoncohesiveCells(0),
    numCorners(NULL),
    getMatId(NULL),
    numGroups(0),
    groupSizes(NULL),
    groupNames(NULL),
    groupTypes(NULL),
    failureExpected(false) {}


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::faults::TestAdjustTopology_Data::~TestAdjustTopology_Data(void) {}


// End of file
