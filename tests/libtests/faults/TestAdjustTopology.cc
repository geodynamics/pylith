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

#include "TestAdjustTopology.hh"

#include "pylith/testing/FaultCohesiveStub.hh" // USES FaultsCohesiveStub

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps::nondimensionalize()
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor

#include "pylith/utils/array.hh" // USES int_array, scalar_array
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart

// ------------------------------------------------------------------------------------------------
// Setup testing _data->
void
pylith::faults::TestAdjustTopology::setUp(void) {
    PYLITH_METHOD_BEGIN;

    _data = new TestAdjustTopology_Data;CPPUNIT_ASSERT(_data);
    _mesh = NULL;

    PYLITH_METHOD_END;
} // setUp


// ------------------------------------------------------------------------------------------------
// Tear down testing _data->
void
pylith::faults::TestAdjustTopology::tearDown(void) {
    PYLITH_METHOD_BEGIN;

    delete _data;_data = NULL;
    delete _mesh;_mesh = NULL;

    PYLITH_METHOD_END;
} // tearDown


// ------------------------------------------------------------------------------------------------
// Test reorder().
void
pylith::faults::TestAdjustTopology::testAdjustTopology(void) {
    PYLITH_METHOD_BEGIN;

    _initialize();
    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_data);

    for (size_t i = 0; i < _data->numFaults; ++i) {
        FaultCohesiveStub fault;

        CPPUNIT_ASSERT(_data->interfaceIds);
        CPPUNIT_ASSERT(_data->faultSurfaceLabels);
        CPPUNIT_ASSERT(_data->faultEdgeLabels);

        fault.setInterfaceId(_data->interfaceIds[i]);
        fault.setSurfaceMarkerLabel(_data->faultSurfaceLabels[i]);
        if (_data->faultEdgeLabels[i]) {
            fault.setBuriedEdgesMarkerLabel(_data->faultEdgeLabels[i]);
        } // if
        if (!_data->failureExpected) {
            fault.adjustTopology(_mesh);
        } else {
            CPPUNIT_ASSERT_THROW(fault.adjustTopology(_mesh), std::runtime_error);
            return;
        } // if/else
    } // for

#if 0 // DEBUGGING
    PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_INFO_DETAIL);
    DMView(_mesh->getDM(), PETSC_VIEWER_STDOUT_WORLD);
    PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD);
#endif

    CPPUNIT_ASSERT_EQUAL(_data->cellDim, size_t(_mesh->getDimension()));
    PetscDM dmMesh = _mesh->getDM();CPPUNIT_ASSERT(dmMesh);

    // Check vertices
    topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
    const PetscInt vStart = verticesStratum.begin();
    const PetscInt vEnd = verticesStratum.end();
    CPPUNIT_ASSERT_EQUAL(_data->numVertices, size_t(verticesStratum.size()));

    topology::CoordsVisitor coordsVisitor(dmMesh);
    const PetscInt spaceDim = _data->spaceDim;

    for (PetscInt v = vStart; v < vEnd; ++v) {
        CPPUNIT_ASSERT_EQUAL(spaceDim, coordsVisitor.sectionDof(v));
    } // for

    // check cells
    CPPUNIT_ASSERT(_data->numCorners);
    PetscErrorCode err = 0;
    pylith::topology::Stratum cellsStratum(dmMesh, pylith::topology::Stratum::HEIGHT, 0);
    const PetscInt cStart = cellsStratum.begin();
    const PetscInt cEnd = cellsStratum.end();
    CPPUNIT_ASSERT_EQUAL(_data->numCells, size_t(cellsStratum.size()));
    for (PetscInt c = cStart, cell = 0; c < cEnd; ++c, ++cell) {
        PetscInt coneSize = 0;
        err = DMPlexGetConeSize(dmMesh, c, &coneSize);PYLITH_CHECK_ERROR(err);
        CPPUNIT_ASSERT_EQUAL(_data->numCorners[cell], coneSize);
    } // for

    // check materials
    CPPUNIT_ASSERT(_data->materialIds);
    PetscDMLabel labelMaterials = NULL;
    const char* const cellsLabelName = pylith::topology::Mesh::getCellsLabelName();
    err = DMGetLabel(dmMesh, cellsLabelName, &labelMaterials);PYLITH_CHECK_ERROR(err);
    CPPUNIT_ASSERT(labelMaterials);
    const PetscInt idDefault = -999;
    for (PetscInt c = cStart, cell = 0; c < cEnd; ++c, ++cell) {
        PetscInt value;

        err = DMLabelGetValue(labelMaterials, c, &value);PYLITH_CHECK_ERROR(err);
        if (value == -1) {
            value = idDefault;
        } // if
        std::ostringstream msg;
        msg << "Mismatch in '"<<cellsLabelName<<"' for cell "<<cell<<".";
        CPPUNIT_ASSERT_EQUAL_MESSAGE(msg.str().c_str(), _data->materialIds[cell], value);
    } // for

    // Check groups
    CPPUNIT_ASSERT(_data->groupSizes);
    CPPUNIT_ASSERT(_data->groupNames);
    CPPUNIT_ASSERT(_data->groupTypes);
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
        ignoreLabels.insert("material-id");
        ignoreLabels.insert("vtk");
        ignoreLabels.insert("ghost");
        ignoreLabels.insert("dim");
        ignoreLabels.insert("celltype");

        err = DMGetLabelName(dmMesh, iLabel, &labelName);PYLITH_CHECK_ERROR(err);
        if (ignoreLabels.count(labelName) > 0) { continue; }
        err = DMGetLabel(dmMesh, labelName, &label);PYLITH_CHECK_ERROR(err);CPPUNIT_ASSERT(label);
        err = DMLabelGetStratumIS(label, 1, &pointIS);PYLITH_CHECK_ERROR(err);
        err = ISGetLocalSize(pointIS, &numPoints);PYLITH_CHECK_ERROR(err);
        err = ISGetIndices(pointIS, &points);PYLITH_CHECK_ERROR(err);
        err = DMGetLabelValue(dmMesh, "depth", points[0], &depth);PYLITH_CHECK_ERROR(err);
        std::string groupType = depth ? "cell" : "vertex";

        bool foundGroup = false;
        for (size_t i = 0; i < _data->numGroups; ++i) {
            if (std::string(_data->groupNames[i]) == std::string(labelName)) {
                std::ostringstream msg;
                msg << "Mismatch for group '" <<labelName<<"'.";
                CPPUNIT_ASSERT_EQUAL_MESSAGE(msg.str().c_str(), std::string(_data->groupTypes[i]), groupType);
                CPPUNIT_ASSERT_EQUAL_MESSAGE(msg.str().c_str(), _data->groupSizes[i], numPoints);
                foundGroup = true;
                break;
            } // if
        } // for
        err = ISRestoreIndices(pointIS, &points);PYLITH_CHECK_ERROR(err);
        err = ISDestroy(&pointIS);PYLITH_CHECK_ERROR(err);
        std::ostringstream msg;
        msg << "Could not find group '" << labelName << "'.";
        CPPUNIT_ASSERT_MESSAGE(msg.str().c_str(), foundGroup);
    } // for

} // testAdjustTopology


// ------------------------------------------------------------------------------------------------
void
pylith::faults::TestAdjustTopology::_initialize(void) {
    PYLITH_METHOD_BEGIN;
    CPPUNIT_ASSERT(_data);

    delete _mesh;_mesh = new pylith::topology::Mesh;CPPUNIT_ASSERT(_mesh);

    pylith::meshio::MeshIOAscii iohandler;
    iohandler.filename(_data->filename);
    iohandler.read(_mesh);
    CPPUNIT_ASSERT(pylith::topology::MeshOps::getNumCells(*_mesh) > 0);
    CPPUNIT_ASSERT(pylith::topology::MeshOps::getNumVertices(*_mesh) > 0);

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
    numCorners(NULL),
    materialIds(NULL),
    numGroups(0),
    groupSizes(NULL),
    groupNames(NULL),
    groupTypes(NULL),
    failureExpected(false) {}


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::faults::TestAdjustTopology_Data::~TestAdjustTopology_Data(void) {}


// End of file
