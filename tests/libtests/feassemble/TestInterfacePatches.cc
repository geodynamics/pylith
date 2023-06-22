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

#include "TestInterfacePatches.hh" // Implementation of class methods

#include "pylith/feassemble/InterfacePatches.hh" // Test subject
#include "tests/src/FaultCohesiveStub.hh" // USES FaultCohesiveStub

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

#include "catch2/catch_test_macros.hpp"

// ------------------------------------------------------------------------------------------------
// Constructor.
pylith::feassemble::TestInterfacePatches::TestInterfacePatches(TestInterfacePatches_Data* data) :
    _data(data) {
    PYLITH_METHOD_BEGIN;

    assert(_data);
    _mesh = NULL;
    _fault = NULL;

    PYLITH_METHOD_END;
} // setUp


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::feassemble::TestInterfacePatches::~TestInterfacePatches(void) {
    PYLITH_METHOD_BEGIN;

    delete _data;_data = NULL;
    delete _mesh;_mesh = NULL;
    delete _fault;_fault = NULL;

    PYLITH_METHOD_END;
} // tearDown


// ------------------------------------------------------------------------------------------------
// Test getLabelName().
void
pylith::feassemble::TestInterfacePatches::testAccessors(void) {
    PYLITH_METHOD_BEGIN;

    InterfacePatches patches;

    const std::string& defaultName = pylith::topology::Mesh::cells_label_name;
    CHECK(defaultName == patches._labelName);

    const std::string& name = "fault patches";
    patches._labelName = name;
    CHECK(name == std::string(patches.getLabelName()));

    PYLITH_METHOD_END;
} // testAccessors


// ------------------------------------------------------------------------------------------------
// Test createMaterialPairs().
void
pylith::feassemble::TestInterfacePatches::testCreateMaterialPairs(void) {
    PYLITH_METHOD_BEGIN;
    assert(_data);

    _initialize();
    assert(_fault);
    InterfacePatches* patches = InterfacePatches::createMaterialPairs(_fault, _mesh->getDM());
    assert(patches);

    const std::string& labelName = std::string(_data->faultLabel) + std::string("-integration-patches");
    CHECK(labelName == std::string(patches->getLabelName()));

    const size_t numPatches = _data->numPatches;
    REQUIRE(numPatches == patches->_keys.size());

    const std::string& cellsLabelName = pylith::topology::Mesh::cells_label_name;
    assert(_data->patchKeys);
    assert(_data->patchNumCells);
    assert(_data->patchCells);
    for (InterfacePatches::keysmap_t::iterator iter = patches->_keys.begin(); iter != patches->_keys.end(); ++iter) {
        const PetscInt labelValue = iter->first;
        const InterfacePatches::WeakFormKeys weakFormKeys = iter->second;

        const int matIdNegative = weakFormKeys.negative._value;
        const int matIdPositive = weakFormKeys.positive._value;
        size_t patchIndex = -1;
        for (size_t i = 0; i < _data->numPatches; ++i) {
            if ((matIdNegative == _data->patchKeys[i].negative_value) &&
                (matIdPositive == _data->patchKeys[i].positive_value) ) {
                patchIndex = i;
                break;
            } // if
        } // for
        assert(patchIndex >= 0);

        CHECK(labelName == weakFormKeys.cohesive._name);

        CHECK(cellsLabelName == weakFormKeys.negative._name);
        CHECK(_data->patchKeys[patchIndex].negative_value == matIdNegative);

        CHECK(cellsLabelName == weakFormKeys.positive._name);
        CHECK(_data->patchKeys[patchIndex].positive_value == matIdPositive);

        // Check labels
        PetscErrorCode err = 0;
        PetscDMLabel label = NULL;
        PetscIS pointsIS = NULL;
        PetscInt numPoints = 0;
        const PetscInt* points = NULL;
        err = DMGetLabel(_mesh->getDM(), labelName.c_str(), &label);assert(!err);
        err = DMLabelGetStratumIS(label, labelValue, &pointsIS);assert(!err);
        err = ISGetSize(pointsIS, &numPoints);assert(!err);
        INFO("Checking integration patch for materials ("<<matIdNegative<<", "<<matIdPositive<<").");
        REQUIRE(_data->patchNumCells[patchIndex] == numPoints);

        err = ISGetIndices(pointsIS, &points);assert(!err);
        assert(_data->patchCells[patchIndex]);
        for (PetscInt iPoint = 0; iPoint < numPoints; ++iPoint) {
            CHECK(_data->patchCells[patchIndex][iPoint] == points[iPoint]);
        } // for
        err = ISRestoreIndices(pointsIS, &points);PYLITH_CHECK_ERROR(err);
        err = ISDestroy(&pointsIS);PYLITH_CHECK_ERROR(err);
    } // for
    delete patches;patches = NULL;

    PYLITH_METHOD_END;
} // testCreateSinglees


// ------------------------------------------------------------------------------------------------
void
pylith::feassemble::TestInterfacePatches::_initialize() {
    PYLITH_METHOD_BEGIN;
    assert(_data);

    delete _mesh;_mesh = new pylith::topology::Mesh;assert(_mesh);

    pylith::meshio::MeshIOAscii iohandler;
    iohandler.setFilename(_data->filename);
    iohandler.read(_mesh);
    assert(pylith::topology::MeshOps::getNumCells(*_mesh) > 0);
    assert(pylith::topology::MeshOps::getNumVertices(*_mesh) > 0);

    assert(_data->faultLabel);
    _fault = new pylith::faults::FaultCohesiveStub();assert(_fault);
    _fault->setCohesiveLabelName(pylith::topology::Mesh::cells_label_name);
    _fault->setCohesiveLabelValue(101);
    _fault->setSurfaceLabelName(_data->faultLabel);
    _fault->setSurfaceLabelValue(1);
    if (_data->edgeLabel) {
        _fault->setBuriedEdgesLabelName(_data->edgeLabel);
        _fault->setBuriedEdgesLabelValue(1);
    } // if
    _fault->adjustTopology(_mesh);

    PYLITH_METHOD_END;
} // _initialize


// ------------------------------------------------------------------------------------------------
// Constructor
pylith::feassemble::TestInterfacePatches_Data::TestInterfacePatches_Data(void) :
    filename(NULL),
    faultLabel(NULL),
    edgeLabel(NULL),
    numPatches(0),
    patchKeys(NULL),
    patchNumCells(NULL),
    patchCells(NULL) {}


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::feassemble::TestInterfacePatches_Data::~TestInterfacePatches_Data(void) {}


// End of file
