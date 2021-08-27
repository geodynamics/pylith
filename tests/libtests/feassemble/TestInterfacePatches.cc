// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestInterfacePatches.hh" // Implementation of class methods

#include "pylith/feassemble/InterfacePatches.hh" // Test subject

#include "pylith/testing/FaultCohesiveStub.hh" // USES FaultCohesiveStub

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

// ------------------------------------------------------------------------------------------------
// Setup testing data.
void
pylith::feassemble::TestInterfacePatches::setUp(void) {
    PYLITH_METHOD_BEGIN;

    _data = new TestInterfacePatches_Data;CPPUNIT_ASSERT(_data);
    _mesh = NULL;
    _fault = NULL;

    PYLITH_METHOD_END;
} // setUp


// ------------------------------------------------------------------------------------------------
// Tear down testing data.
void
pylith::feassemble::TestInterfacePatches::tearDown(void) {
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

    const std::string& defaultName = pylith::topology::Mesh::getCellsLabelName();
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Mismatch in default label name.", defaultName, patches._labelName);

    const std::string& name = "fault patches";
    patches._labelName = name;
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Mismatch in label name.", name, std::string(patches.getLabelName()));

    PYLITH_METHOD_END;
} // testAccessors


// ------------------------------------------------------------------------------------------------
// Test createMaterialPairs().
void
pylith::feassemble::TestInterfacePatches::testCreateMaterialPairs(void) {
    PYLITH_METHOD_BEGIN;
    CPPUNIT_ASSERT(_data);

    _initialize();
    CPPUNIT_ASSERT(_fault);
    InterfacePatches* patches = InterfacePatches::createMaterialPairs(_fault, _mesh->getDM());
    CPPUNIT_ASSERT(patches);

    const std::string& labelName = std::string(_data->faultLabel) + std::string("-integration-patches");
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Mismatch in label name.", labelName, std::string(patches->getLabelName()));

    const size_t numPatches = _data->numPatches;
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Mismatch in number of integration patches",
                                 numPatches, patches->_keys.size());

    const std::string& fieldName = "lagrange_multiplier_fault";
    const std::string& cellsLabelName = pylith::topology::Mesh::getCellsLabelName();
    CPPUNIT_ASSERT(_data->patchKeys);
    CPPUNIT_ASSERT(_data->patchNumCells);
    CPPUNIT_ASSERT(_data->patchCells);
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
        CPPUNIT_ASSERT_MESSAGE("Could not find patch.", patchIndex >= 0);

        CPPUNIT_ASSERT_EQUAL(labelName, weakFormKeys.cohesive._name);
        CPPUNIT_ASSERT_EQUAL(fieldName, weakFormKeys.cohesive._field);

        CPPUNIT_ASSERT_EQUAL(cellsLabelName, weakFormKeys.negative._name);
        CPPUNIT_ASSERT_EQUAL(_data->patchKeys[patchIndex].negative_value, matIdNegative);
        CPPUNIT_ASSERT_EQUAL(fieldName, weakFormKeys.negative._field);

        CPPUNIT_ASSERT_EQUAL(cellsLabelName, weakFormKeys.positive._name);
        CPPUNIT_ASSERT_EQUAL(_data->patchKeys[patchIndex].positive_value, matIdPositive);
        CPPUNIT_ASSERT_EQUAL(fieldName, weakFormKeys.positive._field);

        // Check labels
        PetscErrorCode err = 0;
        PetscDMLabel label = NULL;
        PetscIS pointsIS = NULL;
        PetscInt numPoints = 0;
        const PetscInt* points = NULL;
        err = DMGetLabel(_mesh->getDM(), labelName.c_str(), &label);CPPUNIT_ASSERT(!err);
        err = DMLabelGetStratumIS(label, labelValue, &pointsIS);CPPUNIT_ASSERT(!err);
        err = ISGetSize(pointsIS, &numPoints);CPPUNIT_ASSERT(!err);
        std::ostringstream msg;
        msg << "Mismatch in number of cells for integration patch for materials ("<<matIdNegative<<", "<<matIdPositive<<").";
        CPPUNIT_ASSERT_EQUAL_MESSAGE(msg.str().c_str(), _data->patchNumCells[patchIndex], numPoints);

        err = ISGetIndices(pointsIS, &points);CPPUNIT_ASSERT(!err);
        CPPUNIT_ASSERT(_data->patchCells[patchIndex]);
        for (PetscInt iPoint = 0; iPoint < numPoints; ++iPoint) {
            std::ostringstream msg;
            msg << "Mismatch in cells in integration patch for materials ("<<matIdNegative<<", "<<matIdPositive<<").";
            CPPUNIT_ASSERT_EQUAL_MESSAGE(msg.str().c_str(), _data->patchCells[patchIndex][iPoint], points[iPoint]);
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
    CPPUNIT_ASSERT(_data);

    delete _mesh;_mesh = new pylith::topology::Mesh;CPPUNIT_ASSERT(_mesh);

    pylith::meshio::MeshIOAscii iohandler;
    iohandler.filename(_data->filename);
    iohandler.read(_mesh);
    CPPUNIT_ASSERT(pylith::topology::MeshOps::getNumCells(*_mesh) > 0);
    CPPUNIT_ASSERT(pylith::topology::MeshOps::getNumVertices(*_mesh) > 0);

    CPPUNIT_ASSERT(_data->faultLabel);
    _fault = new pylith::faults::FaultCohesiveStub();CPPUNIT_ASSERT(_fault);
    _fault->setInterfaceId(101);
    _fault->setSurfaceMarkerLabel(_data->faultLabel);
    if (_data->edgeLabel) {
        _fault->setBuriedEdgesMarkerLabel(_data->edgeLabel);
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
