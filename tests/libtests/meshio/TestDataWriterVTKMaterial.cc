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

#include <portinfo>

#include "TestDataWriterVTKMaterial.hh" // Implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/meshio/DataWriterVTK.hh" // USES DataWriterVTK
#include "pylith/meshio/OutputSubfield.hh" // USES OutputSubfield
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*

// ------------------------------------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterVTKMaterial::setUp(void) {
    PYLITH_METHOD_BEGIN;

    TestDataWriterMaterial::setUp();
    _data = NULL;

    PYLITH_METHOD_END;
} // setUp


// ------------------------------------------------------------------------------------------------
// Tear down testing data.
void
pylith::meshio::TestDataWriterVTKMaterial::tearDown(void) {
    PYLITH_METHOD_BEGIN;

    TestDataWriterMaterial::tearDown();
    delete _data;_data = NULL;

    PYLITH_METHOD_END;
} // tearDown


// ------------------------------------------------------------------------------------------------
// Test openTimeStep() and closeTimeStep()
void
pylith::meshio::TestDataWriterVTKMaterial::testTimeStep(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_materialMesh);
    CPPUNIT_ASSERT(_data);

    DataWriterVTK writer;

    writer.filename(_data->timestepFilename);
    writer.timeFormat(_data->timeFormat);

    CPPUNIT_ASSERT_EQUAL(false, writer._wroteVertexHeader);
    CPPUNIT_ASSERT_EQUAL(false, writer._wroteCellHeader);

    const PylithScalar t = _data->time;
    const bool isInfo = false;
    writer.open(*_materialMesh, isInfo);
    writer.openTimeStep(t, *_materialMesh);

    CPPUNIT_ASSERT_EQUAL(false, writer._wroteVertexHeader);
    CPPUNIT_ASSERT_EQUAL(false, writer._wroteCellHeader);

    writer.closeTimeStep();
    writer.close();

    CPPUNIT_ASSERT_EQUAL(false, writer._wroteVertexHeader);
    CPPUNIT_ASSERT_EQUAL(false, writer._wroteCellHeader);

    // Nothing to check. We do not create VTK files without fields anymore.

    PYLITH_METHOD_END;
} // testTimeStep


// ------------------------------------------------------------------------------------------------
// Test writeVertexField.
void
pylith::meshio::TestDataWriterVTKMaterial::testWriteVertexField(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_domainMesh);
    CPPUNIT_ASSERT(_materialMesh);
    CPPUNIT_ASSERT(_data);

    pylith::topology::Field vertexField(*_domainMesh);
    _createVertexField(&vertexField);

    DataWriterVTK writer;
    writer.filename(_data->vertexFilename);
    writer.timeFormat(_data->timeFormat);

    const PylithScalar t = _data->time;
    const bool isInfo = false;
    writer.open(*_materialMesh, isInfo);
    writer.openTimeStep(t, *_materialMesh);

    const pylith::string_vector& subfieldNames = vertexField.getSubfieldNames();
    const size_t numFields = subfieldNames.size();
    for (size_t i = 0; i < numFields; ++i) {
        OutputSubfield* subfield = OutputSubfield::create(vertexField, *_materialMesh, subfieldNames[i].c_str(), 1);
        CPPUNIT_ASSERT(subfield);
        subfield->project(vertexField.getOutputVector());
        writer.writeVertexField(t, *subfield);
        CPPUNIT_ASSERT(writer._wroteVertexHeader);
        CPPUNIT_ASSERT_EQUAL(false, writer._wroteCellHeader);
        delete subfield;subfield = NULL;
    } // for
    writer.closeTimeStep();
    writer.close();

    CPPUNIT_ASSERT_EQUAL(false, writer._wroteVertexHeader);
    CPPUNIT_ASSERT_EQUAL(false, writer._wroteCellHeader);

    checkFile(_data->vertexFilename, t, _data->timeFormat);

    PYLITH_METHOD_END;
} // testWriteVertexField


// ------------------------------------------------------------------------------------------------
// Test writeCellField.
void
pylith::meshio::TestDataWriterVTKMaterial::testWriteCellField(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_materialMesh);
    CPPUNIT_ASSERT(_data);

    pylith::topology::Field cellField(*_materialMesh);
    _createCellField(&cellField);

    DataWriterVTK writer;
    writer.filename(_data->cellFilename);
    writer.timeFormat(_data->timeFormat);

    const PylithScalar t = _data->time;
    const bool isInfo = false;
    writer.open(*_materialMesh, isInfo);
    writer.openTimeStep(t, *_materialMesh);

    const pylith::string_vector& subfieldNames = cellField.getSubfieldNames();
    const size_t numFields = subfieldNames.size();
    for (size_t i = 0; i < numFields; ++i) {
        OutputSubfield* subfield = OutputSubfield::create(cellField, *_materialMesh, subfieldNames[i].c_str(), 0);
        CPPUNIT_ASSERT(subfield);
        subfield->project(cellField.getOutputVector());
        writer.writeCellField(t, *subfield);
        CPPUNIT_ASSERT_EQUAL(false, writer._wroteVertexHeader);
        CPPUNIT_ASSERT(writer._wroteCellHeader);
        delete subfield;subfield = NULL;
    } // for
    writer.closeTimeStep();
    writer.close();
    CPPUNIT_ASSERT_EQUAL(false, writer._wroteVertexHeader);
    CPPUNIT_ASSERT_EQUAL(false, writer._wroteCellHeader);

    checkFile(_data->cellFilename, t, _data->timeFormat);

    PYLITH_METHOD_END;
} // testWriteCellField


// ------------------------------------------------------------------------------------------------
// Get test data.
pylith::meshio::TestDataWriterMaterial_Data*
pylith::meshio::TestDataWriterVTKMaterial::_getData(void) {
    return _data;
} // _getData


// End of file
