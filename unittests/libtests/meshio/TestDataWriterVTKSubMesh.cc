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

#include "TestDataWriterVTKSubMesh.hh" // Implementation of class methods

#include "data/DataWriterData.hh" // USES DataWriterData

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/meshio/DataWriterVTK.hh" // USES DataWriterVTK

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterVTKSubMesh::setUp(void)
{ // setUp
    PYLITH_METHOD_BEGIN;

    TestDataWriterSubMesh::setUp();
    _data = NULL;

    PYLITH_METHOD_END;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::meshio::TestDataWriterVTKSubMesh::tearDown(void)
{ // tearDown
    PYLITH_METHOD_BEGIN;

    TestDataWriterSubMesh::tearDown();
    delete _data; _data = NULL;

    PYLITH_METHOD_END;
} // tearDown

// ----------------------------------------------------------------------
// Test openTimeStep() and closeTimeStep()
void
pylith::meshio::TestDataWriterVTKSubMesh::testTimeStep(void)
{ // testTimeStep
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_submesh);
    CPPUNIT_ASSERT(_data);

    DataWriterVTK writer;

    writer.filename(_data->timestepFilename);
    writer.timeFormat(_data->timeFormat);

    CPPUNIT_ASSERT_EQUAL(false, writer._wroteVertexHeader);
    CPPUNIT_ASSERT_EQUAL(false, writer._wroteCellHeader);

    const PylithScalar t = _data->time;
    const bool isInfo = false;
    writer.open(*_submesh, isInfo);
    writer.openTimeStep(t, *_submesh);

    CPPUNIT_ASSERT_EQUAL(false, writer._wroteVertexHeader);
    CPPUNIT_ASSERT_EQUAL(false, writer._wroteCellHeader);

    writer.closeTimeStep();
    writer.close();

    CPPUNIT_ASSERT_EQUAL(false, writer._wroteVertexHeader);
    CPPUNIT_ASSERT_EQUAL(false, writer._wroteCellHeader);

    // Nothing to check. We do not create VTK files without fields anymore.

    PYLITH_METHOD_END;
} // testTimeStep

// ----------------------------------------------------------------------
// Test writeVertexField.
void
pylith::meshio::TestDataWriterVTKSubMesh::testWriteVertexField(void)
{ // testWriteVertexField
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_submesh);
    CPPUNIT_ASSERT(_data);

    DataWriterVTK writer;

    pylith::topology::Fields vertexFields(*_mesh);
    _createVertexFields(&vertexFields);

    writer.filename(_data->vertexFilename);
    writer.timeFormat(_data->timeFormat);

    const PylithScalar t = _data->time;
    const bool isInfo = false;
    writer.open(*_submesh, isInfo);
    writer.openTimeStep(t, *_submesh);

    const int numFields = 4;
    const char* fieldNames[4] = {"scalar", "vector", "tensor", "other"};
    for (int i = 0; i < numFields; ++i) {
        pylith::topology::Field& field = vertexFields.get(fieldNames[i]);
        writer.writeVertexField(t, field, *_submesh);
        CPPUNIT_ASSERT(writer._wroteVertexHeader);
        CPPUNIT_ASSERT_EQUAL(false, writer._wroteCellHeader);
    } // for

    writer.closeTimeStep();
    writer.close();

    CPPUNIT_ASSERT_EQUAL(false, writer._wroteVertexHeader);
    CPPUNIT_ASSERT_EQUAL(false, writer._wroteCellHeader);

    checkFile(_data->vertexFilename, t, _data->timeFormat);

    PYLITH_METHOD_END;
} // testWriteVertexField

// ----------------------------------------------------------------------
// Test writeCellField.
void
pylith::meshio::TestDataWriterVTKSubMesh::testWriteCellField(void)
{ // testWriteCellField
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_submesh);
    CPPUNIT_ASSERT(_data);

    DataWriterVTK writer;

    pylith::topology::Fields cellFields(*_submesh);
    _createCellFields(&cellFields);

    writer.filename(_data->cellFilename);
    writer.timeFormat(_data->timeFormat);

    const PylithScalar t = _data->time;
    const bool isInfo = false;
    writer.open(*_submesh, isInfo);
    writer.openTimeStep(t, *_submesh);

    const int numFields = 4;
    const char* fieldNames[4] = {"scalar", "vector", "tensor", "other"};
    for (int i = 0; i < numFields; ++i) {
        pylith::topology::Field& field = cellFields.get(fieldNames[i]);
        writer.writeCellField(t, field);
        CPPUNIT_ASSERT_EQUAL(false, writer._wroteVertexHeader);
        CPPUNIT_ASSERT(writer._wroteCellHeader);
    } // for

    writer.closeTimeStep();
    writer.close();

    CPPUNIT_ASSERT_EQUAL(false, writer._wroteCellHeader);
    CPPUNIT_ASSERT_EQUAL(false, writer._wroteCellHeader);

    checkFile(_data->cellFilename, t, _data->timeFormat);

    PYLITH_METHOD_END;
} // testWriteCellField


// ----------------------------------------------------------------------
// Get test data.
pylith::meshio::TestDataWriterSubMesh_Data*
pylith::meshio::TestDataWriterVTKSubMesh::_getData(void)
{ // _getData
    return _data;
} // _getData


// End of file
