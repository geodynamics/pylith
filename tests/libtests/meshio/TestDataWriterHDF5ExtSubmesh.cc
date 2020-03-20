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

#include "TestDataWriterHDF5ExtSubmesh.hh" // Implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/meshio/DataWriterHDF5Ext.hh" // USES DataWriterHDF5Ext

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5ExtSubmesh::setUp(void) {
    PYLITH_METHOD_BEGIN;

    TestDataWriterSubmesh::setUp();
    _data = NULL;

    PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::meshio::TestDataWriterHDF5ExtSubmesh::tearDown(void) {
    PYLITH_METHOD_BEGIN;

    TestDataWriterSubmesh::tearDown();
    delete _data;_data = NULL;

    PYLITH_METHOD_END;
} // tearDown


// ----------------------------------------------------------------------
// Test openTimeStep() and closeTimeStep()
void
pylith::meshio::TestDataWriterHDF5ExtSubmesh::testOpenClose(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_submesh);
    CPPUNIT_ASSERT(_data);

    DataWriterHDF5Ext writer;

    writer.filename(_data->opencloseFilename);

    const bool isInfo = false;
    writer.open(*_submesh, isInfo);
    writer.close();

    checkFile(_data->opencloseFilename);

    PYLITH_METHOD_END;
} // testOpenClose


// ----------------------------------------------------------------------
// Test writeVertexField.
void
pylith::meshio::TestDataWriterHDF5ExtSubmesh::testWriteVertexField(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_submesh);
    CPPUNIT_ASSERT(_data);

    DataWriterHDF5Ext writer;

    topology::Fields vertexFields(*_mesh);
    _createVertexFields(&vertexFields);

    writer.filename(_data->vertexFilename);

    const PylithScalar timeScale = 4.0;
    writer.setTimeScale(timeScale);
    const PylithScalar t = _data->time / timeScale;

    const bool isInfo = false;
    writer.open(*_submesh, isInfo);
    writer.openTimeStep(t, *_submesh);

    const int numFields = 4;
    const char* fieldNames[4] = {"scalar", "vector", "tensor", "other"};
    for (int i = 0; i < numFields; ++i) {
        topology::Field& field = vertexFields.get(fieldNames[i]);
        writer.writeVertexField(t, field, *_submesh);
    } // for
    writer.closeTimeStep();
    writer.close();

    checkFile(_data->vertexFilename);

    PYLITH_METHOD_END;
} // testWriteVertexField


// ----------------------------------------------------------------------
// Test writeCellField.
void
pylith::meshio::TestDataWriterHDF5ExtSubmesh::testWriteCellField(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_submesh);
    CPPUNIT_ASSERT(_data);

    DataWriterHDF5Ext writer;

    topology::Fields cellFields(*_submesh);
    _createCellFields(&cellFields);

    writer.filename(_data->cellFilename);

    const PylithScalar timeScale = 4.0;
    writer.setTimeScale(timeScale);
    const PylithScalar t = _data->time / timeScale;

    const bool isInfo = false;
    writer.open(*_submesh, isInfo);
    writer.openTimeStep(t, *_submesh);

    const int numFields = 4;
    const char* fieldNames[4] = {"scalar", "vector", "tensor", "other"};
    for (int i = 0; i < numFields; ++i) {
        topology::Field& field = cellFields.get(fieldNames[i]);
        writer.writeCellField(t, field);
    } // for
    writer.closeTimeStep();
    writer.close();

    checkFile(_data->cellFilename);

    PYLITH_METHOD_END;
} // testWriteCellField


// ----------------------------------------------------------------------
// Get test data.
pylith::meshio::TestDataWriterSubmesh_Data*
pylith::meshio::TestDataWriterHDF5ExtSubmesh::_getData(void) {
    return _data;
} // _getData


// End of file
