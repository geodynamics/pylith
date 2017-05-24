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

#include "TestDataWriterHDF5ExtMaterial.hh" // Implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/meshio/DataWriterHDF5Ext.hh" // USES DataWriterHDF5Ext

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5ExtMaterial::setUp(void)
{ // setUp
    PYLITH_METHOD_BEGIN;

    TestDataWriterMaterial::setUp();
    _data = NULL;

    PYLITH_METHOD_END;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::meshio::TestDataWriterHDF5ExtMaterial::tearDown(void)
{ // tearDown
    PYLITH_METHOD_BEGIN;

    TestDataWriterMaterial::tearDown();
    delete _data; _data = NULL;

    PYLITH_METHOD_END;
} // tearDown

// ----------------------------------------------------------------------
// Test open() and close()
void
pylith::meshio::TestDataWriterHDF5ExtMaterial::testOpenClose(void)
{ // testOpenClose
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_data);

    DataWriterHDF5Ext writer;

    writer.filename(_data->opencloseFilename);

    const bool isInfo = false;
    writer.open(*_mesh, isInfo, "material-id", _data->materialId);
    writer.close();

    checkFile(_data->opencloseFilename);

    PYLITH_METHOD_END;
} // testOpenClose

// ----------------------------------------------------------------------
// Test writeVertexField.
void
pylith::meshio::TestDataWriterHDF5ExtMaterial::testWriteVertexField(void)
{ // testWriteVertexField
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_data);

    DataWriterHDF5Ext writer;

    pylith::topology::Fields vertexFields(*_mesh);
    _createVertexFields(&vertexFields);

    writer.filename(_data->vertexFilename);

    const PylithScalar timeScale = 4.0;
    writer.timeScale(timeScale);
    const PylithScalar t = _data->time / timeScale;

    const bool isInfo = false;
    writer.open(*_mesh, isInfo, "material-id", _data->materialId);
    writer.openTimeStep(t, *_mesh, "material-id", _data->materialId);

    const int numFields = 4;
    const char* fieldNames[4] = {"scalar", "vector", "tensor", "other"};
    for (int i = 0; i < numFields; ++i) {
        pylith::topology::Field& field = vertexFields.get(fieldNames[i]);
        writer.writeVertexField(t, field, *_mesh);
    } // for

    writer.closeTimeStep();
    writer.close();

    checkFile(_data->vertexFilename);

    PYLITH_METHOD_END;
} // testWriteVertexField

// ----------------------------------------------------------------------
// Test writeCellField.
void
pylith::meshio::TestDataWriterHDF5ExtMaterial::testWriteCellField(void)
{ // testWriteCellField
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_data);

    DataWriterHDF5Ext writer;

    pylith::topology::Fields cellFields(*_mesh);
    _createCellFields(&cellFields);

    writer.filename(_data->cellFilename);

    const PylithScalar timeScale = 4.0;
    writer.timeScale(timeScale);
    const PylithScalar t = _data->time / timeScale;

    const bool isInfo = false;
    writer.open(*_mesh, isInfo, "material-id", _data->materialId);
    writer.openTimeStep(t, *_mesh, "material-id", _data->materialId);

    const int numFields = 4;
    const char* fieldNames[4] = {"scalar", "vector", "tensor", "other"};
    for (int i = 0; i < numFields; ++i) {
        pylith::topology::Field& field = cellFields.get(fieldNames[i]);
        writer.writeCellField(t, field, "material-id", _data->materialId);
    } // for

    writer.closeTimeStep();
    writer.close();

    checkFile(_data->cellFilename);

    PYLITH_METHOD_END;
} // testWriteCellField

// ----------------------------------------------------------------------
// Get test data.
pylith::meshio::TestDataWriterMaterial_Data*
pylith::meshio::TestDataWriterHDF5ExtMaterial::_getData(void)
{ // _getData
    return _data;
} // _getData


// End of file
