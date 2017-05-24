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

#include "TestDataWriterHDF5Mesh.hh" // Implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/meshio/DataWriterHDF5.hh" // USES DataWriterHDF5

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5Mesh::setUp(void)
{ // setUp
    PYLITH_METHOD_BEGIN;

    TestDataWriterMesh::setUp();
    _data = NULL;

    PYLITH_METHOD_END;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::meshio::TestDataWriterHDF5Mesh::tearDown(void)
{ // tearDown
    PYLITH_METHOD_BEGIN;

    TestDataWriterMesh::tearDown();
    delete _data; _data = NULL;

    PYLITH_METHOD_END;
} // tearDown

// ----------------------------------------------------------------------
// Test constructor
void
pylith::meshio::TestDataWriterHDF5Mesh::testConstructor(void)
{ // testConstructor
    PYLITH_METHOD_BEGIN;

    DataWriterHDF5 writer;

    CPPUNIT_ASSERT(!writer._viewer);

    PYLITH_METHOD_END;
} // testConstructor

// ----------------------------------------------------------------------
// Test filename()
void
pylith::meshio::TestDataWriterHDF5Mesh::testFilename(void)
{ // testDebug
    PYLITH_METHOD_BEGIN;

    DataWriterHDF5 writer;

    const char* filename = "data.h5";
    writer.filename(filename);
    CPPUNIT_ASSERT_EQUAL(std::string(filename), writer._filename);

    PYLITH_METHOD_END;
} // testFilename

// ----------------------------------------------------------------------
// Test open() and close()
void
pylith::meshio::TestDataWriterHDF5Mesh::testOpenClose(void)
{ // testOpenClose
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_data);

    DataWriterHDF5 writer;

    writer.filename(_data->timestepFilename);

    const bool isInfo = false;
    writer.open(*_mesh, isInfo);
    writer.close();

    checkFile(_data->timestepFilename);

    PYLITH_METHOD_END;
} // testOpenClose

// ----------------------------------------------------------------------
// Test writeVertexField.
void
pylith::meshio::TestDataWriterHDF5Mesh::testWriteVertexField(void)
{ // testWriteVertexField
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_data);

    DataWriterHDF5 writer;

    pylith::topology::Fields vertexFields(*_mesh);
    _createVertexFields(&vertexFields);

    writer.filename(_data->vertexFilename);

    const PylithScalar timeScale = 4.0;
    writer.timeScale(timeScale);
    const PylithScalar t = _data->time / timeScale;

    const bool isInfo = false;
    writer.open(*_mesh, isInfo);
    writer.openTimeStep(t, *_mesh);

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
pylith::meshio::TestDataWriterHDF5Mesh::testWriteCellField(void)
{ // testWriteCellField
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_data);

    DataWriterHDF5 writer;

    pylith::topology::Fields cellFields(*_mesh);
    _createCellFields(&cellFields);

    writer.filename(_data->cellFilename);

    const PylithScalar timeScale = 4.0;
    writer.timeScale(timeScale);
    const PylithScalar t = _data->time / timeScale;

    const bool isInfo = false;
    writer.open(*_mesh, isInfo);
    writer.openTimeStep(t, *_mesh);

    const int numFields = 4;
    const char* fieldNames[4] = {"scalar", "vector", "tensor", "other"};
    for (int i = 0; i < numFields; ++i) {
        pylith::topology::Field& field = cellFields.get(fieldNames[i]);
        writer.writeCellField(t, field);
    } // for

    writer.closeTimeStep();
    writer.close();

    checkFile(_data->cellFilename);

    PYLITH_METHOD_END;
} // testWriteCellField

// ----------------------------------------------------------------------
// Test _hdf5Filename.
void
pylith::meshio::TestDataWriterHDF5Mesh::testHdf5Filename(void)
{ // testHdf5Filename
    PYLITH_METHOD_BEGIN;

    DataWriterHDF5 writer;

    // Append info to filename if info.
    writer._isInfo = true;
    writer._filename = "output.h5";
    CPPUNIT_ASSERT_EQUAL(std::string("output_info.h5"), writer._hdf5Filename());

    writer._isInfo = false;
    writer._filename = "output_abc.h5";
    CPPUNIT_ASSERT_EQUAL(std::string("output_abc.h5"),
                         writer._hdf5Filename());

    writer._isInfo = false;
    writer._filename = "output_abcd.h5";
    CPPUNIT_ASSERT_EQUAL(std::string("output_abcd.h5"),
                         writer._hdf5Filename());

    PYLITH_METHOD_END;
} // testHdf5Filename


// ----------------------------------------------------------------------
// Get test data.
pylith::meshio::TestDataWriter_Data*
pylith::meshio::TestDataWriterHDF5Mesh::_getData(void)
{ // _getData
    return _data;
} // _getData


// End of file
