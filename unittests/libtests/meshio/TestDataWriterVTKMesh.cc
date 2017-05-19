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

#include "TestDataWriterVTKMesh.hh" // Implementation of class methods

#include "data/DataWriterData.hh" // USES DataWriterData

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/meshio/DataWriterVTK.hh" // USES DataWriterVTK
#include "pylith/faults/FaultCohesiveKin.hh" // USES FaultCohesiveKin

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterVTKMesh::setUp(void)
{ // setUp
    PYLITH_METHOD_BEGIN;

    TestDataWriterMesh::setUp();
    _data = NULL;

    PYLITH_METHOD_END;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::meshio::TestDataWriterVTKMesh::tearDown(void)
{ // tearDown
    PYLITH_METHOD_BEGIN;

    TestDataWriterMesh::tearDown();
    delete _data; _data = NULL;

    PYLITH_METHOD_END;
} // tearDown

// ----------------------------------------------------------------------
// Test constructor
void
pylith::meshio::TestDataWriterVTKMesh::testConstructor(void)
{ // testConstructor
    PYLITH_METHOD_BEGIN;

    DataWriterVTK writer;

    CPPUNIT_ASSERT(!writer._viewer);
    CPPUNIT_ASSERT_EQUAL(false, writer._wroteVertexHeader);
    CPPUNIT_ASSERT_EQUAL(false, writer._wroteCellHeader);

    PYLITH_METHOD_END;
} // testConstructor

// ----------------------------------------------------------------------
// Test filename()
void
pylith::meshio::TestDataWriterVTKMesh::testFilename(void)
{ // testDebug
    PYLITH_METHOD_BEGIN;

    DataWriterVTK writer;

    const char* filename = "data.vtk";
    writer.filename(filename);
    CPPUNIT_ASSERT_EQUAL(std::string(filename), writer._filename);

    PYLITH_METHOD_END;
} // testFilename

// ----------------------------------------------------------------------
// Test timeFormat()
void
pylith::meshio::TestDataWriterVTKMesh::testTimeFormat(void)
{ // testTimeFormat
    PYLITH_METHOD_BEGIN;

    DataWriterVTK writer;

    const char* format = "%4.1f";
    writer.timeFormat(format);
    CPPUNIT_ASSERT_EQUAL(std::string(format), writer._timeFormat);

    PYLITH_METHOD_END;
} // testTimeFormat

// ----------------------------------------------------------------------
// Test timeConstant()
void
pylith::meshio::TestDataWriterVTKMesh::testTimeConstant(void)
{ // testTimeConstant
    PYLITH_METHOD_BEGIN;

    DataWriterVTK writer;

    const PylithScalar value = 4.5;
    writer.timeConstant(value);
    CPPUNIT_ASSERT_EQUAL(value, writer._timeConstant);

    // Verify error with negative time constant.
    CPPUNIT_ASSERT_THROW(writer.timeConstant(-1.0), std::runtime_error);

    PYLITH_METHOD_END;
} // testTimeConstant

// ----------------------------------------------------------------------
// Test precision()
void
pylith::meshio::TestDataWriterVTKMesh::testPrecision(void)
{ // testPrecision
    PYLITH_METHOD_BEGIN;

    DataWriterVTK writer;

    const int value = 4;
    writer.precision(value);
    CPPUNIT_ASSERT_EQUAL(value, writer._precision);

    // Verify error with nonpositive precision.
    CPPUNIT_ASSERT_THROW(writer.precision(0), std::runtime_error);
    CPPUNIT_ASSERT_THROW(writer.precision(-1), std::runtime_error);

    PYLITH_METHOD_END;
} // testPrecision

// ----------------------------------------------------------------------
// Test openTimeStep() and closeTimeStep()
void
pylith::meshio::TestDataWriterVTKMesh::testTimeStep(void)
{ // testTimeStep
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_data);

    DataWriterVTK writer;

    writer.filename(_data->timestepFilename);
    writer.timeFormat(_data->timeFormat);

    CPPUNIT_ASSERT_EQUAL(false, writer._wroteVertexHeader);
    CPPUNIT_ASSERT_EQUAL(false, writer._wroteCellHeader);

    const PylithScalar t = _data->time;
    const bool isInfo = false;
    writer.open(*_mesh, isInfo);
    writer.openTimeStep(t, *_mesh);

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
pylith::meshio::TestDataWriterVTKMesh::testWriteVertexField(void)
{ // testWriteVertexField
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_data);

    DataWriterVTK writer;

    pylith::topology::Fields vertexFields(*_mesh);
    _createVertexFields(&vertexFields);

    writer.filename(_data->vertexFilename);
    writer.timeFormat(_data->timeFormat);

    const PylithScalar t = _data->time;
    const bool isInfo = false;
    writer.open(*_mesh, isInfo);
    writer.openTimeStep(t, *_mesh);

    const int numFields = 4;
    const char* fieldNames[4] = {"scalar", "vector", "tensor", "other"};
    for (int i = 0; i < numFields; ++i) {
        pylith::topology::Field& field = vertexFields.get(fieldNames[i]);
        writer.writeVertexField(t, field, *_mesh);
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
pylith::meshio::TestDataWriterVTKMesh::testWriteCellField(void)
{ // testWriteCellField
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_data);

    DataWriterVTK writer;

    pylith::topology::Fields cellFields(*_mesh);
    _createCellFields(&cellFields);

    writer.filename(_data->cellFilename);
    writer.timeFormat(_data->timeFormat);

    const PylithScalar t = _data->time;
    const bool isInfo = false;
    writer.open(*_mesh, isInfo);
    writer.openTimeStep(t, *_mesh);

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
// Test _vtkFilename.
void
pylith::meshio::TestDataWriterVTKMesh::testVtkFilename(void)
{ // testVtkFilename
    PYLITH_METHOD_BEGIN;

    DataWriterVTK writer;

    writer._isInfo = true;
    writer._filename = "output.vtk";
    CPPUNIT_ASSERT_EQUAL(std::string("output_info.vtk"), writer._vtkFilename(0.0));

    // Use default normalization of 1.0, remove period from time stamp.
    writer._isInfo = false;
    writer._filename = "output.vtk";
    writer.timeFormat("%05.2f");
    CPPUNIT_ASSERT_EQUAL(std::string("output_t0230.vtk"), writer._vtkFilename(2.3));

    // Use normalization of 20.0, remove period from time stamp.
    writer._isInfo = false;
    writer._filename = "output.vtk";
    writer.timeFormat("%05.2f");
    writer.timeConstant(20.0);
    CPPUNIT_ASSERT_EQUAL(std::string("output_t0250.vtk"), writer._vtkFilename(50.0));

    PYLITH_METHOD_END;
} // testVtkFilename


// ----------------------------------------------------------------------
// Get test data.
pylith::meshio::TestDataWriter_Data*
pylith::meshio::TestDataWriterVTKMesh::_getData(void)
{ // _getData
    return _data;
} // _getData


// End of file
