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

#include "TestDataWriterVTKMesh.hh" // Implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/meshio/DataWriterVTK.hh" // USES DataWriterVTK
#include "pylith/meshio/OutputSubfield.hh" // USES OutputSubfield
#include "pylith/utils/error.hh" // USES PYLITH_METHOD*

// ------------------------------------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterVTKMesh::setUp(void) {
    PYLITH_METHOD_BEGIN;

    TestDataWriterMesh::setUp();
    _data = NULL;

    PYLITH_METHOD_END;
} // setUp


// ------------------------------------------------------------------------------------------------
// Tear down testing data.
void
pylith::meshio::TestDataWriterVTKMesh::tearDown(void) {
    PYLITH_METHOD_BEGIN;

    TestDataWriterMesh::tearDown();
    delete _data;_data = NULL;

    PYLITH_METHOD_END;
} // tearDown


// ------------------------------------------------------------------------------------------------
// Test constructor
void
pylith::meshio::TestDataWriterVTKMesh::testConstructor(void) {
    PYLITH_METHOD_BEGIN;

    DataWriterVTK writer;

    CPPUNIT_ASSERT(!writer._viewer);
    CPPUNIT_ASSERT_EQUAL(false, writer._wroteVertexHeader);
    CPPUNIT_ASSERT_EQUAL(false, writer._wroteCellHeader);

    PYLITH_METHOD_END;
} // testConstructor


// ------------------------------------------------------------------------------------------------
// Test filename()
void
pylith::meshio::TestDataWriterVTKMesh::testFilename(void) {
    PYLITH_METHOD_BEGIN;

    DataWriterVTK writer;

    const char* filename = "data.vtk";
    writer.filename(filename);
    CPPUNIT_ASSERT_EQUAL(std::string(filename), writer._filename);

    PYLITH_METHOD_END;
} // testFilename


// ------------------------------------------------------------------------------------------------
// Test timeFormat()
void
pylith::meshio::TestDataWriterVTKMesh::testTimeFormat(void) {
    PYLITH_METHOD_BEGIN;

    DataWriterVTK writer;

    const char* format = "%4.1f";
    writer.timeFormat(format);
    CPPUNIT_ASSERT_EQUAL(std::string(format), writer._timeFormat);

    PYLITH_METHOD_END;
} // testTimeFormat


// ------------------------------------------------------------------------------------------------
// Test timeConstant()
void
pylith::meshio::TestDataWriterVTKMesh::testTimeConstant(void) {
    PYLITH_METHOD_BEGIN;

    DataWriterVTK writer;

    const PylithScalar value = 4.5;
    writer.timeConstant(value);
    CPPUNIT_ASSERT_EQUAL(value, writer._timeConstant);

    // Verify error with negative time constant.
    CPPUNIT_ASSERT_THROW(writer.timeConstant(-1.0), std::runtime_error);

    PYLITH_METHOD_END;
} // testTimeConstant


// ------------------------------------------------------------------------------------------------
// Test precision()
void
pylith::meshio::TestDataWriterVTKMesh::testPrecision(void) {
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


// ------------------------------------------------------------------------------------------------
// Test openTimeStep() and closeTimeStep()
void
pylith::meshio::TestDataWriterVTKMesh::testTimeStep(void) {
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


// ------------------------------------------------------------------------------------------------
// Test writeVertexField.
void
pylith::meshio::TestDataWriterVTKMesh::testWriteVertexField(void) {
    PYLITH_METHOD_BEGIN;
    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_data);

    pylith::topology::Field vertexField(*_mesh);
    _createVertexField(&vertexField);

    DataWriterVTK writer;
    writer.filename(_data->vertexFilename);
    writer.timeFormat(_data->timeFormat);

    const PylithScalar t = _data->time;
    const bool isInfo = false;
    writer.open(*_mesh, isInfo);
    writer.openTimeStep(t, *_mesh);

    const pylith::string_vector& subfieldNames = vertexField.getSubfieldNames();
    const size_t numFields = subfieldNames.size();
    for (size_t i = 0; i < numFields; ++i) {
        OutputSubfield* subfield = OutputSubfield::create(vertexField, *_mesh, subfieldNames[i].c_str(), 1);
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
pylith::meshio::TestDataWriterVTKMesh::testWriteCellField(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_data);

    pylith::topology::Field cellField(*_mesh);
    _createCellField(&cellField);

    DataWriterVTK writer;
    writer.filename(_data->cellFilename);
    writer.timeFormat(_data->timeFormat);

    const PylithScalar t = _data->time;
    const bool isInfo = false;
    writer.open(*_mesh, isInfo);
    writer.openTimeStep(t, *_mesh);

    const pylith::string_vector& subfieldNames = cellField.getSubfieldNames();
    const size_t numFields = subfieldNames.size();
    for (size_t i = 0; i < numFields; ++i) {
        OutputSubfield* subfield = OutputSubfield::create(cellField, *_mesh, subfieldNames[i].c_str(), 0);
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
// Test _vtkFilename.
void
pylith::meshio::TestDataWriterVTKMesh::testVtkFilename(void) {
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


// ------------------------------------------------------------------------------------------------
// Get test data.
pylith::meshio::TestDataWriter_Data*
pylith::meshio::TestDataWriterVTKMesh::_getData(void) {
    return _data;
} // _getData


// End of file
