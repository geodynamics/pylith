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

#include "TestDataWriterVTKMesh.hh" // Implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/meshio/DataWriterVTK.hh" // USES DataWriterVTK
#include "pylith/meshio/OutputSubfield.hh" // USES OutputSubfield
#include "pylith/utils/error.hh" // USES PYLITH_METHOD*

#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers_floating_point.hpp"

// ------------------------------------------------------------------------------------------------
// Setup testing data.
pylith::meshio::TestDataWriterVTKMesh::TestDataWriterVTKMesh(TestDataWriterVTKMesh_Data* data) :
    _data(data) {
    TestDataWriterMesh::_initialize();
}


// ------------------------------------------------------------------------------------------------
// Tear down testing data.
pylith::meshio::TestDataWriterVTKMesh::~TestDataWriterVTKMesh(void) {
    PYLITH_METHOD_BEGIN;

    delete _data;_data = nullptr;

    PYLITH_METHOD_END;
} // tearDown


// ------------------------------------------------------------------------------------------------
// Test accessors.
void
pylith::meshio::TestDataWriterVTKMesh::testAccessors(void) {
    const double tolerance = 1.0e-6;

    DataWriterVTK writer;

    const std::string& filename = "data.vtk";
    writer.filename(filename.c_str());
    CHECK(filename == writer._filename);

    const char* format = "%4.1f";
    writer.timeFormat(format);
    CHECK(std::string(format) == writer._timeFormat);

    const PylithScalar value = 4.5;
    writer.timeConstant(value);
    CHECK_THAT(writer._timeConstant, Catch::Matchers::WithinAbs(value, tolerance));

    // Verify error with negative time constant.
    CHECK_THROWS_AS(writer.timeConstant(-1.0), std::runtime_error);

    const int ivalue = 4;
    writer.precision(ivalue);
    CHECK(ivalue == writer._precision);

    // Verify error with nonpositive precision.
    CHECK_THROWS_AS(writer.precision(0), std::runtime_error);
    CHECK_THROWS_AS(writer.precision(-1), std::runtime_error);
} // testFilename


// ------------------------------------------------------------------------------------------------
// Test _vtkFilename.
void
pylith::meshio::TestDataWriterVTKMesh::testVtkFilename(void) {
    PYLITH_METHOD_BEGIN;

    DataWriterVTK writer;

    writer._isInfo = true;
    writer._filename = "output.vtk";
    CHECK(std::string("output_info.vtk") == writer._vtkFilename(0.0));

    // Use default normalization of 1.0, remove period from time stamp.
    writer._isInfo = false;
    writer._filename = "output.vtk";
    writer.timeFormat("%05.2f");
    CHECK(std::string("output_t0230.vtk") == writer._vtkFilename(2.3));

    // Use normalization of 20.0, remove period from time stamp.
    writer._isInfo = false;
    writer._filename = "output.vtk";
    writer.timeFormat("%05.2f");
    writer.timeConstant(20.0);
    CHECK(std::string("output_t0250.vtk") == writer._vtkFilename(50.0));

    PYLITH_METHOD_END;
} // testVtkFilename


// ------------------------------------------------------------------------------------------------
// Test openTimeStep() and closeTimeStep()
void
pylith::meshio::TestDataWriterVTKMesh::testTimeStep(void) {
    PYLITH_METHOD_BEGIN;

    assert(_mesh);
    assert(_data);

    DataWriterVTK writer;

    writer.filename(_data->timestepFilename);
    writer.timeFormat(_data->timeFormat);

    CHECK(false == writer._wroteVertexHeader);
    CHECK(false == writer._wroteCellHeader);

    const PylithScalar t = _data->time;
    const bool isInfo = false;
    writer.open(*_mesh, isInfo);
    writer.openTimeStep(t, *_mesh);

    CHECK(false == writer._wroteVertexHeader);
    CHECK(false == writer._wroteCellHeader);

    writer.closeTimeStep();
    writer.close();

    CHECK(false == writer._wroteVertexHeader);
    CHECK(false == writer._wroteCellHeader);

    // Nothing to check. We do not create VTK files without fields anymore.

    PYLITH_METHOD_END;
} // testTimeStep


// ------------------------------------------------------------------------------------------------
// Test writeVertexField.
void
pylith::meshio::TestDataWriterVTKMesh::testWriteVertexField(void) {
    PYLITH_METHOD_BEGIN;
    assert(_mesh);
    assert(_data);

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
        assert(subfield);
        subfield->project(vertexField.getOutputVector());
        writer.writeVertexField(t, *subfield);
        assert(writer._wroteVertexHeader);
        CHECK(false == writer._wroteCellHeader);
        delete subfield;subfield = nullptr;
    } // for
    writer.closeTimeStep();
    writer.close();
    CHECK(false == writer._wroteVertexHeader);
    CHECK(false == writer._wroteCellHeader);

    TestDataWriterVTK::checkFile(_data->vertexFilename, t, _data->timeFormat);

    PYLITH_METHOD_END;
} // testWriteVertexField


// ------------------------------------------------------------------------------------------------
// Test writeCellField.
void
pylith::meshio::TestDataWriterVTKMesh::testWriteCellField(void) {
    PYLITH_METHOD_BEGIN;

    assert(_mesh);
    assert(_data);

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
        assert(subfield);
        subfield->project(cellField.getOutputVector());
        writer.writeCellField(t, *subfield);
        CHECK(false == writer._wroteVertexHeader);
        assert(writer._wroteCellHeader);
        delete subfield;subfield = nullptr;
    } // for
    writer.closeTimeStep();
    writer.close();
    CHECK(false == writer._wroteVertexHeader);
    CHECK(false == writer._wroteCellHeader);

    TestDataWriterVTK::checkFile(_data->cellFilename, t, _data->timeFormat);

    PYLITH_METHOD_END;
} // testWriteCellField


// ------------------------------------------------------------------------------------------------
// Get test data.
pylith::meshio::TestDataWriter_Data*
pylith::meshio::TestDataWriterVTKMesh::_getData(void) {
    return _data;
} // _getData


// End of file
