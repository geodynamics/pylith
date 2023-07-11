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

#include "TestDataWriterHDF5ExtMesh.hh" // Implementation of class methods

#include "pylith/utils/types.hh" // HASA PylithScalar

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/meshio/DataWriterHDF5Ext.hh" // USES DataWriterHDF5Ext
#include "pylith/meshio/OutputSubfield.hh" // USES OutputSubfield
#include "pylith/utils/error.hh" // USES PYLITH_METHOD*

#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers_floating_point.hpp"

// ------------------------------------------------------------------------------------------------
// Constructor.
pylith::meshio::TestDataWriterHDF5ExtMesh::TestDataWriterHDF5ExtMesh(TestDataWriterHDF5ExtMesh_Data* data) :
    _data(data) {
    TestDataWriterMesh::_initialize();
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::meshio::TestDataWriterHDF5ExtMesh::~TestDataWriterHDF5ExtMesh(void) {
    PYLITH_METHOD_BEGIN;

    delete _data;_data = nullptr;

    PYLITH_METHOD_END;
} // destructor


// ------------------------------------------------------------------------------------------------
// Test filename()
void
pylith::meshio::TestDataWriterHDF5ExtMesh::testAccessors(void) {
    PYLITH_METHOD_BEGIN;

    DataWriterHDF5Ext writer;

    const char* filename = "data.h5";
    writer.filename(filename);
    CHECK(std::string(filename) == writer._filename);

    PYLITH_METHOD_END;
} // testAccessors


// ------------------------------------------------------------------------------------------------
// Test hdf5Filename().
void
pylith::meshio::TestDataWriterHDF5ExtMesh::testHdf5Filename(void) {
    PYLITH_METHOD_BEGIN;

    DataWriterHDF5Ext writer;

    // Append info to filename if number of time steps is 0.
    writer._isInfo = true;
    writer._filename = "output.h5";
    CHECK(std::string("output_info.h5") == writer.hdf5Filename());

    writer._isInfo = false;
    writer._filename = "output_abc.h5";
    CHECK(std::string("output_abc.h5") == writer.hdf5Filename());

    writer._isInfo = false;
    writer._filename = "output_abcd.h5";
    CHECK(std::string("output_abcd.h5") == writer.hdf5Filename());

    PYLITH_METHOD_END;
} // testHdf5Filename


// ------------------------------------------------------------------------------------------------
// Test _datasetFilename().
void
pylith::meshio::TestDataWriterHDF5ExtMesh::testDatasetFilename(void) {
    PYLITH_METHOD_BEGIN;

    DataWriterHDF5Ext writer;

    // Append info to filename if info.
    writer._isInfo = true;
    writer._filename = "output.h5";
    CHECK(std::string("output_info_ABCD.dat") == writer._datasetFilename("ABCD"));

    writer._isInfo = false;
    writer._filename = "output_abc.h5";
    CHECK(std::string("output_abc_field1.dat") == writer._datasetFilename("field1"));

    writer._isInfo = false;
    writer._filename = "output_abcd.h5";
    CHECK(std::string("output_abcd_field2.dat") == writer._datasetFilename("field2"));

    PYLITH_METHOD_END;
} // testDatasetFilename


// ------------------------------------------------------------------------------------------------
// Test open() and close()
void
pylith::meshio::TestDataWriterHDF5ExtMesh::testOpenClose(void) {
    PYLITH_METHOD_BEGIN;
    assert(_mesh);
    assert(_data);

    DataWriterHDF5Ext writer;

    writer.filename(_data->opencloseFilename);

    const bool isInfo = false;
    writer.open(*_mesh, isInfo);
    writer.close();

    checkFile(_data->opencloseFilename);

    PYLITH_METHOD_END;
} // testOpenClose


// ------------------------------------------------------------------------------------------------
// Test writeVertexField.
void
pylith::meshio::TestDataWriterHDF5ExtMesh::testWriteVertexField(void) {
    PYLITH_METHOD_BEGIN;
    assert(_mesh);
    assert(_data);

    DataWriterHDF5Ext writer;

    topology::Field vertexField(*_mesh);
    _createVertexField(&vertexField);

    writer.filename(_data->vertexFilename);

    const PylithScalar timeScale = 4.0;
    writer.setTimeScale(timeScale);
    const PylithScalar t = _data->time / timeScale;

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
        delete subfield;subfield = NULL;
    } // for
    writer.closeTimeStep();
    writer.close();

    checkFile(_data->vertexFilename);

    PYLITH_METHOD_END;
} // testWriteVertexField


// ------------------------------------------------------------------------------------------------
// Test writeCellField.
void
pylith::meshio::TestDataWriterHDF5ExtMesh::testWriteCellField(void) {
    PYLITH_METHOD_BEGIN;
    assert(_mesh);
    assert(_data);

    DataWriterHDF5Ext writer;

    topology::Field cellField(*_mesh);
    _createCellField(&cellField);

    writer.filename(_data->cellFilename);

    const PylithScalar timeScale = 4.0;
    writer.setTimeScale(timeScale);
    const PylithScalar t = _data->time / timeScale;

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
        delete subfield;subfield = NULL;
    } // for
    writer.closeTimeStep();
    writer.close();

    checkFile(_data->cellFilename);

    PYLITH_METHOD_END;
} // testWriteCellField


// ------------------------------------------------------------------------------------------------
// Get test data.
pylith::meshio::TestDataWriter_Data*
pylith::meshio::TestDataWriterHDF5ExtMesh::_getData(void) {
    return _data;
} // _getData


// End of file
