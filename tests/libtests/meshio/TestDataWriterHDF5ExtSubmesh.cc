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

#include "TestDataWriterHDF5ExtSubmesh.hh" // Implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/meshio/DataWriterHDF5Ext.hh" // USES DataWriterHDF5Ext
#include "pylith/meshio/OutputSubfield.hh" // USES OutputSubfield
#include "pylith/utils/error.hh" // USES PYLITH_METHOD*

// ------------------------------------------------------------------------------------------------
// Constructor.
pylith::meshio::TestDataWriterHDF5ExtSubmesh::TestDataWriterHDF5ExtSubmesh(TestDataWriterHDF5ExtSubmesh_Data* data) :
    _data(data) {
    TestDataWriterSubmesh::_initialize();
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::meshio::TestDataWriterHDF5ExtSubmesh::~TestDataWriterHDF5ExtSubmesh(void) {
    PYLITH_METHOD_BEGIN;

    delete _data;_data = nullptr;

    PYLITH_METHOD_END;
} // destructor


// ------------------------------------------------------------------------------------------------
// Test openTimeStep() and closeTimeStep()
void
pylith::meshio::TestDataWriterHDF5ExtSubmesh::testOpenClose(void) {
    PYLITH_METHOD_BEGIN;
    assert(_submesh);
    assert(_data);

    DataWriterHDF5Ext writer;

    writer.filename(_data->opencloseFilename);

    const bool isInfo = false;
    writer.open(*_submesh, isInfo);
    writer.close();

    checkFile(_data->opencloseFilename);

    PYLITH_METHOD_END;
} // testOpenClose


// ------------------------------------------------------------------------------------------------
// Test writeVertexField.
void
pylith::meshio::TestDataWriterHDF5ExtSubmesh::testWriteVertexField(void) {
    PYLITH_METHOD_BEGIN;
    assert(_mesh);
    assert(_submesh);
    assert(_data);

    DataWriterHDF5Ext writer;

    pylith::topology::Field vertexField(*_mesh);
    _createVertexField(&vertexField);

    writer.filename(_data->vertexFilename);

    const PylithScalar timeScale = 4.0;
    writer.setTimeScale(timeScale);
    const PylithScalar t = _data->time / timeScale;

    const bool isInfo = false;
    writer.open(*_submesh, isInfo);
    writer.openTimeStep(t, *_submesh);

    const pylith::string_vector& subfieldNames = vertexField.getSubfieldNames();
    const size_t numFields = subfieldNames.size();
    for (size_t i = 0; i < numFields; ++i) {
        OutputSubfield* subfield = OutputSubfield::create(vertexField, *_submesh, subfieldNames[i].c_str(), 1);
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
pylith::meshio::TestDataWriterHDF5ExtSubmesh::testWriteCellField(void) {
    PYLITH_METHOD_BEGIN;
    assert(_submesh);
    assert(_data);

    DataWriterHDF5Ext writer;

    pylith::topology::Field cellField(*_submesh);
    _createCellField(&cellField);

    writer.filename(_data->cellFilename);

    const PylithScalar timeScale = 4.0;
    writer.setTimeScale(timeScale);
    const PylithScalar t = _data->time / timeScale;

    const bool isInfo = false;
    writer.open(*_submesh, isInfo);
    writer.openTimeStep(t, *_submesh);

    const pylith::string_vector& subfieldNames = cellField.getSubfieldNames();
    const size_t numFields = subfieldNames.size();
    for (size_t i = 0; i < numFields; ++i) {
        OutputSubfield* subfield = OutputSubfield::create(cellField, *_submesh, subfieldNames[i].c_str(), 0);
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
pylith::meshio::TestDataWriterSubmesh_Data*
pylith::meshio::TestDataWriterHDF5ExtSubmesh::_getData(void) {
    return _data;
} // _getData


// End of file
