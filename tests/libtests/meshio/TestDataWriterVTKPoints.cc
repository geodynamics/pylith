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

#include "TestDataWriterVTKPoints.hh" // Implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/meshio/OutputSolnPoints.hh" // USES OutputSolnPoints
#include "pylith/meshio/DataWriterVTK.hh" // USES DataWriterVTK
#include "pylith/meshio/OutputSubfield.hh" // USES OutputSubfield
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers_floating_point.hpp"

#include <cassert> // USES assert()

// ------------------------------------------------------------------------------------------------
// Constructor.
pylith::meshio::TestDataWriterVTKPoints::TestDataWriterVTKPoints(TestDataWriterVTKPoints_Data* data) :
    _data(data) {
    TestDataWriterPoints::_initialize();
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::meshio::TestDataWriterVTKPoints::~TestDataWriterVTKPoints(void) {
    PYLITH_METHOD_BEGIN;

    delete _data;_data = nullptr;

    PYLITH_METHOD_END;
} // destructor.


// ------------------------------------------------------------------------------------------------
// Test openTimeStep() and closeTimeStep()
void
pylith::meshio::TestDataWriterVTKPoints::testTimeStep(void) {
    PYLITH_METHOD_BEGIN;
    assert(_pointMesh);
    assert(_data);

    DataWriterVTK writer;
    writer.filename(_data->timestepFilename);
    writer.timeFormat(_data->timeFormat);

    CHECK(false == writer._wroteVertexHeader);
    CHECK(false == writer._wroteCellHeader);

    const PylithScalar t = _data->time;
    const bool isInfo = false;
    writer.open(*_pointMesh, isInfo);
    writer.openTimeStep(t, *_pointMesh);

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
pylith::meshio::TestDataWriterVTKPoints::testWriteVertexField(void) {
    PYLITH_METHOD_BEGIN;
    assert(_pointMesh);
    assert(_data);

    pylith::topology::Field vertexField(*_pointMesh);
    _createVertexField(&vertexField);

    DataWriterVTK writer;
    writer.filename(_data->vertexFilename);
    writer.timeFormat(_data->timeFormat);

    const PylithScalar t = _data->time;
    const bool isInfo = false;
    writer.open(*_pointMesh, isInfo);
    writer.openTimeStep(t, *_pointMesh);

    const pylith::string_vector& subfieldNames = vertexField.getSubfieldNames();
    const size_t numFields = subfieldNames.size();
    for (size_t i = 0; i < numFields; ++i) {
        OutputSubfield* subfield = OutputSubfield::create(vertexField, *_pointMesh, subfieldNames[i].c_str());
        assert(subfield);

        const pylith::topology::Field::SubfieldInfo& info = vertexField.getSubfieldInfo(subfieldNames[i].c_str());
        subfield->extractSubfield(vertexField, info.index);

        writer.writeVertexField(t, *subfield);
        assert(writer._wroteVertexHeader);
        CHECK(false == writer._wroteCellHeader);
        delete subfield;subfield = NULL;
    } // for
    writer.closeTimeStep();
    writer.close();

    CHECK(false == writer._wroteVertexHeader);
    CHECK(false == writer._wroteCellHeader);

    TestDataWriterVTK::checkFile(_data->vertexFilename, t, _data->timeFormat);

    PYLITH_METHOD_END;
} // testWriteVertexField


// ------------------------------------------------------------------------------------------------
// Get test data.
pylith::meshio::TestDataWriterPoints_Data*
pylith::meshio::TestDataWriterVTKPoints::_getData(void) {
    return _data;
} // _getData


// End of file
