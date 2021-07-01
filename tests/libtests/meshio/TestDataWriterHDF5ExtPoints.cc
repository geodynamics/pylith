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

#include "TestDataWriterHDF5ExtPoints.hh" // Implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/meshio/OutputSolnPoints.hh" // USES OutputSolnPoints
#include "pylith/meshio/DataWriterHDF5Ext.hh" // USES DataWriterHDF5Ext
#include "pylith/meshio/OutputSubfield.hh" // USES OutputSubfield
#include "pylith/utils/error.hh" // USES PYLITH_METHOD*

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

// ------------------------------------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5ExtPoints::setUp(void) {
    PYLITH_METHOD_BEGIN;

    TestDataWriterPoints::setUp();
    _data = NULL;

    PYLITH_METHOD_END;
} // setUp


// ------------------------------------------------------------------------------------------------
// Tear down testing data.
void
pylith::meshio::TestDataWriterHDF5ExtPoints::tearDown(void) {
    PYLITH_METHOD_BEGIN;

    TestDataWriterPoints::tearDown();
    delete _data;_data = NULL;

    PYLITH_METHOD_END;
} // tearDown


// ------------------------------------------------------------------------------------------------
// Test openTimeStep() and closeTimeStep()
void
pylith::meshio::TestDataWriterHDF5ExtPoints::testTimeStep(void) {
    PYLITH_METHOD_BEGIN;
    CPPUNIT_ASSERT(_pointMesh);
    CPPUNIT_ASSERT(_data);

    DataWriterHDF5Ext writer;
    writer.filename(_data->opencloseFilename);

    const PylithScalar t = _data->time;
    const bool isInfo = false;
    writer.open(*_pointMesh, isInfo);
    writer.writePointNames(_data->names, *_pointMesh);
    writer.openTimeStep(t, *_pointMesh);
    writer.closeTimeStep();
    writer.close();

    checkFile(_data->opencloseFilename);

    PYLITH_METHOD_END;
} // testTimeStep


// ------------------------------------------------------------------------------------------------
// Test writeVertexField.
void
pylith::meshio::TestDataWriterHDF5ExtPoints::testWriteVertexField(void) {
    PYLITH_METHOD_BEGIN;
    CPPUNIT_ASSERT(_pointMesh);
    CPPUNIT_ASSERT(_data);

    pylith::topology::Field vertexField(*_pointMesh);
    _createVertexField(&vertexField);

    DataWriterHDF5Ext writer;
    writer.filename(_data->vertexFilename);

    const PylithScalar t = _data->time;
    const bool isInfo = false;
    writer.open(*_pointMesh, isInfo);
    writer.writePointNames(_data->names, *_pointMesh);
    writer.openTimeStep(t, *_pointMesh);

    const pylith::string_vector& subfieldNames = vertexField.getSubfieldNames();
    const size_t numFields = subfieldNames.size();
    for (size_t i = 0; i < numFields; ++i) {
        OutputSubfield* subfield = OutputSubfield::create(vertexField, *_pointMesh, subfieldNames[i].c_str());
        CPPUNIT_ASSERT(subfield);

        const pylith::topology::Field::SubfieldInfo& info = vertexField.getSubfieldInfo(subfieldNames[i].c_str());
        subfield->extractSubfield(vertexField, info.index);

        writer.writeVertexField(t, *subfield);
        delete subfield;subfield = NULL;
    } // for
    writer.closeTimeStep();
    writer.close();

    checkFile(_data->vertexFilename);

    PYLITH_METHOD_END;
} // testWriteVertexField


// ------------------------------------------------------------------------------------------------
// Get test data.
pylith::meshio::TestDataWriterPoints_Data*
pylith::meshio::TestDataWriterHDF5ExtPoints::_getData(void) {
    return _data;
} // _getData


// End of file
