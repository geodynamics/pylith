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

#include "TestDataWriterSubmesh.hh" // Implementation of class methods

#include "FieldFactory.hh" // USES FieldFactory

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps::nondimensionalize()
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/meshio/DataWriter.hh" // USES DataWriter
#include "pylith/testing/FaultCohesiveStub.hh" // USES FaultCohesiveStub

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <cppunit/extensions/HelperMacros.h>

// ------------------------------------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterSubmesh::setUp(void) {
    PYLITH_METHOD_BEGIN;

    _mesh = NULL;
    _submesh = NULL;

    PYLITH_METHOD_END;
} // setUp


// ------------------------------------------------------------------------------------------------
// Tear down testing data.
void
pylith::meshio::TestDataWriterSubmesh::tearDown(void) {
    PYLITH_METHOD_BEGIN;

    delete _mesh;_mesh = NULL;
    delete _submesh;_submesh = NULL;

    PYLITH_METHOD_END;
} // tearDown


// ------------------------------------------------------------------------------------------------
// Initialize mesh.
void
pylith::meshio::TestDataWriterSubmesh::_initialize(void) {
    PYLITH_METHOD_BEGIN;

    const TestDataWriterSubmesh_Data* data = _getData();CPPUNIT_ASSERT(data);

    delete _mesh;_mesh = new pylith::topology::Mesh;CPPUNIT_ASSERT(_mesh);
    MeshIOAscii iohandler;
    iohandler.filename(data->meshFilename);
    iohandler.read(_mesh);

    spatialdata::geocoords::CSCart cs;
    cs.setSpaceDim(_mesh->getDimension());
    _mesh->setCoordSys(&cs);

    spatialdata::units::Nondimensional normalizer;
    normalizer.setLengthScale(data->lengthScale);
    pylith::topology::MeshOps::nondimensionalize(_mesh, normalizer);

    if (data->faultLabel) {
        pylith::faults::FaultCohesiveStub fault;
        fault.setSurfaceMarkerLabel(data->faultLabel);
        fault.setInterfaceId(data->faultId);
        fault.adjustTopology(_mesh);
    } // if

    CPPUNIT_ASSERT(data->bcLabel);
    delete _submesh;_submesh = pylith::topology::MeshOps::createLowerDimMesh(*_mesh, data->bcLabel);CPPUNIT_ASSERT(_submesh);

    PYLITH_METHOD_END;
} // _initialize


// ------------------------------------------------------------------------------------------------
// Create vertex fields.
void
pylith::meshio::TestDataWriterSubmesh::_createVertexField(pylith::topology::Field* field) {
    PYLITH_METHOD_BEGIN;
    CPPUNIT_ASSERT(field);
    CPPUNIT_ASSERT(_mesh);

    const TestDataWriterSubmesh_Data* data = _getData();CPPUNIT_ASSERT(data);

    FieldFactory factory(*field);
    factory.addScalar(data->vertexDiscretization);
    factory.addVector(data->vertexDiscretization);
    factory.addTensor(data->vertexDiscretization);
    factory.addOther(data->vertexDiscretization);

    field->subfieldsSetup();
    field->createDiscretization();
    field->allocate();

    factory.setValues(data->vertexValues, data->vertexNumPoints, data->vertexNumDOF);

    field->createOutputVector();
    field->scatterLocalToOutput();

    PYLITH_METHOD_END;
} // _createVertexFields


// ------------------------------------------------------------------------------------------------
// Create cell fields.
void
pylith::meshio::TestDataWriterSubmesh::_createCellField(pylith::topology::Field* field) {
    PYLITH_METHOD_BEGIN;
    CPPUNIT_ASSERT(field);

    const TestDataWriterSubmesh_Data* data = _getData();CPPUNIT_ASSERT(data);

    FieldFactory factory(*field);
    factory.addScalar(data->cellDiscretization);
    factory.addVector(data->cellDiscretization);
    factory.addTensor(data->cellDiscretization);
    factory.addOther(data->cellDiscretization);

    field->subfieldsSetup();
    field->createDiscretization();
    field->allocate();

    factory.setValues(data->cellValues, data->cellNumPoints, data->cellNumDOF);

    field->createOutputVector();
    field->scatterLocalToOutput();

    PYLITH_METHOD_END;
} // _createCellFields


// ================================================================================================
void
pylith::meshio::TestDataWriterSubmesh::_setDataTri(void) {
    TestDataWriterSubmesh_Data* data = this->_getData();CPPUNIT_ASSERT(data);

    data->meshFilename = "data/tri3.mesh";
    data->bcLabel = "bc";
    data->faultLabel = "fault";
    data->faultId = 100;
    data->spaceDim = 2;
    data->lengthScale = 10.0;

    data->time = 1.0;
    data->timeFormat = "%3.1f";

    // Vertex fields ------------------------------------------------------------------------------
    static const size_t vertexNumPoints = 6;
    static const size_t vertexNumDOF = 1 + 2 + 3 + 2;
    data->vertexDiscretization = pylith::topology::FieldBase::Discretization(1, 1, 2);

    data->vertexNumPoints = vertexNumPoints;
    data->vertexNumDOF = vertexNumDOF;
    static const PylithScalar vertexValues[vertexNumPoints*vertexNumDOF] = {
        2.1,    1.1,  2.2,   1.1, 1.2, 1.3,    1.2,  2.3,
        3.2,    3.3,  4.4,   2.1, 2.2, 2.3,    3.4,  4.5,
        4.3,    5.5,  6.6,   3.1, 3.2, 3.3,    5.6,  6.7,
        5.4,    7.7,  8.8,   4.1, 4.2, 4.3,    7.8,  8.9,
        6.5,    9.9, 10.0,   5.1, 5.2, 5.3,    9.0, 10.1,
        7.6,   11.1, 12.2,   6.1, 6.2, 6.3,   11.2, 12.3,
    };
    data->vertexValues = const_cast<PylithScalar*>(vertexValues);

    // Cell fields --------------------------------------------------------------------------------
    static const size_t cellNumPoints = 1;
    static const size_t cellNumDOF = 1 + 2 + 3 + 2;
    data->cellDiscretization = pylith::topology::FieldBase::Discretization(0, 0, 1);

    data->cellNumPoints = cellNumPoints;
    data->cellNumDOF = cellNumDOF;
    static const PylithScalar cellValues[cellNumPoints*cellNumDOF] = {
        2.1,   1.1, 2.2,   1.2, 2.3, 3.4,   1.2, 2.3,
    };
    data->cellValues = const_cast<PylithScalar*>(cellValues);
} // setDataTri


// ================================================================================================
void
pylith::meshio::TestDataWriterSubmesh::_setDataQuad(void) {
    TestDataWriterSubmesh_Data* data = this->_getData();CPPUNIT_ASSERT(data);

    // We do not use a fault in this test case.
    data->meshFilename = "data/quad4.mesh";
    data->bcLabel = "bc3";
    data->spaceDim = 2;
    data->lengthScale = 10.0;

    data->time = 1.0;
    data->timeFormat = "%3.1f";

    // Vertex fields ------------------------------------------------------------------------------
    static const size_t vertexNumPoints = 6;
    static const size_t vertexNumDOF = 1 + 2 + 3 + 2;
    data->vertexDiscretization = pylith::topology::FieldBase::Discretization(1, 1, 2);

    data->vertexNumPoints = vertexNumPoints;
    data->vertexNumDOF = vertexNumDOF;
    static const PylithScalar vertexValues[vertexNumPoints*vertexNumDOF] = {
        2.1,    1.1,  2.2,   1.1, 1.2, 1.3,    1.2,  2.3,
        3.2,    3.3,  4.4,   2.1, 2.2, 2.3,    3.4,  4.5,
        4.3,    5.5,  6.6,   3.1, 3.2, 3.3,    5.6,  6.7,
        5.4,    7.7,  8.8,   4.1, 4.2, 4.3,    7.8,  8.9,
        6.5,    9.9, 10.1,   5.1, 5.2, 5.3,    9.1, 10.2,
        7.6,   11.2, 12.3,   6.1, 6.2, 6.3,    6.5,  5.4,
    };
    data->vertexValues = const_cast<PylithScalar*>(vertexValues);

    // Cell fields --------------------------------------------------------------------------------
    static const size_t cellNumPoints = 2;
    static const size_t cellNumDOF = 1 + 2 + 3 + 2;
    data->cellDiscretization = pylith::topology::FieldBase::Discretization(0, 0, 1);

    data->cellNumPoints = cellNumPoints;
    data->cellNumDOF = cellNumDOF;
    static const PylithScalar cellValues[cellNumPoints*cellNumDOF] = {
        2.1,   1.1, 2.2,   1.2, 2.3, 3.4,   1.2, 2.3,
        3.2,   3.3, 4.4,   4.5, 5.6, 6.7,   4.5, 5.6,
    };
    data->cellValues = const_cast<PylithScalar*>(cellValues);
} // setDataQuad


// ================================================================================================
void
pylith::meshio::TestDataWriterSubmesh::_setDataTet(void) {
    TestDataWriterSubmesh_Data* data = this->_getData();CPPUNIT_ASSERT(data);

    data->meshFilename = "data/tet4.mesh";
    data->bcLabel = "boundary";
    data->spaceDim = 3;
    data->lengthScale = 10.0;

    data->time = 1.0;
    data->timeFormat = "%3.1f";

    // Vertex fields ------------------------------------------------------------------------------
    static const size_t vertexNumPoints = 5;
    static const size_t vertexNumDOF = 1 + 3 + 6 + 2;
    data->vertexDiscretization = pylith::topology::FieldBase::Discretization(1, 1, 3);

    data->vertexNumPoints = vertexNumPoints;
    data->vertexNumDOF = vertexNumDOF;
    static const PylithScalar vertexValues[vertexNumPoints*vertexNumDOF] = {
        2.1,    1.1,  2.2,  3.3,   1.1, 1.2, 1.3, 1.4, 1.5, 1.6,    1.2,  2.3,
        3.2,    4.4,  5.5,  6.6,   2.1, 2.2, 2.3, 2.4, 2.5, 2.6,    3.4,  4.5,
        4.3,    7.7,  8.8,  9.9,   3.1, 3.2, 3.3, 3.4, 3.5, 3.6,    5.6,  6.7,
        5.4,   10.0, 11.1, 12.2,   4.1, 4.2, 4.3, 4.4, 4.5, 4.6,    7.8,  8.9,
        6.5,   13.3, 14.4, 15.5,   5.1, 5.2, 5.3, 5.4, 5.5, 5.6,    9.0, 10.1,
    };
    data->vertexValues = const_cast<PylithScalar*>(vertexValues);

    // Cell fields --------------------------------------------------------------------------------
    static const size_t cellNumPoints = 2;
    static const size_t cellNumDOF = 1 + 3 + 6 + 2;
    data->cellDiscretization = pylith::topology::FieldBase::Discretization(0, 0, 2);

    data->cellNumPoints = cellNumPoints;
    data->cellNumDOF = cellNumDOF;
    static const PylithScalar cellValues[cellNumPoints*cellNumDOF] = {
        2.1,   1.1, 2.2, 3.3,   1.2, 2.3, 3.4,  4.5,  5.6,  6.7,   1.2, 2.3,
        3.2,   4.4, 5.5, 6.6,   7.8, 8.9, 9.0, 10.1, 11.2, 12.3,   7.8, 8.9,
    };
    data->cellValues = const_cast<PylithScalar*>(cellValues);
} // setDataTet


// ================================================================================================
void
pylith::meshio::TestDataWriterSubmesh::_setDataHex(void) {
    TestDataWriterSubmesh_Data* data = this->_getData();CPPUNIT_ASSERT(data);

    data->meshFilename = "data/hex8.mesh";
    data->bcLabel = "top";
    data->spaceDim = 2;
    data->lengthScale = 10.0;

    data->time = 1.0;
    data->timeFormat = "%3.1f";

    // Vertex fields ------------------------------------------------------------------------------
    static const size_t vertexNumPoints = 12;
    static const size_t vertexNumDOF = 1 + 3 + 6 + 2;
    data->vertexDiscretization = pylith::topology::FieldBase::Discretization(1, 1, 3);

    data->vertexNumPoints = vertexNumPoints;
    data->vertexNumDOF = vertexNumDOF;
    static const PylithScalar vertexValues[vertexNumPoints*vertexNumDOF] = {
        2.1,     1.1,  2.2,  3.3,    1.1,  1.2,  1.3,  1.4,  1.5,  1.6,   1.2, 2.3,
        3.2,     4.4,  5.5,  6.6,    2.1,  2.2,  2.3,  2.4,  2.5,  2.6,   3.4, 4.5,
        4.3,     7.7,  8.8,  9.9,    3.1,  3.2,  3.3,  3.4,  3.5,  3.6,   5.6, 6.7,
        5.4,    10.1, 11.2, 12.3,    4.1,  4.2,  4.3,  4.4,  4.5,  4.6,   7.8, 8.9,
        6.5,     1.2,  2.3,  3.4,    5.1,  5.2,  5.3,  5.4,  5.5,  5.6,   1.3, 2.4,
        7.6,     4.5,  5.6,  6.7,    6.1,  6.2,  6.3,  6.4,  6.5,  6.6,   3.5, 4.6,
        8.7,     7.8,  8.9,  9.0,    7.1,  7.2,  7.3,  7.4,  7.5,  7.6,   5.7, 6.8,
        9.8,    10.2, 11.3, 12.4,    8.1,  8.2,  8.3,  8.4,  8.5,  8.6,   7.9, 8.0,
        10.9,    1.3,  2.4,  3.5,    9.1,  9.2,  9.3,  9.4,  9.5,  9.6,   8.1, 8.2,
        11.8,    4.6,  5.7,  6.8,   10.1, 10.2, 10.3, 10.4, 10.5, 10.6,   9.2, 9.3,
        12.7,    7.9,  8.1,  9.2,   11.1, 11.2, 11.3, 11.4, 11.5, 11.6,  10.4, 10.5,
        13.6,   10.3, 11.4, 12.5,   12.1, 12.2, 12.3, 12.4, 12.5, 12.6,  11.5, 11.6,
    };
    data->vertexValues = const_cast<PylithScalar*>(vertexValues);

    // Cell fields --------------------------------------------------------------------------------
    static const size_t cellNumPoints = 2;
    static const size_t cellNumDOF = 1 + 3 + 6 + 2;
    data->cellDiscretization = pylith::topology::FieldBase::Discretization(0, 0, 2);

    data->cellNumPoints = cellNumPoints;
    data->cellNumDOF = cellNumDOF;
    static const PylithScalar cellValues[cellNumPoints*cellNumDOF] = {
        2.1,   1.1, 2.2, 3.3,   1.2, 2.3, 3.4,  4.5,  5.6,  6.7,   1.2, 2.3,
        3.2,   4.4, 5.5, 6.6,   7.8, 8.9, 9.0, 10.1, 11.2, 12.3,   7.8, 8.9,
    };
    data->cellValues = const_cast<PylithScalar*>(cellValues);
} // setDataHex


// ------------------------------------------------------------------------------------------------
// Constructor
pylith::meshio::TestDataWriterSubmesh_Data::TestDataWriterSubmesh_Data(void) :
    bcLabel(NULL) {}


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::meshio::TestDataWriterSubmesh_Data::~TestDataWriterSubmesh_Data(void) {}


// End of file
