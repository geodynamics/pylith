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

#include "TestDataWriterMaterial.hh" // Implementation of class methods

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
pylith::meshio::TestDataWriterMaterial::setUp(void) {
    PYLITH_METHOD_BEGIN;

    _domainMesh = NULL;
    _materialMesh = NULL;

    PYLITH_METHOD_END;
} // setUp


// ------------------------------------------------------------------------------------------------
// Tear down testing data.
void
pylith::meshio::TestDataWriterMaterial::tearDown(void) {
    PYLITH_METHOD_BEGIN;

    delete _domainMesh;_domainMesh = NULL;
    delete _materialMesh;_materialMesh = NULL;

    PYLITH_METHOD_END;
} // tearDown


// ------------------------------------------------------------------------------------------------
// Initialize mesh.
void
pylith::meshio::TestDataWriterMaterial::_initialize(void) {
    PYLITH_METHOD_BEGIN;

    const TestDataWriterMaterial_Data* data = _getData();CPPUNIT_ASSERT(data);

    delete _domainMesh;_domainMesh = new topology::Mesh;CPPUNIT_ASSERT(_domainMesh);
    MeshIOAscii iohandler;
    iohandler.filename(data->meshFilename);
    iohandler.read(_domainMesh);

    spatialdata::geocoords::CSCart cs;
    cs.setSpaceDim(_domainMesh->getDimension());
    _domainMesh->setCoordSys(&cs);

    spatialdata::units::Nondimensional normalizer;
    normalizer.setLengthScale(data->lengthScale);
    pylith::topology::MeshOps::nondimensionalize(_domainMesh, normalizer);

    if (data->faultLabel) {
        pylith::faults::FaultCohesiveStub fault;
        fault.setSurfaceMarkerLabel(data->faultLabel);
        fault.setInterfaceId(data->faultId);
        fault.adjustTopology(_domainMesh);
    } // if

    delete _materialMesh;
    _materialMesh = pylith::topology::MeshOps::createSubdomainMesh(*_domainMesh, "material-id", data->materialId, ":UNKNOWN:");
    CPPUNIT_ASSERT(_materialMesh);

    PYLITH_METHOD_END;
} // _initialize


// ------------------------------------------------------------------------------------------------
// Create vertex fields.
void
pylith::meshio::TestDataWriterMaterial::_createVertexField(pylith::topology::Field* field) {
    PYLITH_METHOD_BEGIN;
    CPPUNIT_ASSERT(field);

    const TestDataWriterMaterial_Data* data = _getData();CPPUNIT_ASSERT(data);

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
pylith::meshio::TestDataWriterMaterial::_createCellField(pylith::topology::Field* field) {
    PYLITH_METHOD_BEGIN;
    CPPUNIT_ASSERT(field);

    const TestDataWriterMaterial_Data* data = _getData();CPPUNIT_ASSERT(data);

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
pylith::meshio::TestDataWriterMaterial::_setDataTri(void) {
    TestDataWriterMaterial_Data* data = this->_getData();CPPUNIT_ASSERT(data);

    data->meshFilename = "data/tri3.mesh";
    data->materialId = 0;
    data->faultLabel = "fault";
    data->faultId = 100;
    data->spaceDim = 2;
    data->lengthScale = 10.0;

    data->time = 1.0;
    data->timeFormat = "%3.1f";

    // Vertex fields ------------------------------------------------------------------------------
    static const size_t vertexNumPoints = 6;
    static const size_t vertexNumDOF = 1 + 2 + 3 + 2;
    data->vertexDiscretization = pylith::topology::FieldBase::Discretization(1, 1);

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
    data->cellDiscretization = pylith::topology::FieldBase::Discretization(0, 0);

    data->cellNumPoints = cellNumPoints;
    data->cellNumDOF = cellNumDOF;
    static const PylithScalar cellValues[cellNumPoints*cellNumDOF] = {
        2.1,   1.1, 2.2,   1.2, 2.3, 3.4,   1.2, 2.3,
    };
    data->cellValues = const_cast<PylithScalar*>(cellValues);
} // setDataTri


// ================================================================================================
void
pylith::meshio::TestDataWriterMaterial::_setDataQuad(void) {
    TestDataWriterMaterial_Data* data = this->_getData();CPPUNIT_ASSERT(data);

    // We do not use a fault in this test case.
    data->meshFilename = "data/quad4.mesh";
    data->materialId = 2;
    data->spaceDim = 2;
    data->lengthScale = 10.0;

    data->time = 1.0;
    data->timeFormat = "%3.1f";

    // Vertex fields ------------------------------------------------------------------------------
    static const size_t vertexNumPoints = 6;
    static const size_t vertexNumDOF = 1 + 2 + 3 + 2;
    data->vertexDiscretization = pylith::topology::FieldBase::Discretization(1, 1);

    data->vertexNumPoints = vertexNumPoints;
    data->vertexNumDOF = vertexNumDOF;
    static const PylithScalar vertexValues[vertexNumPoints*vertexNumDOF] = {
        2.1,    1.1,  2.2,   1.1, 1.2, 1.3,    1.2,  2.3,
        3.2,    3.3,  4.4,   2.1, 2.2, 2.3,    3.4,  4.5,
        4.3,    5.5,  6.6,   3.1, 3.2, 3.3,    5.6,  6.7,
        5.4,    7.7,  8.8,   4.1, 4.2, 4.3,    7.8,  8.9,
        6.5,    9.9, 10.1,   5.1, 5.2, 5.3,    9.8,  7.6,
        7.6,   11.2, 12.3,   6.1, 6.2, 6.3,    6.5,  5.4,
    };
    data->vertexValues = const_cast<PylithScalar*>(vertexValues);

    // Cell fields --------------------------------------------------------------------------------
    static const size_t cellNumPoints = 1;
    static const size_t cellNumDOF = 1 + 2 + 3 + 2;
    data->cellDiscretization = pylith::topology::FieldBase::Discretization(0, 0);

    data->cellNumPoints = cellNumPoints;
    data->cellNumDOF = cellNumDOF;
    static const PylithScalar cellValues[cellNumPoints*cellNumDOF] = {
        2.1,   1.1, 2.2,   1.2, 2.3, 3.4,   1.2, 2.3,
    };
    data->cellValues = const_cast<PylithScalar*>(cellValues);
} // setDataQuad


// ================================================================================================
void
pylith::meshio::TestDataWriterMaterial::_setDataTet(void) {
    TestDataWriterMaterial_Data* data = this->_getData();CPPUNIT_ASSERT(data);

    data->meshFilename = "data/tet4.mesh";
    data->materialId = 1;
    data->faultLabel = "fault";
    data->faultId = 100;
    data->spaceDim = 2;
    data->lengthScale = 10.0;

    data->time = 1.0;
    data->timeFormat = "%3.1f";

    // Vertex fields ------------------------------------------------------------------------------
    static const size_t vertexNumPoints = 8;
    static const size_t vertexNumDOF = 1 + 3 + 6 + 2;
    data->vertexDiscretization = pylith::topology::FieldBase::Discretization(1, 1);

    data->vertexNumPoints = vertexNumPoints;
    data->vertexNumDOF = vertexNumDOF;
    static const PylithScalar vertexValues[vertexNumPoints*vertexNumDOF] = {
        2.1,    1.1,  2.2,  3.3,   1.1, 1.2, 1.3, 1.4, 1.5, 1.6,    1.2,  2.3,
        3.2,    4.4,  5.5,  6.6,   2.1, 2.2, 2.3, 2.4, 2.5, 2.6,    3.4,  4.5,
        4.3,    7.7,  8.8,  9.9,   3.1, 3.2, 3.3, 3.4, 3.5, 3.6,    5.6,  6.7,
        5.4,   10.0, 11.1, 12.2,   4.1, 4.2, 4.3, 4.4, 4.5, 4.6,    7.8,  8.9,
        6.5,   13.3, 14.4, 15.5,   5.1, 5.2, 5.3, 5.4, 5.5, 5.6,    9.0, 10.1,
        7.6,   16.6, 17.7, 18.8,   6.1, 6.2, 6.3, 6.4, 6.5, 6.6,   11.2, 12.3,
        8.7,   19.9, 20.0, 21.1,   7.1, 7.2, 7.3, 7.4, 7.5, 7.6,   13.4, 14.5,
        9.8,   22.2, 23.3, 24.4,   8.1, 8.2, 8.3, 8.4, 8.5, 8.6,   15.6, 16.7,
    };
    data->vertexValues = const_cast<PylithScalar*>(vertexValues);

    // Cell fields --------------------------------------------------------------------------------
    static const size_t cellNumPoints = 2;
    static const size_t cellNumDOF = 1 + 3 + 6 + 2;
    data->cellDiscretization = pylith::topology::FieldBase::Discretization(0, 0);

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
pylith::meshio::TestDataWriterMaterial::_setDataHex(void) {
    TestDataWriterMaterial_Data* data = this->_getData();CPPUNIT_ASSERT(data);

    data->meshFilename = "data/hex8.mesh";
    data->materialId = 0;
    data->faultLabel = "fault";
    data->faultId = 100;
    data->spaceDim = 2;
    data->lengthScale = 10.0;

    data->time = 1.0;
    data->timeFormat = "%3.1f";

    // Vertex fields ------------------------------------------------------------------------------
    static const size_t vertexNumPoints = 16;
    static const size_t vertexNumDOF = 1 + 3 + 6 + 2;
    data->vertexDiscretization = pylith::topology::FieldBase::Discretization(1, 1);

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
        10.0,    1.3,  2.4,  3.5,    9.1,  9.2,  9.3,  9.4,  9.5,  9.6,   1.3, 2.4,
        12.1,    4.6,  5.7,  6.8,   10.1, 10.2, 10.3, 10.4, 10.5, 10.6,   3.5, 4.6,
        11.1,    7.9,  8.0,  9.1,   11.1, 11.2, 11.3, 11.4, 11.5, 11.6,   5.7, 6.8,
        13.1,   10.2, 11.3, 12.4,   12.1, 12.2, 12.3, 12.4, 12.5, 12.6,   8.0, 1.4,
        14.1,   13.5, 14.6, 15.7,   13.1, 13.2, 13.3, 13.4, 13.5, 13.6,   2.5, 3.6,
        15.1,   16.8, 17.9, 18.1,   14.1, 14.2, 14.3, 14.4, 14.5, 14.6,   4.8, 1.5,
        16.1,   19.2, 20.3, 21.4,   15.1, 15.2, 15.3, 15.4, 15.5, 15.6,   2.6, 3.7,
        17.1,   22.5, 23.6, 24.7,   16.1, 16.2, 16.3, 16.4, 16.5, 16.6,   4.8, 5.9,
    };
    data->vertexValues = const_cast<PylithScalar*>(vertexValues);

    // Cell fields --------------------------------------------------------------------------------
    static const size_t cellNumPoints = 1;
    static const size_t cellNumDOF = 1 + 3 + 6 + 2;
    data->cellDiscretization = pylith::topology::FieldBase::Discretization(0, 0);

    data->cellNumPoints = cellNumPoints;
    data->cellNumDOF = cellNumDOF;
    static const PylithScalar cellValues[cellNumPoints*cellNumDOF] = {
        2.1,   1.1, 2.2, 3.3,   1.2, 2.3, 3.4, 4.5, 5.6, 6.7,   1.2, 2.3,
    };
    data->cellValues = const_cast<PylithScalar*>(cellValues);
} // setDataHex


// ------------------------------------------------------------------------------------------------
// Constructor
pylith::meshio::TestDataWriterMaterial_Data::TestDataWriterMaterial_Data(void) :
    materialId(0) {}


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::meshio::TestDataWriterMaterial_Data::~TestDataWriterMaterial_Data(void) {}


// End of file
