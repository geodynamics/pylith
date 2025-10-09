// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "TestDataWriterPoints.hh" // Implementation of class methods

#include "FieldFactory.hh" // USES FieldFactory
#include "tests/src/FaultCohesiveStub.hh" // USES FaultCohesiveStub

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps::nondimensionalize()
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/meshio/DataWriter.hh" // USES DataWriter
#include "pylith/meshio/MeshBuilder.hh" // USES MeshBuilder
#include "pylith/utils/error.hh" // USES PYLITH_METHOD*

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/units/Scales.hh" // USES Scales

// ------------------------------------------------------------------------------------------------
// Constructor.
pylith::meshio::TestDataWriterPoints::TestDataWriterPoints(void) :
    _pointMesh(nullptr) {}


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::meshio::TestDataWriterPoints::~TestDataWriterPoints(void) {
    PYLITH_METHOD_BEGIN;

    delete _pointMesh;_pointMesh = nullptr;

    PYLITH_METHOD_END;
} // destructor


// ------------------------------------------------------------------------------------------------
void
pylith::meshio::TestDataWriterPoints::setDataTri(TestDataWriterPoints_Data* data) {
    assert(data);

    data->meshFilename = "data/tri3.mesh";
    data->faultLabel = "fault";
    data->faultId = 100;
    data->spaceDim = 2;
    data->lengthScale = 0.1;

    data->time = 1.0;
    data->timeFormat = "%3.1f";

    static const size_t numPoints = 3;
    data->numPoints = numPoints;
    static const PylithReal points[numPoints*2] = {
        -0.3333333, 0.0,
        0.0000001, 0.0,
        0.9999999, 0.0,
    };
    data->points = const_cast<PylithReal*>(points);
    data->names = pylith::string_vector({"ZZ.A", "ZZ.B", "ZZ.C"});

    // Vertex fields ------------------------------------------------------------------------------
    static const size_t vertexNumPoints = 3;
    static const size_t vertexNumDOF = 1 + 2 + 3 + 2;
    data->vertexDiscretization = pylith::topology::FieldBase::Discretization(1, 1, 2);

    data->vertexNumPoints = vertexNumPoints;
    data->vertexNumDOF = vertexNumDOF;
    static const PylithScalar vertexValues[vertexNumPoints*vertexNumDOF] = {
        3.2,    3.3,  4.4,   2.1, 2.2, 2.3,    3.4,  4.5,
        7.05,  10.5, 11.1,   5.6, 5.7, 5.8,   10.1, 11.2,
        5.4,    7.7,  8.8,   4.1, 4.2, 4.3,    7.8,  8.9,
    };
    data->vertexValues = const_cast<PylithScalar*>(vertexValues);

} // setDataTri


// ------------------------------------------------------------------------------------------------
void
pylith::meshio::TestDataWriterPoints::setDataQuad(TestDataWriterPoints_Data* data) {
    assert(data);

    data->meshFilename = "data/quad4.mesh";
    data->faultId = 100;
    data->spaceDim = 2;
    data->lengthScale = 0.1;

    data->time = 1.0;
    data->timeFormat = "%3.1f";

    static const size_t numPoints = 3;
    data->numPoints = numPoints;
    static const PylithReal points[numPoints*2] = {
        -0.5, 0.0,
        0.00000001, 0.0,
        0.99999999, -0.99999999,
    };
    data->points = const_cast<PylithReal*>(points);
    data->names = pylith::string_vector({"ZZ.A", "ZZ.B", "ZZ.C"});

    // Vertex fields ------------------------------------------------------------------------------
    static const size_t vertexNumPoints = 3;
    static const size_t vertexNumDOF = 1 + 2 + 3 + 2;
    data->vertexDiscretization = pylith::topology::FieldBase::Discretization(1, 1, 2);

    data->vertexNumPoints = vertexNumPoints;
    data->vertexNumDOF = vertexNumDOF;
    static const PylithScalar vertexValues[vertexNumPoints*vertexNumDOF] = {
        3.75,   4.4,  5.5,   2.6, 2.7, 2.8,   4.5,  5.6,
        4.85,   6.6,  7.7,   3.6, 3.7, 3.8,   6.7,  7.8,
        6.50,   9.9, 10.1,   5.1, 5.2, 5.3,   9.8,  7.6,
    };
    data->vertexValues = const_cast<PylithScalar*>(vertexValues);
} // setDataQuad


// ------------------------------------------------------------------------------------------------
void
pylith::meshio::TestDataWriterPoints::setDataTet(TestDataWriterPoints_Data* data) {
    assert(data);

    data->meshFilename = "data/tet4.mesh";
    data->faultId = 100;
    data->spaceDim = 3;
    data->lengthScale = 0.1;

    data->time = 1.0;
    data->timeFormat = "%3.1f";

    static const size_t numPoints = 4;
    data->numPoints = numPoints;
    static const PylithReal points[numPoints*3] = {
        -0.33333333, 0.0, 0.33333333,
        +0.00000001, 0.0, 0.33333333,
        +0.00000001, 0.0, 0.00000001,
        0.0, -0.99999999, 0.00000001,
    };
    data->points = const_cast<PylithReal*>(points);
    data->names = pylith::string_vector({"ZZ.A", "ZZ.B", "ZZ.C", "ZZ.DD"});

    // Vertex fields ------------------------------------------------------------------------------
    static const size_t vertexNumPoints = 4;
    static const size_t vertexNumDOF = 1 + 3 + 6 + 2;
    data->vertexDiscretization = pylith::topology::FieldBase::Discretization(1, 1, 3);

    data->vertexNumPoints = vertexNumPoints;
    data->vertexNumDOF = vertexNumDOF;
    static const PylithScalar vertexValues[vertexNumPoints*vertexNumDOF] = {
        3.566667,    5.333333,  6.433333,  7.533333,     2.433333, 2.533333, 2.633333, 2.733333, 2.833333, 2.933333,    4.133333,  5.233333,
        8.7,        19.566667, 20.333333, 21.433333,     7.1,      7.2,      7.3,      7.4,      7.5,      7.6,        13.4,      14.5,
        8.7,        19.4,      20.5,      21.6,          7.1,      7.2,      7.3,      7.4,      7.5,      7.6,        13.4,      14.5,
        3.2,         4.4,       5.5,       6.6,          2.1,      2.2,      2.3,      2.4,      2.5,      2.6,         3.4,       4.5,
    };
    data->vertexValues = const_cast<PylithScalar*>(vertexValues);
} // setDataTet


// ------------------------------------------------------------------------------------------------
void
pylith::meshio::TestDataWriterPoints::setDataHex(TestDataWriterPoints_Data* data) {
    assert(data);

    data->meshFilename = "data/hex8.mesh";
    data->faultId = 100;
    data->spaceDim = 3;
    data->lengthScale = 0.1;

    data->time = 1.0;
    data->timeFormat = "%3.1f";

    static const size_t numPoints = 4;
    data->numPoints = numPoints;
    static const PylithReal points[numPoints*3] = {
        -0.5, 0.0, 0.5,
        -0.00000001, 0.0, 0.0,
        -0.00000001, 0.0, 0.99999999,
        0.99999999, 0.99999999, -0.99999999,
    };
    data->points = const_cast<PylithReal*>(points);
    data->names = pylith::string_vector({"ZZ.A", "ZZ.B", "ZZ.C", "ZZ.DD"});

    // Vertex fields ------------------------------------------------------------------------------
    static const size_t vertexNumPoints = 4;
    static const size_t vertexNumDOF = 1 + 3 + 6 + 2;
    data->vertexDiscretization = pylith::topology::FieldBase::Discretization(1, 1, 3);

    data->vertexNumPoints = vertexNumPoints;
    data->vertexNumDOF = vertexNumDOF;
    static const PylithScalar vertexValues[vertexNumPoints*vertexNumDOF] = {
        8.55,    5.9375, 7.0375, 7.950,    7.1, 7.2, 7.3, 7.4,  7.5,  7.6,    4.575, 5.4875,
        7.95,    5.9250, 7.0250, 8.125,    6.6, 6.7, 6.8, 6.9,  7.0,  7.1,    4.550, 5.6500,
        11.05,   2.9500, 4.0500, 5.150,    9.6, 9.7, 9.8, 9.9, 10.0, 10.1,    2.400, 3.5000,
        7.60,    4.5000, 5.6000, 6.700,    6.1, 6.2, 6.3, 6.4,  6.5,  6.6,    3.500, 4.6000,

    };
    data->vertexValues = const_cast<PylithScalar*>(vertexValues);
} // setDataHex


// ------------------------------------------------------------------------------------------------
// Initialize mesh.
void
pylith::meshio::TestDataWriterPoints::_initialize(void) {
    PYLITH_METHOD_BEGIN;

    const TestDataWriterPoints_Data* data = _getData();assert(data);

    spatialdata::geocoords::CSCart cs;
    cs.setSpaceDim(data->spaceDim);

    delete _pointMesh;_pointMesh = pylith::topology::MeshOps::createFromPoints(
        data->points, data->numPoints, &cs, data->lengthScale, PETSC_COMM_WORLD, "points");

    spatialdata::units::Scales scales;
    scales.setLengthScale(data->lengthScale);
    pylith::topology::MeshOps::nondimensionalize(_pointMesh, scales);

    PYLITH_METHOD_END;
} // _initialize


// ------------------------------------------------------------------------------------------------
// Create vertex fields.
void
pylith::meshio::TestDataWriterPoints::_createVertexField(pylith::topology::Field* field) {
    PYLITH_METHOD_BEGIN;
    assert(field);

    const TestDataWriterPoints_Data* data = _getData();assert(data);

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
} // _createVertexField


// ------------------------------------------------------------------------------------------------
// Constructor
pylith::meshio::TestDataWriterPoints_Data::TestDataWriterPoints_Data(void) :
    numPoints(0),
    points(NULL) {}


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::meshio::TestDataWriterPoints_Data::~TestDataWriterPoints_Data(void) {}


// End of file
