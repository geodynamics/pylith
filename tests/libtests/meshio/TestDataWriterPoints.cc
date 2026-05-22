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
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/meshio/DataWriter.hh" // USES DataWriter
#include "pylith/meshio/MeshBuilder.hh" // USES MeshBuilder
#include "pylith/utils/error.hh" // USES PYLITH_METHOD*

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "pylith/scales/Scales.hh" // USES Scales

#include "catch2/catch_test_macros.hpp"


namespace pylith {
    namespace meshio {
        class _TestDataWriterPoints {
public:

            static PetscErrorCode analyticalFn(PetscInt spaceDim,
                                               PetscReal t,
                                               const PetscReal x[],
                                               PetscInt numComponents,
                                               PetscScalar* values,
                                               void* context);

        };
    }
}


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
    REQUIRE(data);

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

    data->vertexDiscretization = pylith::topology::FieldBase::Discretization(1, 1, 2);
} // setDataTri


// ------------------------------------------------------------------------------------------------
void
pylith::meshio::TestDataWriterPoints::setDataQuad(TestDataWriterPoints_Data* data) {
    REQUIRE(data);

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

    data->vertexDiscretization = pylith::topology::FieldBase::Discretization(1, 1, 2);
} // setDataQuad


// ------------------------------------------------------------------------------------------------
void
pylith::meshio::TestDataWriterPoints::setDataTet(TestDataWriterPoints_Data* data) {
    REQUIRE(data);

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

    data->vertexDiscretization = pylith::topology::FieldBase::Discretization(1, 1, 3);
} // setDataTet


// ------------------------------------------------------------------------------------------------
void
pylith::meshio::TestDataWriterPoints::setDataHex(TestDataWriterPoints_Data* data) {
    REQUIRE(data);

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

    data->vertexDiscretization = pylith::topology::FieldBase::Discretization(1, 1, 3);
} // setDataHex


// ------------------------------------------------------------------------------------------------
// Initialize mesh.
void
pylith::meshio::TestDataWriterPoints::_initialize(void) {
    PYLITH_METHOD_BEGIN;

    const TestDataWriterPoints_Data* data = _getData();REQUIRE(data);

    spatialdata::geocoords::CSCart cs;
    cs.setSpaceDim(data->spaceDim);

    delete _pointMesh;_pointMesh = pylith::topology::MeshOps::createFromPoints(
        data->points, data->numPoints, &cs, data->lengthScale, PETSC_COMM_WORLD, "points");

    pylith::scales::Scales scales;
    scales.setLengthScale(data->lengthScale);
    pylith::topology::MeshOps::nondimensionalize(_pointMesh, scales);

    PYLITH_METHOD_END;
} // _initialize


// ------------------------------------------------------------------------------------------------
// Create vertex fields.
void
pylith::meshio::TestDataWriterPoints::_createVertexField(pylith::topology::Field* field) {
    PYLITH_METHOD_BEGIN;
    REQUIRE(field);

    const TestDataWriterPoints_Data* data = _getData();REQUIRE(data);

    FieldFactory factory(*field);
    factory.addScalar(data->vertexDiscretization);
    factory.addVector(data->vertexDiscretization);
    factory.addTensor(data->vertexDiscretization);
    factory.addOther(data->vertexDiscretization);

    field->subfieldsSetup();
    field->createDiscretization();
    field->allocate();

    pylith::topology::VecVisitorMesh fieldVisitor(*field);
    PylithScalar* fieldArray = fieldVisitor.localArray();REQUIRE(fieldArray);
    const PylithInt fieldSize = field->getStorageSize();
    for (PylithInt i = 0; i < fieldSize; ++i) {
        fieldArray[i] = 1.1 + PetscScalar(i);
    } // for

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
