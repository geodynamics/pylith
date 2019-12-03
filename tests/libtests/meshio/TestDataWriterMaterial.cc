// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestDataWriterMaterial.hh" // Implementation of class methods

#include "FieldFactory.hh" // USES FieldFactory

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps::nondimensionalize()
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/meshio/DataWriter.hh" // USES DataWriter
#include "pylith/testing/FaultCohesiveStub.hh" // USES FaultCohesiveStub

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <cppunit/extensions/HelperMacros.h>

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterMaterial::setUp(void) {
    PYLITH_METHOD_BEGIN;

    _domainMesh = NULL;
    _materialMesh = NULL;

    PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::meshio::TestDataWriterMaterial::tearDown(void) {
    PYLITH_METHOD_BEGIN;

    delete _domainMesh;_domainMesh = NULL;
    delete _materialMesh;_materialMesh = NULL;

    PYLITH_METHOD_END;
} // tearDown


// ----------------------------------------------------------------------
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
    cs.setSpaceDim(_domainMesh->dimension());
    _domainMesh->coordsys(&cs);

    spatialdata::units::Nondimensional normalizer;
    normalizer.lengthScale(data->lengthScale);
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


// ----------------------------------------------------------------------
// Create vertex fields.
void
pylith::meshio::TestDataWriterMaterial::_createVertexFields(pylith::topology::Fields* fields) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(fields);

    const TestDataWriterMaterial_Data* data = _getData();CPPUNIT_ASSERT(data);

    FieldFactory factory(*fields);
    factory.scalar(data->vertexDiscretization, data->vertexScalarValues, data->vertexNumPoints, data->vertexScalarNumComponents);
    factory.vector(data->vertexDiscretization, data->vertexVectorValues, data->vertexNumPoints, data->vertexVectorNumComponents);
    factory.tensor(data->vertexDiscretization, data->vertexTensorValues, data->vertexNumPoints, data->vertexTensorNumComponents);
    factory.other(data->vertexDiscretization, data->vertexOtherValues, data->vertexNumPoints, data->vertexOtherNumComponents);

    PYLITH_METHOD_END;
} // _createVertexFields


// ----------------------------------------------------------------------
// Create cell fields.
void
pylith::meshio::TestDataWriterMaterial::_createCellFields(pylith::topology::Fields* fields) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(fields);

    const TestDataWriterMaterial_Data* data = _getData();CPPUNIT_ASSERT(data);

    FieldFactory factory(*fields);
    factory.scalar(data->cellDiscretization, data->cellScalarValues, data->cellNumPoints, data->cellScalarNumComponents);
    factory.vector(data->cellDiscretization, data->cellVectorValues, data->cellNumPoints, data->cellVectorNumComponents);
    factory.tensor(data->cellDiscretization, data->cellTensorValues, data->cellNumPoints, data->cellTensorNumComponents);
    factory.other(data->cellDiscretization, data->cellOtherValues, data->cellNumPoints, data->cellOtherNumComponents);

    PYLITH_METHOD_END;
} // _createCellFields


// ----------------------------------------------------------------------
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

    // Vertex fields ------------------------------
    data->vertexNumPoints = 6;
    data->vertexDiscretization.basisOrder = 1;
    data->vertexDiscretization.quadOrder = 1;
    data->vertexDiscretization.isBasisContinuous = true;
    data->vertexDiscretization.feSpace = pylith::topology::FieldBase::POLYNOMIAL_SPACE;

    // Scalar
    data->vertexScalarNumComponents = 1;
    static const PylithScalar vertexScalarValues[6*1] = {
        2.1, 3.2, 4.3, 5.4, 6.5, 7.6,
    };
    data->vertexScalarValues = const_cast<PylithScalar*>(vertexScalarValues);

    // Vector
    data->vertexVectorNumComponents = 2;
    static const PylithScalar vertexVectorValues[6*2] = {
        1.1, 2.2,
        3.3, 4.4,
        5.5, 6.6,
        7.7, 8.8,
        9.9, 10.0,
        11.1, 12.2,
    };
    data->vertexVectorValues = const_cast<PylithScalar*>(vertexVectorValues);

    // Tensor
    data->vertexTensorNumComponents = 3;
    static const PylithScalar vertexTensorValues[6*3] = {
        1.1, 1.2, 1.3,
        2.1, 2.2, 2.3,
        3.1, 3.2, 3.3,
        4.1, 4.2, 4.3,
        5.1, 5.2, 5.3,
        6.1, 6.2, 6.3,
    };
    data->vertexTensorValues = const_cast<PylithScalar*>(vertexTensorValues);

    // Other
    data->vertexOtherNumComponents = 2;
    static const PylithScalar vertexOtherValues[6*2] = {
        1.2, 2.3,
        3.4, 4.5,
        5.6, 6.7,
        7.8, 8.9,
        9.0, 10.1,
        11.2, 12.3,
    };
    data->vertexOtherValues = const_cast<PylithScalar*>(vertexOtherValues);

    // Cell fields ------------------------------
    data->cellNumPoints = 1;
    data->cellDiscretization.basisOrder = 0;
    data->cellDiscretization.quadOrder = 0;
    data->cellDiscretization.isBasisContinuous = true;
    data->cellDiscretization.feSpace = pylith::topology::FieldBase::POLYNOMIAL_SPACE;

    // Scalar
    data->cellScalarNumComponents = 1;
    static const PylithScalar cellScalarValues[1*1] = {
        2.1,
    };
    data->cellScalarValues = const_cast<PylithScalar*>(cellScalarValues);

    // Vector
    data->cellVectorNumComponents = 2;
    static const PylithScalar cellVectorValues[1*2] = {
        1.1, 2.2,
    };
    data->cellVectorValues = const_cast<PylithScalar*>(cellVectorValues);

    // Tensor
    data->cellTensorNumComponents = 3;
    static const PylithScalar cellTensorValues[1*3] = {
        1.2, 2.3, 3.4,
    };
    data->cellTensorValues = const_cast<PylithScalar*>(cellTensorValues);

    // Other
    data->cellOtherNumComponents = 2;
    static const PylithScalar cellOtherValues[1*2] = {
        1.2, 2.3,
    };
    data->cellOtherValues = const_cast<PylithScalar*>(cellOtherValues);

} // setDataTri


// ----------------------------------------------------------------------
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

    // Vertex fields ------------------------------
    data->vertexNumPoints = 6;
    data->vertexDiscretization.basisOrder = 1;
    data->vertexDiscretization.quadOrder = 1;
    data->vertexDiscretization.isBasisContinuous = true;
    data->vertexDiscretization.feSpace = pylith::topology::FieldBase::POLYNOMIAL_SPACE;

    // Scalar
    data->vertexScalarNumComponents = 1;
    static const PylithScalar vertexScalarValues[6*1] = {
        2.1, 3.2, 4.3, 5.4, 6.5, 7.6
    };
    data->vertexScalarValues = const_cast<PylithScalar*>(vertexScalarValues);

    // Vector
    data->vertexVectorNumComponents = 2;
    static const PylithScalar vertexVectorValues[6*2] = {
        1.1, 2.2,
        3.3, 4.4,
        5.5, 6.6,
        7.7, 8.8,
        9.9, 10.1,
        11.2, 12.3,
    };
    data->vertexVectorValues = const_cast<PylithScalar*>(vertexVectorValues);

    // Tensor
    data->vertexTensorNumComponents = 3;
    static const PylithScalar vertexTensorValues[6*3] = {
        1.1, 1.2, 1.3,
        2.1, 2.2, 2.3,
        3.1, 3.2, 3.3,
        4.1, 4.2, 4.3,
        5.1, 5.2, 5.3,
        6.1, 6.2, 6.3,
    };
    data->vertexTensorValues = const_cast<PylithScalar*>(vertexTensorValues);

    // Other
    data->vertexOtherNumComponents = 2;
    static const PylithScalar vertexOtherValues[6*2] = {
        1.2, 2.3,
        3.4, 4.5,
        5.6, 6.7,
        7.8, 8.9,
        9.8, 7.6,
        6.5, 5.4
    };
    data->vertexOtherValues = const_cast<PylithScalar*>(vertexOtherValues);

    // Cell fields ------------------------------
    data->cellNumPoints = 1;
    data->cellDiscretization.basisOrder = 0;
    data->cellDiscretization.quadOrder = 0;
    data->cellDiscretization.isBasisContinuous = true;
    data->cellDiscretization.feSpace = pylith::topology::FieldBase::POLYNOMIAL_SPACE;

    // Scalar
    data->cellScalarNumComponents = 1;
    static const PylithScalar cellScalarValues[1*1] = {
        2.1,
    };
    data->cellScalarValues = const_cast<PylithScalar*>(cellScalarValues);

    // Vector
    data->cellVectorNumComponents = 2;
    static const PylithScalar cellVectorValues[1*2] = {
        1.1, 2.2,
    };
    data->cellVectorValues = const_cast<PylithScalar*>(cellVectorValues);

    // Tensor
    data->cellTensorNumComponents = 3;
    static const PylithScalar cellTensorValues[1*3] = {
        1.2, 2.3, 3.4,
    };
    data->cellTensorValues = const_cast<PylithScalar*>(cellTensorValues);

    // Other
    data->cellOtherNumComponents = 2;
    static const PylithScalar cellOtherValues[1*2] = {
        1.2, 2.3,
    };
    data->cellOtherValues = const_cast<PylithScalar*>(cellOtherValues);

} // setDataQuad


// ----------------------------------------------------------------------
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

    // Vertex fields ------------------------------
    data->vertexNumPoints = 8;
    data->vertexDiscretization.basisOrder = 1;
    data->vertexDiscretization.quadOrder = 1;
    data->vertexDiscretization.isBasisContinuous = true;
    data->vertexDiscretization.feSpace = pylith::topology::FieldBase::POLYNOMIAL_SPACE;

    // Scalar
    data->vertexScalarNumComponents = 1;
    static const PylithScalar vertexScalarValues[8*1] = {
        2.1, 3.2, 4.3, 5.4, 6.5, 7.6, 8.7, 9.8,
    };
    data->vertexScalarValues = const_cast<PylithScalar*>(vertexScalarValues);

    // Vector
    data->vertexVectorNumComponents = 3;
    static const PylithScalar vertexVectorValues[8*3] = {
        1.1, 2.2, 3.3,
        4.4, 5.5, 6.6,
        7.7, 8.8, 9.9,
        10.0, 11.1, 12.2,
        13.3, 14.4, 15.5,
        16.6, 17.7, 18.8,
        19.9, 20.0, 21.1,
        22.2, 23.3, 24.4,
    };
    data->vertexVectorValues = const_cast<PylithScalar*>(vertexVectorValues);

    // Tensor
    data->vertexTensorNumComponents = 6;
    static const PylithScalar vertexTensorValues[8*6] = {
        1.1, 1.2, 1.3, 1.4, 1.5, 1.6,
        2.1, 2.2, 2.3, 2.4, 2.5, 2.6,
        3.1, 3.2, 3.3, 3.4, 3.5, 3.6,
        4.1, 4.2, 4.3, 4.4, 4.5, 4.6,
        5.1, 5.2, 5.3, 5.4, 5.5, 5.6,
        6.1, 6.2, 6.3, 6.4, 6.5, 6.6,
        7.1, 7.2, 7.3, 7.4, 7.5, 7.6,
        8.1, 8.2, 8.3, 8.4, 8.5, 8.6,
    };
    data->vertexTensorValues = const_cast<PylithScalar*>(vertexTensorValues);

    // Other
    data->vertexOtherNumComponents = 2;
    static const PylithScalar vertexOtherValues[8*2] = {
        1.2, 2.3,
        3.4, 4.5,
        5.6, 6.7,
        7.8, 8.9,
        9.0, 10.1,
        11.2, 12.3,
        13.4, 14.5,
        15.6, 16.7,
    };
    data->vertexOtherValues = const_cast<PylithScalar*>(vertexOtherValues);

    // Cell fields ------------------------------
    data->cellNumPoints = 2;
    data->cellDiscretization.basisOrder = 0;
    data->cellDiscretization.quadOrder = 0;
    data->cellDiscretization.isBasisContinuous = true;
    data->cellDiscretization.feSpace = pylith::topology::FieldBase::POLYNOMIAL_SPACE;

    // Scalar
    data->cellScalarNumComponents = 1;
    static const PylithScalar cellScalarValues[2*1] = {
        2.1, 3.2,
    };
    data->cellScalarValues = const_cast<PylithScalar*>(cellScalarValues);

    // Vector
    data->cellVectorNumComponents = 3;
    static const PylithScalar cellVectorValues[2*3] = {
        1.1, 2.2, 3.3,
        4.4, 5.5, 6.6,
    };
    data->cellVectorValues = const_cast<PylithScalar*>(cellVectorValues);

    // Tensor
    data->cellTensorNumComponents = 6;
    static const PylithScalar cellTensorValues[2*6] = {
        1.2, 2.3, 3.4, 4.5, 5.6, 6.7,
        7.8, 8.9, 9.0, 10.1, 11.2, 12.3,
    };
    data->cellTensorValues = const_cast<PylithScalar*>(cellTensorValues);

    // Other
    data->cellOtherNumComponents = 2;
    static const PylithScalar cellOtherValues[2*2] = {
        1.2, 2.3,
        7.8, 8.9,
    };
    data->cellOtherValues = const_cast<PylithScalar*>(cellOtherValues);

} // setDataTet


// ----------------------------------------------------------------------
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

    // Vertex fields ------------------------------
    data->vertexNumPoints = 16;
    data->vertexDiscretization.basisOrder = 1;
    data->vertexDiscretization.quadOrder = 1;
    data->vertexDiscretization.isBasisContinuous = true;
    data->vertexDiscretization.feSpace = pylith::topology::FieldBase::POLYNOMIAL_SPACE;

    // Scalar
    data->vertexScalarNumComponents = 1;
    static const PylithScalar vertexScalarValues[16*1] = {
        2.1, 3.2, 4.3, 5.4, 6.5, 7.6, 8.7, 9.8,
        10.0, 12.1, 11.1, 13.1, 14.1, 15.1, 16.1, 17.1,
    };
    data->vertexScalarValues = const_cast<PylithScalar*>(vertexScalarValues);

    // Vector
    data->vertexVectorNumComponents = 3;
    static const PylithScalar vertexVectorValues[16*3] = {
        1.1, 2.2, 3.3,
        4.4, 5.5, 6.6,
        7.7, 8.8, 9.9,
        10.1, 11.2, 12.3,
        1.2, 2.3, 3.4,
        4.5, 5.6, 6.7,
        7.8, 8.9, 9.0,
        10.2, 11.3, 12.4,
        1.3, 2.4, 3.5,
        4.6, 5.7, 6.8,
        7.9, 8.0, 9.1,
        10.2, 11.3, 12.4,
        13.5, 14.6, 15.7,
        16.8, 17.9, 18.1,
        19.2, 20.3, 21.4,
        22.5, 23.6, 24.7,
    };
    data->vertexVectorValues = const_cast<PylithScalar*>(vertexVectorValues);

    // Tensor
    data->vertexTensorNumComponents = 6;
    static const PylithScalar vertexTensorValues[16*6] = {
        1.1, 1.2, 1.3, 1.4, 1.5, 1.6,
        2.1, 2.2, 2.3, 2.4, 2.5, 2.6,
        3.1, 3.2, 3.3, 3.4, 3.5, 3.6,
        4.1, 4.2, 4.3, 4.4, 4.5, 4.6,
        5.1, 5.2, 5.3, 5.4, 5.5, 5.6,
        6.1, 6.2, 6.3, 6.4, 6.5, 6.6,
        7.1, 7.2, 7.3, 7.4, 7.5, 7.6,
        8.1, 8.2, 8.3, 8.4, 8.5, 8.6,
        9.1, 9.2, 9.3, 9.4, 9.5, 9.6,
        10.1, 10.2, 10.3, 10.4, 10.5, 10.6,
        11.1, 11.2, 11.3, 11.4, 11.5, 11.6,
        12.1, 12.2, 12.3, 12.4, 12.5, 12.6,
        13.1, 13.2, 13.3, 13.4, 13.5, 13.6,
        14.1, 14.2, 14.3, 14.4, 14.5, 14.6,
        15.1, 15.2, 15.3, 15.4, 15.5, 15.6,
        16.1, 16.2, 16.3, 16.4, 16.5, 16.6,
    };
    data->vertexTensorValues = const_cast<PylithScalar*>(vertexTensorValues);

    // Other
    data->vertexOtherNumComponents = 2;
    static const PylithScalar vertexOtherValues[16*2] = {
        1.2, 2.3,
        3.4, 4.5,
        5.6, 6.7,
        7.8, 8.9,
        1.3, 2.4,
        3.5, 4.6,
        5.7, 6.8,
        7.9, 8.0,
        1.3, 2.4,
        3.5, 4.6,
        5.7, 6.8,
        8.0, 1.4,
        2.5, 3.6,
        4.8, 1.5,
        2.6, 3.7,
        4.8, 5.9,
    };
    data->vertexOtherValues = const_cast<PylithScalar*>(vertexOtherValues);

    // Cell fields ------------------------------
    data->cellNumPoints = 1;
    data->cellDiscretization.basisOrder = 0;
    data->cellDiscretization.quadOrder = 0;
    data->cellDiscretization.isBasisContinuous = true;
    data->cellDiscretization.feSpace = pylith::topology::FieldBase::POLYNOMIAL_SPACE;

    // Scalar
    data->cellScalarNumComponents = 1;
    static const PylithScalar cellScalarValues[1*1] = {
        2.1,
    };
    data->cellScalarValues = const_cast<PylithScalar*>(cellScalarValues);

    // Vector
    data->cellVectorNumComponents = 3;
    static const PylithScalar cellVectorValues[1*3] = {
        1.1, 2.2, 3.3,
    };
    data->cellVectorValues = const_cast<PylithScalar*>(cellVectorValues);

    // Tensor
    data->cellTensorNumComponents = 6;
    static const PylithScalar cellTensorValues[1*6] = {
        1.2, 2.3, 3.4, 4.5, 5.6, 6.7,
    };
    data->cellTensorValues = const_cast<PylithScalar*>(cellTensorValues);

    // Other
    data->cellOtherNumComponents = 2;
    static const PylithScalar cellOtherValues[1*2] = {
        1.2, 2.3,
    };
    data->cellOtherValues = const_cast<PylithScalar*>(cellOtherValues);

} // setDataHex


// ----------------------------------------------------------------------
// Constructor
pylith::meshio::TestDataWriterMaterial_Data::TestDataWriterMaterial_Data(void) :
    materialId(0) {}


// ----------------------------------------------------------------------
// Destructor
pylith::meshio::TestDataWriterMaterial_Data::~TestDataWriterMaterial_Data(void) {}


// End of file
