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

#include "TestDataWriterMesh.hh" // Implementation of class methods

#include "FieldFactory.hh" // USES FieldFactory

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps::nondimensionalize()
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/meshio/DataWriter.hh" // USES DataWriter
#include "pylith/faults/FaultCohesiveKin.hh" // USES FaultCohesiveKin

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <cppunit/extensions/HelperMacros.h>


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterMesh::setUp(void)
{ // setUp
    PYLITH_METHOD_BEGIN;

    _mesh = NULL;

    PYLITH_METHOD_END;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::meshio::TestDataWriterMesh::tearDown(void)
{ // tearDown
    PYLITH_METHOD_BEGIN;

    delete _mesh; _mesh = NULL;

    PYLITH_METHOD_END;
} // tearDown

// ----------------------------------------------------------------------
// Initialize mesh.
void
pylith::meshio::TestDataWriterMesh::_initialize(void)
{ // _initialize
    PYLITH_METHOD_BEGIN;

    TestDataWriter_Data* data = _getData();CPPUNIT_ASSERT(data);

    delete _mesh; _mesh = new topology::Mesh;CPPUNIT_ASSERT(_mesh);
    MeshIOAscii iohandler;
    iohandler.filename(data->meshFilename);
    iohandler.read(_mesh);

    spatialdata::geocoords::CSCart cs;
    cs.setSpaceDim(_mesh->dimension());
    _mesh->coordsys(&cs);

    spatialdata::units::Nondimensional normalizer;
    normalizer.lengthScale(data->lengthScale);
    pylith::topology::MeshOps::nondimensionalize(_mesh, normalizer);

    if (data->faultLabel) {
        faults::FaultCohesiveKin fault;
        const bool useLagrangeConstraints = true;
        PetscInt firstFaultVertex = 0;
        PetscInt firstLagrangeVertex = 0, firstFaultCell = 0;
        PetscErrorCode err = DMGetStratumSize(_mesh->dmMesh(), data->faultLabel, 1, &firstLagrangeVertex);PYLITH_CHECK_ERROR(err);
        firstFaultCell = firstLagrangeVertex;
        if (useLagrangeConstraints) {
            firstFaultCell += firstLagrangeVertex;
        } // if
        fault.label(data->faultLabel);
        fault.id(data->faultId);
        fault.adjustTopology(_mesh, &firstFaultVertex, &firstLagrangeVertex, &firstFaultCell);
    } // if

    PYLITH_METHOD_END;
} // _initialize

// ----------------------------------------------------------------------
// Create vertex fields.
void
pylith::meshio::TestDataWriterMesh::_createVertexFields(pylith::topology::Fields* fields)
{ // _createVertexFields
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(fields);
    CPPUNIT_ASSERT(_mesh);

    TestDataWriter_Data* data = _getData();CPPUNIT_ASSERT(data);

    PetscDM dmMesh = _mesh->dmMesh();CPPUNIT_ASSERT(dmMesh);

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
pylith::meshio::TestDataWriterMesh::_createCellFields(pylith::topology::Fields* fields)
{ // _createCellFields
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(fields);
    CPPUNIT_ASSERT(_mesh);

    TestDataWriter_Data* data = _getData();CPPUNIT_ASSERT(data);

    PetscDM dmMesh = _mesh->dmMesh();CPPUNIT_ASSERT(dmMesh);

    FieldFactory factory(*fields);
    factory.scalar(data->cellDiscretization, data->cellScalarValues, data->cellNumPoints, data->cellScalarNumComponents);
    factory.vector(data->cellDiscretization, data->cellVectorValues, data->cellNumPoints, data->cellVectorNumComponents);
    factory.tensor(data->cellDiscretization, data->cellTensorValues, data->cellNumPoints, data->cellTensorNumComponents);
    factory.other(data->cellDiscretization, data->cellOtherValues, data->cellNumPoints, data->cellOtherNumComponents);

    PYLITH_METHOD_END;
} // _createCellFields


// ----------------------------------------------------------------------
void
pylith::meshio::TestDataWriterMesh::_setDataTri(void) {
    TestDataWriter_Data* data = this->_getData();CPPUNIT_ASSERT(data);

    data->meshFilename = "data/tri3.mesh";
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
    data->cellNumPoints = 3;
    data->cellDiscretization.basisOrder = 1;
    data->cellDiscretization.quadOrder = 0;
    data->cellDiscretization.isBasisContinuous = true;
    data->cellDiscretization.feSpace = pylith::topology::FieldBase::POINT_SPACE;

    // Scalar
    data->cellScalarNumComponents = 1;
    static const PylithScalar cellScalarValues[3*1] = {
        2.1, 2.2, 2.3
    };
    data->cellScalarValues = const_cast<PylithScalar*>(cellScalarValues);

    // Vector
    data->cellVectorNumComponents = 2;
    static const PylithScalar cellVectorValues[3*2] = {
        1.1, 2.2,
        3.3, 4.4,
        5.5, 6.6,
    };
    data->cellVectorValues = const_cast<PylithScalar*>(cellVectorValues);

    // Tensor
    data->cellTensorNumComponents = 3;
    static const PylithScalar cellTensorValues[3*3] = {
        1.2, 2.3, 3.4,
        4.5, 5.6, 6.7,
        7.8, 8.9, 9.0
    };
    data->cellTensorValues = const_cast<PylithScalar*>(cellTensorValues);

    // Other
    data->cellOtherNumComponents = 2;
    static const PylithScalar cellOtherValues[3*2] = {
        1.2, 2.3,
        4.5, 5.6,
        7.8, 8.9,
    };
    data->cellOtherValues = const_cast<PylithScalar*>(cellOtherValues);


} // setUp


// End of file
