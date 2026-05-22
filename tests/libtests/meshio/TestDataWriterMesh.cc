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

#include "TestDataWriterMesh.hh" // Implementation of class methods

#include "FieldFactory.hh" // USES FieldFactory
#include "tests/src/FaultCohesiveStub.hh" // USES FaultCohesiveStub

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps::nondimensionalize()
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/meshio/DataWriter.hh" // USES DataWriter

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "pylith/scales/Scales.hh" // USES Scales
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include "catch2/catch_test_macros.hpp"


namespace pylith {
    namespace meshio {
        class _TestDataWriterMesh {
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
pylith::meshio::TestDataWriterMesh::TestDataWriterMesh(void) :
    _mesh(nullptr) {}


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::meshio::TestDataWriterMesh::~TestDataWriterMesh(void) {
    PYLITH_METHOD_BEGIN;

    delete _mesh;_mesh = NULL;

    PYLITH_METHOD_END;
} // destructor


// ------------------------------------------------------------------------------------------------
void
pylith::meshio::TestDataWriterMesh::setDataTri(TestDataWriter_Data* data) {
    REQUIRE(data);

    data->meshFilename = "data/tri3.mesh";
    data->faultLabel = "fault";
    data->faultId = 100;
    data->spaceDim = 2;
    data->lengthScale = 0.1;

    data->time = 1.0;
    data->timeFormat = "%3.1f";

    data->fieldFn = _TestDataWriterMesh::analyticalFn;
    data->vertexDiscretization = pylith::topology::FieldBase::Discretization(1, 1);
    data->cellDiscretization = pylith::topology::FieldBase::Discretization(0, 0);
} // setDataTri


// ------------------------------------------------------------------------------------------------
void
pylith::meshio::TestDataWriterMesh::setDataQuad(TestDataWriter_Data* data) {
    REQUIRE(data);

    // We do not use a fault in this test case.
    data->meshFilename = "data/quad4.mesh";
    data->spaceDim = 2;
    data->lengthScale = 0.1;

    data->time = 1.0;
    data->timeFormat = "%3.1f";

    data->fieldFn = _TestDataWriterMesh::analyticalFn;
    data->vertexDiscretization = pylith::topology::FieldBase::Discretization(1, 1);
    data->cellDiscretization = pylith::topology::FieldBase::Discretization(0, 0);
} // setDataQuad


// ------------------------------------------------------------------------------------------------
void
pylith::meshio::TestDataWriterMesh::setDataTet(TestDataWriter_Data* data) {
    REQUIRE(data);

    data->meshFilename = "data/tet4.mesh";
    data->faultLabel = "fault";
    data->faultId = 100;
    data->spaceDim = 2;
    data->lengthScale = 0.1;

    data->time = 1.0;
    data->timeFormat = "%3.1f";

    data->fieldFn = _TestDataWriterMesh::analyticalFn;
    data->vertexDiscretization = pylith::topology::FieldBase::Discretization(1, 1);
    data->cellDiscretization = pylith::topology::FieldBase::Discretization(0, 0);
} // setDataTet


// ------------------------------------------------------------------------------------------------
void
pylith::meshio::TestDataWriterMesh::setDataHex(TestDataWriter_Data* data) {
    REQUIRE(data);

    data->meshFilename = "data/hex8.mesh";
    data->faultLabel = "fault";
    data->faultId = 100;
    data->spaceDim = 2;
    data->lengthScale = 0.1;

    data->time = 1.0;
    data->timeFormat = "%3.1f";

    data->fieldFn = _TestDataWriterMesh::analyticalFn;
    data->vertexDiscretization = pylith::topology::FieldBase::Discretization(1, 1);
    data->cellDiscretization = pylith::topology::FieldBase::Discretization(0, 0);
} // setDataHex


// ------------------------------------------------------------------------------------------------
// Initialize mesh.
void
pylith::meshio::TestDataWriterMesh::_initialize(void) {
    PYLITH_METHOD_BEGIN;

    const TestDataWriter_Data* data = _getData();REQUIRE(data);

    delete _mesh;_mesh = new topology::Mesh;REQUIRE(_mesh);
    MeshIOAscii iohandler;
    iohandler.setFilename(data->meshFilename);
    iohandler.read(_mesh);

    spatialdata::geocoords::CSCart cs;
    cs.setSpaceDim(_mesh->getDimension());
    _mesh->setCoordSys(&cs);

    pylith::scales::Scales scales;
    scales.setLengthScale(data->lengthScale);
    pylith::topology::MeshOps::nondimensionalize(_mesh, scales);

    if (data->faultLabel) {
        pylith::faults::FaultCohesiveStub fault;
        fault.setCohesiveLabelValue(data->faultId);
        fault.setSurfaceLabelName(data->faultLabel);
        pylith::topology::Mesh* meshNew = fault.transformTopology(_mesh);
        delete _mesh;_mesh = meshNew;
    } // if

    PYLITH_METHOD_END;
} // _initialize


// ------------------------------------------------------------------------------------------------
// Create vertex fields.
void
pylith::meshio::TestDataWriterMesh::_createVertexField(pylith::topology::Field* field) {
    PYLITH_METHOD_BEGIN;

    REQUIRE(field);

    const TestDataWriter_Data* data = _getData();REQUIRE(data);

    FieldFactory factory(*field);
    factory.addScalar(data->vertexDiscretization);
    factory.addVector(data->vertexDiscretization);
    factory.addTensor(data->vertexDiscretization);
    factory.addOther(data->vertexDiscretization);

    field->subfieldsSetup();
    field->createDiscretization();
    field->allocate();

    factory.setValues(data->fieldFn);

    field->createOutputVector();
    field->scatterLocalToOutput();

    PYLITH_METHOD_END;
} // _createVertexField


// ------------------------------------------------------------------------------------------------
// Create cell fields.
void
pylith::meshio::TestDataWriterMesh::_createCellField(pylith::topology::Field* field) {
    PYLITH_METHOD_BEGIN;

    REQUIRE(field);

    const TestDataWriter_Data* data = _getData();REQUIRE(data);

    FieldFactory factory(*field);
    factory.addScalar(data->cellDiscretization);
    factory.addVector(data->cellDiscretization);
    factory.addTensor(data->cellDiscretization);
    factory.addOther(data->cellDiscretization);

    field->subfieldsSetup();
    field->createDiscretization();
    field->allocate();

    factory.setValues(data->fieldFn);

    field->createOutputVector();
    field->scatterLocalToOutput();

    PYLITH_METHOD_END;
} // _createCellField


// ------------------------------------------------------------------------------------------------
PetscErrorCode
pylith::meshio::_TestDataWriterMesh::analyticalFn(PetscInt spaceDim,
                                                  PetscReal t,
                                                  const PetscReal x[],
                                                  PetscInt numComponents,
                                                  PetscScalar* values,
                                                  void* context) {
    PYLITH_METHOD_BEGIN;
    REQUIRE(x);
    REQUIRE(values);
    REQUIRE(context);

    const PetscReal X_FAULT = 0.0;

    PetscErrorCode err = PETSC_SUCCESS;
    PetscDM dmMesh = PetscDM(context);
    PetscInt cell = -1;
    err = DMPlexGetActivePoint(dmMesh, &cell);
    if (err) { PYLITH_METHOD_RETURN(err); }

    PetscReal centroid[3] = { 0.0, 0.0, 0.0 };
    err = DMPlexComputeCellGeometryFVM(dmMesh, cell, nullptr, centroid, nullptr);
    if (err) { PYLITH_METHOD_RETURN(err); }

    const PetscScalar constants[3] = { 11.0, 13.0, 17.0 };
    const PetscScalar faultOffset = (centroid[0] < X_FAULT) ? -1.1 : +1.1;
    const PetscScalar fieldBase = 2.0*numComponents/7.0;
    for (PetscInt iC = 0; iC < numComponents; ++iC) {
        values[iC] = faultOffset + fieldBase + iC;
        for (PetscInt iDim = 0; iDim < spaceDim; ++iDim) {
            values[iC] += x[iDim] / constants[iDim];
        } // for
    } // for

    PYLITH_METHOD_RETURN(err);
} // analyticalFn


// End of file
