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

#include "TestDataWriterSubmesh.hh" // Implementation of class methods

#include "FieldFactory.hh" // USES FieldFactory
#include "tests/src/FaultCohesiveStub.hh" // USES FaultCohesiveStub

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps::nondimensionalize()
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/meshio/DataWriter.hh" // USES DataWriter

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "pylith/scales/Scales.hh" // USES Scales

#include "catch2/catch_test_macros.hpp"


namespace pylith {
    namespace meshio {
        class _TestDataWriterSubmesh {
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
pylith::meshio::TestDataWriterSubmesh::TestDataWriterSubmesh(void) :
    _mesh(nullptr),
    _submesh(nullptr) {}


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::meshio::TestDataWriterSubmesh::~TestDataWriterSubmesh(void) {
    PYLITH_METHOD_BEGIN;

    delete _mesh;_mesh = nullptr;
    delete _submesh;_submesh = nullptr;

    PYLITH_METHOD_END;
} // destructor


// ------------------------------------------------------------------------------------------------
void
pylith::meshio::TestDataWriterSubmesh::setDataTri(TestDataWriterSubmesh_Data* data) {
    REQUIRE(data);

    data->meshFilename = "data/tri3.mesh";
    data->bcLabel = "bc";
    data->faultLabel = "fault";
    data->faultId = 100;
    data->spaceDim = 2;
    data->lengthScale = 0.1;

    data->time = 1.0;
    data->timeFormat = "%3.1f";

    data->fieldFn = _TestDataWriterSubmesh::analyticalFn;
    data->vertexDiscretization = pylith::topology::FieldBase::Discretization(1, 1, 2);
    data->cellDiscretization = pylith::topology::FieldBase::Discretization(0, 0, 1);
} // setDataTri


// ------------------------------------------------------------------------------------------------
void
pylith::meshio::TestDataWriterSubmesh::setDataQuad(TestDataWriterSubmesh_Data* data) {
    REQUIRE(data);

    // We do not use a fault in this test case.
    data->meshFilename = "data/quad4.mesh";
    data->bcLabel = "bc3";
    data->spaceDim = 2;
    data->lengthScale = 0.1;

    data->time = 1.0;
    data->timeFormat = "%3.1f";

    data->fieldFn = _TestDataWriterSubmesh::analyticalFn;
    data->vertexDiscretization = pylith::topology::FieldBase::Discretization(1, 1, 2);
    data->cellDiscretization = pylith::topology::FieldBase::Discretization(0, 0, 1);
} // setDataQuad


// ------------------------------------------------------------------------------------------------
void
pylith::meshio::TestDataWriterSubmesh::setDataTet(TestDataWriterSubmesh_Data* data) {
    REQUIRE(data);

    data->meshFilename = "data/tet4.mesh";
    data->bcLabel = "boundary";
    data->spaceDim = 3;
    data->lengthScale = 0.1;

    data->time = 1.0;
    data->timeFormat = "%3.1f";

    data->fieldFn = _TestDataWriterSubmesh::analyticalFn;
    data->vertexDiscretization = pylith::topology::FieldBase::Discretization(1, 1, 3);
    data->cellDiscretization = pylith::topology::FieldBase::Discretization(0, 0, 2);
} // setDataTet


// ------------------------------------------------------------------------------------------------
void
pylith::meshio::TestDataWriterSubmesh::setDataHex(TestDataWriterSubmesh_Data* data) {
    REQUIRE(data);

    data->meshFilename = "data/hex8.mesh";
    data->bcLabel = "top";
    data->spaceDim = 2;
    data->lengthScale = 0.1;

    data->time = 1.0;
    data->timeFormat = "%3.1f";

    data->fieldFn = _TestDataWriterSubmesh::analyticalFn;
    data->vertexDiscretization = pylith::topology::FieldBase::Discretization(1, 1, 3);
    data->cellDiscretization = pylith::topology::FieldBase::Discretization(0, 0, 2);
} // setDataHex


// ------------------------------------------------------------------------------------------------
// Initialize mesh.
void
pylith::meshio::TestDataWriterSubmesh::_initialize(void) {
    PYLITH_METHOD_BEGIN;

    const TestDataWriterSubmesh_Data* data = _getData();REQUIRE(data);

    delete _mesh;_mesh = new pylith::topology::Mesh;REQUIRE(_mesh);
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
        fault.setSurfaceLabelName(data->faultLabel);
        fault.setCohesiveLabelValue(data->faultId);
        pylith::topology::Mesh* meshNew = fault.transformTopology(_mesh);
        delete _mesh;_mesh = meshNew;
    } // if

    REQUIRE(data->bcLabel);
    const int labelValue = 1;
    delete _submesh;_submesh = pylith::topology::MeshOps::createLowerDimMesh(*_mesh, data->bcLabel, labelValue, "subdomain_bc");REQUIRE(_submesh);

    PYLITH_METHOD_END;
} // _initialize


// ------------------------------------------------------------------------------------------------
// Create vertex fields.
void
pylith::meshio::TestDataWriterSubmesh::_createVertexField(pylith::topology::Field* field) {
    PYLITH_METHOD_BEGIN;
    REQUIRE(field);
    REQUIRE(_mesh);

    const TestDataWriterSubmesh_Data* data = _getData();REQUIRE(data);

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
} // _createVertexFields


// ------------------------------------------------------------------------------------------------
// Create cell fields.
void
pylith::meshio::TestDataWriterSubmesh::_createCellField(pylith::topology::Field* field) {
    PYLITH_METHOD_BEGIN;
    REQUIRE(field);

    const TestDataWriterSubmesh_Data* data = _getData();REQUIRE(data);

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
} // _createCellFields


// ------------------------------------------------------------------------------------------------
// Constructor
pylith::meshio::TestDataWriterSubmesh_Data::TestDataWriterSubmesh_Data(void) :
    bcLabel(NULL) {}


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::meshio::TestDataWriterSubmesh_Data::~TestDataWriterSubmesh_Data(void) {}


// ------------------------------------------------------------------------------------------------
PetscErrorCode
pylith::meshio::_TestDataWriterSubmesh::analyticalFn(PetscInt spaceDim,
                                                     PetscReal t,
                                                     const PetscReal x[],
                                                     PetscInt numComponents,
                                                     PetscScalar* values,
                                                     void* context) {
    PYLITH_METHOD_BEGIN;
    REQUIRE(x);
    REQUIRE(values);

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

    PYLITH_METHOD_RETURN(PETSC_SUCCESS);
} // analyticalFn


// End of file
