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

#include "TestDataWriterMaterial.hh" // Implementation of class methods

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
        class _TestDataWriterMaterial {
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
pylith::meshio::TestDataWriterMaterial::TestDataWriterMaterial(void) :
    _domainMesh(nullptr),
    _materialMesh(nullptr) {}


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::meshio::TestDataWriterMaterial::~TestDataWriterMaterial(void) {
    PYLITH_METHOD_BEGIN;

    delete _domainMesh;_domainMesh = nullptr;
    delete _materialMesh;_materialMesh = nullptr;

    PYLITH_METHOD_END;
} // destructor


// ------------------------------------------------------------------------------------------------
void
pylith::meshio::TestDataWriterMaterial::setDataTri(TestDataWriterMaterial_Data* data) {
    REQUIRE(data);

    data->meshFilename = "data/tri3.mesh";
    data->materialId = 0;
    data->faultLabel = "fault";
    data->faultId = 100;
    data->spaceDim = 2;
    data->lengthScale = 0.1;

    data->time = 1.0;
    data->timeFormat = "%3.1f";

    data->fieldFn = _TestDataWriterMaterial::analyticalFn;
    data->vertexDiscretization = pylith::topology::FieldBase::Discretization(1, 1);
    data->cellDiscretization = pylith::topology::FieldBase::Discretization(0, 0);
} // setDataTri


// ------------------------------------------------------------------------------------------------
void
pylith::meshio::TestDataWriterMaterial::setDataQuad(TestDataWriterMaterial_Data* data) {
    REQUIRE(data);

    // We do not use a fault in this test case.
    data->meshFilename = "data/quad4.mesh";
    data->materialId = 2;
    data->spaceDim = 2;
    data->lengthScale = 0.1;

    data->time = 1.0;
    data->timeFormat = "%3.1f";

    data->fieldFn = _TestDataWriterMaterial::analyticalFn;
    data->vertexDiscretization = pylith::topology::FieldBase::Discretization(1, 1);
    data->cellDiscretization = pylith::topology::FieldBase::Discretization(0, 0);
} // setDataQuad


// ------------------------------------------------------------------------------------------------
void
pylith::meshio::TestDataWriterMaterial::setDataTet(TestDataWriterMaterial_Data* data) {
    REQUIRE(data);

    data->meshFilename = "data/tet4.mesh";
    data->materialId = 1;
    data->faultLabel = "fault";
    data->faultId = 100;
    data->spaceDim = 2;
    data->lengthScale = 0.1;

    data->time = 1.0;
    data->timeFormat = "%3.1f";

    data->fieldFn = _TestDataWriterMaterial::analyticalFn;
    data->vertexDiscretization = pylith::topology::FieldBase::Discretization(1, 1);
    data->cellDiscretization = pylith::topology::FieldBase::Discretization(0, 0);
} // setDataTet


// ------------------------------------------------------------------------------------------------
void
pylith::meshio::TestDataWriterMaterial::setDataHex(TestDataWriterMaterial_Data* data) {
    REQUIRE(data);

    data->meshFilename = "data/hex8.mesh";
    data->materialId = 0;
    data->faultLabel = "fault";
    data->faultId = 100;
    data->spaceDim = 2;
    data->lengthScale = 0.1;

    data->time = 1.0;
    data->timeFormat = "%3.1f";

    data->fieldFn = _TestDataWriterMaterial::analyticalFn;
    data->vertexDiscretization = pylith::topology::FieldBase::Discretization(1, 1);
    data->cellDiscretization = pylith::topology::FieldBase::Discretization(0, 0);
} // setDataHex


// ------------------------------------------------------------------------------------------------
// Initialize mesh.
void
pylith::meshio::TestDataWriterMaterial::_initialize(void) {
    PYLITH_METHOD_BEGIN;

    const TestDataWriterMaterial_Data* data = _getData();REQUIRE(data);

    delete _domainMesh;_domainMesh = new topology::Mesh;REQUIRE(_domainMesh);
    MeshIOAscii iohandler;
    iohandler.setFilename(data->meshFilename);
    iohandler.read(_domainMesh);

    spatialdata::geocoords::CSCart cs;
    cs.setSpaceDim(_domainMesh->getDimension());
    _domainMesh->setCoordSys(&cs);

    pylith::scales::Scales scales;
    scales.setLengthScale(data->lengthScale);
    pylith::topology::MeshOps::nondimensionalize(_domainMesh, scales);

    if (data->faultLabel) {
        pylith::faults::FaultCohesiveStub fault;
        fault.setSurfaceLabelName(data->faultLabel);
        fault.setCohesiveLabelValue(data->faultId);
        fault.adjustTopology(_domainMesh);
    } // if

    delete _materialMesh;
    _materialMesh = pylith::topology::MeshOps::createSubdomainMesh(*_domainMesh, pylith::topology::Mesh::cells_label_name, data->materialId, "material");
    REQUIRE(_materialMesh);

    PYLITH_METHOD_END;
} // _initialize


// ------------------------------------------------------------------------------------------------
// Create vertex fields.
void
pylith::meshio::TestDataWriterMaterial::_createVertexField(pylith::topology::Field* field) {
    PYLITH_METHOD_BEGIN;
    REQUIRE(field);

    const TestDataWriterMaterial_Data* data = _getData();REQUIRE(data);

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
pylith::meshio::TestDataWriterMaterial::_createCellField(pylith::topology::Field* field) {
    PYLITH_METHOD_BEGIN;
    REQUIRE(field);

    const TestDataWriterMaterial_Data* data = _getData();REQUIRE(data);

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
pylith::meshio::TestDataWriterMaterial_Data::TestDataWriterMaterial_Data(void) :
    materialId(0) {}


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::meshio::TestDataWriterMaterial_Data::~TestDataWriterMaterial_Data(void) {}


// ------------------------------------------------------------------------------------------------
PetscErrorCode
pylith::meshio::_TestDataWriterMaterial::analyticalFn(PetscInt spaceDim,
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
