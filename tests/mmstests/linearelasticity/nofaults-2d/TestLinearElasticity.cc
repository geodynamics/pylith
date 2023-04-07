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
// Copyright (c) 2010-2022 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestLinearElasticity.hh" // Implementation of class methods

#include "pylith/problems/TimeDependent.hh" // USES TimeDependent

#include "pylith/materials/Query.hh" // USES Query

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps::nondimensionalize()
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/topology/FieldQuery.hh" // USES FieldQuery
#include "pylith/feassemble/AuxiliaryFactory.hh" // USES AuxiliaryFactory
#include "pylith/problems/SolutionFactory.hh" // USES SolutionFactory
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/meshio/MeshIOPetsc.hh" // USES MeshIOPetsc
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // pythia::journal

#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField

// ---------------------------------------------------------------------------------------------------------------------
// Setup testing data.
void
pylith::mmstests::TestLinearElasticity::setUp(void) {
    MMSTest::setUp();

    _data = new TestLinearElasticity_Data();CPPUNIT_ASSERT(_data);
} // setUp


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate testing data.
void
pylith::mmstests::TestLinearElasticity::tearDown(void) {
    delete _data;_data = NULL;

    MMSTest::tearDown();
} // tearDown


// ---------------------------------------------------------------------------------------------------------------------
// Initialize objects for test.
void
pylith::mmstests::TestLinearElasticity::_initialize(void) {
    PYLITH_METHOD_BEGIN;
    CPPUNIT_ASSERT(_mesh);

    PetscErrorCode err = 0;

    if (_data->useAsciiMesh) {
        pylith::meshio::MeshIOAscii iohandler;
        iohandler.setFilename(_data->meshFilename);
        iohandler.read(_mesh);CPPUNIT_ASSERT(_mesh);
    } else {
        if (_data->meshOptions) {
            err = PetscOptionsInsertString(NULL, _data->meshOptions);PYLITH_CHECK_ERROR(err);
        } // if
        pylith::meshio::MeshIOPetsc iohandler;
        iohandler.setFilename(_data->meshFilename);
        iohandler.read(_mesh);CPPUNIT_ASSERT(_mesh);
    } // if/else

    CPPUNIT_ASSERT_MESSAGE("Test mesh does not contain any cells.",
                           pylith::topology::MeshOps::getNumCells(*_mesh) > 0);
    CPPUNIT_ASSERT_MESSAGE("Test mesh does not contain any vertices.",
                           pylith::topology::MeshOps::getNumVertices(*_mesh) > 0);

    // Set up coordinates.
    _mesh->setCoordSys(&_data->cs);
    pylith::topology::MeshOps::nondimensionalize(_mesh, _data->normalizer);

    // Set up material
    _data->material.setBulkRheology(&_data->rheology);
    _data->material.setAuxiliaryFieldDB(&_data->auxDB);

    for (int i = 0; i < _data->numAuxSubfields; ++i) {
        const pylith::topology::FieldBase::Discretization& info = _data->auxDiscretizations[i];
        _data->material.setAuxiliarySubfieldDiscretization(_data->auxSubfields[i], info.basisOrder, info.quadOrder,
                                                           _data->spaceDim, pylith::topology::FieldBase::DEFAULT_BASIS,
                                                           info.feSpace, info.isBasisContinuous);
    } // for

    // Set up problem.
    CPPUNIT_ASSERT(_problem);
    _problem->setNormalizer(_data->normalizer);
    _problem->setGravityField(_data->gravityField);
    pylith::materials::Material* materials[1] = { &_data->material };
    _problem->setMaterials(materials, 1);
    _problem->setBoundaryConditions(_data->bcs.data(), _data->bcs.size());
    _problem->setStartTime(_data->t);
    _problem->setEndTime(_data->t+_data->dt);
    _problem->setInitialTimeStep(_data->dt);
    _problem->setFormulation(_data->formulation);

    // Set up solution field.
    CPPUNIT_ASSERT(!_solution);
    _solution = new pylith::topology::Field(*_mesh);CPPUNIT_ASSERT(_solution);
    _solution->setLabel("solution");
    pylith::problems::SolutionFactory factory(*_solution, _data->normalizer);
    factory.addDisplacement(_data->solnDiscretizations[0]);
    if (pylith::problems::Physics::QUASISTATIC == _data->formulation) {
        CPPUNIT_ASSERT_EQUAL(1, _data->numSolnSubfields);
    } else {
        CPPUNIT_ASSERT_EQUAL(pylith::problems::Physics::DYNAMIC, _data->formulation);
        CPPUNIT_ASSERT_EQUAL(2, _data->numSolnSubfields);
        factory.addVelocity(_data->solnDiscretizations[1]);
    } // if/else
    _problem->setSolution(_solution);

    pylith::testing::MMSTest::_initialize();

    PYLITH_METHOD_END;
} // _initialize


// ---------------------------------------------------------------------------------------------------------------------
// Constructor
pylith::mmstests::TestLinearElasticity_Data::TestLinearElasticity_Data(void) :
    spaceDim(2),
    meshFilename(NULL),
    meshOptions(NULL),
    boundaryLabel(NULL),
    useAsciiMesh(true),

    t(0.0),
    dt(0.05),
    formulation(pylith::problems::Physics::QUASISTATIC),

    gravityField(NULL),

    numSolnSubfields(0),
    solnDiscretizations(NULL),

    numAuxSubfields(0),
    auxSubfields(NULL),
    auxDiscretizations(NULL) {
    auxDB.setDescription("material auxiliary field spatial database");
    cs.setSpaceDim(spaceDim);
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor
pylith::mmstests::TestLinearElasticity_Data::~TestLinearElasticity_Data(void) {
    for (size_t i = 0; i < bcs.size(); ++i) {
        delete bcs[i];bcs[i] = NULL;
    } // for
    delete gravityField;gravityField = NULL;

} // destructor


// End of file
