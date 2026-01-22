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

#include "TestMomentTensorSource.hh" // Implementation of class methods

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

#include "pylith/scales/ElasticityScales.hh" // USES ElasticityScales

// ------------------------------------------------------------------------------------------------
// Constuctor.
pylith::TestMomentTensorSource::TestMomentTensorSource(TestMomentTensorSource_Data* data) :
    _data(data) {
    assert(_data);

    GenericComponent::setName(_data->journalName);
    _jacobianConvergenceRate = _data->jacobianConvergenceRate;
    _tolerance = _data->tolerance;
    _isJacobianLinear = _data->isJacobianLinear;
    _allowZeroResidual = _data->allowZeroResidual;
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::TestMomentTensorSource::~TestMomentTensorSource(void) {
    delete _data;_data = NULL;
} // destructor


// ------------------------------------------------------------------------------------------------
// Initialize objects for test.
void
pylith::TestMomentTensorSource::_initialize(void) {
    PYLITH_METHOD_BEGIN;
    assert(_mesh);
    assert(_data);

    PetscErrorCode err = PETSC_SUCCESS;

    if (_data->useAsciiMesh) {
        pylith::meshio::MeshIOAscii iohandler;
        iohandler.setFilename(_data->meshFilename);
        iohandler.read(_mesh);assert(_mesh);
    } else {
        if (_data->meshOptions) {
            err = PetscOptionsInsertString(NULL, _data->meshOptions);PYLITH_CHECK_ERROR(err);
        } // if
        pylith::meshio::MeshIOPetsc iohandler;
        iohandler.setFilename(_data->meshFilename);
        iohandler.read(_mesh);assert(_mesh);
    } // if/else

    assert(pylith::topology::MeshOps::getNumCells(*_mesh) > 0);
    assert(pylith::topology::MeshOps::getNumVertices(*_mesh) > 0);

    // Set up coordinates.
    _mesh->setCoordSys(&_data->cs);
    pylith::topology::MeshOps::nondimensionalize(_mesh, _data->scales);

    // Set up material
    _data->material.setBulkRheology(&_data->rheology);
    _data->material.setAuxiliaryFieldDB(&_data->auxDB);

    for (size_t i = 0; i < _data->numAuxSubfields; ++i) {
        const pylith::topology::FieldBase::Discretization& info = _data->auxDiscretizations[i];
        _data->material.setAuxiliarySubfieldDiscretization(_data->auxSubfields[i], info.basisOrder, info.quadOrder,
                                                           _data->spaceDim, pylith::topology::FieldBase::DEFAULT_BASIS,
                                                           info.feSpace, info.isBasisContinuous);
    } // for

    // Set up source
    _data->source.setSourceTimeFunction(&_data->sourceTimeFunction);
    _data->source.setAuxiliaryFieldDB(&_data->sourceAuxDB);

    // Set source point coordinates
    if (_data->numSourcePoints > 0 && _data->sourcePointCoords && _data->sourcePointNames) {
        _data->source.setPoints(_data->sourcePointCoords, _data->numSourcePoints, _data->spaceDim,
                                _data->sourcePointNames, _data->numSourcePoints);
    } // if

    for (size_t i = 0; i < _data->numSourceAuxSubfields; ++i) {
        const pylith::topology::FieldBase::Discretization& info = _data->sourceAuxDiscretizations[i];
        _data->source.setAuxiliarySubfieldDiscretization(_data->sourceAuxSubfields[i], info.basisOrder, info.quadOrder,
                                                         _data->spaceDim, pylith::topology::FieldBase::DEFAULT_BASIS,
                                                         info.feSpace, info.isBasisContinuous);
    } // for

    // Set up problem.
    assert(_problem);
    _problem->setScales(_data->scales);
    pylith::materials::Material* materials[1] = { &_data->material };
    _problem->setMaterials(materials, 1);
    _problem->setBoundaryConditions(_data->bcs.data(), _data->bcs.size());
    pylith::sources::Source* sources[1] = { &_data->source };
    _problem->setSources(sources, 1);
    _problem->setStartTime(_data->t);
    _problem->setEndTime(_data->t+_data->dt);
    _problem->setInitialTimeStep(_data->dt);
    _problem->setFormulation(_data->formulation);

    // Set up solution field.
    assert(!_solution);
    _solution = new pylith::topology::Field(*_mesh);assert(_solution);
    _solution->setLabel("solution");
    pylith::problems::SolutionFactory factory(*_solution, _data->scales);
    factory.addDisplacement(_data->solnDiscretizations[0]);
    assert(pylith::problems::Physics::DYNAMIC == _data->formulation);
    assert(2 == _data->numSolnSubfields);
    factory.addVelocity(_data->solnDiscretizations[1]);
    _problem->setSolution(_solution);

    pylith::testing::MMSTest::_initialize();

    PYLITH_METHOD_END;
} // _initialize


// ------------------------------------------------------------------------------------------------
// Set functions for computing the exact solution and its time derivative.
void
pylith::TestMomentTensorSource::_setExactSolution(void) {
    assert(_data->exactSolnFns);

    const pylith::topology::Field* solution = _problem->getSolution();assert(solution);

    PetscErrorCode err = PETSC_SUCCESS;
    PetscDS ds = NULL;
    err = DMGetDS(solution->getDM(), &ds);PYLITH_CHECK_ERROR(err);
    for (size_t i = 0; i < _data->numSolnSubfields; ++i) {
        err = PetscDSSetExactSolution(ds, i, _data->exactSolnFns[i], NULL);PYLITH_CHECK_ERROR(err);
        if (_data->exactSolnDotFns) {
            err = PetscDSSetExactSolutionTimeDerivative(ds, i, _data->exactSolnDotFns[i], NULL);PYLITH_CHECK_ERROR(err);
        } // if
    } // for
} // _setExactSolution


// ------------------------------------------------------------------------------------------------
// Constructor
pylith::TestMomentTensorSource_Data::TestMomentTensorSource_Data(void) :
    spaceDim(2),
    meshFilename(NULL),
    meshOptions(NULL),
    boundaryLabel(NULL),
    useAsciiMesh(true),

    jacobianConvergenceRate(1.0),
    tolerance(1.0e-9),
    isJacobianLinear(true),
    allowZeroResidual(false),

    t(0.0),
    dt(0.05),
    formulation(pylith::problems::Physics::DYNAMIC),

    numSolnSubfields(0),
    solnDiscretizations(NULL),

    numAuxSubfields(0),
    auxSubfields(NULL),
    auxDiscretizations(NULL),

    numSourceAuxSubfields(0),
    sourceAuxSubfields(NULL),
    sourceAuxDiscretizations(NULL),

    numSourcePoints(0),
    sourcePointCoords(NULL),
    sourcePointNames(NULL) {
    auxDB.setDescription("material auxiliary field spatial database");
    sourceAuxDB.setDescription("source auxiliary field spatial database");
    cs.setSpaceDim(spaceDim);

    const double lengthScale = 8.0e+3;
    const double velocityScale = 3000.0;
    pylith::scales::ElasticityScales::setDynamicElasticity(&scales, lengthScale, velocityScale);
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::TestMomentTensorSource_Data::~TestMomentTensorSource_Data(void) {
    for (size_t i = 0; i < bcs.size(); ++i) {
        delete bcs[i];bcs[i] = NULL;
    } // for

} // destructor


// End of file
