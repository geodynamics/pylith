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

#include "TestFaultKin.hh" // Implementation of class methods

#include "pylith/bc/BoundaryCondition.hh" // HASA BoundaryCondition
#include "pylith/problems/TimeDependent.hh" // USES TimeDependent
#include "pylith/faults/FaultCohesiveKin.hh" // USES FaultCohesiveKin
#include "pylith/faults/KinSrc.hh" // USES KinSrc
#include "pylith/materials/Elasticity.hh" // USES Elasticity
#include "pylith/materials/RheologyElasticity.hh" // USES RheologyElasticity
#include "pylith/materials/Query.hh" // USES Query

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps::nondimensionalize()
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // USES FieldQuery
#include "pylith/feassemble/AuxiliaryFactory.hh" // USES AuxiliaryFactory
#include "pylith/problems/SolutionFactory.hh" // USES SolutionFactory
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/meshio/MeshIOPetsc.hh" // USES MeshIOPetsc
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // pythia::journal

#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

// ------------------------------------------------------------------------------------------------
// Constuctor.
pylith::TestFaultKin::TestFaultKin(TestFaultKin_Data* data) :
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
pylith::TestFaultKin::~TestFaultKin(void) {
    delete _data;_data = nullptr;
} // destructor


// ------------------------------------------------------------------------------------------------
// Initialize objects for test.
void
pylith::TestFaultKin::_initialize(void) {
    PYLITH_METHOD_BEGIN;
    assert(_mesh);
    assert(_data);

    PetscErrorCode err = 0;

    if (_data->useAsciiMesh) {
        pylith::meshio::MeshIOAscii iohandler;
        iohandler.setFilename(_data->meshFilename);
        iohandler.read(_mesh);assert(_mesh);
    } else {
        if (_data->meshOptions) {
            err = PetscOptionsInsertString(nullptr, _data->meshOptions);PYLITH_CHECK_ERROR(err);
        } // if
        pylith::meshio::MeshIOPetsc iohandler;
        iohandler.setFilename(_data->meshFilename);
        iohandler.read(_mesh);assert(_mesh);
    } // if/else

    assert(pylith::topology::MeshOps::getNumCells(*_mesh) > 0);
    assert(pylith::topology::MeshOps::getNumVertices(*_mesh) > 0);

    // Set up coordinates.
    _mesh->setCoordSys(&_data->cs);
    pylith::topology::MeshOps::nondimensionalize(_mesh, _data->normalizer);

    // Set up materials
    for (size_t iMat = 0; iMat < _data->materials.size(); ++iMat) {
        assert(_data->materials[iMat]);
        _data->materials[iMat]->setAuxiliaryFieldDB(&_data->matAuxDB);

        for (size_t iSubfield = 0; iSubfield < _data->matNumAuxSubfields; ++iSubfield) {
            const pylith::topology::FieldBase::Discretization& info = _data->matAuxDiscretizations[iSubfield];
            _data->materials[iMat]->setAuxiliarySubfieldDiscretization(_data->matAuxSubfields[iSubfield], info.basisOrder, info.quadOrder,
                                                                       _data->spaceDim, info.cellBasis, info.feSpace, info.isBasisContinuous);
        } // for
    } // for

    // Set up faults
    for (size_t iFault = 0; iFault < _data->faults.size(); ++iFault) {
        assert(_data->faults[iFault]);
        _data->faults[iFault]->adjustTopology(_mesh);

        assert(_data->kinSrc);
        _data->kinSrc->auxFieldDB(&_data->faultAuxDB);
        for (size_t i = 0; i < _data->faultNumAuxSubfields; ++i) {
            const pylith::topology::FieldBase::Discretization& info = _data->faultAuxDiscretizations[i];
            _data->faults[iFault]->setAuxiliarySubfieldDiscretization(_data->faultAuxSubfields[i], info.basisOrder, info.quadOrder,
                                                                      _data->spaceDim-1, info.cellBasis, info.feSpace, info.isBasisContinuous);
        } // for
    } // for

    // Set up problem.
    assert(_problem);
    _problem->setNormalizer(_data->normalizer);
    _problem->setGravityField(_data->gravityField);
    _problem->setMaterials(_data->materials.data(), _data->materials.size());
    _problem->setInterfaces(_data->faults.data(), _data->faults.size());
    _problem->setBoundaryConditions(_data->bcs.data(), _data->bcs.size());
    _problem->setStartTime(_data->t);
    _problem->setEndTime(_data->t+_data->dt);
    _problem->setInitialTimeStep(_data->dt);
    _problem->setFormulation(_data->formulation);

    // Set up solution field.
    assert(!_solution);
    _solution = new pylith::topology::Field(*_mesh);assert(_solution);
    _solution->setLabel("solution");
    pylith::problems::SolutionFactory factory(*_solution, _data->normalizer);
    int iField = 0;
    factory.addDisplacement(_data->solnDiscretizations[iField++]);
    if (pylith::problems::Physics::QUASISTATIC == _data->formulation) {
        assert(1 == _data->numSolnSubfieldsDomain);
    } else {
        assert(pylith::problems::Physics::DYNAMIC == _data->formulation);
        assert(2 == _data->numSolnSubfieldsDomain);
        factory.addVelocity(_data->solnDiscretizations[1]);
    } // if/else
    factory.addLagrangeMultiplierFault(_data->solnDiscretizations[iField++]);
    _problem->setSolution(_solution);

    pylith::testing::MMSTest::_initialize();

    PYLITH_METHOD_END;
} // _initialize


// ------------------------------------------------------------------------------------------------
// Set functions for computing the exact solution and its time derivative.
void
pylith::TestFaultKin::_setExactSolution(void) {
    assert(_data->exactSolnFns);

    const pylith::topology::Field* solution = _problem->getSolution();assert(solution);

    PetscErrorCode err = PETSC_SUCCESS;
    PetscDM dm = solution->getDM();
    PetscDS ds = nullptr;
    err = DMGetDS(dm, &ds);PYLITH_CHECK_ERROR(err);
    for (size_t i = 0; i < _data->numSolnSubfieldsDomain; ++i) {
        err = PetscDSSetExactSolution(ds, i, _data->exactSolnFns[i], dm);PYLITH_CHECK_ERROR(err);
        if (_data->exactSolnDotFns) {
            err = PetscDSSetExactSolutionTimeDerivative(ds, i, _data->exactSolnDotFns[i], dm);PYLITH_CHECK_ERROR(err);
        } // if
    } // for

    PetscDMLabel label = nullptr;
    PetscIS is = nullptr;
    PetscInt cohesiveCell = -1;
    err = DMGetLabel(dm, pylith::topology::Mesh::cells_label_name, &label);PYLITH_CHECK_ERROR(err);
    err = DMLabelGetStratumIS(label, _data->faults[0]->getCohesiveLabelValue(), &is);PYLITH_CHECK_ERROR(err);
    err = ISGetMinMax(is, &cohesiveCell, nullptr);PYLITH_CHECK_ERROR(err);
    err = ISDestroy(&is);PYLITH_CHECK_ERROR(err);
    err = DMGetCellDS(dm, cohesiveCell, &ds, nullptr);PYLITH_CHECK_ERROR(err);
    for (size_t i = 0; i < _data->numSolnSubfieldsDomain; ++i) {
        err = PetscDSSetExactSolution(ds, i, _data->exactSolnFns[i], nullptr);PYLITH_CHECK_ERROR(err);
        if (_data->exactSolnDotFns) {
            err = PetscDSSetExactSolutionTimeDerivative(ds, i, _data->exactSolnDotFns[i], nullptr);PYLITH_CHECK_ERROR(err);
        } // if
    } // for
    for (size_t i = 0; i < _data->numSolnSubfieldsFault; ++i) {
        const size_t iSoln = _data->numSolnSubfieldsDomain + i;
        err = PetscDSSetExactSolution(ds, iSoln, _data->exactSolnFns[iSoln], nullptr);PYLITH_CHECK_ERROR(err);
    } // for
} // _setExactSolution


// ------------------------------------------------------------------------------------------------
// Constructor
pylith::TestFaultKin_Data::TestFaultKin_Data(void) :
    spaceDim(2),
    meshFilename(nullptr),
    meshOptions(nullptr),
    boundaryLabel(nullptr),
    useAsciiMesh(true),

    jacobianConvergenceRate(1.0),
    tolerance(1.0e-9),
    isJacobianLinear(true),
    allowZeroResidual(false),

    t(0.0),
    dt(0.05),
    formulation(pylith::problems::Physics::QUASISTATIC),

    gravityField(nullptr),

    numSolnSubfieldsDomain(0),
    numSolnSubfieldsFault(0),
    solnDiscretizations(nullptr),

    matNumAuxSubfields(0),
    matAuxSubfields(nullptr),
    matAuxDiscretizations(nullptr),

    faultNumAuxSubfields(0),
    faultAuxSubfields(nullptr),
    faultAuxDiscretizations(nullptr),
    kinSrc(nullptr) {
    matAuxDB.setDescription("material auxiliary field spatial database");
    faultAuxDB.setDescription("fault auxiliary field spatial database");

    cs.setSpaceDim(spaceDim);
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::TestFaultKin_Data::~TestFaultKin_Data(void) {
    for (size_t i = 0; i < materials.size(); ++i) {
        delete materials[i];materials[i] = nullptr;
    } // for
    for (size_t i = 0; i < faults.size(); ++i) {
        delete faults[i];faults[i] = nullptr;
    } // for
    for (size_t i = 0; i < bcs.size(); ++i) {
        delete bcs[i];bcs[i] = nullptr;
    } // for
    delete kinSrc;kinSrc = nullptr;

    delete gravityField;gravityField = nullptr;
} // destructor


// End of file
