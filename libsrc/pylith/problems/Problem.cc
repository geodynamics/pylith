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

#include "pylith/problems/Problem.hh" // implementation of class methods

#include "pylith/feassemble/IntegrationData.hh" // HOLDSA IntegrationData
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // HASA Field
#include "pylith/topology/FieldOps.hh" // USES FieldOps

#include "pylith/materials/Material.hh" // USES Material
#include "pylith/faults/FaultCohesive.hh" // USES FaultCohesive
#include "pylith/bc/BoundaryCondition.hh" // USES BoundaryCondition
#include "pylith/feassemble/Integrator.hh" // USES Integrator
#include "pylith/feassemble/IntegratorDomain.hh" // USES IntegratorDomain
#include "pylith/feassemble/IntegratorInterface.hh" // USES IntegratorInterface
#include "pylith/feassemble/Constraint.hh" // USES Constraint
#include "pylith/problems/ObserversSoln.hh" // USES ObserversSoln
#include "pylith/topology/MeshOps.hh" // USES MeshOps
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor::optimizeClosure()

#include "pylith/scales/Scales.hh" // USES Scales
#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField

#include "pylith/utils/EventLogger.hh" // HASA EventLogger
#include "pylith/utils/error.hh" // USES PylithCallPetsc
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
#include "pylith/utils/Exceptions.hh" // USES Exception

#include <cassert> // USES assert()
#include <typeinfo> // USES typeid()

// ------------------------------------------------------------------------------------------------
namespace pylith {
    namespace problems {
        class _Problem {
public:

            /** Create null space for solution subfield.
             *
             * @param[inout] solution Solution field.
             * @param[in] subfieldName Name of solution subfield with null space.
             */
            static
            void createNullSpace(const pylith::topology::Field* solution,
                                 const char* subfieldName);

            /** Set data needed to integrate domain faces on interior interface.
             *
             * @param[inout] solution Solution field.
             * @param[in] integrators Array of integrators for problem.
             */
            static
            void setInterfaceData(const pylith::topology::Field* solution,
                                  std::vector<pylith::feassemble::Integrator*>& integrators);

            /** Get subset of integrators matching template type T.
             *
             * @param[in] Array of integrators for problem
             * @returns Array of integrators of type T.
             */
            template<class T> static std::vector<T*> subset(const std::vector<pylith::feassemble::Integrator*>& integrators);

            // Logging events
            class Events {
public:

                static
                void init(void);

                static pylith::utils::EventLogger logger;
                static PylithInt setSolution;
                static PylithInt preinitialize;
                static PylithInt verifyConfiguration;
                static PylithInt initialize;
                static PylithInt checkMaterials;
                static PylithInt createIntegrators;
                static PylithInt createConstraints;
                static PylithInt setupSolution;
            };

        };
    }
}
pylith::utils::EventLogger pylith::problems::_Problem::Events::logger;
PylithInt pylith::problems::_Problem::Events::setSolution;
PylithInt pylith::problems::_Problem::Events::preinitialize;
PylithInt pylith::problems::_Problem::Events::verifyConfiguration;
PylithInt pylith::problems::_Problem::Events::initialize;
PylithInt pylith::problems::_Problem::Events::checkMaterials;
PylithInt pylith::problems::_Problem::Events::createIntegrators;
PylithInt pylith::problems::_Problem::Events::createConstraints;
PylithInt pylith::problems::_Problem::Events::setupSolution;

// ------------------------------------------------------------------------------------------------
void
pylith::problems::_Problem::Events::init(void) {
    logger.setClassName("Problem");
    logger.initialize();
    setSolution = logger.registerEvent("PL:Problem:setSolution");
    preinitialize = logger.registerEvent("PL:Problem:preinitialize");
    verifyConfiguration = logger.registerEvent("PL:Problem:verifyConfiguration");
    initialize = logger.registerEvent("PL:Problem:initialize");
    checkMaterials = logger.registerEvent("PL:Problem:checkMaterials");
    createIntegrators = logger.registerEvent("PL:Problem:createIntegrators");
    createConstraints = logger.registerEvent("PL:Problem:createConstraints");
    setupSolution = logger.registerEvent("PL:Problem:setupSolution");
} // init


// ------------------------------------------------------------------------------------------------
// Constructor
pylith::problems::Problem::Problem() :
    _integrationData(new pylith::feassemble::IntegrationData),
    _scales(NULL),
    _gravityField(NULL),
    _observers(new pylith::problems::ObserversSoln),
    _formulation(pylith::problems::Physics::QUASISTATIC),
    _solverType(LINEAR),
    _petscDefaults(pylith::utils::PetscDefaults::SOLVER | pylith::utils::PetscDefaults::TESTING) {
    _Problem::Events::init();
}


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::problems::Problem::~Problem(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::problems::Problem::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    const size_t numIntegrators = _integrators.size();
    for (size_t i = 0; i < numIntegrators; ++i) {
        delete _integrators[i];_integrators[i] = NULL;
    } // for

    const size_t numConstraints = _constraints.size();
    for (size_t i = 0; i < numConstraints; ++i) {
        delete _constraints[i];_constraints[i] = NULL;
    } // for

    delete _integrationData;_integrationData = NULL;
    delete _scales;_scales = NULL;
    _gravityField = NULL; // Held by Python. :KLUDGE: :TODO: Use shared pointer.
    delete _observers;_observers = NULL;

    pylith::topology::FieldOps::deallocate();

    PYLITH_METHOD_END;
} // deallocate


// ------------------------------------------------------------------------------------------------
// Set formulation for solving equation.
void
pylith::problems::Problem::setFormulation(const pylith::problems::Physics::FormulationEnum value) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG(pylith::journal::debug_parameters, "formulation="<<value);

    _formulation = value;

    PYLITH_METHOD_END;
} // setFormulation


// ------------------------------------------------------------------------------------------------
// Get formulation for solving equation.
pylith::problems::Physics::FormulationEnum
pylith::problems::Problem::getFormulation(void) const {
    return _formulation;
} // getFormulation


// ------------------------------------------------------------------------------------------------
// Set problem type.
void
pylith::problems::Problem::setSolverType(const SolverTypeEnum value) {
    PYLITH_COMPONENT_DEBUG(pylith::journal::debug_parameters, "solver="<<value);

    _solverType = value;
} // setSolverType


// ------------------------------------------------------------------------------------------------
// Get problem type.
pylith::problems::Problem::SolverTypeEnum
pylith::problems::Problem::getSolverType(void) const {
    return _solverType;
} // getSolverType


// ------------------------------------------------------------------------------------------------
// Specify whether to set defaults for PETSc solver appropriate for problem.
void
pylith::problems::Problem::setPetscDefaults(const int flags) {
    _petscDefaults = flags;
    const bool hasFault = _interfaces.size() > 0;
    pylith::utils::PetscDefaults::set(_materials[0], hasFault, _petscDefaults);
} // setPetscDefaults


// ------------------------------------------------------------------------------------------------
// Set manager of scales used to nondimensionalize problem.
void
pylith::problems::Problem::setScales(const pylith::scales::Scales& dim) {
    PYLITH_COMPONENT_DEBUG(pylith::journal::application_flow, "Setting scales for nondimensionalization.");

    if (!_scales) {
        _scales = new pylith::scales::Scales(scales);
    } else {
        *_scales = scales;
    } // if/else
} // setScales


// ------------------------------------------------------------------------------------------------
// Get manager of scales used to nondimensionalize problem.
const pylith::scales::Scales&
pylith::problems::Problem::getScales(void) const {
    assert(_scales);
    return *_scales;
}


// ------------------------------------------------------------------------------------------------
// Set gravity field.
void
pylith::problems::Problem::setGravityField(spatialdata::spatialdb::GravityField* const g) {
    PYLITH_COMPONENT_DEBUG(pylith::journal::application_flow, "Setting gravity field.");

    _gravityField = g;
} // setGravityField


// ----------------------------------------------------------------------
// Register observer to receive notifications.
void
pylith::problems::Problem::registerObserver(pylith::problems::ObserverSoln* observer) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG(pylith::journal::application_flow, "Register observer.");

    assert(_observers);
    assert(_scales);
    _observers->registerObserver(observer);
    _observers->setTimeScale(_scales->getTimeScale());

    PYLITH_METHOD_END;
} // registerObserver


// ----------------------------------------------------------------------
// Remove observer from receiving notifications.
void
pylith::problems::Problem::removeObserver(pylith::problems::ObserverSoln* observer) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG(pylith::journal::application_flow, "Removing observerr.");

    assert(_observers);
    _observers->removeObserver(observer);

    PYLITH_METHOD_END;
} // removeObserver


// ------------------------------------------------------------------------------------------------
// Set solution field.
void
pylith::problems::Problem::setSolution(pylith::topology::Field* field) {
    PYLITH_COMPONENT_DEBUG(pylith::journal::application_flow, "Setting solution.");
    _Problem::Events::logger.eventBegin(_Problem::Events::setSolution);

    assert(_integrationData);
    _integrationData->setField(pylith::feassemble::IntegrationData::solution, field);

    _Problem::Events::logger.eventEnd(_Problem::Events::setSolution);
} // setSolution


// ------------------------------------------------------------------------------------------------
// Get solution field.
const pylith::topology::Field*
pylith::problems::Problem::getSolution(void) const {
    PYLITH_METHOD_BEGIN;

    assert(_integrationData);
    pylith::topology::Field* solution = NULL;
    if (_integrationData->hasField(pylith::feassemble::IntegrationData::solution)) {
        solution = _integrationData->getField(pylith::feassemble::IntegrationData::solution);
    } // if

    PYLITH_METHOD_RETURN(solution);
} // getSolution


// ------------------------------------------------------------------------------------------------
// Get time derivative solution field.
const pylith::topology::Field*
pylith::problems::Problem::getSolutionDot(void) const {
    assert(_integrationData);
    pylith::topology::Field* solutionDot = NULL;
    if (_integrationData->hasField(pylith::feassemble::IntegrationData::solution_dot)) {
        solutionDot = _integrationData->getField(pylith::feassemble::IntegrationData::solution_dot);
    } // if

    PYLITH_METHOD_RETURN(solutionDot);
} // getSolutionDot


// ------------------------------------------------------------------------------------------------
// Set materials.
void
pylith::problems::Problem::setMaterials(pylith::materials::Material* materials[],
                                        const int numMaterials) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG(pylith::journal::debug_parameters, "Setting "<<numMaterials<<" materials.");

    assert( (!materials && 0 == numMaterials) || (materials && 0 < numMaterials) );

    _materials.resize(numMaterials);
    for (int i = 0; i < numMaterials; ++i) {
        _materials[i] = materials[i];
    } // for

    PYLITH_METHOD_END;
} // setMaterials


// ------------------------------------------------------------------------------------------------
// Get materials.
const std::vector<pylith::materials::Material*>&
pylith::problems::Problem::getMaterials(void) const {
    return _materials;
} // getMaterials


// ------------------------------------------------------------------------------------------------
// Set boundary conditions.
void
pylith::problems::Problem::setBoundaryConditions(pylith::bc::BoundaryCondition* bc[],
                                                 const int numBC) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG(pylith::journal::debug_parameters, "Setting "<<numBC<<" boundary conditions.");

    assert( (!bc && 0 == numBC) || (bc && 0 < numBC) );

    _bc.resize(numBC);
    for (int i = 0; i < numBC; ++i) {
        _bc[i] = bc[i];
    } // for

    PYLITH_METHOD_END;
} // setBoundaryConditions


// ------------------------------------------------------------------------------------------------
// Get boundary conditions.
const std::vector<pylith::bc::BoundaryCondition*>&
pylith::problems::Problem::getBoundaryConditions(void) const {
    return _bc;
} // getBoundaryConditions


// ------------------------------------------------------------------------------------------------
// Set fault interfaces.
void
pylith::problems::Problem::setInterfaces(pylith::faults::FaultCohesive* interfaces[],
                                         const int numInterfaces) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG(pylith::journal::debug_parameters, "Setting "<<numInterfaces<<" interior interfaces.");

    assert( (!interfaces && 0 == numInterfaces) || (interfaces && 0 < numInterfaces) );

    _interfaces.resize(numInterfaces);
    for (int i = 0; i < numInterfaces; ++i) {
        _interfaces[i] = interfaces[i];
    } // for

    PYLITH_METHOD_END;
} // setInterfaces


// ------------------------------------------------------------------------------------------------
// Get fault interfaces.
const std::vector<pylith::faults::FaultCohesive*>&
pylith::problems::Problem::getInterfaces(void) const {
    return _interfaces;
} // getInterfaces


// ----------------------------------------------------------------------
// Do minimal initialization.
void
pylith::problems::Problem::preinitialize(const pylith::topology::Mesh& mesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG(pylith::journal::application_flow, "Preinitializing problem.");
    _Problem::Events::logger.eventBegin(_Problem::Events::preinitialize);

    assert(_scales);

    const size_t numMaterials = _materials.size();
    for (size_t i = 0; i < numMaterials; ++i) {
        assert(_materials[i]);
        _materials[i]->setScales(*_scales);
        _materials[i]->setGravityField(_gravityField);
        _materials[i]->setFormulation(_formulation);
    } // for

    const size_t numInterfaces = _interfaces.size();
    for (size_t i = 0; i < numInterfaces; ++i) {
        assert(_interfaces[i]);
        _interfaces[i]->setScales(*_scales);
        _interfaces[i]->setFormulation(_formulation);
    } // for

    const size_t numBC = _bc.size();
    for (size_t i = 0; i < numBC; ++i) {
        assert(_bc[i]);
        _bc[i]->setScales(*_scales);
        _bc[i]->setFormulation(_formulation);
    } // for

    _Problem::Events::logger.eventEnd(_Problem::Events::preinitialize);
    PYLITH_METHOD_END;
} // preinitialize


// ----------------------------------------------------------------------
// Verify configuration.
void
pylith::problems::Problem::verifyConfiguration(void) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG(pylith::journal::application_flow, "Verifying configuration.");
    _Problem::Events::logger.eventBegin(_Problem::Events::verifyConfiguration);

    assert(_integrationData);
    const pylith::topology::Field* solution = _integrationData->getField("solution");
    assert(solution);

    // Check to make sure materials are compatible with the solution.
    const size_t numMaterials = _materials.size();
    for (size_t i = 0; i < numMaterials; ++i) {
        assert(_materials[i]);
        _materials[i]->verifyConfiguration(*solution);
    } // for

    // Check to make sure interfaces are compatible with the solution.
    const size_t numInterfaces = _interfaces.size();
    for (size_t i = 0; i < numInterfaces; ++i) {
        assert(_interfaces[i]);
        _interfaces[i]->verifyConfiguration(*solution);
    } // for

    // Check to make sure boundary conditions are compatible with the solution.
    const size_t numBC = _bc.size();
    for (size_t i = 0; i < numBC; ++i) {
        assert(_bc[i]);
        _bc[i]->verifyConfiguration(*solution);
    } // for

    _checkMaterialLabels();

    assert(_observers);
    _observers->verifyObservers(*solution);

    _Problem::Events::logger.eventEnd(_Problem::Events::verifyConfiguration);
    PYLITH_METHOD_END;
} // verifyConfiguration


// ----------------------------------------------------------------------
// Initialize.
void
pylith::problems::Problem::initialize(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG(pylith::journal::application_flow, "Initializing problem.");
    _Problem::Events::logger.eventBegin(_Problem::Events::initialize);

    assert(_integrationData);
    pylith::topology::Field* solution = _integrationData->getField("solution");
    assert(solution);

    // Initialize solution field.
    _setupSolution();
    pylith::topology::CoordsVisitor::optimizeClosure(solution->getDM());

    // Initialize integrators.
    _createIntegrators();
    const size_t numIntegrators = _integrators.size();
    for (size_t i = 0; i < numIntegrators; ++i) {
        assert(_integrators[i]);
        _integrators[i]->initialize(*solution);
    } // for

    // Initialize constraints.
    _createConstraints();
    const size_t numConstraints = _constraints.size();
    for (size_t i = 0; i < numConstraints; ++i) {
        assert(_constraints[i]);
        _constraints[i]->initialize(*solution);
    } // for

    solution->allocate();
    solution->createGlobalVector();
    solution->createOutputVector();

    switch (_formulation) {
    case pylith::problems::Physics::DYNAMIC:
    case pylith::problems::Physics::DYNAMIC_IMEX:
        break;
    case pylith::problems::Physics::QUASISTATIC:
        _Problem::createNullSpace(solution, "displacement");
        break;
    default:
        PYLITH_COMPONENT_FIREWALL(pylith::InternalLogicError, pylith::journal::logic, "Unknown formulation '"<<_formulation<<".");
    } // switch
    _Problem::setInterfaceData(solution, _integrators);

    pythia::journal::debug_t debug(pylith::journal::solution);
    if (debug.state()) {
        PYLITH_COMPONENT_DEBUG(pylith::journal::solution, "Displaying solution field layout");
        solution->view("Solution field", pylith::topology::Field::VIEW_LAYOUT);
    } // if

    _Problem::Events::logger.eventEnd(_Problem::Events::initialize);
    PYLITH_METHOD_END;
} // initialize


// ------------------------------------------------------------------------------------------------
// Check material and interface ids.
void
pylith::problems::Problem::_checkMaterialLabels(void) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG(pylith::journal::application_flow, "Checking material labels.");
    _Problem::Events::logger.eventBegin(_Problem::Events::checkMaterials);

    const size_t numMaterials = _materials.size();
    const size_t numInterfaces = _interfaces.size();

    pylith::int_array labelValues(numMaterials + numInterfaces);
    size_t count = 0;
    for (size_t i = 0; i < numMaterials; ++i) {
        assert(_materials[i]);
        labelValues[count++] = _materials[i]->getLabelValue();
    } // for
    for (size_t i = 0; i < numInterfaces; ++i) {
        assert(_interfaces[i]);
        labelValues[count++] = _interfaces[i]->getCohesiveLabelValue();
    } // for

    assert(_integrationData);
    const pylith::topology::Field* solution = _integrationData->getField("solution");
    assert(solution);
    pylith::topology::MeshOps::checkMaterialLabels(solution->getMesh(), labelValues);

    _Problem::Events::logger.eventEnd(_Problem::Events::checkMaterials);
    PYLITH_METHOD_END;
} // _checkMaterialLabels


// ------------------------------------------------------------------------------------------------
// Create array of integrators from materials, interfaces, and boundary conditions.
void
pylith::problems::Problem::_createIntegrators(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG(pylith::journal::application_flow, "Creating integrators.");
    _Problem::Events::logger.eventBegin(_Problem::Events::createIntegrators);

    const size_t numMaterials = _materials.size();
    const size_t numInterfaces = _interfaces.size();
    const size_t numBC = _bc.size();

    const size_t maxSize = numMaterials + numInterfaces + numBC;
    _integrators.resize(maxSize);
    size_t count = 0;

    assert(_integrationData);
    const pylith::topology::Field* solution = _integrationData->getField("solution");
    assert(solution);

    for (size_t i = 0; i < numMaterials; ++i) {
        assert(_materials[i]);
        pylith::feassemble::Integrator* integrator = _materials[i]->createIntegrator(*solution);
        assert(count < maxSize);
        if (integrator) { _integrators[count++] = integrator;}
    } // for

    for (size_t i = 0; i < numInterfaces; ++i) {
        assert(_interfaces[i]);
        pylith::feassemble::Integrator* integrator = _interfaces[i]->createIntegrator(*solution, _materials);
        assert(count < maxSize);
        if (integrator) { _integrators[count++] = integrator;}
    } // for

    // Check to make sure boundary conditions are compatible with the solution.
    for (size_t i = 0; i < numBC; ++i) {
        assert(_bc[i]);
        pylith::feassemble::Integrator* integrator = _bc[i]->createIntegrator(*solution);
        assert(count < maxSize);
        if (integrator) { _integrators[count++] = integrator;}
    } // for

    _integrators.resize(count);

    _Problem::Events::logger.eventEnd(_Problem::Events::createIntegrators);
    PYLITH_METHOD_END;
} // _createIntegrators


// ------------------------------------------------------------------------------------------------
// Create array of constraints from materials, interfaces, and boundary conditions.
void
pylith::problems::Problem::_createConstraints(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG(pylith::journal::application_flow, "Creating constraints.");
    _Problem::Events::logger.eventBegin(_Problem::Events::createConstraints);

    const size_t numMaterials = _materials.size();
    const size_t numInterfaces = _interfaces.size();
    const size_t numBC = _bc.size();

    assert(_integrationData);
    const pylith::topology::Field* solution = _integrationData->getField("solution");
    assert(solution);

    _constraints.resize(0); // insure we start with an empty array.

    for (size_t i = 0; i < numMaterials; ++i) {
        assert(_materials[i]);
        std::vector<pylith::feassemble::Constraint*> constraints = _materials[i]->createConstraints(*solution);
        _constraints.insert(_constraints.end(), constraints.begin(), constraints.end());

    } // for

    for (size_t i = 0; i < numInterfaces; ++i) {
        assert(_interfaces[i]);
        std::vector<pylith::feassemble::Constraint*> constraints = _interfaces[i]->createConstraints(*solution);
        _constraints.insert(_constraints.end(), constraints.begin(), constraints.end());

    } // for

    for (size_t i = 0; i < numBC; ++i) {
        assert(_bc[i]);
        std::vector<pylith::feassemble::Constraint*> constraints = _bc[i]->createConstraints(*solution);
        _constraints.insert(_constraints.end(), constraints.begin(), constraints.end());

    } // for

    _Problem::Events::logger.eventEnd(_Problem::Events::createConstraints);
    PYLITH_METHOD_END;
} // _createConstraints


// ------------------------------------------------------------------------------------------------
// Setup solution subfields and discretization.
void
pylith::problems::Problem::_setupSolution(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG(pylith::journal::application_flow, "Setting up solution field.");
    _Problem::Events::logger.eventBegin(_Problem::Events::createConstraints);

    assert(_integrationData);
    pylith::topology::Field* solution = _integrationData->getField("solution");
    assert(solution);
    solution->subfieldsSetup();
    solution->createDiscretization();

    // Mark fault fields as implicit.
    const pylith::string_vector& subfieldNames = solution->getSubfieldNames();
    for (size_t i = 0; i < subfieldNames.size(); ++i) {
        const pylith::topology::Field::SubfieldInfo& subfieldInfo = solution->getSubfieldInfo(subfieldNames[i].c_str());
        if (subfieldInfo.fe.isFaultOnly) {
            PetscDS ds = NULL;
            PetscInt cStart = 0, cEnd = 0;
            PetscDM dmSoln = solution->getDM();assert(dmSoln);
            PylithCallPetsc(DMPlexGetHeightStratum(dmSoln, 0, &cStart, &cEnd));
            PetscInt cell = cStart;
            bool found = false;
            for (; cell < cEnd; ++cell) {
                if (pylith::topology::MeshOps::isCohesiveCell(dmSoln, cell)) {
                    found = true;
                    break;
                } // if
            } // for
            if (!found) {
                continue;
            } // if
            PylithCallPetsc(DMGetCellDS(dmSoln, cell, &ds, NULL));
            assert(ds);
            PylithCallPetsc(PetscDSSetImplicit(ds, subfieldInfo.index, PETSC_TRUE));
        } // if
    } // for

    _Problem::Events::logger.eventEnd(_Problem::Events::setupSolution);
    PYLITH_METHOD_END;
} // _setupSolution


// ------------------------------------------------------------------------------------------------
// Create null space for solution subfield.
void
pylith::problems::_Problem::createNullSpace(const pylith::topology::Field* solution,
                                            const char* subfieldName) {
    PYLITH_METHOD_BEGIN;
    assert(solution);

    const int spaceDim = solution->getSpaceDim();
    const PetscInt m = (spaceDim * (spaceDim + 1)) / 2;assert(m > 0 && m <= 6);

    PetscInt numDofUnconstrained = 0;
    PylithCallPetsc(PetscSectionGetConstrainedStorageSize(solution->getLocalSection(), &numDofUnconstrained));
    if (m > numDofUnconstrained) {
        PYLITH_METHOD_END;
    } // if

    const PetscDM dmSoln = solution->getDM();
    const pylith::topology::Field::SubfieldInfo info = solution->getSubfieldInfo(subfieldName);
    MatNullSpace nullSpace = NULL;
    PylithCallPetsc(DMPlexCreateRigidBody(dmSoln, info.index, &nullSpace));

    PetscObject field = NULL;
    PylithCallPetsc(DMGetField(dmSoln, info.index, NULL, &field));
    PylithCallPetsc(PetscObjectCompose(field, "nearnullspace", (PetscObject) nullSpace));
    PylithCallPetsc(MatNullSpaceDestroy(&nullSpace));

    PYLITH_METHOD_END;
} // createNullSpace


// ------------------------------------------------------------------------------------------------
// Set data needed to integrate domain faces on interior interface.
void
pylith::problems::_Problem::setInterfaceData(const pylith::topology::Field* solution,
                                             std::vector<pylith::feassemble::Integrator*>& integrators) {
    PYLITH_METHOD_BEGIN;

    const std::vector<pylith::feassemble::IntegratorDomain*>& integratorsDomain =
        subset<pylith::feassemble::IntegratorDomain>(integrators);
    const std::vector<pylith::feassemble::IntegratorInterface*>& integratorsInterface =
        subset<pylith::feassemble::IntegratorInterface>(integrators);

    for (size_t i = 0; i < integratorsDomain.size(); ++i) {
        integratorsDomain[i]->setInterfaceData(solution, integratorsInterface);
    } // for

    PYLITH_METHOD_END;
} // setInterfaceData


// ------------------------------------------------------------------------------------------------
// Get subset of integrators matching template type T.
template<class T>
std::vector<T*>
pylith::problems::_Problem::subset(const std::vector<pylith::feassemble::Integrator*>& integrators) {
    // Count number of matches
    size_t numFound = 0;
    for (size_t i = 0; i < integrators.size(); ++i) {
        if (dynamic_cast<T*>(integrators[i])) {
            numFound++;
        } // if/else
    } // for

    // Collect matches
    std::vector<T*> matches(numFound);
    size_t index = 0;
    for (size_t i = 0; i < integrators.size(); ++i) {
        T* integrator = dynamic_cast<T*>(integrators[i]);
        if (integrator) {
            matches[index++] = integrator;
        } // if/else
    } // for

    return matches;
} // subset


// End of file
