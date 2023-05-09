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

#include "FaultCohesiveDyn.hh" // implementation of object methods

#include "pylith/faults/TractionPerturbation.hh" // USES TractionPerturbation
#include "pylith/faults/AuxiliaryFactoryDynamic.hh" // USES AuxiliaryFactoryDynamic
#include "pylith/feassemble/IntegratorInterface.hh" // USES IntegratorInterface
#include "pylith/feassemble/InterfacePatches.hh" // USES InterfacePatches
#include "pylith/feassemble/ConstraintSimple.hh" // USES ConstraintSimple

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldOps.hh" // USES FieldOps::checkDiscretization()
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh

// #include "pylith/fekernels/FaultCohesiveDyn.hh" // USES FaultCohesiveDyn

#include "pylith/utils/EventLogger.hh" // USES EventLogger
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensionalizer
#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB

#include <cmath> // USES pow(), sqrt()
#include <strings.h> // USES strcasecmp()
#include <cstring> // USES strlen()
#include <cstdlib> // USES atoi()
#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error
#include <typeinfo> // USES typeid()

// ------------------------------------------------------------------------------------------------
typedef pylith::feassemble::IntegratorInterface::ResidualKernels ResidualKernels;
typedef pylith::feassemble::IntegratorInterface::JacobianKernels JacobianKernels;

// ------------------------------------------------------------------------------------------------
namespace pylith {
    namespace faults {
        class _FaultCohesiveDyn {
            // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////
public:

            static const char* pyreComponent;

        };
        const char* _FaultCohesiveDyn::pyreComponent = "faultcohesivedyn";

        // _FaultCohesiveDyn

    } // faults
} // pylith

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::faults::FaultCohesiveDyn::FaultCohesiveDyn(void) :
    _auxiliaryFactory(new pylith::faults::AuxiliaryFactoryDynamic),
    _rheology(NULL),
    _tractionVecPerturbation(NULL),
    _tractionVecTotal(NULL) {
    pylith::utils::PyreComponent::setName(_FaultCohesiveDyn::pyreComponent);
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::faults::FaultCohesiveDyn::~FaultCohesiveDyn(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::faults::FaultCohesiveDyn::deallocate(void) {
    FaultCohesive::deallocate();

    PetscErrorCode err = VecDestroy(&_tractionVecPerturbation);PYLITH_CHECK_ERROR(err);
    err = VecDestroy(&_tractionVecTotal);PYLITH_CHECK_ERROR(err);
    delete _auxiliaryFactory;_auxiliaryFactory = NULL;
    delete _rheology;_rheology = NULL;
    _perturbations.clear(); // :TODO: Use shared pointers for earthquake ruptures
} // deallocate


// ------------------------------------------------------------------------------------------------
// Set fault rheology.
void
pylith::faults::FaultCohesiveDyn::setFaultRheology(pylith::faults::FaultRheology* const rheology) {}


// ------------------------------------------------------------------------------------------------
// Get fault rheology.
pylith::faults::FaultRheology*
pylith::faults::FaultCohesiveDyn::getFaultRheology(void) const {}


// ------------------------------------------------------------------------------------------------
// Set kinematic earthquake ruptures.
void
pylith::faults::FaultCohesiveDyn::setTractionPerturbations(const char* const * names,
                                                           const int numNames,
                                                           TractionPerturbation** perturbations,
                                                           const int numPerturbations) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("setTractionPerturbations(names="<<names<<", numNames="<<numNames<<", perturbations="<<perturbations<<", numPerturbations="<<numPerturbations<<")");

    assert(numNames == numPerturbations);

    // :TODO: Use shared pointers for traction perturbations
    _perturbations.clear();
    for (int i = 0; i < numPerturbations; ++i) {
        if (!perturbations[i]) {
            std::ostringstream msg;
            msg << "Null traction perturbation object for fault traction perturbation '" << names[i] << "'.";
            throw std::runtime_error(msg.str());
        } // if
        _perturbations[std::string(names[i])] = perturbations[i];
    } // for

    PYLITH_METHOD_END;
} // setTractionPerturbations


// ------------------------------------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::faults::FaultCohesiveDyn::verifyConfiguration(const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("verifyConfiguration(solution="<<solution.getLabel()<<")");

    if (!solution.hasSubfield("lagrange_multiplier_fault")) {
        std::ostringstream msg;
        msg << "Cannot find 'lagrange_multiplier_fault' subfield in solution field for fault implementation in component '"
            << PyreComponent::getIdentifier() << "'.";
        throw std::runtime_error(msg.str());
    } // if

    switch (_formulation) {
    case QUASISTATIC:
        if (!solution.hasSubfield("displacement")) {
            std::ostringstream msg;
            msg << "Cannot find 'displacement' subfield in solution field for fault implementation in component '"
                << PyreComponent::getIdentifier() << "'.";
            throw std::runtime_error(msg.str());
        } // if
        break;
    case DYNAMIC_IMEX:
        if (!solution.hasSubfield("velocity")) {
            std::ostringstream msg;
            msg << "Cannot find 'velocity' subfield in solution field for fault implementation in component '"
                << PyreComponent::getIdentifier() << "'.";
            throw std::runtime_error(msg.str());
        } // if
        break;
    case DYNAMIC:
        PYLITH_COMPONENT_LOGICERROR("Fault implementation is incompatible with 'dynamic' formulation. Use 'dynamic_imex'.");
        break;
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown formulation for equations (" << _formulation << ").");
    } // switch

    PYLITH_METHOD_END;
} // verifyConfiguration


// ------------------------------------------------------------------------------------------------
// Create integrator and set kernels.
pylith::feassemble::Integrator*
pylith::faults::FaultCohesiveDyn::createIntegrator(const pylith::topology::Field& solution,
                                                   const std::vector<pylith::materials::Material*>& materials) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createIntegrator(solution="<<solution.getLabel()<<")");

    pylith::feassemble::IntegratorInterface* integrator = new pylith::feassemble::IntegratorInterface(this);assert(integrator);
    integrator->setLabelName(getCohesiveLabelName());
    integrator->setLabelValue(getCohesiveLabelValue());
    integrator->setSurfaceLabelName(getSurfaceLabelName());

    pylith::feassemble::InterfacePatches* patches =
        pylith::feassemble::InterfacePatches::createMaterialPairs(this, solution.getDM());
    integrator->setIntegrationPatches(patches);

    _setKernelsResidual(integrator, solution, materials);
    _setKernelsJacobian(integrator, solution, materials);
    _setKernelsUpdateStateVars(integrator, solution);

    PYLITH_METHOD_RETURN(integrator);
} // createIntegrator


// ------------------------------------------------------------------------------------------------
// Create auxiliary field.
pylith::topology::Field*
pylith::faults::FaultCohesiveDyn::createAuxiliaryField(const pylith::topology::Field& solution,
                                                       const pylith::topology::Mesh& domainMesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createAuxiliaryField(solution="<<solution.getLabel()<<", domainMesh=)"<<typeid(domainMesh).name()<<")");

    assert(_normalizer);

    pylith::topology::Field* auxiliaryField = new pylith::topology::Field(domainMesh);assert(auxiliaryField);
    auxiliaryField->setLabel("FaultCohesiveDyn auxiliary field");

    // Set default discretization of auxiliary subfields to match lagrange_multiplier_fault subfield in solution.
    assert(_auxiliaryFactory);
    const pylith::topology::FieldBase::Discretization& discretization = solution.getSubfieldInfo("lagrange_multiplier_fault").fe;
    const PylithInt dimension = -1;
    const bool isFaultOnly = false;
    _auxiliaryFactory->setSubfieldDiscretization("default", discretization.basisOrder, discretization.quadOrder, dimension,
                                                 isFaultOnly, discretization.cellBasis, discretization.feSpace,
                                                 discretization.isBasisContinuous);

    assert(_auxiliaryFactory);
    assert(_normalizer);
    _auxiliaryFactory->initialize(auxiliaryField, *_normalizer, solution.getSpaceDim());

    // :ATTENTION: The order for adding subfields must match the order of the auxiliary fields in the FE kernels.
    _auxiliaryFactory->addTractionRheology(); // 0
    _auxiliaryFactory->addTractionPerturbation(); // 1

    auxiliaryField->subfieldsSetup();
    auxiliaryField->createDiscretization();
    pylith::topology::FieldOps::checkDiscretization(solution, *auxiliaryField);
    auxiliaryField->allocate();
    auxiliaryField->createOutputVector();

    // We don't populate the auxiliary field via a spatial database, because they will be set from
    // the traction perturbations.

    // Initialize auxiliary fields for dynamic (spontaneous) ruptures.
    assert(auxiliaryField);
    const perturbations_type::const_iterator perturbationsEnd = _perturbations.end();
    for (perturbations_type::iterator p_iter = _perturbations.begin(); p_iter != perturbationsEnd; ++p_iter) {
        TractionPerturbation* perturbation = p_iter->second;
        assert(perturbation);
        perturbation->initialize(*auxiliaryField, *_normalizer, solution.getMesh().getCoordSys());
    } // for

    // Create local PETSc vector to hold current traction perturbation.
    PetscErrorCode err = 0;
    err = DMCreateLocalVector(auxiliaryField->getDM(), &_tractionVecPerturbation);PYLITH_CHECK_ERROR(err);
    err = DMCreateLocalVector(auxiliaryField->getDM(), &_tractionVecTotal);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_RETURN(auxiliaryField);
} // createAuxiliaryField


// ------------------------------------------------------------------------------------------------
// Update auxiliary fields at beginning of time step.
void
pylith::faults::FaultCohesiveDyn::updateAuxiliaryField(pylith::topology::Field* auxiliaryField,
                                                       const double t) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("updateAuxiliaryField(auxiliaryField="<<auxiliaryField<<", t="<<t<<")");

    assert(auxiliaryField);
    assert(_normalizer);

    this->_updateTractionPerturbation(auxiliaryField, t);

    PYLITH_METHOD_END;
} // updateAuxiliaryField


// ------------------------------------------------------------------------------------------------
// Get auxiliary factory associated with physics.
pylith::feassemble::AuxiliaryFactory*
pylith::faults::FaultCohesiveDyn::_getAuxiliaryFactory(void) {
    return _auxiliaryFactory;
} // _getAuxiliaryFactory


// ------------------------------------------------------------------------------------------------
// Update traction subfield in auxiliary field at beginning of time step.
void
pylith::faults::FaultCohesiveDyn::_updateTractionPerturbation(pylith::topology::Field* auxiliaryField,
                                                              const double t) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("updateTractionPerturbation(auxiliaryField="<<auxiliaryField<<", t="<<t<<")");

    assert(auxiliaryField);
    assert(_normalizer);

    // Update traction subfield at current time step
    PetscErrorCode err = VecSet(_tractionVecTotal, 0.0);PYLITH_CHECK_ERROR(err);
    const perturbations_type::const_iterator perturbationsEnd = _perturbations.end();
    for (perturbations_type::iterator p_iter = _perturbations.begin(); p_iter != perturbationsEnd; ++p_iter) {
        err = VecSet(_tractionVecPerturbation, 0.0);PYLITH_CHECK_ERROR(err);

        TractionPerturbation* perturbation = p_iter->second;assert(perturbation);
        perturbation->getTraction(_tractionVecPerturbation, auxiliaryField, t, _normalizer->getTimeScale());
        err = VecAYPX(_tractionVecTotal, 1.0, _tractionVecPerturbation);
    } // for

    // Transfer traction values from local PETSc vector to fault auxiliary field.
    PetscInt pStart = 0, pEnd = 0;
    err = PetscSectionGetChart(auxiliaryField->getLocalSection(), &pStart, &pEnd);PYLITH_CHECK_ERROR(err);
    pylith::topology::VecVisitorMesh auxiliaryVisitor(*auxiliaryField);
    PylithScalar* auxiliaryArray = auxiliaryVisitor.localArray();

    const PylithScalar* tractionArray = NULL;
    err = VecGetArrayRead(_tractionVecTotal, &tractionArray);PYLITH_CHECK_ERROR(err);

    for (PetscInt point = pStart, iSlip = 0; point < pEnd; ++point) {
        for (PetscInt iSubfield = 0; iSubfield < numSubfields; ++iSubfield) {
            const PetscInt subfieldIndex = subfieldIndices[iSubfield];
            const PetscInt subfieldDof = auxiliaryVisitor.sectionSubfieldDof(subfieldIndex, point);
            const PetscInt subfieldOff = auxiliaryVisitor.sectionSubfieldOffset(subfieldIndex, point);
            for (PetscInt iDof = 0; iDof < subfieldDof; ++iDof, ++iSlip) {
                auxiliaryArray[subfieldOff+iDof] = tractionArray[iSlip];
            } // for
        } // for
    } // for
    err = VecRestoreArrayRead(_tractionVecTotal, &tractionArray);PYLITH_CHECK_ERROR(err);

    pythia::journal::debug_t debug(pylith::utils::PyreComponent::getName());
    if (debug.state()) {
        auxiliaryField->view("Fault auxiliary field after updating fault tractions.");
    } // if

    PYLITH_METHOD_END;
} // _updateTractionPerturbation


// ------------------------------------------------------------------------------------------------
// Set kernels for residual.
void
pylith::faults::FaultCohesiveDyn::_setKernelsResidual(pylith::feassemble::IntegratorInterface* integrator,
                                                      const pylith::topology::Field& solution,
                                                      const std::vector<pylith::materials::Material*>& materials) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelsResidual (integrator = "<<integrator<<", solution = "<<solution.getLabel()<<") ");
    typedef pylith::feassemble::IntegratorInterface integrator_t;

    std::vector<ResidualKernels> kernels;
    switch (_formulation) {
    case pylith::problems::Physics::QUASISTATIC: {
        // Elasticity equation (displacement) for negative side of the fault.
        const PetscBdPointFunc f0u_neg = pylith::fekernels::FaultCohesiveDyn::f0u_neg;
        const PetscBdPointFunc f1u_neg = NULL;

        // Elasticity equation (displacement) for positive side of the fault.
        const PetscBdPointFunc f0u_pos = pylith::fekernels::FaultCohesiveDyn::f0u_pos;
        const PetscBdPointFunc f1u_pos = NULL;

        // Fault slip constraint equation.
        const PetscBdPointFunc f0l = pylith::fekernels::FaultCohesiveDyn::f0l_slip;
        const PetscBdPointFunc f1l = NULL;

        kernels.resize(3);
        kernels[0] = ResidualKernels("displacement ", integrator_t::LHS, integrator_t::NEGATIVE_FACE,
                                     f0u_neg, f1u_neg);
        kernels[1] = ResidualKernels("displacement ", integrator_t::LHS, integrator_t::POSITIVE_FACE,
                                     f0u_pos, f1u_pos);
        kernels[2] = ResidualKernels("lagrange_multiplier_fault ", integrator_t::LHS, integrator_t::FAULT_FACE,
                                     f0l, f1l);

        break;
    } // QUASISTATIC
    case pylith::problems::Physics::DYNAMIC_IMEX: {
        // Elasticity equation (displacement) for negative side of the fault.
        const PetscBdPointFunc g0v_neg = pylith::fekernels::FaultCohesiveDyn::f0u_neg;
        const PetscBdPointFunc g1v_neg = NULL;

        // Elasticity equation (displacement) for positive side of the fault.
        const PetscBdPointFunc g0v_pos = pylith::fekernels::FaultCohesiveDyn::f0u_pos;
        const PetscBdPointFunc g1v_pos = NULL;

        // Fault slip constraint equation.
        const PetscBdPointFunc f0l_slip = pylith::fekernels::FaultCohesiveDyn::f0l_slip;
        const PetscBdPointFunc f1l_slip = NULL;

        // Fault DAE equation.
        const PetscBdPointFunc f0l_dae = pylith::fekernels::FaultCohesiveDyn::f0l_slipAcc;
        const PetscBdPointFunc f1l_dae = NULL;

        kernels.resize(4);
        kernels[0] = ResidualKernels("velocity ", integrator_t::LHS, integrator_t::NEGATIVE_FACE, g0v_neg, g1v_neg);
        kernels[1] = ResidualKernels("velocity ", integrator_t::LHS, integrator_t::POSITIVE_FACE, g0v_pos, g1v_pos);
        kernels[2] = ResidualKernels("lagrange_multiplier_fault ", integrator_t::LHS, integrator_t::FAULT_FACE, f0l_slip, f1l_slip);
        kernels[3] = ResidualKernels("lagrange_multiplier_fault ", integrator_t::LHS_WEIGHTED, integrator_t::FAULT_FACE, f0l_dae, f1l_dae);
        break;
    } // DYNAMIC_IMEX
    case pylith::problems::Physics::DYNAMIC:
        PYLITH_COMPONENT_LOGICERROR("Fault implementation is incompatible with 'dynamic' formulation.Use 'dynamic_imex'.");
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown formulation '"<<_formulation<<"'.");
    } // switch

    assert(integrator);
    integrator->setKernels(kernels, solution, materials);

    PYLITH_METHOD_END;
} // _setKernelsResidual


// ------------------------------------------------------------------------------------------------
// Set kernels for Jacobian.
void
pylith::faults::FaultCohesiveDyn::_setKernelsJacobian(pylith::feassemble::IntegratorInterface* integrator,
                                                      const pylith::topology::Field& solution,
                                                      const std::vector<pylith::materials::Material*>& materials) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelsJacobian (integrator = "<<integrator<<", solution = "<<solution.getLabel()<<") ");
    typedef pylith::feassemble::IntegratorInterface integrator_t;

    std::vector<JacobianKernels> kernels;
    switch (_formulation) {
    case QUASISTATIC: {
        const PetscBdPointJac Jf0ul_neg = pylith::fekernels::FaultCohesiveDyn::Jf0ul_neg;
        const PetscBdPointJac Jf1ul_neg = NULL;
        const PetscBdPointJac Jf2ul_neg = NULL;
        const PetscBdPointJac Jf3ul_neg = NULL;

        const PetscBdPointJac Jf0ul_pos = pylith::fekernels::FaultCohesiveDyn::Jf0ul_pos;
        const PetscBdPointJac Jf1ul_pos = NULL;
        const PetscBdPointJac Jf2ul_pos = NULL;
        const PetscBdPointJac Jf3ul_pos = NULL;

        const PetscBdPointJac Jf0lu = pylith::fekernels::FaultCohesiveDyn::Jf0lu;
        const PetscBdPointJac Jf1lu = NULL;
        const PetscBdPointJac Jf2lu = NULL;
        const PetscBdPointJac Jf3lu = NULL;

        kernels.resize(3);
        const char* nameDisplacement = "displacement ";
        const char* nameLagrangeMultiplier = "lagrange_multiplier_fault ";
        kernels[0] = JacobianKernels(nameDisplacement, nameLagrangeMultiplier, integrator_t::LHS,
                                     integrator_t::NEGATIVE_FACE, Jf0ul_neg, Jf1ul_neg, Jf2ul_neg, Jf3ul_neg);
        kernels[1] = JacobianKernels(nameDisplacement, nameLagrangeMultiplier, integrator_t::LHS,
                                     integrator_t::POSITIVE_FACE, Jf0ul_pos, Jf1ul_pos, Jf2ul_pos, Jf3ul_pos);
        kernels[2] = JacobianKernels(nameLagrangeMultiplier, nameDisplacement, integrator_t::LHS,
                                     integrator_t::FAULT_FACE, Jf0lu, Jf1lu, Jf2lu, Jf3lu);
        break;
    } // QUASISTATIC
    case pylith::problems::Physics::DYNAMIC_IMEX: {
        const PetscBdPointJac Jf0lu = pylith::fekernels::FaultCohesiveDyn::Jf0lu;
        const PetscBdPointJac Jf1lu = NULL;
        const PetscBdPointJac Jf2lu = NULL;
        const PetscBdPointJac Jf3lu = NULL;

        kernels.resize(1);
        const char* nameDisplacement = "displacement ";
        const char* nameLagrangeMultiplier = "lagrange_multiplier_fault ";
        kernels[0] = JacobianKernels(nameLagrangeMultiplier, nameDisplacement, integrator_t::LHS,
                                     integrator_t::FAULT_FACE, Jf0lu, Jf1lu, Jf2lu, Jf3lu);
        break;
    } // DYNAMIC_IMEX
    case pylith::problems::Physics::DYNAMIC:
        PYLITH_COMPONENT_LOGICERROR("Fault implementation is incompatible with 'dynamic' formulation.Use 'dynamic_imex'.");
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown formulation '"<<_formulation<<"'.");
    } // switch

    assert(integrator);
    integrator->setKernels(kernels, solution, materials);

    PYLITH_METHOD_END;
} // _setKernelsJacobian


// End of file
