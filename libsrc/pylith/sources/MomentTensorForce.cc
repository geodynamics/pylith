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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "pylith/sources/MomentTensorForce.hh" // implementation of object methods

#include "pylith/sources/SourceTimeFunctionMomentTensorForce.hh" // HASA SourceTimeFunctionMomentTensorSource
#include "pylith/sources/AuxiliaryFactoryMomentTensorForce.hh" // USES AuxiliaryFactoryMomentTensorForce
#include "pylith/sources/DerivedFactoryMomentTensorForce.hh" // USES DerivedFactoryMomentTensorForce
#include "pylith/feassemble/IntegratorDomain.hh" // USES IntegratorDomain
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field::SubfieldInfo
#include "pylith/topology/FieldOps.hh" // USES FieldOps
#include "pylith/fekernels/DispVel.hh" // USES DispVel kernels

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <typeinfo> // USES typeid()

// ---------------------------------------------------------------------------------------------------------------------
typedef pylith::feassemble::IntegratorDomain::ResidualKernels ResidualKernels;
typedef pylith::feassemble::IntegratorDomain::JacobianKernels JacobianKernels;
typedef pylith::feassemble::IntegratorDomain::ProjectKernels ProjectKernels;
typedef pylith::feassemble::Integrator::EquationPart EquationPart;

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::sources::MomentTensorForce::MomentTensorForce(void) :
    _sourceTimeFunction(NULL),
    _derivedFactory(new pylith::sources::DerivedFactoryMomentTensorForce) {
    pylith::utils::PyreComponent::setName("momenttensorforce");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::sources::MomentTensorForce::~MomentTensorForce(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::sources::MomentTensorForce::deallocate(void) {
    Source::deallocate();

    delete _derivedFactory;_derivedFactory = NULL;
    _sourceTimeFunction = NULL;
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Set source time function.
void
pylith::sources::MomentTensorForce::setSourceTimeFunction(pylith::sources::SourceTimeFunctionMomentTensorForce* const sourceTimeFunction) {
    _sourceTimeFunction = sourceTimeFunction;
} // setSourceTimeFunction


// ---------------------------------------------------------------------------------------------------------------------
// Get bulk source time function.
pylith::sources::SourceTimeFunctionMomentTensorForce*
pylith::sources::MomentTensorForce::getSourceTimeFunction(void) const {
    return _sourceTimeFunction;
} // getSourceTimeFunction


// ---------------------------------------------------------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::sources::MomentTensorForce::verifyConfiguration(const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("verifyConfiguration(solution="<<solution.getLabel()<<")");

    // Verify solution contains expected fields.
    if (!solution.hasSubfield("velocity")) {
        throw std::runtime_error("Cannot find 'velocity' field in solution; required for 'MomentTensorForce'.");
    } // if

    PYLITH_METHOD_END;
} // verifyConfiguration


// ---------------------------------------------------------------------------------------------------------------------
// Create integrator and set kernels.
pylith::feassemble::Integrator*
pylith::sources::MomentTensorForce::createIntegrator(const pylith::topology::Field &solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createIntegrator(solution="<<solution.getLabel()<<")");

    pylith::sources::Source::locateSource(solution);

    pylith::feassemble::IntegratorDomain* integrator = new pylith::feassemble::IntegratorDomain(this);assert(integrator);
    integrator->setLabelName(getLabelName());
    integrator->setLabelValue(getLabelValue());

    _setKernelsResidual(integrator, solution);
    // _setKernelsJacobian(integrator, *this, solution);

    _setKernelsUpdateStateVars(integrator, solution);
    _setKernelsDerivedField(integrator, solution);

    PYLITH_METHOD_RETURN(integrator);
} // createIntegrator


// ---------------------------------------------------------------------------------------------------------------------
// Create auxiliary field.
pylith::topology::Field*
pylith::sources::MomentTensorForce::createAuxiliaryField(const pylith::topology::Field& solution,
                                                         const pylith::topology::Mesh& domainMesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createAuxiliaryField(solution="<<solution.getLabel()<<", domainMesh="<<typeid(domainMesh).name()<<")");

    pylith::topology::Field* auxiliaryField = new pylith::topology::Field(domainMesh);assert(auxiliaryField);
    auxiliaryField->setLabel("MomentTensorForce auxiliary field");

    assert(_sourceTimeFunction);
    pylith::sources::AuxiliaryFactoryMomentTensorForce* auxiliaryFactory = _sourceTimeFunction->getAuxiliaryFactory();assert(auxiliaryFactory);

    assert(_normalizer);
    auxiliaryFactory->initialize(auxiliaryField, *_normalizer, domainMesh.getDimension());

    // :ATTENTION: The order for adding subfields must match the order of the auxiliary fields in the FE kernels.

    auxiliaryFactory->addMomentTensor(); // 0
    auxiliaryFactory->addTimeDelay(); // 1

    _sourceTimeFunction->addAuxiliarySubfields();

    auxiliaryField->subfieldsSetup();
    auxiliaryField->createDiscretization();
    pylith::topology::FieldOps::checkDiscretization(solution, *auxiliaryField);
    auxiliaryField->allocate();
    auxiliaryField->createOutputVector();

    assert(auxiliaryFactory);
    auxiliaryFactory->setValuesFromDB();

    // Debug option
    auxiliaryField->view("MomentTensor auxiliary field.");

    PYLITH_METHOD_RETURN(auxiliaryField);
} // createAuxiliaryField


// ---------------------------------------------------------------------------------------------------------------------
// Create derived field.
pylith::topology::Field*
pylith::sources::MomentTensorForce::createDerivedField(const pylith::topology::Field& solution,
                                                       const pylith::topology::Mesh& domainMesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createDerivedField(solution="<<solution.getLabel()<<", domainMesh=)"<<typeid(domainMesh).name()<<") empty method");

    PYLITH_METHOD_RETURN(NULL);
} // createDerivedField


// ---------------------------------------------------------------------------------------------------------------------
// Get auxiliary factory associated with physics.
pylith::feassemble::AuxiliaryFactory*
pylith::sources::MomentTensorForce::_getAuxiliaryFactory(void) {
    assert(_sourceTimeFunction);
    return _sourceTimeFunction->getAuxiliaryFactory();
} // _getAuxiliaryFactory


// ---------------------------------------------------------------------------------------------------------------------
// Update kernel constants.
void
pylith::sources::MomentTensorForce::_updateKernelConstants(const PylithReal dt) {
    assert(_sourceTimeFunction);
    _sourceTimeFunction->updateKernelConstants(&_kernelConstants, dt);
} // _updateKernelConstants


// ---------------------------------------------------------------------------------------------------------------------
// Get derived factory associated with physics.
pylith::topology::FieldFactory*
pylith::sources::MomentTensorForce::_getDerivedFactory(void) {
    return _derivedFactory;
} // _getDerivedFactory


// ---------------------------------------------------------------------------------------------------------------------
// Set kernels for RHS residual G(t,s).
void
pylith::sources::MomentTensorForce::_setKernelsResidual(pylith::feassemble::IntegratorDomain* integrator,
                                                        const topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelsResidual(integrator="<<integrator<<", solution="<<solution.getLabel()<<")");

    const spatialdata::geocoords::CoordSys* coordsys = solution.getMesh().getCoordSys();

    std::vector<ResidualKernels> kernels(1);
    switch (_formulation) {
    case QUASISTATIC: {
        break;
    } // QUASISTATIC
    case DYNAMIC_IMEX:
    case DYNAMIC: {
        // Velocity
        const PetscPointFunc g0v = NULL;
        const PetscPointFunc g1v = _sourceTimeFunction->getKernelg1v_explicit(coordsys);
        kernels[0] = ResidualKernels(getSubfieldName(),  pylith::feassemble::Integrator::RHS, g0v, g1v);
        // kernels[0] = ResidualKernels("velocity",  pylith::feassemble::Integrator::RHS, g0v, g1v);
        break;
    } // DYNAMIC
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown formulation for equations (" << _formulation << ").");
    } // switch

    assert(integrator);
    integrator->setKernelsResidual(kernels, solution);

    PYLITH_METHOD_END;
} // _setKernelsResidual


// //
// ---------------------------------------------------------------------------------------------------------------------
// // Set kernels for RHS Jacobian G(t,s).
// void
// pylith::sources::MomentTensorForce::_setKernelsJacobian(pylith::feassemble::IntegratorDomain* integrator,
//                                                         const pylith::sources::MomentTensorForce& source,
//                                                         const topology::Field& solution) const {
//     PYLITH_METHOD_BEGIN;
//     PYLITH_COMPONENT_DEBUG("_setKernelsJacobian(integrator="<<integrator<<", solution="<<solution.getLabel()<<")");
//     assert(integrator);

//     const spatialdata::geocoords::CoordSys* coordsys = solution.getMesh().getCoordSys();

//     std::vector<JacobianKernels> kernels;

//     switch (_formulation) {
//     case QUASISTATIC: {
//         break;
//     } // QUASISTATIC
//     case DYNAMIC_IMEX:
//     case DYNAMIC: {
//         typedef pylith::feassemble::JacobianValues::JacobianKernel ValueKernel;
//         std::vector<ValueKernel> valueKernelsJacobian(2);
//         std::vector<ValueKernel> valueKernelsPrecond;
//         valueKernelsJacobian[0] = ValueKernel("displacement", "displacement",
// pylith::feassemble::JacobianValues::blockDiag_tshift);
//         valueKernelsJacobian[1] = ValueKernel("velocity", "velocity",
// pylith::feassemble::JacobianValues::blockDiag_tshift);
//         integrator->setKernelsJacobian(valueKernelsJacobian, valueKernelsPrecond);

//         const PetscPointJac Jf0uu = pylith::fekernels::DispVel::Jg0uv;
//         const PetscPointJac Jf1uu = NULL;
//         const PetscPointJac Jf2uu = NULL;
//         const PetscPointJac Jf3uu = NULL;

//         const PetscPointJac Jf0vv = pylith::fekernels::DispVel::Jg0uv;
//         const PetscPointJac Jf1vv = NULL;
//         const PetscPointJac Jf2vv = NULL;
//         const PetscPointJac Jf3vv = NULL;

//         integrator->setLHSJacobianTriggers(pylith::feassemble::Integrator::NEW_JACOBIAN_TIME_STEP_CHANGE);

//         kernels.resize(2);
//         const EquationPart equationPart = pylith::feassemble::Integrator::LHS_LUMPED_INV;
//         kernels[0] = JacobianKernels("displacement", "displacement", equationPart, Jf0uu, Jf1uu, Jf2uu, Jf3uu);
//         kernels[1] = JacobianKernels("velocity", "velocity", equationPart, Jf0vv, Jf1vv, Jf2vv, Jf3vv);
//         break;
//     } // DYNAMIC
//     default:
//         PYLITH_COMPONENT_LOGICERROR("Unknown formulation for equations (" << _formulation << ").");
//     } // switch

//     integrator->setKernelsJacobian(kernels, solution);

//     PYLITH_METHOD_END;

// } // _setKernelsJacobian

// ---------------------------------------------------------------------------------------------------------------------
// Set kernels for computing updated state variables in auxiliary field.
void
pylith::sources::MomentTensorForce::_setKernelsUpdateStateVars(pylith::feassemble::IntegratorDomain* integrator,
                                                               const topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelsUpdateStateVars(integrator="<<integrator<<", solution="<<solution.getLabel()<<")");

    const spatialdata::geocoords::CoordSys* coordsys = solution.getMesh().getCoordSys();
    assert(coordsys);

    std::vector<ProjectKernels> kernels;
    _sourceTimeFunction->addKernelsUpdateStateVars(&kernels, coordsys);

    integrator->setKernelsUpdateStateVars(kernels);

    PYLITH_METHOD_END;
} // _setKernelsUpdateStateVars


// ---------------------------------------------------------------------------------------------------------------------
// Set kernels for computing derived field.
void
pylith::sources::MomentTensorForce::_setKernelsDerivedField(pylith::feassemble::IntegratorDomain* integrator,
                                                            const topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelsDerivedField(integrator="<<integrator<<", solution="<<solution.getLabel()<<")");

    // Default is to do nothing.

    PYLITH_METHOD_END;
} // _setKernelsDerivedField


// End of file
